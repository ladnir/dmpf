# AES and reusable DPF implementation audit

2026-09-07. Scope: the concrete DPF tree generators, correction encoding,
dense-to-sparse boundary, and cached leaf conversion. This is an implementation
transfer audit, not a change to the approved chosen-output-share functionality.

Inspected worktree: `C:/Users/peter/repo/libOTe-codex-rev-analysis`.
The original audit was read-only for C++ and TeX. Peter then approved fixing
the correction-bit disclosure. The implementation now erases that bit before
transmission, repacks corrected tree tags, and exports cached seed prefixes
with low bit zero. Its immediate-payload API hashes those prefixes before
using them as full-block masks. No extra OTs or cached-expansion AES calls
are required. The dense generator/export and reusable key schedule remain
unchanged. See AUD-213 and `analysis/sparse_dpf_correction_encoding.md` in
the worktree. No paper TeX or functionality change was needed.
All 17 focused WSL/GCC tests and nine exact Python encoding tests passed
for that correction-bit fix. The later leaf-rekey pass is recorded below.

**Latest update: internal-node rekeying (2026-09-07).** Peter subsequently
approved a counter that rekeys after 1,024 parent expansions. This is now
implemented in RegularDPF and SparseDPF. Sparse traversal advances the counter
across trees and levels, including dense-to-sparse boundary expansions.
Dense recomputation restores each level's counter range, so replay and
evaluation select the same keys despite different traversal orders.
Purpose-separated roots reuse DMPF public setup coins; no DMPF messages or
OTs are added. Standalone calls without supplied coins exchange 16 bytes in
each direction. RegularDPF serialized keys grow by 16 bytes per key object
and are incompatible with the old format. Both peers must update together.
See `analysis/internal_node_rekey.md` in the implementation worktree for
the precise schedule and counter convention. All 23 focused WSL/GCC tests
passed with default settings and with `-domain 4097`; the aggregate test
registration object also compiled. This is implementation hardening, not a
proof of the concrete generator. Earlier schedule descriptions below record
the design history and are superseded by the implementation note.
Two sequential ABBA timing passes checked 432 full expansions and unchanged
communication. Setup variation prevents a reliable small-overhead estimate.
Repeated expansion was roughly flat for Waterfall 4N, but 3--7% slower for
Waterfall 3N and about 5% slower for Rev-Cuckoo 3N. The leaf source kernel is
unchanged; the timing increases remain unexplained. The implementation note
records both passes without claiming zero overhead.

**Timing follow-up:** Those larger slowdowns do not reproduce with guest-CPU
affinity and explicit warmup: the same frozen binaries reverse the sign of
the difference. A one-executable generator comparison leaves 1--2%
differences. The cached-leaf instruction sequences are unchanged after address
normalization. Execution conditions are the leading explanation, without an
isolated physical cause or a claimed speedup. See
`analysis/internal_node_timing_investigation.md` and its raw data in the
implementation worktree. No production changes were needed.

**Previous update:** Peter approved rekeying cached-leaf conversion first,
without changing internal-node expansion. Rev-Cuckoo and both Waterfall
cached-leaf paths now derive a key for each 1,024 leaves in public flattened
row order. Purpose roots come from the existing public setup seed and advance
per expansion call. No OTs, messages, or rounds are added. See
`analysis/cached_leaf_rekey.md` in the implementation worktree for the exact
schedule, scope, and validation. Earlier statements that the leaf schedule
is unchanged describe the correction-bit fix, not this later implementation.
Final validation: all 18 focused WSL/GCC tests passed with default parameters
and with `-domain 1024`. Sequential synthetic kernel comparisons found no
remaining slowdown after explicit eight-block AES batching. See the
implementation note for timings; no end-to-end performance claim is implied.

Sections 1--2 below record the pre-fix evidence and repair rationale;
source line numbers refer to that snapshot, not the subsequently edited file.

## 1. Confirmed sparse correction disclosure

`libOTe/Dpf/SparseDpf.h:831` constructs the tag corrections. Lines 843--846
select and open a full 128-bit seed correction. The overload at line 955
sends all 128 bits together with both tag corrections; it does not erase
or replace the seed correction's low bit.

Fix one processed level. Let a in {0,1} be the reconstructed index bit,
and let z0,z1 in {0,1}^128 be the reconstructed child sums after activity
masking. Write b0=lsb(z0) and b1=lsb(z1). The public opening is

```text
sigma = z[1-a]
tau0  = b0 xor a xor 1
tau1  = b1 xor a
```

Consequently,

```text
a=0 implies lsb(sigma)=tau1
a=1 implies lsb(sigma)=tau0.
```

Whenever tau0 differs from tau1, the public opening determines a exactly.
This implication holds for every generator and every choice of masks.
It is not an AES cryptanalysis claim.

For the directly sampled root, fix either corrupt party's root shares.
The honest root shares make z0,z1 independent uniform blocks under the
project's ideal-randomness model. Then tau0 differs from tau1 with probability
1/2. On the other half of the openings the two possible index bits have
identical opening distributions. Thus the opening distributions for a=0
and a=1 have statistical distance 1/2. With an equiprobable challenge bit,
the optimal classifier succeeds with probability 3/4 (advantage 1/4 over
guessing). This already gives a non-negligible setup distinguisher for a
two-point sparse set with no dense prefix.

The same exact law holds for an inactive level with fresh independent
activity masks. The mask fix therefore does not remove this encoding leak.
For production dense-prefix configurations, the identity still holds at
every opened sparse level; this pass does not claim a measured frequency
or an end-to-end production support-recovery attack.

Evidence: `experiments/test_sparse_dpf_correction_encoding.py` exhaustively
checks block sizes 1--4, the statistical distance, and the existing paper
encoding. All nine tests pass. These are exact checks of the inspected
opening equations, not C++ socket tests or AES experiments.

## 2. Representation repair required before theorem transfer

The paper opens only a 127-bit prefix and two tag corrections:

```text
prefix = msbs(z[1-a])
CW[0]  = prefix || tau0
CW[1]  = prefix || tau1.
```

For independent uniform child sums, this opening is uniform and independent
of a. `RegularDpf.h:628` clears the selection input's low bit and packs the
tag correction into that position. Its opening at lines 640--656 therefore
already has the paper's format. The full-low-bit disclosure is specific to
the inspected sparse path, not the dense correction encoding.

The sparse repair must also change how corrected seeds are represented.
Clearing the transmitted bit alone is not sufficient: the current sparse
generator consumes a raw corrected block, while its logical tag is computed
separately from the uncorrected child and tau.

The least ambiguous representation is a packed corrected seed s whose low
bit is its logical tag. Apply CW[j] using the parent's logical tag, then
expand s. If separate storage is preferred, store the same 127-bit prefix
and logical tag separately and repack them before tree expansion. Cached
leaf conversion should use the hidden prefix specified by the abstract
proof, rather than silently treating the public logical tag as entropy.

The exact test also checks that the packed and explicitly split
representations agree and that off-path children become equal in both seed
and tag. This repair needs scalar, eight-lane, leaf, and dense-boundary
coverage. It requires no extra OT or interaction; actual runtime is not
benchmarked here. Do not deploy a wire-format change to only one peer.

## 3. Sparse fixed-key AES generator

`SparseDpf.h:472` and lines 508--509 implement

```text
H_k(x) = E_k(x) xor x
G(s)   = (H_k(s), H_k(s xor 1)),
```

where k is the public fixed AES key and 1 is the block with only its low bit
set. `cryptoTools/Crypto/AES.h:720` implements the feed-forward XOR.

There is a direct ideal-permutation starting point. For a fixed input s,
the two cipher inputs are distinct. Before any cipher queries or auxiliary
exposure, the pair G(s) is uniform over ordered block pairs (u,v) satisfying
u xor v != 1. Its statistical distance from two independent uniform blocks
is exactly 2^-128. This distribution is the same for either fixed input tag.
Unlike the custom dense evaluator, the sparse generator really does make
two distinct full-cipher calls.

This local observation is not yet the full public-oracle reduction. That
reduction must charge forward queries at hidden inputs, inverse queries at
hidden outputs, and overlaps among parent input pairs across the complete
execution. Adjacent inputs share a generator pair in reversed order:
G(s xor 1) is the swap of G(s). The proof must therefore count collisions
of 127-bit prefixes, not just exact 128-bit seed collisions. Establishing
the canonical seed/tag representation is a prerequisite.

## 4. Dense generator and dense-to-sparse export are distinct obligations

The dense tree computes, for a corrected seed s,

```text
t       = E_k(s)
child0  = AESRound(t, s)
child1  = add64x2(t, s).
```

Here AESRound is one public AES round including its round-key XOR, and
add64x2 adds the two 64-bit lanes separately modulo 2^64. See
`RegularDpf.h:452` and the eight-lane equivalent at line 530. These children
are not two independent ideal-cipher evaluations. Replacing E_k by an ideal
cipher does not license replacing this output pair by independent blocks;
the postprocessing still needs its own argument. No efficient attack on this
dense generator is established by this audit.

There is another transformation at the boundary. Even with an empty payload
vector, `RegularDpf::expand` calls the ordinary evaluator (line 380).
The evaluator exports

```text
value = AESRound(correctedLeaf, correctedLeaf)
tag   = lsb(correctedLeaf),
```

using `CoeffCtxInteger` with block output in the sparse caller. See
`RegularDpf.h:988,1077` and `SparseDpf.h:703`. The caller saves this value
as a seed, and either uses it directly as a cached leaf or expands it into
the sparse suffix. It is not the corrected tree seed itself. This is not
an off-path correctness objection: equal inputs still produce equal values.
It is a missing entropy/distribution argument at the proof boundary.

Two possible directions require Peter's choice before implementation:

- Preserve these optimizations and prove the precise custom generator and
  export properties. Do not assume an AES round is an independent oracle.
- Use the two-call fixed-key generator for dense levels too, and provide
  an explicit corrected-seed export for sparse/cached use. Preserve the
  current batched memory layout. This changes dense setup work and needs
  measurement before making performance claims.

## 5. Cached leaf conversion: public keys are not themselves a problem

`CachedDpfExpansion.h:46` constructs AES under a public state k_q, advances
that state using a fixed public block c, and hashes saved leaf inputs:

```text
k_(q+1) = E_(k_q)(c) xor c
L(q,s)  = convert_G(E_(k_q)(s) xor s).
```

Waterfall's row-to-column path uses the same pattern with separate initial
key and step constant (`WaterfallDmpf.h:769`). Column expansion uses the
shared cached helper (line 894). Initial keys and constants are fixed and
public. Different states can use the same round key; their saved inputs
must then be handled jointly rather than assigning them independent oracles.

The hidden saved input, not the public AES key, supplies the secrecy.
Thus the key schedule's publicity alone does not imply a reuse attack.
However, this is not ordinary AES secret-key PRF security. The appropriate
ideal-cipher proof must expose the public schedule and all known-input
evaluations, and handle both forward and inverse queries.

A complete proof must count repeated schedule keys, overlap with the fixed
tree key or other key streams, saved-input collisions, and queries that hit
hidden cipher points. It must include schedule evaluation at c under the
same key used for leaf hashing. Rounds with repeated keys cannot simply be
declared fresh. Reuse across multiple trees and states belongs in the same
experiment. No numeric bound or full instantiation theorem is claimed here.

For G=Z_(2^64), converting a uniform block by taking its low 64 bits is
exactly uniform. Other coefficient contexts need their actual conversion
law; a generic `fromBlock` name does not establish uniformity. The present
audit does not certify every supported context.

## Next action

The sparse correction encoding and canonical seed/tag repair are implemented,
with a transcript regression in the C++ path. Next settle the dense generator
and export choice before completing a joint ideal-cipher proof of setup and
cached expansion. The approved functionality and the paper's abstract
correction format do not need changing for this defect.

## Literature follow-up and proposed outer key schedule

Peter identified the XOR/addition precedent in *Correlated Pseudorandomness
from Expand-Accumulate Codes*, https://eprint.iacr.org/2022/1014.pdf.
The inspected full version has 59 pages; the relevant discussion is Section
6.2, pages 47--48, Conjecture 6.4. It proposes children P(x) xor x and
P(x)+x modulo 2^lambda, where P is a public invertible random permutation.
The authors give partial analysis and conjecture unpredictability of the
punctured function. They explicitly identify carry-based relations that
prevent their stronger unpredictability property. Section 6.3 supplies a
separate hashing transformation for converting UPFs to PPRFs; the
random-permutation instantiation remains conditional on the conjecture.

This supports the motivation for Peter's dense generator, but does not prove
the DPF simulation requirement. The implementation adds an AES round to the
XOR branch and uses two independent 64-bit additions rather than one 128-bit
addition. The extra round disrupts the specific direct bit/carry relations
of the unmodified candidate. This is a hardening rationale, not a reduction
showing that every attack transfers to the original construction. In
particular, DPF setup openings need their own joint distribution argument;
hashing final leaves cannot erase an earlier setup disclosure.

Peter also proposed an outer key schedule to limit collision effects.
Treat the following as a design proposal, not implemented or proved:

- Agree one public setup identifier/seed during existing setup messages.
- Derive public AES keys indexed by setup, purpose, tree batch, and level.
  Leaf conversion uses its own purpose and the expansion counter instead
  of the tree level. Both parties use the same context at corresponding
  nodes; no label depends on the hidden active path.
- Keep a key constant across the existing SIMD batch. Dense evaluation
  interleaves trees, so per-tree keys could force a different multi-key
  kernel. Batch-level keys retain the current batching structure.
- If needed, refine a scope by public subtree/chunk boundaries. Rekeying
  every individual node is not the initial performance-conscious choice.
- A PRP applied to injectively encoded local context counters can make
  derived keys distinct within one setup. Cross-setup collisions and overlap
  with the outer derivation key still belong in the ideal-cipher accounting.

For intuition only, suppose a birthday-style term counts n_c independent
b-bit seed values within key scope c. Partitioning the evaluations changes
the relevant same-scope collision sum from binom(sum_c n_c,2)/2^b to
sum_c binom(n_c,2)/2^b. Equal-sized scopes remove the cross-scope pairs.
The complete DPF bound must also count oracle hits and other dependencies;
this formula is not a claimed bound for the custom dense generator.

Key separation limits where equal seeds are evaluated identically; it does
not increase their entropy. Collisions within one scope remain. In
particular, a hidden active prefix matching the corrupt prefix at the same
node cannot be separated by public context labels while preserving off-path
cancellation. Changing only leaf-round keys cannot repair collisions that
already merged states during setup. The schedule also does not remove
algebraic relations within the two children of one generator call.

Recommended next step: retain the current dense candidate while analyzing
its specific one-round hardening, and specify a public batch/level key
schedule alongside that analysis. No generator or schedule code was changed
in this literature follow-up.

Peter subsequently approved describing the dense optimization as conjectural
and citing the EA paper as related work. That note is now in Appendix C and
the `RegularDpf` class comment. He selected a rekeying interval of **1,024
node expansions**, explicitly not 1,024 complete trees. This schedule is not
implemented yet. Derive key selection from stable public node coordinates,
not an execution counter: compact key generation recomputes earlier nodes,
and evaluation can traverse them in a different order. A key schedule should
be prepared outside each 1,024-node inner loop and reused by the SIMD kernel.
Context derivation must also distinguish setups; resetting the node counter
with the same global keys would not isolate separate executions. Treat
near-zero overhead as an implementation target until measured.
