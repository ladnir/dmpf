# Sol-to-Astra handoff: distributed multi-point functions

**Project:** *Fully Distributed Multi-Point Functions for PCGs and Beyond*  
**Authors:** Amit Agarwal, Srinivasan Raghuraman, and Peter Rindal  
**Snapshot date:** 2026-09-05  
**Primary paper repository:** `C:\Users\peter\repo\dmpf`  
**Packaged implementation worktree:** `C:\Users\peter\repo\libOTe-codex-rev-analysis`  

This is an internal operational handoff. It records the design history,
motivation, current evidence, implementation state, unresolved questions, and
audit targets that are not all present in the paper. It is not a substitute
for reading the current PDF and the cited experiment notes.

Reviewed and corrected on 2026-09-05. The subsequent hash-correctness pass
also updated the paper and added exact analysis/tests in the package worktree;
it did not change the C++ implementation or parameters.
The roadmap is a recommendation. Peter asked to be consulted before larger
or questionable changes, including protocol and parameter changes.

## Latest checkpoint: approved CWC corrections

**2026-09-07 prose deduplication:** Shortened repeated roadmaps in the
introduction, overview, and DMPF section; compressed repeated DPF, Waterfall,
PCG, and Reverse-Cuckoo scope explanations. Historical benchmark limitations
are collected at the start of Evaluation, with local qualifications retained
where they limit individual claims. This was a prose-only CWC pass: all 29
figure, table, theorem/lemma/assumption, and proof blocks in the eight edited
files were compared with their pre-pass text and are unchanged. No protocol,
parameter, performance number, or implementation changed. The 72-page PDF
builds with stable references and no overfull boxes; affected transitions,
the Waterfall protocol and cost table, and Evaluation were rendered and
checked. Figure 10 is now on page 23 and Figure 11 on page 26. Next: review
the remaining proof/assumption boundaries, especially the application-level
correlation-distribution argument; do not treat prose cleanup as closing it.

**2026-09-07 named links:** All 116 raw named-reference occurrences across
the root TeX sources now use the shared named-reference macros, including
the retained inactive Reverse-Cuckoo source. Names and numbers are inside one
link; equation references preserve parentheses via `namedeqref`/`equationref`.
Loading hyperref and cleveref after the float packages removed the duplicate
figure/table destination warnings. All 28 figure/table label entries agree
with the PDF destination pages; full-phrase links for Figures 10 and 14 were
checked in exported PDF link data. There are 617 unique named destinations.
The 72-page build has no overfull boxes or undefined references. Rendered
pages 23--25, 28--29, and 67 were checked. Float pagination shifted: Figure 10
is now on page 24 and Figure 11 on page 27. Next: resume the substantive
paper review; the duplicate-destination follow-up is closed.

**2026-09-07 classical-cuckoo cleanup:** The older worked example and rank
formula are repaired; all four displayed matrices and their sparse sets were
checked against the TeX. The section now defines objects before using the
figure and separates sparse expansion, address translation, and placement.
It distinguishes at-most-wN general leaf work from exact wN for distinct
candidates. The protected insertion paragraph is unchanged. Matrix spacing
also removes the old 11.8429pt overfull box. The 72-page PDF builds with stable
references; pages 19--22 were rendered and checked. No implementation or
protocol changes. Next: the duplicate PDF figure/table destination warnings.

**2026-09-07 narrative pass:** The front matter now introduces the reusable
setup/expand task before routing machinery, explains placement-count leakage
in plain language, and qualifies the historical measurements for AES rekeying
as well as mask refresh and support filtering. The overview now introduces
Reverse Cuckoo after its Ring-LPN motivation. The DMPF roadmap is shorter;
the old baseline's "no satisfactory MPC insertion" conclusion now points to
Waterfall's repair circuit. Global optimality wording for 4N and the misleading
K=1-only 2N rationale were removed. No protocols, theorem statements, numerical
results, or implementation changed. The 72-page build has stable references;
12 changed/transition pages were rendered and checked, with no new overfull
boxes. Existing duplicate PDF figure/table destinations remain. Next: repair
the classical-cuckoo section's older worked-example and address-translation
indexing errors; details are in `PAPER_NARRATIVE_AUDIT.md`. Ask before moving
the long subroutine and auxiliary-construction sections to appendices.

**2026-09-07 attribution clarification:** Peter clarified that the dense AES
evaluator is an older libOTe implementation choice, not a new generator
design introduced by this paper. Appendix C.4 now says so, cites libOTe,
and removes the misleading "our extra AES round" / "our variant" wording.
The security caveat remains; the public rekeying schedules retain the
existing child formulas. Keep this provenance clear in future revisions.

**2026-09-07 paper synchronization:** Appendix C.4 now describes the concrete
dense and sparse AES evaluators, dense seed export, and both 1,024-position
rekeying schedules. It defines the public root/counter encodings and explains
dense replay, sparse traversal, and deterministic per-expansion root updates.
The concrete leaf evaluator also uses public leaf-position context; the
abstract figures and Theorem C.3 remain unchanged. Rekeying is a collision-
mitigation rationale, not a proved instantiation of that theorem. Section 4
links to C.4. Section 13 now explicitly says its original timing tables
predate rekeying; no numerical table entries or performance claims were
replaced with the recent local diagnostics. The PDF builds at 72 pages with
stable references and only the existing 11.8429pt overfull box. Pages 8, 49,
50, 71, and 72 were rendered and inspected. All 53 targeted correction,
output-boundary, and sparse-schedule tests pass. Copies are at `main.pdf`
and `output/pdf/main.pdf`. Next: a final narrative/claim consistency pass,
keeping Waterfall generic and Reverse Cuckoo explicitly LPN-specialized.

**2026-09-07 timing investigation:** The apparent 5--7% online slowdown after
internal rekeying does not reproduce under tighter controls. The same frozen
binaries, with guest-CPU affinity and explicit warmup, measured slightly
faster after rekeying for all three profiles; these are not claimed speedups.
A one-executable fixed/rekeyed generator comparison leaves 1--2% differences.
The cached-leaf instruction sequences match after address normalization;
correction/scatter accounts for roughly 60--70% of online time. Identical-cache
buffer-offset replay also shows few-percent timing variation. Execution
conditions are the leading explanation, but frequency, scheduling, and layout
were not individually isolated. See `analysis/internal_node_timing_investigation.md`
and its raw CSVs in the package worktree. No production changes were made.
Recommendation: retain rekeying, use controlled benchmark conditions, and
return to synchronizing the paper's concrete AES description. The earlier
slowdown figures below are historical observations, not established overhead.

**2026-09-07 internal-node rekey implementation:** Peter approved the simple
counter policy. `DpfTreeHash.h` now derives and caches one public AES key per
1,024 parent expansions. Sparse traversal carries the counter across trees,
levels, and dense-to-sparse boundary roots. Dense kernels restore a padded
per-tree, per-level counter range on replay; the cleartext generator maps its
active prefix into the same physical buffer order. The child formulas and
eight-lane kernels are retained. AES schedules remain outside coroutine
frames. DMPF setup supplies purpose-separated seeds from existing public
coins, adding no messages or OTs. Standalone DPF calls without supplied seeds
exchange 16 bytes each way. RegularDPF key serialization gains a 16-byte
public seed per key object and is not backward compatible. All 23 focused
WSL/GCC tests pass with defaults and `-domain 4097`; the aggregate test
registration object compiles. The audit's stale empty-set rejection fixture
was replaced by a valid zero-output check. See
`analysis/internal_node_rekey.md` in the package worktree. The cached-leaf
execution path is unchanged in this pass. No TeX, functionality, or abstract
security theorem was changed. Rekeying remains hardening, not a generator
security proof. Next: synchronize the paper's concrete implementation note
with both schedules, without overstating the AES assumption.
Two sequential ABBA benchmark passes checked 432 complete expansions.
Setup timings were too variable to isolate a small overhead. Repeated
expansion was approximately flat for Waterfall 4N, 3--7% slower for Waterfall
3N, and about 5% slower for Rev-Cuckoo 3N. The leaf source kernel is unchanged,
but the recurring increases have not been explained; do not call them proven
noise or claim zero end-to-end overhead. Communication was unchanged. Raw
timings and fixture limitations are in `analysis/internal_node_rekey.md`
and `analysis/internal_node_e2e_results.csv`. A same-cache, stage-level
comparison is the next performance diagnostic if needed.

**2026-09-07 cached-leaf rekey implementation:** Peter approved doing only
leaf rekeying first. The common cached-leaf kernel now rekeys after each
1,024 leaves across public row boundaries. Existing public setup coins derive
purpose-separated roots for Rev-Cuckoo, Waterfall columns, and Waterfall
scatter. Roots advance per expansion call. Waterfall scatter uses the same
kernel. The explicit eight-lane loops remain; AES schedules live in a
non-inlined non-coroutine helper. Internal-node hashing and the immediate
SparseDPF API are unchanged. No OTs/messages/rounds are added. Both peers must
update together. See `analysis/cached_leaf_rekey.md` and its test/benchmark
drivers in the package worktree. The older "not implemented" checkpoint below
is historical. Rekeying is not a proof of the concrete AES construction.
All 18 focused WSL/GCC tests passed after the final edit, with both defaults
and `-domain 1024`. A sequential local microbenchmark initially found a
sign-negation regression; explicit `hashBlocks<8>` batching removed it.
The final combined change reduced synthetic leaf-kernel time by about
7--35% versus the old kernel, depending on party/layout. This is not an
end-to-end speed claim. Exact settings/results are in the implementation
note. Next: reflect the public leaf-position context in the paper's concrete
implementation note; internal-node rekeying remains deferred.

**2026-09-07 conjecture note and rekey granularity:** Peter approved a concise
paper/code comment that the dense AES optimization is conjecturally secure,
with EA (2022/1014, Section 6.2, Conjecture 6.4) as related work, not a proof.
The note is in Appendix C and `RegularDpf.h`. He clarified that rekeying
should occur after **1,024 node expansions**, not after 1,024 complete DPF
trees. This is the selected next implementation policy, not implemented in
this documentation pass. Use stable public node coordinates so recomputation
and different traversal orders select the same keys; amortize schedule setup
outside SIMD chunks. Separate setup instances too. Do not claim measured
near-zero overhead before benchmarking. No evaluator arithmetic changed.
Validation: `latexmk` stabilized at 71 pages with no undefined references
and no new overfull boxes (the existing 11.8429pt box remains). Rendered
and inspected the citation, proof transition, and final conjecture note;
the disclaimer fits on the final page.

**2026-09-07 literature follow-up:** Peter recalled the XOR/addition PCG
candidate and proposed an outer key schedule. The reference is the full
version of *Correlated Pseudorandomness from Expand-Accumulate Codes*,
ePrint 2022/1014, Section 6.2, pages 47--48, Conjecture 6.4. Its candidate
has a partial analysis and an unpredictability conjecture, not a full PRG
theorem. Peter's extra AES round disrupts the cited simple carry relations,
but the exact dense DPF generator still needs its own argument. The new
section in `AES_DPF_IMPLEMENTATION_AUDIT.md` records the comparison and a
public setup/purpose/batch/level key-schedule proposal, including SIMD
constraints and residual same-scope collisions. No schedule or generator
implementation changes were made in this follow-up.

**2026-09-07 approved correction-bit fix:** Peter approved masking out the
selected correction's low bit and then discussing the AES issues. Sparse
setup now does that before serialization, applies packed tag corrections in
both kernels, and exports zero-low-bit cached prefixes. The immediate-payload
API hashes prefixes to avoid an unmasked payload bit. The singleton cached
path was also corrected to avoid unallocated payload buffers. No extra OTs
or cached-expansion AES calls were added. See AUD-213 and
`analysis/sparse_dpf_correction_encoding.md` in the worktree. The custom dense
generator, transformed dense export, and reusable public key schedule are
unchanged and remain the next proof topics. The paper's abstract protocol
already had the corrected encoding; no TeX/functionality change was needed.
All 17 focused WSL/GCC tests and nine exact Python encoding checks passed.
No new benchmark or MSVC build was run. Changes remain uncommitted alongside
the earlier activity-mask package changes; preserve unrelated hash edits.

**2026-09-07 AES/DPF audit:** Peter requested the AES and DPF implementation
details only. `AES_DPF_IMPLEMENTATION_AUDIT.md` documents a newly confirmed
sparse correction-encoding leak, independent of AES: opening the full selected
block plus both tag corrections reveals the index bit whenever the tag
corrections differ. At a directly sampled root the reveal probability is
1/2 and the two bit-conditioned opening laws have statistical distance 1/2.
The activity mask does not remove this defect. The paper and dense correction
path already use the safe packed format (127-bit prefix plus two tags).
Nine exact equation-level tests pass; no production C++ or TeX was changed.
Next obtain approval to repair sparse encoding and canonical seed/tag use,
then decide whether to justify or replace the custom dense generator and
single-round transformed seed export. The note separately records the actual
public-key cached-leaf family and its remaining ideal-cipher obligations.

**2026-09-07 paper cost reconciliation:** the Waterfall tables and executable
cost model now include two OTs per shared DPF correction selection, the
secret-activity masks below the implementation's dense prefix, and cached
sign-conversion correlations for non-characteristic-two output groups.
For t=16, N=2^20, kappa=128, and G=Z_(2^64), the default 4N profile costs
39,052 OTs and 942,888 logical payload bits; the 3N profile costs 55,358 OTs
and 3,623,410 bits. These are initial-correlation reservation bounds, not
runtime measurements or total wire traffic. The same profiles remain selected
by the cost sweeps; their 43.37/41.13 batch-correctness bits and repeated
expansion costs are unchanged. Appendix alternatives are updated too.

The full paper experiment suite passes: 207 tests and 9 subtests. The PDF
has 71 pages, stable references, and no new overfull boxes; changed pages,
tables, and adjacent protocol figures were rendered and inspected. No C++
or benchmark changes were made in this reconciliation pass. Historical
cost figures below are superseded. Next return to the paper proof audit:
concrete generator/representation transfer and the full OLE correlation law
remain open, independently of this accounting correction.

**2026-09-07 approved implementation fix:** Peter asked to fix the sparse
activity-mask discrepancy. The package worktree now tracks public
`mHasSplit` separately from secret `mActivity`, accumulates corrected parent
tags once per split, and masks both child sums with one shared
bit-times-256-bit-string multiplication. The reservation changes from Dn
to (d+2(D-d))n OTs per direction, where D is domain depth, d is dense depth,
and n is the number of parallel instances. Entirely unused public levels
consume no OTs. Cached expansion is unchanged; both peers need updated OT
provisioning. The new masking call also uses hardened aligned storage in
the generic DpfMult path. No threads or benchmarks were added.

All 16 tests in `analysis/run_sparse_dpf_mask_tests.sh` passed under WSL/GCC,
including the new regression, existing sparse-DPF tests, and Waterfall and
Reverse-Cuckoo integration tests. See AUD-212 and
`analysis/sparse_dpf_activity_mask.md` in the package worktree. No MSVC test
was run. Next reconcile paper setup-cost accounting with this operation;
do not quote old implementation setup totals as current without checking
whether their model already included the mask. Concrete generator-security
and representation arguments remain open. No TeX/PDF changes in this pass.

**2026-09-07 generator-reduction pass:** Appendix C now gives Theorem C.3
for the abstract standard/sparse DPF protocols: tag-conditioned tree-PRG
security plus a joint leaf-PRF assumption imply online chosen-share
simulation. The reduction costs a_S tree challenges and one Q-query leaf
PRF challenge. It does not introduce fresh public randomness or change the
protocol. Concrete AES evaluators still need separate justification.

Read `REUSABLE_DPF_GENERATOR_AUDIT.md` before transferring that theorem to
the implementation. The read-only C++ audit found that sparse setup masks
publicly empty levels, not secret-inactive levels. A publicly present
off-path split therefore opens zero reconstructed child sums without the
paper's secret-activity mask. The support {0,1,4}, alpha=4, no dense prefix,
gives an exact algebraic witness: raw sigma=0 and tau=(1,0) at bit 1.
The new toy-model regression checks this witness while reconstruction
remains correct. It is not an end-to-end C++ test or a production leakage
estimate. No implementation fix has been made; discuss it with Peter.

The audit also inventories the dense custom two-child AES evaluator,
sparse raw-seed/separate-tag representation, and public-key leaf schedule.
Their equivalence to the abstract primitive requirements is not proved;
no attack on the dense evaluator is claimed. Next confirm the inactive
mask witness in a targeted implementation test and agree on the fix, then
finish the concrete evaluator arguments. The full OLE correlation law
remains a separate obligation. The PDF has 70 pages; all 70 targeted tests
pass. Revised pages were rendered and inspected. No undefined references
or new overfull boxes were introduced.

**2026-09-07 approved text and proof pass:** Peter approved the updated
functionality and asked to update the text and work on the proofs. The
preliminaries, overview, PCG scope, and Waterfall theorem/appendix now follow
corrupt-party output-share selection. Waterfall no longer conditions column
maps on a prescribed final output: it runs the corrupt share-selection
strategies, sums their maps, and submits that sum to the top functionality.
The proof preserves joint honest outputs and supports adaptive calls; its
routing charge remains once per state. Concrete realization errors remain
in the separately assumed subprotocol bound.

`ReusableDpfAppendix.tex` (Appendix C) now proves two exact statements:
the binary/wide correction-word laws in the independent-child experiment,
including masked inactive sparse levels, and the online leaf-update
simulator in the fresh-active-value experiment. The latter includes both
corruption choices, non-characteristic-two sign correction, and the honest
output complement. It does not yet prove the real-to-hybrid generator
transitions. Close those transitions under the intended generator model
next: account for tagged active inputs, joint reuse with setup auxiliary
information, collisions, and any ideal-cipher forward/inverse queries.
Do not call the concrete reusable-DPF or full OLE proof complete.

The targeted suite has 69 passing tests (14 new correction/opening tests).
Historical uniform-output-fiber tests remain labeled as historical algebra.
The PDF has 69 pages, no new overfull boxes, and no undefined references.
Changed text and appendix pages were rendered and inspected. No C++,
protocol, parameter, or performance-number changes were made in this pass.

**2026-09-07 functionality edit for Peter's review:** Peter approved the
simple ideal-world rule: receive the corrupt party's output share and
output the correct complement to the honest party. `DPF.tex` and `DMPF.tex`
now use that rule, with matching adjacent simulation definitions. Setup,
cached-seed protocols, and implementation are unchanged. The PDF builds to
66 pages; Figures 1 and 6 are on pages 7 and 14. Both figures and their
continuation pages were rendered and inspected. No new overfull boxes or
undefined references were introduced. The downstream Waterfall proof still
uses the former prescribed-output argument and must be revised next;
neither the concrete reusable-DPF simulation nor the OLE output-distribution
proof is completed by this definition edit. Wait for Peter's review before
expanding the paper changes.

**2026-09-07 contract investigation:** Peter questioned adding fresh public
coins merely to satisfy stronger output-sharing semantics than applications
need. `REUSABLE_DPF_SHARING_CONTRACT.md` proposes an online simulator that
produces the corrupt share, with the ideal functionality completing the
honest share. It retains cached seeds and adaptive payloads. SHARK's IFSS
definition supplies a precedent for this output convention, not a theorem
for our reusable shared-input API. The note gives conditional leaf-update,
Waterfall, and OLE reverse-sampling arguments. Concrete setup/hidden-active-pad
simulation and the all-honest OLE output-distribution hybrid remain to prove.
All 55 targeted algebra/hybrid tests pass, including 18 new checks. This is
an investigation only: no TeX, functionality, protocol, or C++ change was
authorized or made. Discuss the proposed contract before adopting it.

**Newest proof checkpoint:** the concrete DPF realization audit found a
same-tag output-difference check against the fresh-uniform-share ideal
functionality under online simulation. See `REUSABLE_DPF_SIMULATION_AUDIT.md`
for the chronology, exact 1/|G| ideal acceptance probability, and scope.
The hidden index is not recovered, and Waterfall's ideal-subprotocol theorem
is not refuted. All 37 targeted tests pass. No TeX, functionality, or
protocol changes were made during this audit. Discuss the output-sharing
contract with Peter before attempting to complete the realization proof.

**Subsequent implementation fix:** Peter confirmed that the public Ring-LPN
masks must be resampled on every expansion. The local `codex/dmpf-package`
worktree now exchanges fresh 16-byte contributions from both parties and
regenerates the masks and their products in `RingLpnTriple::expand`.
Supports and DPF state remain cached. The patch and benchmark are committed
as `b2e9b7b` on `codex/dmpf-package`. All eight
Ring-LPN test groups pass in WSL, plus ten-round P=4 stationary checks with
both Reverse-Cuckoo and Sum-DMPF over F12289, Fp31, and Goldilocks. These are
correctness tests; the full harness uses mock OTs. The older mismatch finding
below is now fixed in this worktree, not necessarily in other branches.
Historical paper timings still predate mask refresh; remeasure under their
original configuration before replacing those numbers. A quick WSL/O2 check
at P=4, t=16, N=2^20 found 3,761.20 ms per repeated expansion and 322.084 ms
for isolated mask refresh, both summed over the two parties on one thread.
The 8.56% ratio is a cost share, not a measured before/after slowdown.
The committed `analysis/ringlpn_mask_bench.cpp` and dated JSON preserve the
fixture and measurements. All six full expansions passed output checks.

The ten approved CWC corrections have been applied to the paper; see the
current checkpoint in `PAPER_NARRATIVE_AUDIT.md` for the complete list.
The PDF now has 66 pages. The test suite passes 153 tests plus nine subtests.
That CWC pass did not change C++ implementation, benchmark data, or Waterfall
profile totals. The subsequent mask fix above is a separate implementation change.

Most importantly, the Reverse-Cuckoo application now states a joint
**decisional** assumption (`ass:rc-decisional-lpn`), not merely support-search
hardness. It includes the concrete descriptors and a polynomial sequence of
samples. All P noise polynomials enter one ring sample per party; r_0=1 and
the other public ring masks weight both the inputs and their cross products.
The stationary model refreshes those masks across expansions.

The now-resolved mismatch was discovered by reading `RingLpnTriple.h`:
`init` called `sampleA` with a fixed seed, and later expansions reused that vector.
The evaluation now distinguishes these prototype timings from the paper's
fresh-mask reusable model. Fixed supports plus fixed masks confine complete
samples to a subspace of dimension at most Pt; the significance for the
protocol view depends on the target correlation-security definition. This
motivated the approved fresh-mask fix. It does not reopen that settled schedule.

DPF cleanup adds valid-input preconditions, singleton identity resharing,
w-ary/truncated-tree correctness, and explicit digit order. Group-By count
width and overflow behavior are repaired. Shuffle is an ideal functionality.
The generic Waterfall hybrid proof remains intact. Concrete reusable-DPF
simulation and full PCG composition are still separate proof obligations.

**Recommended next step:** the mask schedule is settled. The paper now
distinguishes silent PCG expansion from interactive reusable batches, explains
sampling in NTT coordinates, and labels the historical measurements in the
front matter and evaluation. Continue with the formal correlation interface
and DPF/PCG realization proofs. Do not
resume hash-leakage refinements or broaden the construction hierarchy merely
because those proof obligations require care. Consult Peter before changing
implementation behavior, functionality semantics, or cryptographic parameters.

The older snapshot below records history; this checkpoint takes precedence
where descriptions of the paper or recommended next work differ.

## 1. Executive state

The project now has a deliberately split construction hierarchy.

1. **Waterfall Cuckoo is the generic construction.** It samples public hashes
   before routing and obtains an input-independent candidate matrix. A fixed
   MPC circuit routes the rows. Its security claim is an ordinary simulation
   statement with a bounded statistical correctness/simulation loss. It uses
   no LPN or syndrome-decoding assumption.
2. **Reverse Cuckoo is a Ring-LPN specialization.** It hides a random placement
   first and then programs public hashes around that placement. It reduces the
   repeated local expansion to $2N$, but the public descriptor is genuinely
   support dependent. It is not a generic DMPF. Its use requires an explicit
   leakage-robust Ring-LPN assumption.

The most recent paper audit aligned the Reverse-Cuckoo hashing analysis with
the Ring-LPN caller's unconditioned uniform support sampling. That analysis
gives 41.935 bits of programmed-system correctness at ring degree
$N_{\mathsf R}=2^{20}$. At the same evaluator widths, the bound gives
34.213 bits at $2^{18}$ and 26.218 bits at $2^{16}$. These are upper bounds on
failure probability, not measured failure rates or demonstrated limits on
achievable correctness.

The proposed weak-factor rejection changes the support distribution used by
that proof. The blanket conditioning bound gave only 39.946 bits. This gap is
now closed: a polynomial-pair marginal bound gives at least 40.940 bits for
independent local support-only rejection; see Section 8.6. No C++ or parameter
change was needed. The prototype rejection sampler remains unimplemented.

The immediate high-value work is an audit, not another design expansion:

- verify that the Waterfall simulation proof matches the C++ transcript and
  failure behavior exactly;
- preserve the completed rejection-conditioned hashing proof and its exact
  regression tests;
- formalize the full leakage experiment for the two Ring-LPN instances and all
  product-list descriptors, rather than relying on prose alone;
- prepare the weak-factor rejection needed to match the cited Ring-LPN
  parameter recipe, together with its effect on the hashing proof;
- rerun and preserve one-at-a-time Waterfall and Reverse-Cuckoo benchmarks;
- decide whether the paper needs secure $2^{16}$ or $2^{18}$ profiles. If
  so, evaluator widths alone cannot fix the current affine-cube obstruction.

## 2. Peter's motivation and decision criteria

The main goal is a practical, reusable, fully distributed DMPF for PCGs. The
important interface is setup/expand separation:

- setup receives secret-shared indices and establishes all address-dependent
  state;
- each later expansion receives only a new secret-shared payload vector;
- no new routing correlations are needed after setup;
- repeated expansion should be almost entirely local.

For stationary constructions, the same support is reused many times. The
coefficient multiplying $N$ in every local expansion is therefore a central
performance target. This explains the persistent interest in $2N$, $3N$,
and $4N$ profiles even when a lower-expansion profile has a more expensive
one-time setup.

The cost model should report two coordinates:

- **OT count**, treating one shared-bit AND as two OTs and one shared bit times
  a shared string as two wide OTs;
- **OT payload bits**, because two OTs on 128-bit strings and two OTs on bits
  have the same OT count but very different communication.

Round complexity and circuit depth are secondary. Total nonlinear work and
payload are the main setup metrics. Bin-solver work must be counted when a
construction uses it.

Peter wants the generic construction to have a conventional simulation proof.
He is willing to use an additional decoding assumption for a clearly isolated
Ring-LPN optimization when the $2N$ benefit is real. He does not accept
removing a useful protocol path and then calling the remaining system secure.
When a proposed change is only a proof device, say so; do not describe it as a
protocol efficiency improvement.

The preferred exposition is concrete and notation disciplined:

- define domains at first use, for example $i\in[t]$ and $j\in[m]$;
- preserve paper notation instead of inventing synonyms;
- use $B$ for payload values;
- use $P$ for the position map and $S_j$ for sparse sets;
- use one placement matrix $O\in\{0,1\}^{t\times m}$, not separate ad hoc
  placement flags;
- say “honest party” and “semi-honest corrupt party,” never “insider”;
- protocol boxes must be code-only, must define every invoked local helper or
  cite it as a previously defined standard primitive, and must give typed
  signatures;
- setup and expand belong in the same box and should mirror the functionality;
- persistent state is stored by the protocol/functionality, not returned as an
  explicit output;
- function declarations must be left aligned and followed by a new line;
- render protocol figures after editing. Do not trust source inspection alone.

For C++ work, performance is a first-class constraint. Preserve unrolled and
fixed-width hot loops, batching, data-oriented layouts, and explicit ownership.
Do not add threads yet. Never run two benchmarks concurrently.

## 3. How the design reached the current split

### 3.1 The original Reverse-Cuckoo idea

The original construction padded $t$ real rows to $m$, secretly permuted
the padded rows, interpreted the resulting rows as physical bins, and solved
for public hash descriptors consistent with those placements. This is
attractive in MPC because it replaces iterative cuckoo insertion with small
binary systems.

The initial hope was that the hidden permutation and dummy padding made the
published descriptors essentially uniform. The core issue became clear after
comparing two idealized games for a fixed support $A$:

1. sample a uniform hidden injection $\tau$, then sample a uniform hash tuple
   $h$ consistent with $\tau$;
2. sample a uniform hash tuple conditioned only on the existence of some valid
   injection for $A$.

Let $Z_A(h)$ be the number of valid injective cuckoo placements of $A$
supported by $h$. Game 1 gives

\[
  \Pr[H=h\mid A]
  =\Pr_U[H=h]\frac{Z_A(h)}{\mathbb E_U[Z_A(H)]},
\]

whereas Game 2 gives equal weight to every hash tuple satisfying
$Z_A(h)>0$. The hidden permutation has already been summed into $Z_A(h)$;
it does not remove the multiplicity bias. This was the decisive obstruction.

Dummy rows provide slack and can smooth the distribution, but they impose no
hash equations and do not make $Z_A(h)$ constant. Increasing the number of
choices $w$ helps ordinary cuckoo feasibility but does not by itself restore
input-independent Reverse-Cuckoo descriptors.

Normal cuckoo executed inside MPC can also leak its success event if that event
is opened. The relevant distinction is quantitative and distributional. With
standard parameters, ordinary cuckoo failure can be negligible for each fixed
input sampled before the public hash. Reverse Cuckoo changes the probability
of every feasible descriptor according to its placement multiplicity; the
difference is not merely a rare abort bit.

### 3.2 Discarded attempts to repair Reverse Cuckoo generically

Several structural fixes were considered:

- program extra collisions to compensate for the missing collision mass;
- plant higher-degree collision structures;
- increase $w$ and combine several programmed structures;
- generate multiple placement-first proposals and select among them with an
  inverse-likelihood rule.

The collision-planting ideas tended to create an effective $w=5$ object or a
sequence of patches whose distribution was harder to state than the original
problem. They were not adopted as the generic construction.

The multiple-proposal idea remains a useful LPN-specific heuristic. It changes
the soft likelihood ratios, but every finite number $K$ of feasible proposals
retains the hard support condition $Z_A(h)>0$. It therefore cannot give a
reduction to ordinary DMPF security by itself.

### 3.3 Waterfall Cuckoo

Waterfall was introduced to obtain an input-independent public hash
distribution without running data-dependent cuckoo insertion. The parties
sample the hash candidates directly and process the partitions in public
order. Colliding rows flow to the next partition. A bounded canonical
augmenting-path circuit repairs the small final overflow.

The routing design evolved through the following experiments:

- a basic ordered waterfall scan;
- bounded random evictions;
- biased first evictions that avoided the final partition;
- targeted conditional simulation of rare overflow cases;
- shallow breadth-first repair;
- a complete multi-source reachability search.

Random eviction gave useful intuition but was not a clean proof strategy at
the desired tail probabilities. The current main construction uses complete
repair for a bounded number of roots. The shallow one-repair circuit remains a
compact appendix alternative.

### 3.4 Scatter and permutation history

The scatter step originally used an aggressive Waksman-like permutation whose
distribution was not exactly uniform. That was intentionally challenged and
replaced.

The current design uses two serial Waksman passes. Party $p$ privately
controls pass $p$. Their composition is exactly uniform conditioned on either
party's view, and the same sampled network supports the required inverse pass.
This is deliberately the “boring” exact solution. The permutation is not the
dominant cost, so do not reintroduce an approximate network without a measured
reason and a complete distribution proof.

The implementation in `C:\Users\peter\repo\secure-join` contains the heavier
construction associated with ePrint 2024/547. It motivated the
permute-then-reveal scatter pattern but was not imported because it was too
heavy for this use case.

Once the permutation scatter was adopted, the earlier initial shuffle in the
top-level Waterfall protocol became redundant and was removed.

### 3.5 Reusable DPF cleanup

The paper initially fused address setup and payload expansion. Peter required
the full reusable interface instead of a halfway state-passing formulation.
The intended semantics are simple:

1. `Setup` receives the index and expands the DPF tree to leaf seeds $s$.
2. The protocol stores those seeds internally.
3. During expansion $q$, a round-indexed PRF/PRG
   $\mathsf L_q(s):=\mathsf L(q,s)$ maps each saved seed to a fresh group
   element.
4. `updateLeaves` applies the new payload $B$.

Only `updateLeaves` needs $B$. Setup does not receive payloads. The
functionality and protocol boxes now follow this API, but this area still
deserves a semantic audit because older explanatory prose remains in `DPF.tex`.

### 3.6 Implementation and optimization history

Waterfall was implemented end to end in libOTe. Early performance comparisons
were misleading because the profiles did not consistently use $t=16$ and
because Waterfall used a less efficient post-routing expansion path than
Reverse Cuckoo. The implementation was corrected to use the same expansion
organization after routing.

The sparse-set builder was optimized without threads. Sparse-DPF execution was
optimized while preserving the protocol. The key execution idea was to expand
eight active subtrees through one batched AES kernel and share accumulation
state, similar in spirit to the memory-layout work under
`C:\Users\peter\repo\ecdsa-vole\ecdsaVole\Tree`. The explicit eight-lane loops
are intentional hot-path code.

A sparse-bucket overflow was addressed in commit `84cba33` by adding batch
slack to the estimated capacity. The packaged code subsequently replaced that
capacity estimate with exact counting and prefix offsets. It now allocates one
contiguous buffer of the exact required size; see Section 10.4.

A separate crash exposed a C++ coroutine-frame alignment bug. Over-aligned SIMD
scratch locals crossed a suspension point and were placed in a coroutine frame
whose allocation was not sufficiently aligned. Commit `d6ce62f` moves the
complete level-expansion kernel into a non-coroutine, non-inlined helper. This
preserves eight-way AES batching while keeping aligned scratch off the
coroutine frame.

## 4. Repository and branch map

### 4.1 Paper repository

`C:\Users\peter\repo\dmpf` is on `main` at base commit `b55bcc4`, tracking
`origin/main`. The working tree is intentionally very dirty. The current paper
rewrite is not committed.

Important new, currently untracked source files are:

- `WaterfallCuckoo.tex`
- `WaterfallCuckooAppendix.tex`
- `RevCuckooSpecialized.tex`
- `RevCuckooHashAppendix.tex`

The active `main.tex` inputs those files. It no longer inputs the old
`RevCuckoo.tex`, although that legacy file still contains a large modified
draft. Treat `RevCuckooSpecialized.tex` as authoritative. Do not accidentally
edit or cite the inactive `RevCuckoo.tex` as the current construction.

The working tree also contains many untracked experiment scripts and CSVs.
They are research evidence, not disposable build products. Do not clean them
with a blanket command. `experiments/__pycache__`, `tmp`, LaTeX auxiliaries,
and rendered output are disposable, but verify exact paths before removing
anything.

The latest validated PDF is `C:\Users\peter\repo\dmpf\main.pdf`. It has 69
pages. `latexmk` succeeds, references stabilize, and the rewrite introduced no
new overfull boxes. Five older overfull warnings remain elsewhere in the
document.

### 4.2 Packaged libOTe implementation

The isolated implementation worktree is
`C:\Users\peter\repo\libOTe-codex-rev-analysis` on
`codex/dmpf-package` at `d6ce62f`.
It now has uncommitted hash-analysis updates: the conditioning calculator and
tests, the updated analysis report, and a corrected compression-proof
docstring. The C++ implementation is unchanged.

The package sequence is:

- `952ddca` — analyze and harden Reverse-Cuckoo hash widths;
- `e7a8355` — implement Waterfall DMPF and the optimized reusable sparse-DPF
  path;
- `6d2d999` — document the bundled DMPF constructions;
- `d6ce62f` — keep aligned SIMD scratch out of the sparse-DPF coroutine.

This worktree is the best place for isolated inspection and testing.

The ordinary checkout `C:\Users\peter\repo\libOTe` is on
`codex/foleage-trace` at merge commit `d0e6161`, which merged `d6ce62f`. That
checkout currently has unrelated uncommitted changes from another party in
Foleage, RingLPN, AnyField, tests, and a new `LpnParameters.h`. Do not use it for
isolated cleanup and do not overwrite those changes.

An older standalone Waterfall worktree exists at
`C:\Users\peter\.codex\worktrees\waterfall-dmpf\libOTe` on
`codex/waterfall-dmpf` at `c2a3aa9`. The packaged branch above supersedes it.

### 4.3 Historical AnyField scope

The package branch has AnyField Reverse-Cuckoo work in its ancestry. Earlier
in the project Peter asked to avoid carrying unrelated AnyFieldOLE changes
while isolating the hash analysis. The resulting package nevertheless descends
from the branch containing those commits. If a minimal upstream patch series
is needed, audit the commit ancestry rather than assuming the package contains
only Waterfall and hash-width changes.

More importantly, the paper's Reverse-Cuckoo correctness proof is specific to
the Ring-LPN caller. The AnyField input distribution has not received the same
relation-counting analysis. Do not infer AnyField security from the Ring-LPN
proposition.

### 4.4 Source material already used

The main external and local references behind the current direction are:

- ePrint 2022/1035, for the earlier DMPF/cuckoo context that motivated this
  project;
- `C:\Users\peter\Downloads\2023-1568.pdf`, for regular syndrome-decoding
  asymptotics and the heuristic that structured sparse noise can remain hard;
- `C:\Users\peter\Downloads\2024-547.pdf`, for the secure permutation and
  reveal-mapping pattern used when redesigning Waterfall scatter;
- ePrint 2025/295 and the corresponding dissertation material, for Stationary
  Syndrome Decoding, linear-test evidence, and cross-expansion caveats;
- ePrint 2023/176 and 2025/415, for algebraic attacks and the limits of
  interpreting XL-style soundness evidence;
- `C:\Users\peter\repo\secure-join`, which contains Peter's heavier
  implementation related to the 2024/547 permutation construction;
- `C:\Users\peter\repo\ecdsa-vole\ecdsaVole\Tree`, whose half-tree execution
  and memory-layout ideas motivated the sparse-DPF optimization;
- ePrint 2025/2294, the public version of the current DMPF paper referenced by
  the standalone libOTe hash-analysis report.

When citing or changing claims from these works, re-open the primary paper.
The experiment notes contain paraphrases and design conclusions, not a
replacement for source verification.

## 5. Current reusable functionality contract

The paper now targets a reusable DMPF functionality:

- `Setup(A)` stores the address vector and initializes the expansion counter;
- `Expand(B)` returns shares of the multi-point vector for the stored addresses
  and the new payloads, then advances the counter;
- polynomially many expansions are allowed;
- setup cannot be invoked twice without reinitializing the instance;
- payloads can be chosen adaptively across expansions;
- routing correlations and address-dependent DPF OTs are setup-only.

Duplicate addresses are merged in setup. On each expansion, the payloads of
all duplicates are summed into the retained first occurrence. The public dummy
address is $N\notin[0,N)$ and always carries zero payload.

The DPF protocol stores leaf-seed shares, not output values. Expansion uses
the public round number $q$ to select $\mathsf L_q$, derives fresh group
shares from the saved seeds, and calls `updateLeaves` with $B$. This design is
important for stationary reuse and must remain consistent across the DPF,
Waterfall, and Reverse-Cuckoo protocol boxes.

## 6. Waterfall Cuckoo: technical state

### 6.1 Objects and distribution

Let $t$ be the maximum number of input rows. Let the $w$ partitions have
sizes $\mathbf d=(d_0,\ldots,d_{w-1})$, let $c$ be the stash capacity, and
let

\[
  m=\sum_{s\in[w]}d_s+c.
\]

For every row $i\in[t]$ and partition $s\in[w]$, the candidate value
$C_{i,s}\in[d_s]$ selects one physical main-bin position. The logical
position map $P\in\{0,1\}^{[0,N)\times m}$ contains every admissible position
for every real address. The actual placement
$O\in\{0,1\}^{t\times m}$ has row and column weights at most one. A real row
may select only a column admitted by its address in $P$. A dummy row instead
uses its independently sampled candidates in $C$. On successful routing,
every one of the $t$ rows, including dummies, has weight one. An unmatched row
has weight zero on a failed routing. Physical placement and real-row activity
are separate properties.

This notation replaced earlier separate output/flag records. Preserve it.

For arbitrary $N$, the paper embeds real addresses into
$\mathbb F_{2^\ell}$, where
$\ell=\max\{1,\lceil\log_2N\rceil\}$. A degree-$(t-1)$ polynomial with
uniform field coefficients gives a $t$-wise independent hash. Each partition
uses an independent seed. Real deduplicated addresses obtain their candidates
from these hashes. Dummy rows obtain independent uniform candidates directly.

For every fixed deduplicated input vector, the complete matrix

\[
  C\in\prod_{s\in[w]}[d_s]^t
\]

is exactly uniform. This lemma is the foundation of the generic simulation
argument. Do not replace it with an informal “hashes look random” claim.

### 6.2 Placement algorithm

The initial waterfall scan processes all rows through partition 0, forwards
colliding live rows to partition 1, and continues in order. Each partition has
a secret occupancy array. Public row order provides deterministic tie breaking.

After the scan, a complete repair round treats all unmatched rows as roots of
one alternating-path reachability search. It expands the first discovered but
unexpanded row in public order and scans partitions in public order. The first
free candidate determines a unique parent path, which is flipped. A round
repairs one row whenever any augmenting path exists. The generator runs a
public fixed number $r$ of rounds and then compacts any remaining rows into
the $c$ stash columns.

If $L$ rows remain after the initial waterfall stage and the candidate graph
has a size-$t$ matching, $r$ complete rounds plus a stash of size $c$
succeed whenever $L\le r+c$. Therefore

\[
  p_\Gamma
  \le \Pr[L>r+c]
     +\Pr[G_C\text{ has no matching of size }t].
\]

The first term is computed by the exact occupancy recurrence. The second uses
a Hall-witness union bound. The detailed circuit, recurrence, and proof are in
`WaterfallCuckooAppendix.tex` and the Python reference models.

The implementation calls the traversal a multi-source reachability search.
It does not process strict distance layers and should not be described as a
shortest-path BFS.

### 6.3 Scatter

`WaterfallScatter` hides the physical column identities before opening $O$.
Party $p$ samples a uniform permutation $\sigma_p$ on $[m]$ and privately
programs one Waksman pass. Their composition $\rho=\sigma_1\circ\sigma_0$
permutes the column dimension:

\[
 \widehat O_{i,\rho(j)}=O_{i,j}.
\]

The $t$ logical row labels remain fixed. In C++, each column is stored as a
$t$-bit record, so these are the rows of the transposed storage matrix
`hiddenPlacement`. They must not be confused with the rows of mathematical
$O$. The opening gives a uniform labeled injection into the permuted columns
on routing success, in the ideal permutation-subprotocol view.

Public wiring uses $\widehat O$ to scatter shared activity bits and masked
addresses into those permuted columns. The inverse passes then restore the
original physical columns. The hidden row destinations are cached separately
for reusable payload scatter.

Every dummy row must also have a physical placement on success; otherwise the
opening reveals the number of real rows. On failure, zero rows of
$\widehat O$ reveal exactly which logical row labels were unmatched. The paper
charges the entire failed-state probability $p_\Gamma$ in the simulation
loss. The transcript-level audit must verify that this replacement is valid.

### 6.4 Sparse sets and expansion

For each physical column $j\in[m]$, define

\[
  S_j=\{x\in[0,N):P_{x,j}=1\}.
\]

For a main partition, its columns partition the domain. A stash column admits
the full domain. The total number of sparse leaves is

\[
  \sum_j |S_j|=(w+c)N.
\]

Setup creates one reusable sparse DPF for each physical column and caches a
second family of row-to-column point DPFs for scatter. Each later expansion
maps the new $B$ values into their hidden columns, updates the cached column
leaves, and accumulates them at their original domain positions.

An all-zero/empty sparse column is valid. It contributes no real address and
needs no explicit failure condition.

### 6.5 Main profiles

Both main profiles use $t=16$, $c=0$, and the same complete-repair
generator.

| profile | $\mathbf d$ | $r$ | $m$ | local expansion | OTs per state | OT payload per state | batch bits for 256 states |
|---|---:|---:|---:|---:|---:|---:|---:|
| default compact $4N$ | $(16,16,16,16)$ | 3 | 64 | $4N$ | 35,980 | 363,304 bits | 43.37 |
| lower-expansion $3N$ | $(128,128,64)$ | 2 | 320 | $3N$ | 42,382 | 1,268,210 bits | 41.13 |

These setup counts include direct hashing, placement, scatter, and reusable
DPF initialization. The count is per cached state. The batch size 256 affects
the displayed correctness bits, not the per-state OT count.

The compact $4N$ profile is the recommended default. The $3N$ profile is
worth using only when enough expansions amortize the larger setup, payload,
and 320-column state.

### 6.6 Alternatives that should remain secondary

- The asymmetric complete-repair profile
  $\mathbf d=(32,32,16,8),c=0,r=2$ has $m=88$, $4N$ expansion,
  41.20 batch bits, 34,142 OTs, and 414,166 payload bits. It saves 1,838 OTs
  relative to the default but increases state and payload.
- The proved shallow one-repair profile
  $(t,w,d,c)=(16,4,16,2)$ has $m=66$, $6N$ expansion, 45.2277 batch bits,
  at most 27,932 OTs, and 263,784 payload bits. It is attractive for minimal
  setup but not for repeated expansion.
- The initial scan with $r=0$ is a parameter choice, not a separately named
  construction.
- Random-eviction variants and earlier bounded-walk generators are experiment
  history, not paper-level schemes.

### 6.7 C++ implementation state

The implementation is under `libOTe/Dpf/Waterfall/` and is orchestrated by
`libOTe/Dpf/WaterfallDmpf.h`.

- `WaterfallConfig.h` supplies only the two $c=0,t=16$ main profiles.
- `WaterfallHash.h` implements the exact $t$-wise field-polynomial hash. At
  20 field bits it uses a generated $\mathbb F_{16}^5$ tower multiplier with
  81 shared bit products. Address powers are prepared once and reused across
  partitions and sets.
- `WaterfallBasic.h` implements the ordered initial scan.
- `WaterfallReachability.h` implements bounded complete repair.
- `WaterfallScatter.h` and `SerialWaksmanPermute.h` implement exact scatter.
- `CachedDpfExpansion.h` and `SparseDpf.h` contain the optimized repeated
  expansion path.
- `Waterfall_Tests.cpp` contains component and end-to-end tests.

Current implementation restrictions are narrower than the abstract paper:

- the domain must be a power of two and at most $2^{32}$;
- $t$ must be a power of two in the placement MPC;
- partition sizes must be powers of two;
- built-in profiles require $t=16$;
- `WaterfallConfig` does not expose a stash, because both implemented main
  profiles have $c=0$.

The paper can state the general construction, but the implementation section
must not imply that every general profile is implemented.

## 7. Reverse Cuckoo: technical state

### 7.1 Current $K=1$ construction

The main specialization has $w=2$, $d=16$, $m=32$, and $t=16$. It
pads with (m-t) dummies, applies a uniform secret shuffle, and interprets the
post-shuffle rows as a uniform real-row injection into the 32 physical bins.

Each partition uses a public degree-two evaluator $F_s$. For partition
$s$, the protocol forms the evaluator matrix on its $d$ assigned rows and
uses `bin-solver` to find a descriptor mapping every independent active row to
its public row label. Dummy rows evaluate to zero and impose no constraint.
The solver also skips a row that reduces to zero with a nonzero right-hand
side; any resulting inconsistency is a bounded correctness error, not an
explicit abort.

The public position map contains one position per partition for every real
domain address. Hence the total sparse-leaf count is exactly $2N$.

### 7.2 Exact ideal-hash leakage channel

For a fixed support $A$, let $Z_A(h)$ be the number of injective placements
compatible with public hash tuple $h$. The idealized $K=1$ descriptor law is

\[
  \Pr[H=h\mid A]_{\mathsf{RC}_1}
  =\Pr_U[H=h]\frac{Z_A(h)}{\mathbb E_U[Z_A(H)]}.
\]

This is a size-biased cuckoo channel. It is not uniform conditioned on
feasibility. Dummies change the placement space and can reduce visible bias,
but do not change the multiplicity weighting.

For $w=2$, each support point is an edge in a bipartite multigraph. A valid
placement orients every edge toward a distinct endpoint. A tree component with
$e$ edges contributes $e+1$ orientations, a unicyclic component contributes
two, and an overfull component contributes zero. The global matching count is
the product of component contributions.

This graph view is useful intuition but not a complete attack model. Pair
collisions are only a surrogate for a higher-order matching statistic.

### 7.3 Finite-$K$ debiasing

The explored heuristic draws $K$ independent raw placement-first proposals
and selects among them using a public score derived from endpoint occupancy.
The exact inverse-multiplicity selector would weight proposal $i$ by
(1/Z_i). The cheaper proposed selector uses the public table

\[
  a(D)=1/\mathbb E_U[Z\mid D],
\]

where $D$ is the total endpoint occupancy deficit. A finite-precision
exponential race can select a proposal with probability proportional to
$a(D_i)$ without secret division.

For every finite $K$, all score weights remain positive. Therefore

\[
  \Pr_K[H=h\mid A]>0
  \quad\Longleftrightarrow\quad
  Z_A(h)>0.
\]

Finite $K$ suppresses soft multiplicity leakage but never removes the hard
feasibility language. This is why $K=4$ cannot be presented as a generic
security repair.

Historical cost modeling for $t=16,(d_0,d_1)=(16,32),K=4$ estimated
9,927,168 OTs over ten expansions, about 2.76 times the original symmetric
$K=1$ reference. The score/selector note reports 2,371,584 proposal-scoring
OTs, 155,136 selector OTs, and 3,499,008 OTs for the three extra proposals
across 256 descriptors. Treat these as model outputs, not current benchmark
claims.

`experiments/rev_cuckoo_k_tradeoff.md` says that a selector prototype existed,
but the present `codex/dmpf-package` tree contains no `OccupancySelect` or
`AsymmetricPlacement` symbol. This status is stale or the prototype lived on a
different branch. Do not claim a shipped finite-$K$ implementation without
locating it.

The paper now leads with $K=1$. This is intentional: it preserves the full
$2N$ motivation, avoids a large heuristic setup, and analyzes the strongest
raw leakage rather than choosing $K=4$ because it “looks like a knee.”

## 8. Reverse-Cuckoo programmed-hash correctness

### 8.1 Concrete evaluator

Let $n=\lceil\log_2(N+1)\rceil$, encode the dummy $N$ as $\nu$, and
sample ideal uniform matrices

\[
 A_0,A_1\in\mathbb F_2^{n\times q'},
 \qquad
 A_2\in\mathbb F_2^{(n+q')\times q}.
\]

The evaluator computes two random linear images, multiplies them
coordinate-wise, anchors both the linear and quadratic parts at $N$, and
applies the final random compression. Thus $F(N)=0$.

The implementation retains `linearSecParam` as an independent parameter and
uses the fixed construction margin

\[
 q'=q=d+\lambda+8
\]

up to byte alignment. Do not rename the setting to a 48-bit linear security
parameter. The user explicitly chose a hard-coded eight-bit evaluator margin.

The ideal-cipher model is assumed for the AES/PRG expansion. Within that model,
the matrices can be treated as truly random. Do not spend time attacking PRG
pseudorandomness unless the model changes.

### 8.2 Why independent-row reasoning is invalid

Fix the real rows assigned to one partition. A nonzero row selector
$\alpha$ is harmful when it annihilates the evaluator rows but not their
public labels. If its anchored linear moment is nonzero, the relation survives
final compression with probability $2^{-q}$. If the linear moment vanishes
and the alternating quadratic moment has rank $\rho$, the exact survival
probability is

\[
 \left(\frac12+2^{-\rho-1}\right)^{q'}
 +\left(1-
   \left(\frac12+2^{-\rho-1}\right)^{q'}\right)2^{-q}.
\]

Affine planes have quadratic rank two and survive the nonlinear lift with
$(5/8)^{q'}$. Affine three-cubes have rank zero and are deterministic
degree-two identities. Consequently, sampling the hash after fixing the rows
is necessary but does not justify treating the evaluator rows as independent
uniform strings.

### 8.3 Exact Ring-LPN caller mapping

Write $N_{\mathsf R}$ for the ring degree and $L=N_{\mathsf R}/t$. Each
party $p\in\{0,1\}$ has $P=4$ regular sparse polynomials. Polynomial $a$
has one offset $X_{p,a,j}\in[0,L)$ in each block $j\in[t]$, with $t=16$.

For a polynomial pair $(a,b)$ and product block $k$, the DMPF list is

\[
 A_{a,b,k}
 =\left\{
 X_{0,a,j}+X_{1,b,(k-j)\bmod t}:j\in[t]
 \right\}\subseteq[0,2L).
\]

There are $P^2t=256$ lists. Each has at most 16 points after deduplication.
Under unconditioned support sampling, the pre-dedup points within one list
are independent samples from the
triangular sum of two uniform offsets. Different lists share base offsets and
are correlated. The proof uses an expectation/union bound and does not assume
independent list failures. Rejection of a complete support changes this law;
the product-list formula still holds, but independence must not be carried
over without proof.

At $N_{\mathsf R}=2^{20}$, $L=2^{16}$, so each DPF list has domain
$[0,2^{17})$. Against one semi-honest corruption, the corrupt party's four
base support lists are known and the honest party's $Pt=64$ block offsets
remain hidden. Before considering the public syndromes and descriptors, those
offsets range over $L^{Pt}=2^{1024}$ possibilities before weak-factor
conditioning.

### 8.4 Current bounds

For unconditioned uniform supports at
$w=2,d=16,P=4,t=16,N_{\mathsf R}=2^{20}$:

| $q'=q$ | structural term | final compression | complete batch bound |
|---:|---:|---:|---:|
| 56 | $2^{-39.261}$ | $2^{-38.288}$ | $2^{-37.694}$ |
| 64 | $2^{-42.007}$ | $2^{-46.288}$ | $2^{-41.935}$ |

At the current width 64, the ring-size dependence is:

| $N_{\mathsf R}$ | $L$ | DPF domain | proved batch correctness bits |
|---:|---:|---:|---:|
| $2^{16}$ | $2^{12}$ | $2^{13}$ | 26.218 |
| $2^{18}$ | $2^{14}$ | $2^{15}$ | 34.213 |
| $2^{20}$ | $2^{16}$ | $2^{17}$ | 41.935 |

The bounds for the smaller two are dominated by the upper-bound contribution
from deterministic affine-cube relations. Increasing $q$ or $q'$ cannot
remove those relations, although a sharper count could lower the bound.
If those ring sizes
must become claimed secure profiles, the candidate fixes are a higher-degree
feature layer, input/rank rejection, or a sharper caller-specific argument.
Do not simply add more output bits.

For comparison, the appendix reports that $w=3,d=8$ gives 32.577 bits at
width 48 and 40.194 bits at the fixed-margin width 56. This is a proof
parameter comparison, not the main implemented $2N$ profile.

### 8.5 Scope limitations

The correctness analysis assumes:

- an ideal uniform secret shuffle;
- independent uniform Ring-LPN block offsets before weak-factor rejection;
- ideal random linear maps derived from the evaluator seed;
- semi-honest joint sampling, with no malicious bias;
- the exact Ring-LPN caller distribution above.

It bounds programmed-system inconsistency. It is separate from the descriptor
leakage assumption. It does not prove AnyField correctness and does not prove
the concrete descriptor family pseudorandom for arbitrary chosen supports.

### 8.6 Resolved: weak-factor rejection and hash correctness

The initial blanket transfer divided the batch bound by both parties'
acceptance probability and gave only 39.946 bits. A sharper proof now closes
this gap without changing the design.

Let $D$ count the occupied residues in one polynomial's 16 positions modulo
128. Since 128 divides $L$, this is the occupancy law of 16 independent uniform
draws into 128 bins. With independent copies $D_1,\ldots,D_4$, put

\[
 \alpha:=\Pr[D_1+\cdots+D_4\ge61]\approx0.5019144583,
 \qquad
 \beta:=\Pr[D_1+D_2+D_3\ge45]\approx0.7174085984.
\]

Fix one polynomial's full support. Its folded weight is at most 16, so
acceptance over the other three supports has probability at most $\beta$.
Its marginal density after rejection is therefore at most $\beta/\alpha$
times its unfiltered density. Each product list depends on one polynomial
from each party. Independent local rejection multiplies each list's expected
structural union contribution by at most
$(\beta/\alpha)^2\approx2.043025290$. Sum these contributions; do not apply
this factor directly to an arbitrary event involving all eight polynomials.

The compression bound is pointwise in the supports and needs no multiplier:
for fixed lifted rows, union over at most $2^r-1$ nonzero lifted relations,
each killed by the independent final matrix with probability $2^{-q}$.
This does not require independent lifted rows.

Consequently, at the stated $2^{20}$ ring degree,

\[
 \epsilon_{\rm filtered}
 \le(\beta/\alpha)^2\epsilon_{\rm str}+\epsilon_{\rm cmp}
 =4.738726728814093\ldots\cdot10^{-13}<2^{-40.940}.
\]

The paper now contains this argument in Appendix B.4 and extends the existing
correctness proposition. The package adds
`analysis/rev_cuckoo_hash_conditioning.py` and its regression tests.
Integer occupancy counts, an integer Walsh transform, and rational arithmetic
verify the comparison with $2^{-40}$ exactly; this is not a Monte Carlo claim.
The original unfiltered exponent is 41.934779..., rounded to 41.935 in tables.
The literal proposition now uses the conservative exponent 41.934.

Scope: independent local support-only rejection, residue counts separate per
polynomial, and independent ideal shuffle/evaluator coins. Coefficient-based
or joint rejection needs another argument. This closes the correctness
dependency, not the leakage-assumption or attack-evidence questions.

## 9. Ring-LPN leakage and attack evidence

### 9.1 The security question

Reverse-Cuckoo leakage does not make the solution space dense. It removes or
downweights candidate supports while preserving regular sparsity. The question
is whether this nonlinear posterior reweighting makes the planted sparse
support materially easier to find.

The analogy to regular syndrome decoding and Stationary Syndrome Decoding is
the main intuition:

- regular supports occupy a tiny subset of all supports, yet the leading known
  regular-ISD exponent remains comparable to ordinary ISD;
- SSD gives strong evidence against its fixed linear-test model even when a
  support is reused;
- fresh coefficients remain hidden and cannot be directly cancelled from a
  placement descriptor.

This analogy is not a reduction. Leakage-adaptive dual-word selection,
block-circulant Ring-LPN structure, higher-order Hall information, and the
shared concrete feature root are outside the cited no-leakage theorems.

### 9.2 Full observer model that must not be simplified away

The intended analysis is not one isolated random support list. The adversary
sees the ordinary public Ring-LPN samples for the complete PCG instance and all
Reverse-Cuckoo descriptors generated for the product lists. One semi-honest
corrupt party also knows its own four sparse polynomials. The hidden variables
are the honest party's four regular supports, their coefficients/secrets, and
the derived product-list placements.

There are two coupled Ring-LPN sides/instances in the PCG, and the leaked
descriptors concern sums/products of their support positions. Any formal game
must include this coupling. The current paper assumption names the public
Ring-LPN samples and every concrete descriptor, but it remains prose rather
than a fully procedural experiment. This is a priority audit item.

In the stationary setting, setup publishes one descriptor for a support and
later expansions change only coefficients. Repeated expansion does not reveal
fresh independent placements, but the adversary may amortize arbitrary
preprocessing of the fixed descriptor. The assumption intentionally permits
such preprocessing.

### 9.3 Why direct linear attacks do not immediately consume the leakage

The descriptor exposes collision, Hall, and placement-multiplicity facts about
candidate support images. It does not reveal error coefficients or an explicit
linear equation in the support positions.

For a fixed dual word, its syndrome-test bias is controlled by the posterior
probability that its support avoids the hidden error. Leakage can therefore
help only by guiding the choice of dual word. The SSD theorem gives the
no-leakage analogy, but does not cover this adaptive conditional-dual-avoidance
problem. Do not claim that “linear attacks cannot use the leakage.” The correct
claim is that no direct cancellation equation appears, and the remaining
adaptive problem has been isolated but not fully solved.

### 9.4 Attack portfolio already studied

The paper now gives a compact coverage table. The detailed methods and raw
results are in `experiments/rev_cuckoo_lpn_attack_study.md`.

Keep the parameter tuple and evidence type attached to every number. The
measured leakage selectors below mostly use the unconditioned ideal-hash
experiment at $(P,t,L)=(4,16,4096)$, hence $N_{\mathsf R}=2^{16}$. They do
not establish the same gain at $N_{\mathsf R}=2^{20}$ or after weak-factor
rejection. The older notes often call the smaller point “production.”

- **Linear, Fourier, and BKW tests:** no direct coefficient cancellation;
  leakage-adaptive dual avoidance remains open.
- **Stationary output-span and module tests:** reuse permits preprocessing but
  no fresh descriptor samples; independent ring masks are critical.
- **Likelihood-ordered support search:** exact toy models confirm real leakage.
- **Quadratic system solving:** the hidden-placement formulation is low degree,
  with about $P N_{\mathsf R}$ support variables. At $P=4,t=d=16$ and
  $N_{\mathsf R}=2^{20}$, the direct formulation has 131,072 placement
  variables and around $10^{13}$ quadratic support monomials. This is a
  combinatorial size calculation, not a solving-time measurement.
- **Regular Prange and posterior-guided ISD:** marginal, pairwise, sequential,
  and local matching-aware selectors were tested. At
  $(P,t,L,q_{\mathrm{field}},K)=(4,16,4096,65537,1)$, a six-instance
  selector screen combined with an optimistic attack-cost model gives
  $2^{127.123}$ versus its $2^{128}$ baseline. Polynomial elimination costs
  are omitted and the selected-matrix rank caveat remains. This is an
  estimated work exponent for the tested selector.
- **Equation-coupled lists and blockwise beam search:** clear advantages exist
  on small enumerated instances. Enumerating $L^t$ supports per polynomial
  requires $2^{192}$ candidates at $(t,L)=(16,4096)$ and $2^{256}$ at
  $(16,65536)$. These are list sizes; the tested bounded beams do not give
  a practical decoder at $L=4096$.
- **Syndrome projection and Wagner-style merging:** the tested projected attack
  loses to the flat merge; guessed intermediate shares repay the filtering
  gain in repetitions.
- **Regular Stern/MMT/BJMM:** at $(P,t,L,q_{\mathrm{field}})=(4,16,4096,65537)$,
  the tested $K=1$ Stern model gives $2^{146.526}$. The screened optimistic
  MMT model remains at least $2^{145.55}$ after granting the measured leakage
  gain. These estimates use the random-matrix/list-cost model; they do not
  bound every attack on the structured ring.
- **Sparse-factor statistical decoding:** the main factor reduction gives a
  serious ordinary Ring-LPN tail issue, but the tested leakage-aware projection
  produced no additional gain.
- **Concrete feature and structured-rank diagnostics:** reduced experiments did
  not find stable extra leakage or a guided-rank advantage. These are not
  production-rank proofs.

Exact toy experiments find multi-bit improvements. This matters: the leakage
is exploitable and must not be dismissed merely because it is nonlinear.
At the tested $P=4,t=16,L=4096$ point, the recorded scalable selectors give
roughly sub-bit to one-bit improvements, with small samples and unresolved
particle-degeneracy and rank issues. This evidence supplies neither a bound
on all decoders nor a measured gain for the $2^{20}$ ring.

Two earlier estimates were explicitly withdrawn in the experiment notes:
30--50-bit importance-sampling gains had effective sample fraction below
$10^{-3}$, and an 8.6-bit four-polynomial extrapolation was retracted after
unstable SMC normalizers changed with additional replicates. Do not revive
these figures from historical notes or CSVs as attack findings.

### 9.5 Weak-factor rejection and the two decoding-cost calculations

Reducing modulo $X^{128}+a$ folds the support. At the current regular-noise
parameters $(P,t,n_{\mathrm{factor}})=(4,16,128)$, the expected total folded
support weight is 60.383. The analysis uses two distinct cost calculations:

| calculation | unconditioned support law | conditioned on folded weight $\ge61$ |
|---|---:|---:|
| optimistic independent-coordinate expression from the parameter analysis | 122.482 bits | 130.572 bits |
| explicit uniform dual-check containment, using sampling without replacement | 136.431 bits | 150.526 bits |

These are logarithms of averaged model costs over the folded-weight law.
The often quoted 128.450-bit value instead evaluates the first expression
at the mean weight. The 122.482-bit figure is not a demonstrated attack cost;
the explicit dual-check calculation has the larger exponent in the second
row. Neither row certifies security against every decoder.

The rejection is required to match the cited Ring-LPN parameter recipe. Its
motivation predates the Reverse-Cuckoo leakage question, but its effect on
the product-support distribution must be included in the hash analysis.

For each party separately, define the folded support of polynomial $a$ by

\[
 T_{p,a}:=\{(jL+X_{p,a,j})\bmod128:j\in[t]\},
 \qquad W_p:=\sum_{a\in[P]}|T_{p,a}|.
\]

The proposed local sampler resamples that party's complete set of four
supports until $W_p\ge61$. Equal residues in different polynomials count as
different coordinates. The acceptance probability is about 0.501914 per
party when residues are uniform, as they are for the displayed block lengths.
This counts support positions before coefficient cancellation; the cited
model treats cancellation as negligible for the large fields.

The packaged `RingLpnTriple.h` contains a TODO for this sampler and does not
implement it. The factor degree, supported field/ring, and threshold must be
checked together before implementation. Section 8.6 records the resulting
completed proof for programmed-hash correctness.

### 9.6 Finite-$K$ diagnostics

The independent-list exhaustive-ranking diagnostic at $t=d=16$ reports the
following gains over an irrelevant $L^{64}$ exhaustive baseline:

| channel | modeled gain across 256 ideal lists |
|---|---:|
| raw $K=1$ | 166.8 bits |
| occupancy $K=2$ | 115.4 bits |
| occupancy $K=4$ | 85.3 bits |
| occupancy $K=8$ | 69.5 bits |
| feasibility only | 23.7 bits |

These numbers are not losses from a 128-bit decoder. The underlying uniform
exhaustive exponent is roughly 767 bits at $N=2^{16}$. They measure channel
strength only.

Increasing the regular weight while retaining $d=t$ and $2N$ expansion did
not look attractive. From $t=16$ to $t=18$, the model gains about 16 nominal
regular-Prange bits but roughly 26 additional $K=4$ independent-list leakage
bits. This is a heuristic warning, not an attack.

The original instinct that $K=4$ was a natural compromise was superseded by
the full attack study. $K$ cannot be chosen without $(P,t,L)$, and finite
$K$ cannot remove hard feasibility leakage. The paper therefore sensibly
leads with $K=1$ and states the assumption directly.

### 9.7 Remaining attack surface

If further cryptanalysis is commissioned, prioritize:

1. a decoder using the block-circulant Ring-LPN operator or
   leakage-adaptive dual codewords;
2. a global Hall-aware or non-rectangle optimizer beyond local matching-factor
   messages;
3. conditional dual-avoidance experiments on reduced random-code instances;
4. BP/MCMC on the pair-potential planted CSP at $t=16$;
5. tiny exact algebraic systems across multiple stationary expansions;
6. attacks using the actual shared Goldreich feature family, not only ideal
   independent hashes;
7. high likelihood-ratio quantiles or smooth max-information, rather than
   average KL alone.

## 10. Implementation details that are easy to lose

### 10.1 Sparse-DPF expansion

The optimized expansion executes one sparse-DPF protocol; it does not invent a
new multipoint DPF. After routing and correction mapping, it expands cached
leaf seeds and applies the update correction. The eight-lane AES loops in
`CachedDpfExpansion.h` and `SparseDpf.h` reduce register pressure and share
accumulators across active subtrees. Preserve their fixed batching unless a
benchmark justifies a change.

### 10.2 No threads

Peter explicitly rejected multithreading at this stage. Optimize data layout,
batching, allocation, and traversal first. Do not hide performance problems by
adding worker threads.

### 10.3 Coroutine lifetime rule

Treat a coroutine as a heap-allocated state machine. Locals that cross
`co_await` live in the coroutine frame even when their lexical block looks
small. References and spans do not extend referent lifetime. Keep over-aligned
SIMD scratch, large temporary arrays, and locks outside coroutine frames.
Use an owned object or a non-coroutine helper whose lifetime is explicit.

The current sparse expansion helper is intentionally marked non-inlined so the
compiler cannot fold aligned scratch back into the coroutine frame.

### 10.4 Sparse-bucket storage

The packaged RevCuckoo sparse-set builder makes two passes over the public
mapping. The first counts the exact number of addresses in every bucket.
Prefix sums give disjoint offsets into one contiguous buffer, allocated to
the exact total size. The second pass fills those intervals, and each sparse
set is a span over its interval.

There is no probabilistic bucket-capacity bound in this path. The historical
batch-slack fix was superseded by exact sizing. Preserve the contiguous layout,
checked bucket indices, and count/fill agreement. This implementation is at
`libOTe/Dpf/RevCuckooDmpf.h`, starting around line 406 in commit `d6ce62f`.

### 10.5 Benchmark status

Waterfall and Reverse Cuckoo were benchmarked on Peach during development. The
important qualitative result was that Waterfall's post-routing expansion can
match the Reverse-Cuckoo organization; the initial slower $3N$ result was not
inherent and came from parameterization and expansion-layout differences.
After aligning $t=16$, copying the post-routing expansion organization, and
optimizing sparse-set construction, Peter considered performance fast enough
to stop micro-optimizing.

The exact Peach output is not preserved in either repository. The paper
therefore still says that its displayed performance table measures the
Reverse-Cuckoo Ring-LPN path and does not include optimized Waterfall numbers.
Before adding any Waterfall benchmark claim, rerun the benchmark and save the
command, host/compiler configuration, commit, complete stdout, and a machine-
readable result. Run only one benchmark process at a time.

## 11. Paper state and known editorial debt

### 11.1 Current structure

The active paper order is approximately:

1. reusable DPF and DMPF definitions;
2. core MPC subroutines;
3. cuckoo background;
4. Waterfall Cuckoo as the generic construction;
5. complementary bucketing/big-state material;
6. Ring-LPN and stationary PCG discussion;
7. Reverse Cuckoo near the end as the Ring-LPN specialization;
8. evaluation;
9. Waterfall and Reverse-Cuckoo appendices.

The old section titled “Reverse-Cuckoo Leakage and Ring-LPN Analysis” was
hidden/replaced. The current active section is `RevCuckooSpecialized.tex`.

### 11.2 Waterfall rewrite already completed

The main Waterfall section now presents one scheme and one named generator,
`WaterfallGen`. The $4N$ and $3N$ points are profiles, not versions or
separate schemes. Gate-level complete repair, exact recurrences, Hall
calculations, multiplier details, and component OT formulas are in the
appendix. The main body ends with one combined parameter/cost table.

### 11.3 Claims recently corrected

- The abstract and introduction now promote Waterfall as the generic
  construction and Reverse Cuckoo as an LPN specialization.
- The paper's unfiltered correctness exponent is about 41.935 at
  $N_{\mathsf R}=2^{20}$. The proposition now also covers the proposed local
  support-only rejection, with a conservative 40.940-bit bound.
- Evaluation rows $2^{16}$ and $2^{18}$ are explicitly labeled as scaling
  data.
- The Ring-LPN support mapping now says four base support lists per party,
  eight total, not four total.
- The reduced attack point $L=4096$ is identified as
  $N_{\mathsf R}=2^{16}$, not as the $2^{20}$ production correctness point.
- “Insider” terminology was removed from active prose.

### 11.4 Remaining editorial and semantic debt

The paper is not submission-clean. An initial search still finds:

- live Amit/Peter discussion comments in `DPF.tex`;
- older awkward or incorrect prose in `DPF.tex`, including sparse-tree
  explanation and grammar;
- obvious typos such as “seceret” and “come” in `Cuckoo.tex`, “conider” in
  `eval.tex`, and “suprisingly” in `DMPF.tex`;
- inactive but heavily modified `RevCuckoo.tex`, which can confuse reviewers
  and future edits;
- broad changes in legacy sections that have not received a final consistency
  review.

The reusable DPF section deserves a focused pass. Verify that every setup box
actually expands and stores leaf seeds, every expand box uses $\mathsf L_q$
and `updateLeaves`, and no prose still describes returning protocol state.

The Waterfall hybrid simulator is now written explicitly in the main theorem
and Appendix A.8. It samples unconditioned public hashes, simulates the hidden
column opening, preserves local share dependencies, and samples ideal column
DPF outputs conditioned on their prescribed final sum. Eight finite tests
check the key distributional identities. The concrete C++/subprotocol
realizations still need their own transcript checks. In particular, verify:

- what is opened on a failed placement;
- whether all failed transcripts can be charged wholesale to $p_\Gamma$;
- the exact independence argument for the two serial Waksman passes;
- that the candidate generator implemented in `WaterfallHash.h` matches the
  field embedding and coefficient sampling in the paper;
- that dummy-row candidates are independently uniform in the implementation;
- that all named subprotocol simulators provide exactly the outputs used in
  the composed theorem;
- that every cached state has the required marginal failure bound. The union
  bound $B_0p_\Gamma$ itself needs no independence between failures; any
  stronger independence assumption used elsewhere must be identified.

Peter clarified the proof boundary after the narrative audit. In the
semi-honest hybrid proof, use the ideal subprotocols and simulate the honest
messages from the corrupted party's inputs and outputs. Cached implementation
seeds inside an ideal subprotocol are not independently simulated by Waterfall.
The initial audit's blanket caching objection and proposed functionality
changes were withdrawn. Keep the functionality and cached protocol unchanged;
check the reusable DPF realization separately.

The 2026-09-05 follow-up fixed the sparse-DPF paper box and partition helper.
Corrections are generated after expanding the nodes at their actual split
bit; queued children retain that split bit until they consume the correction.
The root correction runs even with empty expansion queues. The pass also
fixed the activity-bit initialization, child numbering, inclusive partition
boundary, leaf-address lookup, correction helper's argument indexing, and
the explicit integer floor in the non-characteristic-two leaf update.
No C++ code, functionality, caching policy, parameter, or cost-table entry
changed.

Proposition `prop:sparse-leaf-invariant` gives an algebraic correctness proof,
including arbitrary-group sign conversion, for valid inputs with |S|>1.
Eleven tests in `experiments/test_sparse_dpf_schedule.py` exhaust small
supports/all active indices, check delayed correction dependencies and exact
tree/correction counts, and exercise repeated payloads over six cyclic groups.
The tests use a toy deterministic PRG, not a cryptographic primitive.

The separate DPF simulator is not completed. Do not claim correctness tests
establish simulation. The concrete expansion relation is
y_p=u_p+(1-2p)gamma*t_p, where u_p is the locally generated provisional map.
The realization proof must handle this relation jointly with the prescribed
ideal output and the actual corrupted-party view. Identify ideal-subprotocol
boundaries before replacing any PRG values: known local PRG inputs do not
permit an independent fresh-value replacement. The prior Waterfall-level
caching objection remains withdrawn; no functionality change is authorized
by this observation. Singleton handling and promised invalid-index rejection
also remain open interface issues. Details are in `PAPER_NARRATIVE_AUDIT.md`.

The Reverse-Cuckoo leakage assumption should eventually become a procedural
game with explicit actors, inputs, public Ring-LPN samples, concrete descriptor
generation, persistent stationary state, preprocessing, and success event.

## 12. Validation and reproduction

### 12.1 Paper

From `C:\Users\peter\repo\dmpf`, build the paper with:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The subsequent Waterfall-simulator pass rebuilt the PDF successfully: 71 pages, stable
references, no undefined references, and the same five pre-existing overfull
boxes. The changed main-text and appendix pages were rendered for review.

The experiment suite mixes `unittest.TestCase` methods with plain pytest
test functions. The previous command
`python -m unittest discover -s experiments -p 'test_*.py'` ran only 25
tests and skipped the 100 top-level test functions.

For complete collection, use a Python environment with pytest and the
experiment dependencies, including NumPy:

```powershell
python -m pytest -q experiments
```

During this review, the system Python was `C:\Python314\python.exe` and had
no pytest. A local environment was created using
`python -m venv --system-site-packages tmp/handoff-review-venv`; pytest was
installed there, and the existing system/user experiment dependencies were
inherited. The actual validation command was:

```powershell
& '.\tmp\handoff-review-venv\Scripts\python.exe' -m pytest -q experiments
```

Initial review: 125 tests passed plus 9 subtests. After adding the eight hybrid
simulation checks and eleven sparse-DPF correctness checks, the full suite
passed **144 tests plus 9 subtests in 5.09 seconds** on 2026-09-05.
The local environment is disposable. A fresh machine must install the
dependencies in its chosen environment before using the generic command.
These tests validate the recorded deterministic models; they do not replace
a cryptographic proof or extend the attack evidence to untested parameters.

Relevant experiment notes are:

- `experiments/waterfall_evictions.md`
- `experiments/waterfall_reachability.md`
- `experiments/rev_cuckoo_leakage.md`
- `experiments/rev_cuckoo_k_tradeoff.md`
- `experiments/rev_cuckoo_lpn_attack_study.md`
- `experiments/rev_cuckoo_regular_isd.md`

### 12.2 Reverse-Cuckoo hash analyzer

From `C:\Users\peter\repo\libOTe-codex-rev-analysis`:

```powershell
python analysis/rev_cuckoo_hash_relations.py --trials 1000 --ring-size 1048576
python -m unittest -v analysis/test_rev_cuckoo_hash_relations.py
```

The report `REV_CUCKOO_HASH_ANALYSIS.md` describes the probability experiment,
relation formula, bounded analyzer, and scope limitations. The
`--width-slack` argument is an analysis override, not a protocol security
parameter.

The rejection-conditioned calculation and combined test suite are:

```powershell
python -m analysis.rev_cuckoo_hash_conditioning
python -m unittest -v analysis.test_rev_cuckoo_hash_conditioning analysis.test_rev_cuckoo_hash_relations
```

All 20 tests passed on 2026-09-05. The new script does not sample supports or
implement the filter; it computes an exact rational correctness bound.

### 12.3 Waterfall C++ tests and benchmark

The frontend accepts the Waterfall benchmark flags
`--waterfallDmpf --profile 3|4`; both built-in profiles require
`--numPoints 16`. Inspect the locally built frontend for the exact executable
path and remaining flags before running. Preserve stdout and run only one
profile at a time.

`libOTe_Tests/Waterfall_Tests.cpp` covers configuration validation, empty
sparse columns, clear placement, hashing, candidates, the basic scan, complete
repair, scatter, and end-to-end expansion for multiple coefficient types.

## 13. Priority audit plan for Astra

### Priority 0: do not lose the current work

1. Read this handoff, `main.pdf`, and the six experiment notes above.
2. Record `git status` in both repositories before editing.
3. Work in the isolated package worktree; preserve its new analysis updates.
4. Do not clean or reset the dirty paper or `codex/foleage-trace` checkouts.

### Priority 1: check the subprotocol realizations

The hashing dependency is resolved in Section 8.6. Return to the paper's main
narrative and correctness task. The Waterfall hybrid simulator is explicit;
the sparse-DPF schedule and reconstructed-map proof are now fixed and tested.
Next settle the remaining DPF interface cases and its separate realization
proof, checking with Peter before changing the model or construction.
Trace one complete Waterfall setup and two expansions from the C++
implementation into the functionality and theorem. Produce a short audit table
with:

- object or message;
- paper name and type;
- C++ representation;
- visibility to each party;
- persistence across expansion;
- failure behavior.

Resolve semantic mismatches before polishing prose. For Waterfall, verify the
column permutation, placement of dummy rows, and replacement of the complete
failed-state transcript in the simulation.

### Priority 2: full Reverse-Cuckoo assumption game

Replace the prose assumption with a complete experiment that includes:

- both parties' $P=4$ regular supports and coefficients;
- the public Ring-LPN samples used by the PCG;
- all $P^2t=256$ product-list descriptors;
- the concrete degree-two evaluator and its shared seeds;
- the semi-honest corrupt party's input and view;
- arbitrary preprocessing of one stationary descriptor set;
- the precise recovery/distinguishing success event;
- weak-factor rejection.

Then ensure every attack experiment is described as evidence about a stated
projection or relaxation of that full game.

### Priority 3: prepare the weak-factor rejection change

The recorded proposal is a local check by each support holder before product
positions or DPF state are derived. Specify that change together with the
conditioned support law and updated hash-correctness accounting, then discuss
the protocol change with Peter before implementing it.

Once agreed, test weights 60 and 61, count residues separately per polynomial,
and check the predicted acceptance probability. Verify applicability to the
field and ring in use. A joint check or a different threshold would require
its own distribution analysis.

### Priority 4: preserve performance evidence

Rerun Reverse-Cuckoo and both Waterfall profiles on Peach, one benchmark at a
time. Use identical $t=16$, domain, set count, coefficient type, compiler,
and expansion count. Record setup time, expansion throughput, setup bytes,
expand bytes, and sparse-expansion subprofile. Only after that comparison
should the paper claim that Waterfall is at least as efficient as Reverse
Cuckoo.

### Priority 5: decide the smaller-ring story

If $2^{16}$ and $2^{18}$ are only scaling rows, the current caveat is
adequate. If they must be deployable 40-bit profiles, evaluate in this order:

1. caller-side rejection of dangerous affine-cube/rank patterns;
2. a minimal cubic feature layer;
3. a sharper joint caller-distribution count.

Increasing $q$ or $q'$ alone cannot suppress deterministic degree-two
identities.

### Priority 6: submission cleanup

After the semantic audits, remove author comments and typos, decide what to do
with inactive `RevCuckoo.tex`, update the implementation/evaluation narrative,
rebuild to stable references, and render every touched protocol/table page.

## 14. Stable terminology and claims

Use these phrases consistently:

- **Waterfall Cuckoo:** generic reusable DMPF with ordinary simulation.
- **Reverse Cuckoo:** $2N$ Ring-LPN specialization with descriptor leakage.
- **profile:** a public Waterfall parameter choice; avoid “version.”
- **candidate matrix $C$:** one candidate bin per row and partition.
- **position map $P$:** all admissible address-to-column positions.
- **placement $O$:** selected row-to-column submatrix.
- **sparse set $S_j$:** domain addresses belonging to physical column $j$.
- **semi-honest corrupt party:** observer terminology.
- **programmed-system correctness:** rank/consistency property of the concrete
  Reverse-Cuckoo evaluator.
- **descriptor leakage:** size-biased matching distribution; separate from
  programmed-system correctness.
- **leakage-robust Ring-LPN assumption:** application-specific assumption, not
  a generic DMPF theorem.

Do not claim:

- that hidden permutation makes Reverse-Cuckoo descriptors uniform;
- that finite $K$ removes feasibility leakage;
- that nonlinear leakage is unusable by linear attacks;
- that 41.935 bits applies to $N_{\mathsf R}=2^{16}$ or $2^{18}$;
- that the 41.935-bit proof automatically survives weak-factor rejection;
- that the 122.482-bit conservative expression is a demonstrated attack;
- that the measured $L=4096$ selector gains have been established at $L=65536$;
- that the AnyField caller is covered by the Ring-LPN relation analysis;
- that weak-factor rejection is already implemented;
- that Waterfall benchmark numbers are in the paper;
- that an ideal-cipher assumption resolves low-degree algebraic identities.

## 15. Recommended next turn

Return to the paper-wide narrative and correctness pass. Peter prioritizes a
clear, accurate paper over further leakage refinements to the already leaky
Reverse-Cuckoo specialization. The hash filter now has its required correctness
bound; do not keep treating that gap as open.

Start with the remaining sparse-DPF pseudocode and subprotocol proof checks;
the Waterfall hybrid now has an explicit simulator. Keep Waterfall generic and central,
and present Reverse Cuckoo later as the $2N$ Ring-LPN specialization. Check
with Peter before larger or questionable changes. The later benchmark task
remains a preserved, one-at-a-time Peach comparison.
