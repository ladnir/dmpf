# Reusable DPF: generator reductions and implementation transfer

2026-09-07. The original pass was paper proof work and read-only implementation
inspection. Peter subsequently approved fixing the implementation.

**Approved encoding repair:** the selected correction's low bit is now
cleared before transmission; corrected tree seeds pack the logical tag.
Cached outputs contain only the 127-bit prefix, and the immediate-payload
API hashes that prefix. See AUD-213 and
`analysis/sparse_dpf_correction_encoding.md` in the implementation worktree.
This supersedes the unresolved encoding status in the original audit below.
The dense AES generator/export and reusable key-schedule proofs remain open.

**New implementation-transfer finding:** `AES_DPF_IMPLEMENTATION_AUDIT.md`
records a separate confirmed defect: sparse setup opens the selected seed
correction's low bit as well as both tag corrections. Their algebra reveals
the index bit whenever the tag corrections differ. At a directly sampled
root this occurs with probability 1/2; the two index-bit opening laws have
statistical distance 1/2. Nine new exact tests pass. The activity mask does
not fix this encoding error. No C++ fix was made in this audit. Dense
correction encoding already uses the paper's packed format; its custom AES
generator and transformed seed export still need separate justification.
See the new note for the concrete equations and the cached-leaf audit.

**Paper cost checkpoint:** the main and appendix Waterfall tables now match
the revised reservation model, including two OTs per shared correction
selection, residual sparse activity masks, and cached non-characteristic-two
sign correlations. For the stated 64-bit cyclic group, the default 4N and
3N totals are respectively 39,052 / 55,358 OTs and 942,888 / 3,623,410
logical payload bits. These are setup-correlation bounds, not measured
wire traffic. The full paper suite passes 207 tests and 9 subtests; the
71-page PDF was rebuilt and changed pages visually checked. The concrete
generator arguments below remain open.

**Implementation checkpoint:** the package worktree now separates public
`mHasSplit` from secret `mActivity` and masks both child sums using the
existing shared bit-times-256-bit-string primitive. All 16 focused WSL tests
passed, including the new activity-mask regression and Waterfall/Reverse-Cuckoo
integration tests. See `libOTe-codex-rev-analysis/analysis/sparse_dpf_activity_mask.md`
and AUD-212 in that worktree for details. The additional setup cost is one OT
per direction per instance and scheduled sparse level; cached expansion is
unchanged. The paper cost checkpoint above supersedes historical setup totals.
The concrete generator arguments below remain open.

## What is now proved in the paper

Appendix C (`ReusableDpfAppendix.tex`) gives a conditional realization theorem
for the paper's standard and sparse DPF protocols. It uses the approved
functionality in which the corrupt party selects its output map and the
honest party receives the correct complement.

The sufficient primitive requirements are:

- For each known tag theta, G(K || theta) is pseudorandom for uniform secret
  (kappa-1)-bit K. The tag is selected before the challenge. This is stronger
  than an unqualified claim about a uniformly sampled full input block.
- L(q, K) is a group-valued PRF under uniform secret K, jointly for all
  expansion counters q. The public family remains evaluable under known keys.

The reduction replaces honest active tree expansions in root-to-leaf order.
At each replacement, earlier random-child corrections leave the next honest
active prefix uniform and unexposed. The reduction embeds an unknown PRG
challenge seed there; the parent is not used after its expansion. It knows
the honest inputs to run the surrounding experiment, unlike the final
privacy simulator. A single PRF challenge then covers all expansions at the
final honest active leaf. Both reductions preserve joint honest outputs and
the corrupt view, including adaptive payload selection.

The bound is a_S * epsilon_G + epsilon_L(Q) + epsilon_circuit(Q), where a_S
counts active-parent generator calls after the directly sampled root.
It is at most D-1 for standard DPF and min(D-1, |S|-2) for sparse binary DPF.
The singleton uses only ideal resharing. There is no additional factor Q on
the joint PRF advantage. The paper does not claim a numeric AES bound.

The protocol already allowed the PRF leaf-conversion option. This pass
selects that option for the theorem; it does not change protocol steps,
introduce fresh public coins, or prescribe independently uniform output maps.
An arbitrary round-indexed family needs its own joint-security justification.

## Implementation discrepancy: public occupancy is not secret activity

Inspected worktree: `C:/Users/peter/repo/libOTe-codex-rev-analysis`.
Relevant source: `libOTe/Dpf/SparseDpf.h`.

In `expandSparseLevel`, `tree[level].mC` is set to one when the public queue
has any nodes (`if (size)`, around lines 512-517). The correction loop masks
child sums only when `tree[d].mC == 0` (around lines 776-795). Despite the
comment describing the paper's `(1 xor v_d) * r`, this is not a computation
of the secret active-split indicator v_d.

The paper's indicator is the reconstructed XOR of parent tags. It is zero
when every parent splitting at that level is off-path, even if public nodes
exist. The paper then masks both reconstructed child sums with independent
random blocks. Correctness is unaffected because all those parents have
equal tags in the two parties' states.

### Exact small witness

Take domain [0,8), public support S = {0,1,4}, active address alpha = 4,
and no dense prefix. The root splits at bit 3. Its left subtree {0,1} splits
at bit 1, whereas its right subtree {4} is already a leaf.

At bit 1, a public splitting node exists but it is off-path. Its two corrected
seed copies are equal, so deterministic expansion gives identical children
to both parties. Both reconstructed child sums are therefore exactly zero.
The C++ public-occupancy test omits the mask. The resulting opening has raw
sigma = 0 and tau = (1,0), since alpha's low bit is zero.

This equality holds for every deterministic generator, not just a toy PRG.
In the paper's packed correction representation it is the pair (1,0).
In the independent-child experiment at an active split, that packed pair
has probability 2^(-(kappa+1)); at an inactive split the paper's mask restores
the same distribution. This last probability is a hybrid calculation, not
a proved probability for the custom AES implementation.

The regression in `experiments/test_sparse_dpf_schedule.py` omits the mask
in an optional toy-model branch, checks the deterministic opening above,
and checks that reconstruction remains correct. It is an algebraic witness,
not an end-to-end C++ run. No production leakage rate or support-recovery
attack has been measured. Dense prefixes can remove this particular example,
but public occupancy and secret activity still differ in uneven residual
sparse subtrees. Production reachability and the least-cost fix need auditing.

This was a setup-privacy issue. The approved output-share functionality did
not remove it. The implementation checkpoint above resolves the mask mismatch;
the generator/representation obligations below still limit theorem transfer.

## Concrete evaluator inventory: proof obligations, not attacks

1. **Sparse tree evaluator.** `SparseDpf.h` uses fixed-key AES hashing on
   `seed xor ZeroBlock` and `seed xor OneBlock`, including the feed-forward
   XOR in `hashBlock`. The implementation stores a raw seed block and a
   separate corrected logical tag. The raw seed's low bit is not simply the
   corrected tag used in the paper's combined representation. Establish the
   representation correspondence before applying the conditioned-tag theorem.
2. **Dense tree evaluator.** `RegularDpf.h` computes `tmp = AES_K(seed)` and
   derives children as `AES::roundEnc(tmp, seed)` and `tmp.add_epi64(seed)`.
   This is not two independent ideal-cipher evaluations. Ideal AES alone
   does not establish the required generator property for this postprocessing.
   No distinguishing attack on this evaluator is established by this audit;
   the theorem's sufficient PRG condition need not be necessary for DPF security.
3. **Reusable leaf evaluator.** `CachedDpfExpansion.h` hashes saved raw seeds
   under a public, deterministically advanced AES key. It advances the public
   key using another AES hash on a fixed block. The paper instead expresses
   the sufficient property as a PRF keyed by the hidden saved prefix. An
   ideal-cipher argument must justify the implemented family, query bounds,
   conversion to the output group, key-schedule collisions, and domain
   separation. A PRF notation change alone does not supply that argument.

The public Ring-LPN masks refreshed by the earlier implementation commit
are unrelated to this leaf-conversion key schedule. Do not conflate them.

## Recommended next work

The activity-mask fix and cost reconciliation are complete. The next priority
is the sparse correction encoding defect documented in
`AES_DPF_IMPLEMENTATION_AUDIT.md`, followed by the dense evaluator/export
decision and a joint ideal-cipher argument for setup and cached expansion.
Do not yet claim that the concrete implementation realizes the generic
functionality proved for the abstract protocol.
