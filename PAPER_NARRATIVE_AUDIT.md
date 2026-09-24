# Paper narrative and correctness audit

## Upload presentation pass (2026-09-23)

Peter approved retaining the historical performance numbers. The abstract,
introduction, and evaluation now also identify the later exact-permutation
and per-list evaluator-seed changes. Performance tables and plotted data are
unchanged; no protocol, assumption, or theorem was revised.

Checked the 20 TeX inputs used by the full build and the extracted PDF text:
no TODOs, placeholders, or author comments occur in the included paper.
Historical draft markers remain in `notes.tex`, which is not included.
Disabled author-comment rendering for the upload build, corrected big-state
prose and its DMPF caption, and removed redundant package/font commands.
Explicit plot bounds remove the empty-axis warnings; the small Dedup+ label
is placed above its bar. Allowing the large big-state figure to float nearby
removes a largely blank page. The protected insertion paragraph is unchanged.

The final 71-page PDF has stable references and no LaTeX warnings, undefined
references, overfull boxes, or duplicate destinations. Rendered pages were
reviewed for layout, with close inspection of the revised front matter,
big-state figure, and performance charts. This is presentation validation,
not a new security audit or benchmark. Next: commit and freeze the upload
snapshot, keeping historical working notes out of any submission bundle.

## Classical-cuckoo cleanup (2026-09-07)

Completed the worked-example follow-up below. `Cuckoo.tex` now defines the
position map and sparse sets before the figure, distinguishes zero-based
addresses/bins from one-based input-row/hash-choice labels, and uses
A=(2,0,5), c=(1,2,3). Corrected all sparse-set lists and the transposed rank
formula. Renamed the fourth displayed matrix to D|_S so S consistently denotes
the sparse sets. The leaf count is at most wN for unrestricted candidates,
with equality for distinct choices, including partitioned layouts. Removed
unqualified concrete load ratios and the suggestion that arbitrary linear
hashes suffice; the section now separates sparse expansion, translation, and
placement. The explicitly protected insertion paragraph is byte-for-byte
unchanged (SHA256 D24CCB7A485D0F6AD70185CAE5407865781BB9396BBB56684BCBA7BDBC2E05CB).

Checked all four matrices directly from TeX: dimensions, placement, payloads,
sparse masks, the caption indices, all four sparse sets, the 16-leaf count,
and the example's translated rank. The 72-page build has stable references
and no overfull boxes; reducing the matrix column spacing removed the previous
11.8429pt overflow. Rendered the section and Waterfall transition on pages
19--22. No protocol, implementation, parameter, or performance number changed.
Next: fix duplicate PDF figure/table destinations, which remain a separate
cross-reference-navigation issue.

## Narrative pass (2026-09-07, after AES synchronization)

Revised the abstract, introduction, technical overview, DMPF roadmap, and
construction/evaluation transitions. The narrative now starts with the
setup/expand task and the tN leaf-work baseline, presents Waterfall as the
generic construction, and introduces Reverse Cuckoo after the Ring-LPN
application in the overview. The abstract explains placement-count leakage
without the unexplained phrase "matching-weight channel." Front-matter
measurement qualifications now also mention the later AES rekeying.

Removed the obsolete conclusion that efficient MPC insertion was unresolved;
the baseline now leads into Waterfall's bounded repair circuit. Scoped the
default 4N comparison to the two main profiles rather than claiming global
state optimality. Clarified that K=1 avoids extra descriptor proposals;
larger K does not improve the 2N expansion factor. The idealized distribution
law is distinguished from the concrete evaluator in the section opening.

This is a presentation pass, not a new proof audit. Protocols, functionalities,
theorem statements, parameters, and numerical results are unchanged. The
72-page PDF builds with stable references and no new overfull boxes. Rendered
pages 1, 2, 5--7, 15--16, 21--22, 28, 41, and 49 were inspected. The existing
11.8429pt overfull box and duplicate PDF figure/table destination warnings
remain. No implementation changes or benchmarks were made.

Next focused work:

- Repair the old classical-cuckoo worked example before polishing the section.
  Its declared zero-based address domain does not match the one-based example
  (A1,A2,A3)=(3,1,6). The bold rows are actually (2,0,5). Column 1's listed
  sparse domain also disagrees with the displayed matrix, and the address-
  translation formula transposes P's row/column indices. These older issues
  were not silently changed during the narrative pass. Preserve the earlier
  "test, do not change" annotated insertion passage.
- Consider moving long core-MPC and auxiliary-construction material out of the
  main development, but ask Peter before changing the section hierarchy.
- Resolve duplicate PDF destination identifiers in a separate layout pass.
  This warning is distinct from undefined LaTeX references.

The chronological checkpoints below are historical; later checkpoints
supersede their page counts and descriptions of remaining implementation work.

2026-09-05. The initial pass was an audit. Following discussion with Peter,
the Waterfall hybrid proof has now been expanded in the paper and checked by
finite-distribution tests. A subsequent pass corrected the sparse-DPF
pseudocode and added its reconstructed-output correctness proof.
Those earlier passes did not change functionality, C++ implementation,
parameters, or benchmark results. The approved CWC revision below subsequently
clarified interface domains and added explicit singleton behavior.

Scope: the abstract, introduction, technical overview, reusable DPF/DMPF
interfaces, Waterfall construction and cost analysis, Reverse-Cuckoo claims,
PCG narrative, and evaluation. Selected implementation paths were checked to
resolve representation and reuse questions. This is not a line-by-line proof
certification of every construction or a new survey of the cited literature.

## Approved CWC revisions (2026-09-05)

**Current checkpoint, 2026-09-07:** the approved sparse activity-mask fix
has passed 16 focused implementation tests in WSL. Paper cost accounting now
includes that mask, two OTs per shared correction selection, and cached sign
correlations. The default 4N / lower-expansion 3N profiles now require
39,052 / 55,358 initial OTs and 942,888 / 3,623,410 logical payload bits
for G=Z_(2^64). Parameters, correctness bounds, and reuse costs are unchanged.
The prose distinguishes these modeled reservation bounds from measured
performance and wire traffic. All 207 paper tests and 9 subtests pass;
the PDF has 71 pages, stable references, and no new overfull boxes. Changed
pages and figures were rendered and checked. Concrete generator transfer
and the full OLE correlation-distribution proof remain open. The checkpoint
below predates the implementation fix and cost reconciliation.

**Latest update, 2026-09-07:** the corrupt-output-share functionality is now
approved. The current Waterfall theorem and Appendix A.8 use direct
aggregation of chosen corrupt DPF shares, not the former uniform-output
fiber sampler. Appendix C gives exact correction-word and adaptive
leaf-update simulation arguments, followed by reductions under explicit
tag-conditioned tree-PRG and joint leaf-PRF requirements (Theorem C.3).
Concrete AES instantiation arguments and the OLE correlation-distribution
proof remain open. The read-only implementation audit also found a missing
secret-activity mask at publicly present, inactive sparse levels; see
`REUSABLE_DPF_GENERATOR_AUDIT.md`. No C++ change was made. The PDF has 70
pages; all 70 targeted tests pass, including the new algebraic mask witness.
Earlier descriptions of fresh prescribed DPF shares below are historical.

This checkpoint supersedes stale descriptions of the text in the historical
findings below. Peter approved the following ten corrections.

1. **Index domains:** retain `[n]=[1,n]`, explicitly define
   `[0,n)={0,...,n-1}`, and use zero-based ranges for the new routing,
   hash, polynomial, and protocol arrays. The candidate graph's right domain
   and hash truncation range were also corrected.
2. **Assumption:** `ass:rc-decisional-lpn` now defines a joint distinguishing
   experiment with concrete descriptors, corrupted-party supports and
   coefficients, retained honest supports, and polynomially many samples.
   Adversarial work includes descriptor preprocessing. Attack work estimates
   are not presented as a proof of this assumption.
3. **PCG algebra:** include the public ring-mask weights. For P polynomials,
   x_p=sum_i r_i e_{p,i} and y=sum_{i,j} r_i r_j e_{0,i}e_{1,j}, with r_0=1.
   All P polynomials enter one sample per party, not P/2 independent pairs.
   The regular product lists include their negacyclic block-wrap signs.
   Reuse is a joint decisional statement, not a marginal one.
4. **DPF boundaries:** fix w-ary child numbering and base-w digit order;
   pad the standard tree and return its first n leaves. Add singleton
   identity resharing, valid-input preconditions, and standard-tree
   correctness. The sparse schedule uses the same little-endian digit
   convention. No uncosted input-membership circuit is claimed.
5. **Shuffle:** define the shared shuffle's ideal interface, retained
   permutation on reuse, fresh output sharings, and call-order behavior.
   Define the privately controlled permutation and identity resharing used
   by the protocol boxes.
6. **Group-By:** use ceil(log2(t+1)) count bits, saturate padding at capacity,
   and define overflow truncation. The no-overflow opening is a uniform
   permutation of a fixed multiset of joint (label,selection-bit) records.
   State the output-materialization and grouping costs.
7. **Performance:** scope ratios to the measured Reverse-Cuckoo rows and
   individual baselines. Correct Goldilocks, MB/GB, and throughput claims;
   remove unmeasured multicore and all-PCG claims. Do not claim lower
   Goldilocks expansion bytes than the sum-of-DPFs baseline.
8. **Cost semantics:** distinguish initial selected-message OT correlations,
   logical OT payload, per-expansion masked messages, and wire traffic.
   Explain the cached deduplication controls. Preserve all profile totals.
9. **Narrative:** rewrite the abstract, introduction, overview, and application
   summary around generic Waterfall and the later leaky 2N specialization.
   Bounded complete repair is not described as eviction-free insertion.
10. **Proof scope:** retain the explicit Waterfall ideal-subprotocol simulator;
    separate it from the concrete DPF realization and full PCG composition.
    Correct the repair-box loop nesting, first-parent record type, and stash
    recomputation. The DPF reconstructed-output proof is not a simulator.

### Resolved mask-schedule mismatch; historical timings retained

**Later resolution:** Peter confirmed the intended fresh-mask behavior.
The local `codex/dmpf-package` worktree now refreshes masks in every
`RingLpnTriple::expand` using two fresh seed contributions. All eight
Ring-LPN test groups and the additional P=4 stationary checks pass in WSL.
The tests verify mask agreement/freshness, refreshed mask products, unchanged
supports, and correlation correctness across three fields and both DMPF
backends. The fix and a standalone benchmark are committed as `b2e9b7b`
on `codex/dmpf-package`. Historical paper benchmark numbers have not been
rerun under the same measurement configuration. The following records the
original mismatch, not the current state of that worktree.

A sequential WSL check on the i7-13700H used Goldilocks, P=4, t=16,
N=2^20, GCC 13 / O2, and one thread for both parties. Median repeated
expansion took 3,761.20 ms; two isolated mask-refresh calls took 322.084 ms.
The latter includes all P^2 mask products, with storage already allocated.
Their ratio is 8.56%, not a measured before/after slowdown. Six expansions
passed full output-correlation checks. Locally generated base correlations
were excluded from timing. This does not replace the paper's O3/LTO,
single-thread-per-party benchmark rows. Source and raw timings are in the
package worktree's `analysis/ringlpn_mask_bench.cpp` and
`analysis/ringlpn_mask_bench_2026-09-05.json`.

The follow-up paper pass distinguishes historical timings from the current
mask-refresh implementation in the abstract, introduction, and evaluation.
It also separates the silent PCG expansion interface from the interactive
reusable batch procedure. The stationary section explains mask sampling
directly in NTT coordinates and its per-round work. No protocol, assumption,
parameter, or performance-table number changed in this editorial pass.
The rebuilt PDF remains 66 pages, with stable references and no undefined
references or citations. The only overfull box remains the 11.8429pt cuckoo
matrix. Visual review covered pages 1, 2, 6, 37--40, and 48--51. The PCG
chart's small-segment labels were moved above the bar to remove overlap;
both chart captions were shortened without changing data. The log also
reports duplicate figure/table PDF destinations; hyperlink-anchor cleanup
remains a separate submission task.

Read-only inspection of the package worktree's
`libOTe/Triple/RingLpn/RingLpnTriple.h` found that `sampleA` makes its first
public mask polynomial the identity and samples P-1 other masks. `init`
calls it with a fixed seed; `expand` reuses those masks. The paper's repeated
sample experiment uses fresh public mask vectors across expansions. The
evaluation now states this mismatch rather than claiming that the prototype
implements the reusable security model. No C++ change was made.

This distinction is material: with both masks and supports fixed, the honest
ring samples lie in a fixed subspace generated by at most Pt shifted mask
polynomials. If complete samples are exposed for more than Pt rounds and
N_R>Pt, a rank test can distinguish them from independent uniform ring
samples. Do not claim an attack on the corrupt party's protocol view without
checking which values that view exposes. Resolve the required correlation
security notion and mask schedule before claiming end-to-end reusable PCG
security. The local support-rejection filter remains a separate unimplemented
prototype requirement already documented earlier.

### Validation and remaining work

- `latexmk` builds the 66-page PDF with stable references and no undefined
  references or citations. The sole remaining overfull box is the preexisting
  11.8429pt cuckoo-matrix display; no new overfull boxes remain.
- Rendered the final PDF and inspected the revised DPF, DMPF, Waterfall,
  scatter, Group-By, Reverse-Cuckoo, and appendix protocol figures, the
  combined table, front matter, application/assumption pages, and evaluation
  transitions. Function declarations are left aligned and followed by line
  breaks; the revised figures have no clipping.
- 153 tests and nine subtests pass. New tests check standard w-ary trees,
  non-power domains, singleton resharing, Group-By boundaries and joint
  openings, negacyclic product-list signs, and two-/four-polynomial mask
  weights. They are algebraic checks, not cryptographic security tests.
- The 4N totals remain 35,980 OTs / 363,304 payload bits / 43.37 batch bits.
  The 3N totals remain 42,382 / 1,268,210 / 41.13.
- The concrete reusable-DPF honest-message simulator and full PCG composition
  remain proof work. The published analysis still distinguishes ideal-hash
  attack experiments from the concrete evaluator and larger correctness point.

## Historical assessment

**Subsequent DPF realization audit:** `REUSABLE_DPF_SIMULATION_AUDIT.md`
records an explicit online output-consistency test. For any three leaves,
the corrupt party has two equal local tags. Their output difference equals
the provisional difference computed before the payload arrives, regardless
of the opened scalar correction. An independently sampled uniform ideal
map matches a previously recorded difference with probability 1/|G|.
This concerns the concrete DPF realization boundary, not Waterfall's ideal
subprotocols or a hidden-index attack. The paper and protocol were left
unchanged pending discussion of the output-sharing contract and simulation
chronology. All 37 targeted boundary, scheduling, and hybrid checks pass.

Keep the current construction. Waterfall has a coherent core: public hashes
first, a uniform candidate matrix, bounded deterministic placement, a hidden
column permutation, and reusable sparse DPFs. The two main profiles belong
together. Reverse Cuckoo should remain the later, explicitly leaky $2N$
Ring-LPN specialization.

The Waterfall hybrid proof is now explicit. Cached seeds inside an ideal
subprotocol are not a Waterfall-level objection; the initial finding
conflated the hybrid proof with a separate subprotocol-realization proof.
The sparse-DPF correction-scheduling issue below is now fixed and tested.
Its realization-level simulation remains separate. No attack on the C++
DMPF has been established by this audit.

## 1. Revised: separate the Waterfall hybrid and DPF realization proofs

Sources: [DPF functionality](C:/Users/peter/repo/dmpf/DPF.tex:49),
[DPF simulation definition](C:/Users/peter/repo/dmpf/DPF.tex:70),
[leaf correction](C:/Users/peter/repo/dmpf/DPF.tex:202),
[DMPF functionality](C:/Users/peter/repo/dmpf/DMPF.tex:42),
[Waterfall simulation](C:/Users/peter/repo/dmpf/WaterfallCuckoo.tex:489).

Peter clarified the intended semi-honest proof: use the ideal subprotocol
interfaces, prove correctness, and simulate the honest party's messages from
the corrupted party's inputs and outputs. Protocols and functionalities may
retain internal state. Waterfall does not inspect or independently simulate
the implementation state of an ideal reusable DPF.

The earlier fixed-seed consistency example was not a distinguisher for that
Waterfall hybrid. It does not justify changing the functionality, making
corrupted output shares simulator-selected, or adding fresh vector masks.
Those recommendations are withdrawn.

The expanded proof now handles prescribed ideal DMPF output shares directly.
Column DPFs return independent uniform maps in the hybrid. The simulator
samples their shares uniformly conditioned on their sum equaling the
prescribed output. At each address, one first-partition entry completes the
sum, and all other entries are uniform. No DPF seed is accessed.

The reusable DPF's own realization still requires its honest-message
simulator and correctness proof. That is a separate obligation, not
established by the Waterfall composition argument and not refuted merely by
caching. Finding 2 records the completed pseudocode correction and its
correctness proof.

CWC-04, CWC-08, CWC-11: keep the proof boundary and the simulated interfaces
explicit.

## 2. Resolved: sparse-DPF correction scheduling

Sources: [sparse setup](C:/Users/peter/repo/dmpf/DPF.tex:347),
[empty-level guard](C:/Users/peter/repo/dmpf/DPF.tex:368),
[leaf correction lookup](C:/Users/peter/repo/dmpf/DPF.tex:386).

Trace $S=\{0,1\}$ at depth $D=1$. The root splits at depth 1. Both children
are singleton intervals and are queued in $u_0$. Consequently $u_1$ is
empty. The guard skips depth 1, but leaf processing subsequently reads the
parent correction $\sigma_{1,j}$, which was never computed.

Independently, the correction step assigns an unindexed $\sigma$, while
later steps read the depth-indexed table $\sigma_{\rho,j}$. The box needs
to store the corrections by their producing depth.

The follow-up trace found additional indexing errors: the child sums were
assigned to depth d-1 instead of their split bit d; v started at one instead
of zero outside the root; root child numbers mixed zero- and one-based
indexing; the partition boundary included one extra address; and leaf storage
used a rank/interval instead of the actual address in S.

**Completed correction:** initialize v to zero except at the root split.
At each public split bit d, first process u_d using the already stored
parent corrections. Accumulate new child sums at d, then generate and store
sigma_d. Compute the root correction even when its children are both leaves.
Every child records its parent's split bit and consumes that correction at
its own split, or during final leaf processing.

The partition helper now has explicit domains and correct inclusive
intervals. The correction helper uses one-based child numbering and its
declared alpha argument. The non-characteristic-two update writes the
required integer floor and permits c_p=|S|.

Proposition `prop:sparse-leaf-invariant` proves the leaf seed/tag invariant
and reconstructed-map correctness for every valid alpha in S, |S|>1,
all setup randomness, and each payload. The proof includes the sign
conversion for arbitrary additive groups.

Eleven tests in `experiments/test_sparse_dpf_schedule.py` cover every
support of size at least two in [0,8), every active address, and three
deterministic randomness tapes; deferred corrections; exact split-depth
counts; and repeated expansion over cyclic groups of orders 2,3,4,5,8,17.
They use a toy deterministic PRG and ideal shared circuit operations.
They are correctness tests, not a cryptographic simulation proof.
The full suite now passes 144 tests plus 9 subtests.

**Separate simulation work remains:** at the DPF realization boundary,
identify which local values belong to the corrupted party's view and which
state remains inside ideal subprotocols. Preserve locally computed PRG/PRF
values wherever their inputs are available to that party. For an expansion,
if u_p is its provisional map and t_p its local tag vector, its returned map
is y_p=u_p+(1-2p)gamma*t_p. The honest-message simulator must account for
that equation jointly with its prescribed ideal output, not just show that
the two parties' outputs sum correctly. No blanket fresh-value hybrid or
new functionality has been inserted to bypass this obligation.

CWC-03, CWC-10: an invoked object must exist on every reachable path.
This is a paper pseudocode finding, not evidence that the implementation
has the same bug.

## 3. P1: separate the generic contribution from measured specialized gains

Sources: [abstract](C:/Users/peter/repo/dmpf/abstract.tex:4),
[introduction](C:/Users/peter/repo/dmpf/intro.tex:16),
[evaluation scope](C:/Users/peter/repo/dmpf/eval.tex:98),
[OLE table](C:/Users/peter/repo/dmpf/eval.tex:112),
[PCG extensions](C:/Users/peter/repo/dmpf/pcg.tex:36).

The abstract moves from generic Waterfall to order-of-magnitude performance
and communication gains across fields and rings. The evaluation reports a
particular Reverse-Cuckoo prime-field pipeline. Those are different scopes.
For example, the sum-of-DPFs Goldilocks row has 12.4 MB per expansion versus
12.99 MB for Reverse Cuckoo: it supports a throughput gain, not a communication
reduction against that baseline. The lower-communication comparison is against
a different baseline.

**Correction:** attach each measured claim to its construction, field,
baseline, hardware, and security-qualified parameter point. Distinguish:

- Waterfall's proved generic construction and modeled costs;
- the implemented Waterfall DMPF, whose measurements are not in these tables;
- the measured Reverse-Cuckoo Ring-LPN pipeline;
- possible transfers to other PCGs, without claiming those transfers were
  implemented or their reuse assumptions proved here.

Keep the $2^{16}$ and $2^{18}$ performance-scaling rows separate from the
proved $2^{20}$ hash-correctness point. Remove the unsupported near-linear
multicore-scaling footnote unless an experiment is supplied. Do not run new
benchmarks as part of this editorial pass.

CWC-08, CWC-09: scope of evidence and claims.

## 4. P2: replace the obsolete anti-eviction narrative

Sources: [introduction](C:/Users/peter/repo/dmpf/intro.tex:9),
[overview motivation](C:/Users/peter/repo/dmpf/overview.tex:35),
[Waterfall placement](C:/Users/peter/repo/dmpf/WaterfallCuckoo.tex:23).

The introduction says the construction eliminates iterative insertion; the
overview says the alternatives avoid iterative eviction entirely. Both main
Waterfall profiles perform complete augmentations that move existing rows.
The actual improvement is an efficient fixed-budget circuit with explicit
cost and failure accounting, not the absence of repair.

Use a positive description: an inexpensive forward scan followed by a small
public number of canonical repairs. Explain why that circuit is economical
in total nonlinear operations; do not center the motivation on circuit depth.

The Waterfall intuition also contrasts itself with Reverse Cuckoo before the
latter is defined. Explain hash-first sampling positively here and make the
comparison when Reverse Cuckoo is introduced.

CWC-01, CWC-02, CWC-05: presentation order and stable construction description.

## 5. P2: distinguish reusable correlations from online work and wire traffic

Sources: [Waterfall state/API](C:/Users/peter/repo/dmpf/WaterfallCuckoo.tex:138),
[cost-table scope](C:/Users/peter/repo/dmpf/WaterfallCuckoo.tex:521),
[common costs](C:/Users/peter/repo/dmpf/WaterfallCuckooAppendix.tex:235),
[model value masks](C:/Users/peter/repo/dmpf/experiments/waterfall_mpc_cost.py:152),
[cached implementation sessions](C:/Users/peter/repo/libOTe-codex-rev-analysis/libOTe/Dpf/WaterfallDmpf.h:451).

The table calls its totals one-time setup costs. Its common term includes
$t(t+1)/2$ payload-width products, although payloads arrive during Expand.
This need not invalidate the OT-count model: controls fixed during setup can
have reusable OT correlations. The implementation does cache multiplication
sessions for deduplication and consumes them during expansion.

The paper needs to identify that persistent resource. Its stored equality
matrix alone is not a description of the reusable multiplication mechanism.
Separate the number of initial correlations from online payload operations
and from actual bytes exchanged. Logical OT payload lengths are not total
network traffic, particularly when seeds are expanded across many rounds.

**Correction:** retain the two cost coordinates and existing profile numbers
pending a phase-by-phase reconciliation. State what is cached for duplicate
aggregation, what messages remain per expansion, and which costs the table
excludes. Do not infer fresh OTs per expansion or add a new routing protocol.

CWC-03, CWC-04, CWC-08: state, interface, and cost semantics.

## 6. P2: make the claimed domains complete

Sources: [hash output restriction](C:/Users/peter/repo/dmpf/WaterfallCuckoo.tex:80),
[arbitrary-domain statement](C:/Users/peter/repo/dmpf/WaterfallCuckoo.tex:116),
[sparse-DPF domain](C:/Users/peter/repo/dmpf/DPF.tex:347).

Waterfall uses a field of dimension $\ell=\max(1,\lceil\log_2N\rceil)$
and requires $\log_2d_s\le\ell$. Its generic profile declaration does not
state that restriction. For example, the $d_s=128$ profile is not defined
by this hash family when $N=2$.

The sparse-DPF box requires $|S|>1$, but Waterfall invokes it for every
nonempty sparse set, including singletons. The general DPF API also promises
invalid-index rejection that its boxes do not implement.

**Correction:** state admissible public parameters and explicit singleton
behavior, and choose either valid-input preconditions or implemented
validation. Increasing the embedding field to accommodate every partition
is possible, but would change the generic cost formula and should be discussed.
These issues do not change the displayed large-domain profile calculations.

CWC-03, CWC-04, CWC-10: domains and reachable branches.

## 7. P2: keep evidence qualifications throughout the later sections

Sources: [attack coverage](C:/Users/peter/repo/dmpf/RevCuckooSpecialized.tex:506),
[bucketing complexity](C:/Users/peter/repo/dmpf/Bucketing.tex:10),
[PCG overview](C:/Users/peter/repo/dmpf/overview.tex:100),
[unbounded-reuse wording](C:/Users/peter/repo/dmpf/overview.tex:108).

- The attack table and surrounding prose still use “production” for numbers
  measured at $L=4096$, distinct from the $L=65536$ correctness point.
  Carry over the handoff's parameter tuples and distinguish attacks,
  optimistic model costs, and limited diagnostics. Keep the closed hashing
  result closed; the introduction should now qualify 41.935 as unfiltered
  or simply state the filtered 40-bit guarantee.
- Bucketing claims total computation and output concatenation take $O(t)$,
  despite producing $N$ outputs. State at least the $N$ output term and
  express communication through the invoked DMPF and grouping costs.
  Its prose also reverses quotient and remainder relative to the box.
- The PCG overview confuses the number of polynomials, support size, and
  evaluation work. Replace the purported “$O(Nt)$-sparse structure” by the
  actual product-list decomposition, with separate $P,t,N_{\mathsf R}$ roles.
- Change “arbitrarily many expansions” to the stated polynomial-reuse scope.
  Do not infer stationary Ring-LPN security from a restricted SSD attack model.
- The introduction gives Goldilocks as $2^{64}-2^{32}-1$, whereas the evaluation
  correctly uses $2^{64}-2^{32}+1$. Fix the introduction.

CWC-05, CWC-07, CWC-08: notation, quantities, and evidence scope.

## Waterfall simulation: retain the core argument, write the hybrid explicitly

The column-permutation and physical-dummy-placement description is now right.
On a successful placement, a fresh uniform hidden column permutation makes
the opened labeled injection uniform for each fixed support and hash.

For the proof, keep the hashes unconditioned. Replace the opening by a fresh
uniform injection for every hash sample; couple to the real opening on success
and charge the whole failed transcript once. Do not instead assert that
conditioning on routing success leaves the public hashes uniform.

The phrase “conditioned on either party's local view” should identify the
pre-opening ideal-subprotocol view. The full view includes the opening itself
and cannot leave that same opening uniformly undetermined.

These steps are now explicit in the main theorem and
the appendix subsection labeled `app:waterfall-simulator`. The appendix follows the setup calls,
constructs the honest opening message, preserves local linear dependencies,
and gives the exact conditional sampler for column-DPF outputs. It also
explains the once-per-state failure charge and adaptive expansions.

Eight tests in
`experiments/test_waterfall_hybrid_simulation.py` check the finite
distributional identities: either corrupted permutation, opening-share
completion, exposed unmatched labels, uniform output fibers over groups of
orders 2, 3, and 4, empty columns, and an adaptive continuation.
These are checks of the hybrid argument, not security tests of the concrete
subprotocols.

## Recommended narrative and edit order

**2026-09-07 interface investigation:** see
`REUSABLE_DPF_SHARING_CONTRACT.md` before revising the DPF/DMPF definitions.
The candidate contract simulates the corrupt output rather than prescribing
an independent uniform map. This preserves a substantive input-privacy
requirement and offers a conditional Waterfall composition route without
new runtime messages. The concrete realization and application's marginal
OLE distribution still need proofs. No definition change has been approved.

Main narrative:

1. Fix supports once; provide new payloads to each expansion.
2. Reusable sparse DPF: the leaf-processing building block.
3. A short cuckoo-layout motivation: placement is the distributed bottleneck.
4. Waterfall: normalize, sample candidates, place, scatter/cache; one generator,
   two public profiles, one correctness and simulation argument.
5. PCG integration: what the generic DMPF replaces and what reuse assumes.
6. Reverse Cuckoo: the $2N$ specialization, its exact idealized leakage law,
   concrete-evaluator distinction, leakage assumption, and qualified evidence.
7. Evaluation: name what was measured and keep analytic estimates distinct.

Keep bucketing and big-state as clearly secondary material. Moving their long
derivations to the appendix is a presentation choice to confirm with Peter,
not a prerequisite for repairing the two main constructions.

Next edit order: finish the DPF interface boundaries and check their separate
realization proofs; reconcile costs; then rewrite
the front matter and tighten the application/evaluation claims. Ask before
changing functionality semantics, protocol behavior, parameters, or the
major section hierarchy.
