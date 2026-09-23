# Composing reusable DMPFs with Ring-LPN OLE generation

Standalone proof note, 2026-09-23. This document is not included in the paper's
PDF. It records a conditional composition lemma for the existing protocol;
it proposes no protocol, functionality, or implementation changes.

**Audit result:** the simulator bridge below closes in the paper's
ideal-subprotocol model. A source cross-check at libOTe commit `c2fa1bc`
found distribution mismatches, recorded in Section 7. The user approved
exact-uniform permutations and fresh per-list evaluator seeds; the working
implementation now contains these changes. Validation is recorded below.

## 1. Question and scope

The reusable DMPF functionality lets the simulator choose the corrupt party's
output share. The functionality supplies the honest complement. Does this
contract suffice for the Ring-LPN application?

For security of the honest OLE output given one party's view, the answer is
yes, provided the distributed computation has the simulation property stated
below. The argument uses the existing joint Ring-LPN assumption. Cached seeds
and nonuniform DMPF output shares do not obstruct this argument.

We consider one static semi-honest corruption and a polynomial number of
sequential expansions. Coefficients and masks follow the paper's sampling
rule; this is not a claim for adversarially chosen Ring-LPN noises. Expansion
is interactive, as in the implementation, not communication-free PCG expansion.

## 2. Construction and correctness

Let $R=\mathbb F_q[X]/(X^{N_R}+1)$, where $2N_R\mid(q-1)$, and let
$\mathsf{NTT}:R\to\mathbb F_q^{N_R}$ be the negacyclic NTT. It is a bijection
that maps multiplication to coordinate-wise multiplication, denoted $\odot$.
There are $P\geq2$ noise polynomials per party and $Q$ expansions.

Party $p\in\{0,1\}$ samples its supports once. At expansion $u\in[0,Q)$,
it samples fresh coefficients to obtain $e_{p,i}^{(u)}\in R$, for $i\in[0,P)$.
Supports use the specified regular distribution and, for the analyzed
Reverse-Cuckoo profile, the local support filter. Coefficients are independent
uniform field elements, including zero, as specified by the paper's assumption.
Here supports mean the retained selected positions, not necessarily the
nonzero support of each expansion's polynomial. The support filter tests
those positions and does not reject coefficients.

For each expansion, the parties sample public masks
$r_0^{(u)}:=1$ and $r_i^{(u)}\gets R$ for $i\in[1,P)$.
The masks are fresh across expansions and independent of the noises.
Each party computes

$$
x_p'^{(u)}:=\sum_{i\in[0,P)}r_i^{(u)}e_{p,i}^{(u)}.
$$

For each pair $(i,j)\in[0,P)^2$, the DMPFs and public block-assembly operations
produce ring shares $v_{p,i,j}^{(u)}\in R$ satisfying

$$
v_{0,i,j}^{(u)}+v_{1,i,j}^{(u)}=e_{0,i}^{(u)}e_{1,j}^{(u)}.
$$

This equation includes duplicate-address aggregation, block shifts, and
negacyclic signs. Define

$$
z_p'^{(u)}:=\sum_{i,j\in[0,P)}r_i^{(u)}r_j^{(u)}v_{p,i,j}^{(u)},
\qquad
(x_p^{(u)},z_p^{(u)}):=
(\mathsf{NTT}(x_p'^{(u)}),\mathsf{NTT}(z_p'^{(u)})).
$$

Distributivity gives $z_0'^{(u)}+z_1'^{(u)}=x_0'^{(u)}x_1'^{(u)}$.
Consequently, the parties' outputs satisfy the packed OLE relation

$$
z_0^{(u)}+z_1^{(u)}=x_0^{(u)}\odot x_1^{(u)}.
$$

No uniformity of the product shares is used in this calculation. Let
$\delta$ bound the probability of any correctness failure across this
setup and its $Q$ expansions. A bad placement or programmed hash is a
setup event; reusing it does not multiply its probability by $Q$.

## 3. Security target and simulation premise

Fix the corrupt party $p$ and put $h:=1-p$. Let $V_p$ be its complete view:
private coins, retained state, public values, messages, and outputs throughout
the execution. In particular, $V_p$ determines every $(x_p^{(u)},z_p^{(u)})$.

For a view $V$ with these output fields and ring elements
$a=(a^{(u)})_{u\in[0,Q)}\in R^Q$, define its honest completion by

$$
\mathsf{Complete}_p(V,a):=
\left(V,
\left(\mathsf{NTT}(a^{(u)}),\;
x_p^{(u)}\odot\mathsf{NTT}(a^{(u)})-z_p^{(u)}\right)_{u\in[0,Q)}
\right).
$$

The security target keeps the corrupt view and replaces the honest output by
this completion with independent $U^{(u)}\gets R$. These replacements are
independent of the view. The honest multiplicative inputs are therefore
uniform; the honest product shares are determined by correctness.

Here are the two premises needed for the composition.

**Distributed simulation.** Let $W_p$ contain the corrupt supports and
coefficients, the public masks, and the permitted setup leakage. There is a
stateful simulator $\mathcal S_p$ whose generated view
$\widetilde V_p:=\mathcal S_p(W_p)$ satisfies

$$
(V_p,(x_h'^{(u)})_{u<Q})
\approx_{\epsilon_{\rm sim}}
(\widetilde V_p,(x_h'^{(u)})_{u<Q}).
$$

The simulator does not receive the honest noises or ring samples. It processes
setup and expansions in order, retaining its own state. The displayed
comparison includes the honest samples as auxiliary outputs; simulation of
the marginal transcript alone would not suffice. The notation
$\approx_\epsilon$ bounds distinguishing advantage by $\epsilon$ at the
considered computational budget.

For Reverse Cuckoo, the leakage in $W_p$ is the tuple $\mathsf H$ of public
evaluator seeds and programmed descriptors, jointly across all product lists.
For Waterfall there is no support-dependent descriptor leakage; independent
public setup randomness can be generated by the simulator.

**Joint Ring-LPN.** For this same sampling distribution and leakage, assume

$$
(W_p,(x_h'^{(u)})_{u<Q})
\approx_{\epsilon_{\rm LPN}}
(W_p,(U^{(u)})_{u<Q}),
\qquad U^{(u)}\gets R.
$$

The replacements are mutually independent and independent of $W_p$.
For Reverse Cuckoo this is the paper's leakage-robust decisional assumption,
not an assumption that permits arbitrary additional transcript leakage.
It covers the entire stationary sequence; a one-expansion assumption is
not substituted for it.

## 4. Conditional composition lemma

Under these two premises and the correctness bound $\delta$, the real joint
distribution of the corrupt view and honest outputs is indistinguishable
from

$$
\mathsf{Complete}_p(\widetilde V_p,U)
$$

with advantage at most
$\delta+\epsilon_{\rm sim}+\epsilon_{\rm LPN}$.
If the target retains the *real* corrupt view instead, the bound is
$\delta+2\epsilon_{\rm sim}+\epsilon_{\rm LPN}$.
The premises hold for either choice of $p$.

**Proof.** Couple the real execution with its correctness-based completion
$\mathsf{Complete}_p(V_p,x'_h)$. They agree unless a correctness failure occurs,
which costs at most $\delta$. This step does not condition the support
distribution on success.

Apply distributed simulation to replace $V_p$ by $\widetilde V_p$, retaining
the honest samples. Applying the deterministic completion map preserves the
advantage bound $\epsilon_{\rm sim}$.

To justify the next transition, a Ring-LPN reduction receives $(W_p,a)$.
It runs $\mathcal S_p(W_p)$ with fresh simulation coins, computes
$\mathsf{Complete}_p(\mathcal S_p(W_p),a)$, and returns the distinguisher's bit.
The real challenge gives exactly the preceding hybrid. The uniform challenge
gives $\mathsf{Complete}_p(\widetilde V_p,U)$. This costs at most
$\epsilon_{\rm LPN}$. The reduction uses one joint challenge and neither
rewinds nor guesses an expansion.

For the real-view formulation, distributed simulation also implies marginal
indistinguishability of $V_p$ and $\widetilde V_p$. Append independent $U$ and
apply the completion map once more, at cost $\epsilon_{\rm sim}$.
The triangle inequality gives both bounds. All premise bounds are evaluated
at the distinguisher's budget plus the simulator and completion overhead. The
completion uses $Q$ NTTs and $O(QN_R)$ field operations. $\square$

## 5. Simulator in the paper's hybrid model

Use the interface convention in `prelim.tex`: ideal shared circuits and
permutations return fresh random shares. DPF calls instead use their
stateful chosen-share simulators. Product preprocessing is a private
shared-output computation, instantiated by addition/conversion and tensor
product protocols. Its realization security is included in the subprotocol
loss, not in the Ring-LPN assumption.

The simulator keeps a table of the corrupt party's shares and persistent
subprotocol states. It samples fresh ideal secret-output shares, but computes
public constants, local linear operations, and repeated uses from that table.
It never resamples a retained share merely because another call uses it.

**Product addresses.** Each party's contribution to a product address is its
own block offset. The simulator knows those contributions from $W_p$. It
simulates the arithmetic-to-binary conversion, then supplies the resulting
corrupt shares to the DMPF setup. The block pairing and row order are public;
no reconstructed product address is needed or opened.

**Reverse-Cuckoo setup.** Simulate the equality, activity, deduplication, and
reusable-shuffle calls in order, retaining their corrupt output shares. The
shuffle's hidden permutation stays inside its ideal functionality. Take the
public evaluator seeds from $\mathsf H$. Simulate feature evaluation and the
bin solver as shared circuits with the prescribed fixed schedule. Pivots,
rank, skipped rows, and real/dummy indicators are not opened.

For each descriptor matrix $R_s$ in $\mathsf H$, let $D_{p,s}$ denote the
corrupt solver-output share. Send $R_s\oplus D_{p,s}$ as the honest opening
share. In the fresh-sharing interface, $D_{p,s}$ is uniform even for a fixed
$R_s$. In a gate-level simulation, retain the internal local linear relations
instead of sampling the final share independently. The shared-circuit
simulator handles that internal transcript consistently with its interface.

This construction uses a descriptor already sampled by the leakage
experiment; it does not sample a descriptor independently of the supports.
The simulator does not need a compatible hidden support, permutation, or
solution witness. The public position map, sparse domains, representatives,
and empty-column decisions are deterministic functions of $\mathsf H$.
Feed the simulated corrupt address shares and public domains into the
reusable-DPF setup simulators, and retain their resulting local states.

**Payloads and expansion.** At each expansion, simulate tensor preprocessing
on the corrupt party's new coefficients. In its ideal sharing, each tensor
entry can have a fresh corrupt share; the honest share completes the product.
The public reorderings and negacyclic signs are applied to those same shares.
Feed them into deduplication and the retained shuffle, then into the column
DPFs. Neither the equality state nor the permutation is sampled again.

At each DPF update, use the appendix's online simulator: compute the corrupt
leaf conversions from its retained seeds, sample the opened correction with
its hybrid distribution, and complete the honest opening message. In odd
characteristic, the sign-multiplication call is simulated before that opening.
The corrupt output map is computed from this state and correction. It is not
replaced by a fresh uniform map. Sum the maps and apply the public ring
weights and NTT to obtain the corrupt OLE outputs.

**Joint guarantee.** Fix both parties' sampled noises and a descriptor from
their prescribed joint distribution. Fresh ideal shares have the same law
for every such hidden input. The simulator computes every local dependency
from previously coupled shares; the descriptor opening reconstructs exactly
the given descriptor. The DPF simulators preserve the joint corrupt view and
honest outputs. These facts inductively couple setup and all expansions,
while leaving the honest ring samples unchanged as auxiliary outputs.

Invalid column indices caused by a failed programmed equation are outside
the valid-input DPF guarantee. Charge the entire continuation of that setup
to its existing correctness-failure bound, without conditioning and
resampling the support distribution. Include this charge in
$\epsilon_{\rm sim}$ when applying Section 4; that section's bound may count
the same bad event twice and is conservative.

For Waterfall, its existing generic simulator replaces the Reverse-Cuckoo
setup and expansion steps. Product preprocessing and public ring assembly
are identical. Replacing the ideal subprotocols by their assumed secure
realizations adds their call-counted simulation losses. There are $P^2t$
product-list DMPFs; for Reverse Cuckoo, at most $2dP^2t$ column DPFs each
receive $Q$ expansions. This completes the abstract simulation bridge under
the paper's existing subprotocol guarantees. It does not reprove those
primitives or the conjectured concrete AES instantiation.

## 6. Claim boundary and next step

The lemma establishes honest-output completion security. It does not by itself
prove that the joint all-honest outputs have the standard random-OLE distribution
obtained by sampling independent uniform $x_0,x_1,z_0$ and setting
$z_1=x_0\odot x_1-z_0$. Such a claim also needs an output-distribution argument
for the actual construction; it does not follow solely from the chosen-share
DMPF interface.

Nothing in this composition argument requires fresh DPF seeds, extra
interaction, or a stronger DMPF functionality. The abstract simulator check
is complete at the stated interfaces. Implementation/model alignment is a
separate question, addressed next.

## 7. Source cross-check: what still needs a decision

Checked against paper commit `ee0f67b` and libOTe commit `c2fa1bc` in
`../libOTe-codex-rev-analysis`. This is a call/interface check, not a new
primitive-security proof or end-to-end benchmark.

- `RingLpnTriple.h::genDpf` builds the expected block-pair address lists and
  calls `arithmeticToBinary` before `setPoints`. Its insecure reveal branch
  is disabled by `if (0)`.
- `tensorRecv` and `tensorSend` use noisy VOLE for the coefficient tensor.
  The fresh tensor OTs are consumed and cleared. Expansion clears the
  coefficients and tensor shares; a later expansion obtains a new tensor.
- `expand` samples a fresh joint public mask seed on every invocation,
  retains the address-dependent DPF state, and implements the signs, shifts,
  wraparound, and NTT weighting in Section 2. Public seed expansion remains
  under the project's ideal-cipher convention, not an additional claim
  proved here.
- `RevCuckooDmpf::setPoints/expand` opens the hash seeds, descriptors, and
  DPF corrections discussed above. Diagnostic `mPrint` and `mDebug` paths
  reveal secrets and must remain disabled; they are outside the security claim.

The historical audit found three differences from the then-current paper:

| Detail | Paper model at audit | Implementation at `c2fa1bc` | Consequence |
| --- | --- | --- | --- |
| Hidden permutation | Uniform permutation | `RevCuckooDmpf` uses `WaksmanPermute::init`; its shared controls come from random OT choice bits | This is not the exact-uniform serial sampler used by Waterfall. The uniform-placement model does not automatically apply. |
| Evaluator randomness | Fresh evaluator randomness for each product list | One public root derives a seed per partition, reused across all sets in the batch | The leakage assumption must describe that joint reuse, or the implementation must use the stated per-list sampling. |
| Noise coefficients | Independent nonzero field elements | Both tensor samplers use `mCtx.fromBlock`, permitting zero | In the ideal uniform-field model this difference has a small explicit statistical bound; it need not force a protocol change. |

The permutation distinction is exact, not a measured security loss. A fixed
network driven by $b$ independent fair bits assigns permutation probabilities
in multiples of $2^{-b}$. For 32 positions it cannot be exactly uniform over
$32!$ permutations, since $32!$ has odd factors. This observation gives no
useful bound on its distance from uniform and is not an attack. Nevertheless,
the paper's uniform-placement leakage law and correctness calculation cannot
be transferred to it without an argument. Waterfall's `SerialWaksmanPermute`
is a different call path; its presence does not change Reverse Cuckoo's path.

For coefficients, coupling uniform field sampling to nonzero sampling fails
only if at least one sampled coefficient is zero. Across both parties and
$Q$ expansions, the distance is at most $2PtQ/q$. At $P=4$, $t=16$, and the
Goldilocks prime, this is approximately $Q\,2^{-57}$. Charge a sampling
replacement wherever it is used in a hybrid. This estimate does not include
primitive-realization losses and must not be described as independent of $Q$.

**Approved follow-up.** Reverse Cuckoo now uses `SerialWaksmanPermute`, samples
both parties' private permutations once during setup, and retains them for
payload expansion. Each (product list, partition) receives independent joint
seed contributions. Evaluator caches, dummy anchors, and public sparse-set
construction use that same flattened instance index. Public domain hashing
streams one evaluator at a time, preserving the 32-row inner-product kernel
without allocating all feature tables simultaneously.

These implementation changes address the first two rows without changing
the stated leakage experiment. The coefficient sampler is unchanged; the
paper now specifies uniform field coefficients, including zero.
A zero coefficient is not a correctness failure: it contributes a zero
payload through the tensor, DMPF, and ring arithmetic. No hardness ordering
between uniform and nonzero coefficients has been established. The bound
above is a sampling-distance bound, not an identified attack or a bound on
all OLE distinguishing advantage. It is needed only when comparing to a
nonzero-coefficient model, not for matching the implementation to the assumption.

Validation of the approved changes: the 24-suite focused WSL/GCC runner
passes with its defaults and with domain 4097, 16 points, three sets, four
expansions, and either two or three partitions. The new sparse-set regression
compares per-instance public hashes against a scalar reference, including
empty and maximally loaded buckets and a partial 32-row batch. Repeated
expansion checks that the sampled private permutations remain unchanged.
Seven hash-conditioning checks and 18 sharing-contract checks also pass.
The four-suite Ring-LPN runner also passes with `-trials 2`, covering support
filtering, its audit checks, stationary reuse, and OLE correctness. These
integration tests use small rings; they are not a full-size performance run.
These are correctness/regression checks, not empirical proofs of uniformity
or new performance measurements.

### Local sources

- [DMPF interface and simulation convention](DMPF.tex).
- [Ring-LPN algebra, product lists, and stationary sampling](ringLpnOverview.tex).
- [Generic Waterfall simulation](WaterfallCuckoo.tex), theorem
  `thm:waterfall-simulation`.
- [Reusable DPF simulation](ReusableDpfAppendix.tex), theorem
  `thm:reusable-dpf-simulation`.
- [Reverse-Cuckoo protocol and leakage assumption](RevCuckooSpecialized.tex),
  assumption `ass:rc-decisional-lpn`.
- [Earlier output-sharing investigation](REUSABLE_DPF_SHARING_CONTRACT.md),
  Section 4. Its historical status notes are not current release status.
