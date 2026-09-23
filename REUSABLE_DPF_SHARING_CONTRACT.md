# Reusable DPF: a candidate output-sharing contract

**Date:** 2026-09-07. **Status update:** Peter approved the output-sharing
convention. The DPF/DMPF functionality boxes and adjacent simulation
definitions have now been updated for his review. The corrupt party chooses
its own share; the functionality returns the honest complement. The
investigation below predates that edit. Its conditional proofs and open
obligations remain unchanged; no protocol or C++ implementation was changed.

**Subsequent approved proof pass:** Waterfall's theorem and Appendix A.8 now
use the chosen-share aggregation proof. New Appendix C proves the exact
correction distributions under independent honest active children and the
online update simulator under fresh honest active values. Theorem C.3 now
adds reductions under tag-conditioned tree-PRG and joint leaf-PRF security.
The concrete evaluator arguments and OLE output distribution remain to
prove. `REUSABLE_DPF_GENERATOR_AUDIT.md` records an additional implementation
discrepancy in sparse-level masking; the chosen-share contract does not
remove that setup-privacy issue. The targeted tests now total 70; the
historical investigation below is retained as context, not current status.

## Decision in brief

There is a meaningful alternative to independently uniform ideal output
shares: simulate the corrupt party's view and its locally computed share,
then complete the honest share to the specified function. This removes the
particular output-prescription obstruction in
`REUSABLE_DPF_SIMULATION_AUDIT.md`. It does not establish the concrete
simulator or the application's pseudorandom-correlation guarantee by itself.

Cached seeds and polynomially many adaptive payloads can remain. No fresh
coin-flipping exchange is part of this proposal. Whether the existing
construction realizes the contract remains a proof obligation.

## 1. The contract, with its observer and chronology

Let G be the finite additive output group. Let D be the public output
domain: a public sparse set for sparse DPF, or [0,N) for DMPF. Setup takes
the parties' index shares and stores the reconstructed indices inside the
functionality. Expand takes their current payload shares. Write
f_q in G^D for the specified point or multi-point map on expansion q.
The same dummy and repeated-index conventions as the paper apply.

Fix one static semi-honest corrupt party p in {0,1}. The other party is
h=1-p. Define an online ideal experiment as follows.

1. The functionality receives both parties' input shares. The simulator
   receives only the corrupt party's input shares, public parameters, and
   explicitly permitted public information. For generic Waterfall there is
   no extra index-dependent leakage. Reverse Cuckoo needs its separate
   descriptor-leakage interface.
2. At setup the simulator generates the corrupt local state and honest-party
   messages visible to that party. It keeps its simulation state.
3. For expansion q, the simulator processes the current corrupt input and
   simulates the messages and local computation online. It produces the
   corrupt output map Y_p in G^D as part of that computation.
4. The functionality delivers Y_h = f_q - Y_p to the honest party. The
   simulated corrupt output is Y_p, not a newly sampled prescribed map.

Require one stateful PPT simulator, for either choice of p, that works for
all valid honest inputs and polynomially many adaptively chosen calls. It
does not rewind the caller or obtain future inputs. Compare the joint
external execution, including the corrupt view and honest output maps,
not just the marginal transcript. The caller can use prior outputs in
choosing later inputs; any information it explicitly communicates is then
available in both experiments.

The simulator is not given the hidden reconstructed indices, hidden
payloads, f_q, or current honest output to choose Y_p. Completing the
honest output is the functionality's job. Sending Y_p to the functionality
is an ideal-world operation, not a protocol message or a new application API.

This contract remains substantive: the corrupt view, including its output,
must be generatable without the honest private inputs. A transcript that
reveals an honest index cannot be simulated under that rule. Correctness
separately fixes the sum of the shares, except for the stated failure bound.

It deliberately does not promise that the individual shares are fresh
uniform maps independent of prior local state. A full UC functionality
would additionally need its corruption, scheduling, and all-honest branches
specified. The experiment above is a candidate static, sequential,
stateful contract, not a claim of arbitrary UC composition.

### Existing precedent

SHARK's interactive-FSS definition lets its stateful simulator generate the
corrupt output and completes the honest output to the correct function
value; the comparison includes the honest output and corrupt view, without
rewinding. Its interface is not our distributed shared-input, reusable API.
We can adopt the output convention, but must prove our extension.
See [SHARK, Section 2.4, Definitions 3–4](https://eprint.iacr.org/2025/716.pdf).

## 2. Why the leaf-update obstruction disappears

Keep the actual local computation. For party p, address x in D, and current
round q, let u_p(x) be its signed leaf conversion and t_p(x) its saved tag.
It returns

    Y_p(x) = u_p(x) + (1-2p) gamma t_p(x).

The simulator keeps its generated seeds and the same public conversion
algorithm, then computes this output. Equal-tag differences therefore
remain exactly as predicted from the local state. They no longer have to
agree with an independently prescribed ideal map.

Here is an algebraic route to simulating gamma, conditional on a cryptographic
setup/leaf-conversion lemma. Suppose the active address is alpha in D.
Off alpha, u_0(x)+u_1(x)=0 and t_0(x)=t_1(x). At alpha the tags differ.
Set sigma=t_0(alpha)-t_1(alpha), which is either +1 or -1. Correctness uses

    gamma = sigma [beta-u_0(alpha)-u_1(alpha)].

If the honest active value u_h(alpha) is independent uniform in G,
conditioned on the corrupt view before the current opening and on the
current payload, then gamma is uniform:
negation, multiplication by sigma, and translation are bijections of G.
The simulator can sample gamma without knowing alpha, beta, or sigma,
simulate the ideal shared-arithmetic messages consistently with that opening,
and compute Y_p. Correctness determines the honest output.

For reuse, the required hybrid replaces the sequence of honest active
conversions by fresh hidden uniform elements across distinct round inputs.
It must remain valid with adaptive payloads and all previous messages.
We never replace evaluations under seeds known to the corrupt party.

**Unproved cryptographic obligations:** input-private setup; the joint hidden
active-seed/conversion property with all correction words and off-path keys
as auxiliary information; cross-instance/session domain separation; seed
collisions and oracle guesses; all valid dummy/empty-domain cases; and the
exact shared-arithmetic transcript.
Ordinary one-shot PRG security alone does not discharge these obligations.
The implementation uses AES hashing with public round-dependent keys and
hidden leaf inputs, so an ideal-cipher proof must also account for inverse
queries. Do not describe it without justification as a secret-key PRF.

The small executable checks below establish the conditional algebra, not
these cryptographic claims.

## 3. Waterfall composition

Use subroutine-compatible online simulators with the joint-output guarantee
above. Preserve the existing public-hash, uniform-injection, and ideal
shared-circuit setup simulation. Its placement-failure charge is still
once per setup, not once per expansion.

During an expansion, simulate the row-scatter instances to obtain the
corrupt column payload shares. Feed those shares into the column-DPF
simulators. If S_j is the public sparse domain of column j and Y_{p,j} its
simulated corrupt output, return the local sum

    Y_p(x) = sum_{j : x in S_j} Y_{p,j}(x),  x in [0,N).

The honest sum is f_q-Y_p by subprotocol correctness and the placement
invariant. No assumption that these maps are uniform or independent is
needed for this algebra. The nonlinear shared-circuit calls must handle
arbitrary valid additive input sharings, not only fresh uniform ones.

This is a route to a replacement proof, not an unchanged use of the current
appendix: its prescribed-output pivot/fiber sampler would be removed.
The theorem would refer to the new contract. The existing theorem for
fresh-uniform ideal DPF calls is not refuted by this proposal.
Subroutine security losses require a call-counted hybrid argument; the
output-sum identity alone is not a composition theorem.

The API, caching, parameters, and implementation costs need not change.
No LPN assumption enters generic Waterfall's simulation.

## 4. Ring-LPN / OLE: the application check

The PCG definition distinguishes the expanded pair's target distribution
from security given one party's actual seed. The latter uses reverse
sampling: keep that party's actual output and sample a compatible honest
output. See [Ring-LPN PCGs, Definitions 2.5–2.6](https://eprint.iacr.org/2022/1035.pdf).
Our batches are interactive; this comparison does not turn them into silent
PCG expansion.

For packed OLE over a field F, write a party's output as (x_p,z_p) in
F^N times F^N. Correctness is z_0+z_1=x_0 componentwise-times x_1.
The ideal honest completion of a fixed corrupt output is

    U uniform in F^N,
    (x_h,z_h) = (U, x_p componentwise-times U - z_p).

This matches the proposed sharing convention. More precisely, suppose the
preprocessing and DMPF simulations generate the entire corrupt view V
(including x_p,z_p) from the corrupt noises, public masks, allowed descriptor
leakage, and simulation randomness, without the honest noises. The
replacement must preserve the honest ring samples as auxiliary outputs.
Then the paper's joint leakage-robust Ring-LPN challenge can be passed through
this simulator. Replacing the honest samples by independent U gives exactly
the completion above, round by round. This is efficient postprocessing of
the joint challenge, not a new assumption about arbitrary transcript leakage.

This conditional argument addresses the reverse-sampling/security part.
It still requires the actual simulators and the stated joint Ring-LPN
assumption, including fixed supports, refreshed coefficients/masks, and
descriptor leakage across all rounds.

**Separate requirement:** the all-honest expanded output pair must have the
desired pseudorandom OLE distribution. Sum correctness and private sharing
do not alone assert marginal uniformity of the output masks. For example,
the constant-zero map can be privately shared as two zero maps under the
weak contract; those maps are not independent random sharings of zero.
This is a logical counterexample to an implication, not a proposed DPF
construction or an attack on the current OLE protocol.

Accordingly, the OLE proof needs a separate output-distribution hybrid for
the actual leaf masks and ring arithmetic, or another demonstrated source
of the requisite masking. Do not claim that weakening the DMPF interface
automatically proves the full correlation-generation theorem.

## 5. Recommendation and validation

Prefer this output-completion contract over adding fresh public coins solely
to repair the old output prescription. First establish the concrete
setup/active-leaf simulation lemma and check the OLE marginal-distribution
hybrid. Then, with Peter's approval, update DPF/DMPF together and replace
Waterfall's output-fiber proof with the compositional simulator above.

`experiments/test_reusable_dpf_sharing_contract.py` checks the uniform-pad
lemma, both corruption choices, adaptive two-round algebra, additive
aggregation, and the OLE completion identity. The tests are exact small
finite enumerations and identities, not an AES/PRF simulation proof.

Validation: 18 new checks pass; together with the output-boundary,
sparse-schedule, and Waterfall-hybrid checks, **55 tests pass**. Run:

```powershell
& './tmp/handoff-review-venv/Scripts/python.exe' -m pytest -q experiments/test_reusable_dpf_sharing_contract.py experiments/test_reusable_dpf_output_boundary.py experiments/test_sparse_dpf_schedule.py experiments/test_waterfall_hybrid_simulation.py
```

No performance measurements are required for this contract investigation.
No TeX or implementation files were changed. CWC's interface, observer,
chronology, and claim-scope checks motivate the explicit qualifications here.
