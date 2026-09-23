# Reusable DPF: output-sharing boundary

2026-09-05. This audit does not change the paper, functionality, protocol,
parameters, or implementation. It records why the current concrete
reusable-DPF realization claim cannot yet be completed in the online model.
The issue is the ideal output-sharing law, not seed caching by itself.

## Interfaces being compared

`DPF.tex`, Figure `fig:DpfFunc`, samples an independent uniform output map
for one party at every expansion. The other party also has a uniform map,
since its map is a fixed point map minus the first party's uniform map.
This sampling uses fresh randomness after both payload shares are submitted.

Figures `fig:Dpf` and `fig:sparseDpf` instead let each party locally evaluate
its saved seed shares. The leaf-update box opens one group element and
applies it locally using the saved tags. Private ideal circuit state is
not exposed or used in this argument.

## A check using only the corrupt party's local values

Fix party p, a public ordered set S with at least three elements, and an
expansion counter q. For each x in S, define

\[
t_{p,x}:=\operatorname{lsb}(s_{p,x}),\qquad
u_{p,x}:=(1-2p)L_q(\operatorname{msbs}(s_{p,x})).
\]

These are the party's local tag and provisional group value. They can be
computed after setup, before a payload is supplied. The concrete expansion
returns

\[
y_{p,x}=u_{p,x}+(1-2p)\gamma t_{p,x},
\]

where gamma is the opened leaf correction. The internal sign correction
for groups of characteristic other than two only changes gamma.

Among three tags, two agree. Choose the lexicographically first pair
(a,b) of distinct addresses with equal tags and define

\[
\delta:=u_{p,a}-u_{p,b}.
\]

For every possible gamma, the concrete outputs satisfy

\[
y_{p,a}-y_{p,b}=\delta.
\]

This follows by subtracting the two update equations. The correction terms
cancel because the tags agree. No assumption about the PRG, index,
payload, or other party's seeds is needed.

## Chronology and the resulting test

The conclusion here concerns an online realization of the reusable
functionality, consistent with adaptive calls and composition:

1. Complete setup.
2. The semi-honest corrupt party computes and records (a,b,delta).
   An external observer receives this record before supplying the next payload.
   This is local observation; it does not change any protocol message.
3. Supply the expansion inputs and obtain the corrupt party's output map.
4. Check whether its values at a and b have difference delta.

The real execution always passes. In the ideal execution, (a,b,delta)
is fixed before the functionality samples its next map Y_p. Its coordinates
at a and b are independent uniform group elements, so

\[
\Pr[Y_p(a)-Y_p(b)=\delta]=1/|\mathbb G|.
\]

The distinguishing gap is therefore 1-1/|G|. The test is efficient and
applies to every nontrivial finite additive group. It already uses only
one expansion after setup. A programmable leaf oracle does not bypass
the test if the local evaluations were queried and recorded before the
ideal output was sampled.

Chronology is essential. This calculation is not, by itself, an impossibility
proof against every offline transcript sampler given the entire future
output sequence before simulating setup. Such a sampler has a different
interface; its sufficiency for the paper's adaptive composition would need
to be justified separately. This audit does not silently impose an online
restriction on an otherwise standalone theorem.

## Scope and proposed next decision

**2026-09-07 follow-up:** `REUSABLE_DPF_SHARING_CONTRACT.md` specifies a
candidate online experiment: simulate the corrupt share from its inputs and
view, then complete the honest share to the correct map. It records an
interactive-FSS precedent, conditional leaf-update and Waterfall arguments,
and the separate OLE distribution obligation. This removes the particular
prescribed-output obstruction, but the concrete setup/hidden-active-pad
lemma is not yet proved. No TeX or implementation change has been made.

- There is no attack here on index privacy, payload privacy, or Ring-LPN.
- The calculation does not inspect hidden state of an ideal DPF invoked by
  Waterfall. Its existing ideal-subprotocol theorem is not refuted.
- A PRF value can be pseudorandom to someone without its key and still be
  exactly computable by the party holding that key. Ordinary PRF security
  does not justify replacing those local computations by independent values.
- Merely retaining hidden state inside an ideal functionality is allowed.
  The relevant clause is its prescribed fresh uniform output map, independent
  of the already recorded local prediction.
- Refreshing the public Ring-LPN masks does not address or cause this DPF
  interface issue.

Before changing the paper, discuss the intended output-sharing convention
with Peter. The preferred investigation is an interface/proof correction
that models DPF-generated shares and preserves the Setup/Expand API,
cached seeds, and efficient implementation. That correction must retain an
explicit privacy guarantee and support the actual caller's composition.
It is not sufficient to delete the uniformity clause, let a simulator choose
outputs without restriction, or assume the existing Waterfall output-fiber
proof transfers unchanged.

The follow-up contract restricts the simulator's information and requires
joint online simulation with honest outputs. Those restrictions, not merely
allowing output choice, supply its proposed privacy guarantee.

Adding a final fresh resharing is a different protocol and cost choice.
No such change is proposed for implementation without discussion.

## Executable checks

`experiments/test_reusable_dpf_output_boundary.py` checks:

- cancellation of the scalar correction for all three-bit tag patterns;
- the exact 1/|G| acceptance fraction for independently uniform maps;
- the identity on states generated by the existing sparse-tree model.

Cyclic groups of orders 2, 3, 4, 5, 8, and 17 are covered. These tests
establish finite algebraic identities, not a cryptographic security theorem.
The boundary checks, existing sparse scheduling tests, and Waterfall hybrid
tests pass together: 37 tests.

Run from the paper repository:

```powershell
& './tmp/handoff-review-venv/Scripts/python.exe' -m pytest -q experiments/test_reusable_dpf_output_boundary.py experiments/test_sparse_dpf_schedule.py experiments/test_waterfall_hybrid_simulation.py
```

CWC-04, CWC-08, and CWC-11 require keeping the interface, observer,
chronology, and distributional claim explicit before completing the proof.
