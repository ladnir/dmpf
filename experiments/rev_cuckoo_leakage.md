# Reverse-cuckoo leakage and cuckoo failure

## Model and exact leakage channel

Let `n` labeled active items have one candidate in each of `w` disjoint
partitions of size `d`, and let `H` denote the resulting partitioned cuckoo
graph.  Write `Z_A(H)` for the number of injective placements of support `A`
in `H`.  Under an input-independent uniformly random hash family `U`,

\[
  \mathbb E_U[Z_A(H)] = \frac{(wd)_n}{d^n}.
\]

Reverse cuckoo first samples a uniform injection of the active items into the
`wd` bins and then samples the other `w-1` choices uniformly.  Consequently,
its idealized public-table distribution is not `U` conditioned on
`Z_A(H)>0`.  It is the size-biased distribution

\[
  \Pr[H\mid A]_{\rm RC}
    = \Pr[H]_U\frac{Z_A(H)}{\mathbb E_U Z_A(H)}.
\]

Thus the exact information density relative to the input-independent base is
`log2 Z_A(H) - log2 E[Z_A(H)]`.  For `w=2`, every connected cuckoo component
is a bipartite multigraph.  A placement exists iff every component has
`|E|<=|V|`; a tree component has `|V|` placements and a unicyclic component
has two.  This gives exact linear-time failure and matching-count samples.
For `w=3`, the experiment uses exact augmenting-path matching for failure.

The KL figures below are divergences from the input-independent base `U`, not
claims of equal mutual information about a particular Ring-LPN support prior.
They upper-bound that mutual information because

\[
 \mathbb E_A D(P_{H\mid A}\|U)
 = I(A;H)+D(P_H\|U).
\]

The actual Goldreich feature family and its shared batching structure are
covered by the leakage-robust Ring-LPN assumption in the paper; the experiment
characterizes the standard independent-choice cuckoo channel.

## Current parameters

The implementation default is `w=2,d=n`, i.e. total load `n/(wd)=1/2`, the
two-choice threshold.  Four million exact trials per row give:

| `n` | exact failure | 95% CI | KL from `U` | TV from `U` |
|---:|---:|---:|---:|---:|
| 16 | 0.062030 | [0.061799, 0.062261] | 1.049 bits | 0.476 |
| 64 | 0.100536 | [0.100248, 0.100824] | 4.613 bits | 0.831 |
| 128 | 0.115777 | [0.115471, 0.116083] | 9.370 bits | 0.951 |

The failure probability increases with `n`; it does not decay toward a
security parameter.  The KL/TV values use 1,048,576 direct samples from the
size-biased real distribution.  Increasing `d` lowers size-bias leakage only
slowly: for `n=16,w=2`, measured KL at `d=32,64,128,256` is respectively
0.430, 0.200, 0.0958, and 0.0464 bits.  Therefore a low cuckoo failure
probability by itself does not recover ordinary leakage-free DMPF security.
The current `cuckooSecParam` argument is not used in the partition-size
calculation, so the label “cuckoo sec. 2” does not select or certify a failure
exponent.

If 256 sets are generated, a safe failure target for 40-bit batch security is
`2^-48` per set.  Reusing public feature matrices across the batch does not
invalidate this failure union bound, but it means that multiplying the
single-set leakage figures is not an empirical characterization of the joint
leakage channel.

## Failure extrapolation and 40-bit choices

For each `k`, failure implies a Hall witness: some `k` items collectively use
at most `k-1` bins.  The script `rev_cuckoo_hall.py` sums all such witnesses in
log space.  This is a conservative union bound, not just a fitted line.  In
the observable three-choice tail for `n=16`, measured failure at
`d=9,10,11,12` was `6.35e-5, 1.04e-5, 2.03e-6, 1.67e-6`; the corresponding
Hall bounds were `2^-11.26, 2^-14.36, 2^-16.91, 2^-18.97`.  The bound converges
to the leading four-item witness `binom(n,4)/d^9` as slack grows.

Power-of-two partition sizes sufficient under the full Hall union bound are:

| choices | `n` | `d` for one 40-bit set | `d` for a 256-set 40-bit batch |
|---:|---:|---:|---:|
| 2 | 16 | 8,192 | 32,768 |
| 2 | 64 | 32,768 | 65,536 |
| 2 | 128 | 32,768 | 131,072 |
| 3 | 16 | 64 | 128 |
| 3 | 64 | 128 | 256 |
| 3 | 128 | 256 | 256 |
| 4 | 16 | 16 | 16 |
| 4 | 64 | 32 | 32 |
| 4 | 128 | 64 | 64 |

For the 256-set three-choice rows, the log2 batch failure bounds are
respectively `-44.16`, `-44.71`, and `-40.63`.  With four choices the much
smaller rows above give `-43.83`, `-49.08`, and `-60.01`.  The two-choice
sizes are impractical, but three choices are not globally optimal once the
number of choices is allowed to change.

The simple four-choice schedule `d=2n` suggested during review is also valid
with ample margin: its log2 batch bounds are `-59.90`, `-81.14`, and `-92.02`
for `n=16,64,128`.  Relative to the three-choice rows it reduces total bins
and the `w(d+40)` hash-width proxy at `n=16,64`, but is overprovisioned at
`n=128`.  The Hall-minimal four-choice schedule `d=16,32,64` is smaller still.

This does not settle performance from parameters alone.  Each additional
choice adds another full domain-sized partition pass, and the expanded leaf
work is proportional to `wN`; smaller `d` reduces the hash width, solver size,
and number of bins.  Therefore `w=3` versus `w=4` should be selected with an
end-to-end setup/expansion benchmark, especially for Stationary LPN where the
expansion cost is repeated.  These bounds address placement failure only. The
size-biased support leakage remains and must be covered by the explicit
leakage-robust Ring-LPN/Stationary-LPN assumption.

The production implementation currently supports only `w=2,3` and derives
`d` directly from `n,w`; it has no override for these recommended values.  The
exact experiment harness supports arbitrary `w>=2`, but adopting a four-choice
schedule requires a separate implementation change to generalize initialization
and plumb an explicit partition size or validated failure target, followed by
performance measurement.

## Reproduction

Build `frontend_libOTe`, then run, for example:

```text
frontend_libOTe -invMtx -cuckooLeak -n 16 -w 2 -d 16 -trials 4194304 -seed 10
frontend_libOTe -invMtx -cuckooLeak -real -n 16 -w 2 -d 16 -trials 1048576 -seed 2
python experiments/rev_cuckoo_hall.py --batch 256 --target 40
```

Raw rows, seeds, Wilson intervals, information-density quantiles, and the
three-choice threshold/tail sweeps are in `rev_cuckoo_results.csv`.
