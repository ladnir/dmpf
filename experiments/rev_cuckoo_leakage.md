# Reverse-cuckoo leakage and cuckoo failure

The attack-focused continuation, including the exact finite-`K` channel,
conditional linear-test bias, SSD comparison, and nonlinear attack plan, is in
[`rev_cuckoo_lpn_attack_study.md`](rev_cuckoo_lpn_attack_study.md).

This work concerns a research prototype that is not known to be deployed or
used in production.  The accompanying manuscript is a recent unpublished
ePrint preprint, and the parameter rows below are experimental comparison
points rather than deployment recommendations.

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

## Empirical failure extrapolation with `w` as the primary cost

The protocol performs a full partition/hash pass for each choice, so minimizing
the total expansion `e=wd/n` is not the right primary objective.  We instead
compare only `w=2` and `w=3`, and treat the partition size `d` as an unrestricted
integer (and as a continuous variable when fitting and inverting the curve).

We ran 469,762,048 exact-matching trials over 82 `(w,n,d,seed)` rows.  In the
observable tail, rows with at least 20 failures and `lambda>=8` fit

```
lambda = -log2(p_fail) = a log2(d) + b log2(n) + c.
```

| `w` | points | observed `d` slope `a` | RMS residual | maximum residual |
|---:|---:|---:|---:|---:|
| 2 | 18 | 4.759 | 0.381 bits | 0.959 bits |
| 3 | 21 | 16.065 | 0.915 bits | 1.889 bits |

The observed three-choice slope is attractive but cannot be extrapolated
unchanged to 40 bits.  Deep in the tail, failure is dominated by the smallest
local obstruction: `w+1` items receive exactly the same `w` choices.  Its
probability is

```
p_fail ~ binom(n, w + 1) / d^(w^2).
```

Thus the ultimate `d` slopes are only 4 for `w=2` and 9 for `w=3`.  The
experiments show the crossover.  At `n=16,w=2`, the ratio of measured failure
to the leading obstruction falls from 2.11 at `d=40` to 1.23 at `d=96`
(later points have wider counting error).  For `w=3`, the ratio falls from
40.4 at `d=8`, to 6.16 at `d=10`, 3.63 at `d=11`, and 2.37 at `d=12`.
Consequently, extending the steeper observable slopes would materially
underestimate the required `d`; the continuous extrapolation below uses the
asymptotic local-obstruction slope.

### Sanity check against ordinary three-choice cuckoo

The familiar `e` around 1.27 is a large-enough-set ordinary-cuckoo parameter,
not a set-size-independent constant.  The published tables begin at `n=256`.
Ordinary cryptoTools cuckoo also constructs three distinct candidate bins for
each item.  For four fixed items its smallest obstruction therefore has
probability `1/binom(m,3)^3`.  In the partitioned design it has probability
`1/d^9 = 3^9/m^9`, about 91 times larger than ordinary cuckoo's asymptotic
`6^3/m^9`.  This is a 6.51-bit tail penalty, or a factor of about 1.65 in total
table size when this obstruction dominates.

An ordinary-cuckoo control using exact matching and distinct choices confirms
the finite-size effect.  At `n=128`, ordinary cuckoo at `m=156` (`e=1.219`)
failed 71 of 262,144 trials, while partitioned cuckoo failed 98; at `m` about
165 (`e` about 1.29), both were around a few failures per million, only about
19 bits.  At `n=256,m=300` (`e=1.172`), ordinary and partitioned failure were
respectively `40/65536` and `37/65536`.  Thus partitioning has little threshold
penalty; the important difference is the deep local-obstruction tail.

Inverting just this necessary local obstruction gives the following common
scale.  The ordinary `n=256` single-set row recovers `e=1.277`, the standard
`1.27` sanity anchor.  Both ordinary and partitioned expansion improve as `n`
grows; displaying only absolute `d` in the earlier table obscured this.

| `n` | ordinary `e`, one 40-bit set | partitioned `e`, one 40-bit set | ordinary `e`, 256-set batch | partitioned `e`, 256-set batch |
|---:|---:|---:|---:|---:|
| 16 | 5.813 | 9.400 | 10.625 | 17.407 |
| 64 | 2.750 | 4.504 | 5.078 | 8.341 |
| 128 | 1.875 | 3.081 | 3.469 | 5.706 |
| 256 | 1.277 | 2.102 | 2.363 | 3.892 |

| `w` | `n` | continuous `d` for one 40-bit set | continuous `d` for a 256-set 40-bit batch |
|---:|---:|---:|---:|
| 2 | 16 | 4,981.35 | 19,925.40 |
| 2 | 64 | 14,629.86 | 58,519.44 |
| 2 | 128 | 24,751.88 | 99,007.52 |
| 3 | 16 | 50.14 | 92.84 |
| 3 | 64 | 96.09 | 177.94 |
| 3 | 128 | 131.46 | 243.44 |

The table above isolates the leading obstruction and is a sanity lower scale,
not the final selector.  For `w=3` we fit the full failure curve as the sum of
two mechanisms:

```
p_model(n,d) = 2^-(a_n * (3d/n) + b_n) + L4(n,d) + L5(n,d) + L6(n,d).
```

The first term is the threshold/global failure component.  `L4` is the exact
four-item term `binom(n,4)/d^9`; `L5` counts the `(2,1,1)` five-item
obstructions; and `L6` counts the `(3,1,1)` and `(2,2,1)` six-item
obstructions.  Threshold-dominated observations give:

| `n` | threshold slope `a_n` | intercept `b_n` | fit RMS |
|---:|---:|---:|---:|
| 16 | 18.699 | -17.160 | 0.098 bits |
| 32 | 32.563 | -33.653 | 0.075 bits |
| 64 | 53.228 | -56.932 | 0.117 bits |
| 128 | 84.701 | -91.945 | 0.144 bits |

The `n=16,d=10..14` crossover points were held out.  The model error on them
is respectively `-0.208,-0.249,-0.247,+0.002,+0.430` bits.  This validates
both the handoff between terms and the conservative inclusion of six-item
events.  Solving the summed model gives:

| target | `n` | integer `d` | expansion `3d/n` | modeled exponent |
|---|---:|---:|---:|---:|
| one 40-bit set | 16 | 51 | 9.563 | 40.162 |
|  | 64 | 97 | 4.547 | 40.037 |
|  | 128 | 133 | 3.117 | 40.056 |
| 256-set 40-bit batch | 16 | 93 | 17.438 | 48.005 |
|  | 64 | 179 | 8.391 | 48.053 |
|  | 128 | 244 | 5.719 | 48.003 |
| batch plus 2 fitted bits | 16 | 109 | 20.438 | 50.071 |
|  | 64 | 208 | 9.750 | 50.009 |
|  | 128 | 285 | 6.680 | 50.027 |

At every selected row the fitted threshold term is beyond 160 bits; local
obstructions completely determine the result.  Because the held-out error is
below half a bit and unmodeled larger obstructions decay by additional powers
of `d`, the two-bit-margin rows are candidate benchmark parameters.
They cost only about 17% more `d` than the model-minimal batch rows because
the dominant probability scales as `d^-9`.

These estimates address the probability of a non-matchable cuckoo graph.
They do not remove the inherent support leakage of a successful public hash
descriptor; that leakage is covered by the explicit leakage-robust LPN
assumption.  The experiment already accepts arbitrary integer `d`.  The prototype
currently rounds `d` to a power of two because its hash output is interpreted
directly as a bit-string index; supporting continuous integer `d` requires an
unbiased range-reduction step in that circuit and corresponding sparse-set
changes.

## Full Ring-LPN leakage game

The next experiment must preserve the support correlations in the
implementation.  Write `P=4` for `mNumPolys`, `t=16` for `mPolyWeight`, and
`L=N/t` for the regular-polynomial block length.  Party 0 has offsets
`x[a,u]` and party 1 has offsets `y[b,v]` in `[0,L)`, for polynomial indices
`a,b in [P]` and block indices `u,v in [t]`.  For each polynomial pair and
output block, the DMPF support list is

```
A[a,b,k] = { x[a,u] + y[b,(k-u) mod t] : u in [t] }.
```

The sums lie in the DPF tree domain `[0,2L)`.  There are `P^2*t=256` lists,
but only `2*P*t=128` underlying position variables, and an insider already
knows one party's `P*t=64` positions.  Consequently, multiplying a single-list
divergence by 256 does not give the mutual information of the full game.
Likelihood scores do add under the independent-descriptor idealization below,
but they constrain the same 64 unknown offsets rather than 256 fresh secrets.

For the base two-noise notation, let `X=(e,f)`, `R=r`, `Y=re+f`, and let
`Lambda` be the complete joint descriptor transcript.  The average extra
information for decoding and the incremental KL evidence in the decisional
game are respectively

```
k_dec  = I(X; Lambda | R,Y)
k_dist = I(Y; Lambda | R).
```

They satisfy `I(X;Lambda)=k_dist+k_dec`.  The implementation uses four regular
sparse polynomials and public dense multipliers, so its corresponding syndrome
is `Y=E[0]+sum_j R[j]*E[j]`; the first reduced experiment should implement both
the two-noise game and this implementation-aligned module version.

For an ideal independently represented descriptor batch, a candidate support
has likelihood proportional to

```
W(X; Lambda) = product_i Z_(A_i(X))(h_i).
```

The exact small-parameter experiments enumerate syndrome-compatible candidates
and verify the information decomposition above.  They also expose a limitation
of Shannon leakage in the full game: at larger dimensions the low-weight error
is normally information-theoretically unique, even though finding it remains
computationally hard.  Thus `k_dec` can approach zero while the descriptor still
helps prioritize support candidates.  Our operational exhaustive-search metric
is

```
g_search = -log2(E[rank_Lambda(X) / rank_random(X)]),
```

where candidates are ordered by `W(X;Lambda)` and then checked against the
syndrome.  This is an ideal likelihood-ordering benchmark, not yet an efficient
ISD attack.

`rev_cuckoo_full_lpn.py` implements the first exact structural model over a
small binary cyclic ring.  It conditions on one party's known regular supports,
enumerates the other party's supports, retains all correlated product lists,
and checks the information-density decomposition on every trial.  A sanity run
with `P=2,t=2,L=4,d=2`, 256 candidates, eight descriptor lists, and 1,024 trials
gave total information `0.198 +/- 0.026` bits, split into decisional and decoding
terms `0.079 +/- 0.018` and `0.119 +/- 0.018` bits.  Its expected support-search
gain was `0.492` bits and the unsmoothed posterior min-entropy loss was
`0.510 +/- 0.011` bits.  The prime-field enumerator uses nonzero coefficients
and true negacyclic multiplication; it confirms that coefficients affect the
syndrome ambiguity but not the support-prioritization channel.

### Current-scale ideal rank model

For one support list, let `U` be the independent hash distribution, let `P` be
the planted Reverse-Cuckoo distribution, and write

```
r(h) = P(h)/U(h) = Z_A(h)/E_U[Z_A(H)]
ell(h) = log2 r(h).
```

For a random false support disjoint from the planted support, its descriptor
score has the `U` distribution while the true score has the `P` distribution.
With `B` independent descriptors, the false-before-true rank fraction is

```
q_B = Pr[sum_i ell(U_i) >= sum_i ell(P_i)],
```

with half weight on ties.  Its large-deviation exponent per descriptor is

```
C_pair = -2 log2 E_U[sqrt(r)] = -2 log2 E_P[1/sqrt(r)].
```

The simulator estimates this moment from planted samples.  It also estimates
the variance under the midpoint distribution proportional to `sqrt(PU)` and
uses the corresponding normal saddlepoint prefactor to report
`g_search ~= -log2(2 q_B)`.  For the fully enumerable `n=2,w=2,d=2,B=8` case,
the exact expected-work gain is `0.6650` bits and the estimate is `0.6635`.

For `t=16,B=256`, the ideal independent-descriptor results are:

| `w` | `d` | KL/list | `C_pair`/list | modeled `g_search` |
|---:|---:|---:|---:|---:|
| 2 | 16 | 1.0507 | 0.6342 | 166.46 bits |
| 2 | 64 | 0.1995 | 0.1066 | 30.28 bits |
| 2 | 256 | 0.04644 | 0.02411 | 8.20 bits |
| 2 | 8,192 | 0.001397 | 0.000705 | 0.715 bits |
| 3 | 51 | 0.06747 | 0.03455 | 11.08 bits |
| 3 | 93 | 0.03658 | 0.01870 | 6.65 bits |
| 3 | 109 | 0.03261 | 0.01745 | 6.23 bits |
| 3 | 2,048 | 0.001596 | 0.000798 | 0.759 bits |

Thus the current `w=2,d=t` batch has a very strong ideal support-ranking
signal.  The three-choice `d=93` batch-failure point and candidate
two-bit-margin `d=109` reduce it to about 6--7 modeled rank bits, but negligible
placement failure does not imply negligible size-bias leakage.  Driving this
particular search metric below one bit takes roughly `d=8192` for `w=2` or
`d=2048` for `w=3`, both far outside the intended performance regime.  Even
these are low search-leakage points, not 40-bit statistical closeness.

The `166.46` figure is not a loss from a 128-bit Ring-LPN attack.  It is a gain
over uninformed exhaustive enumeration of regular supports.  In the insider
game, one party's `4*16=64` offsets are known and the other 64 offsets each
range over a block of length `L=N/16`, giving `L^64` candidates.  The modeled
expected-rank exponents are therefore:

| `N` | uniform expected rank | after ideal `w=2,d=16` gain |
|---:|---:|---:|
| `2^16` | `2^767` | about `2^600.5` |
| `2^18` | `2^895` | about `2^728.5` |
| `2^20` | `2^1023` | about `2^856.5` |

On this baseline the remaining search is nowhere near the 128-bit target.
It would be incorrect to subtract 166.46 from a nominal 128-bit Ring-LPN work
factor, since that work factor comes from a syndrome-decoding algorithm rather
than support enumeration.  This observation also does not establish security:
the open question is whether the likelihood factors can accelerate the best
syndrome-decoding attack.

### Comparison with regular noise

Regular noise is a useful control example.  Among all weight-`t` supports in a
domain of size `N`, the probability that a uniform support has exactly one
position in each of `t` equal blocks is

```
Pr[regular] = (N/t)^t / binom(N,t).
```

For `t=16`, this is `2^-19.747` at `N=2^16` and approaches `2^-19.750` as `N`
grows.  For `t=64` and `t=128`, the corresponding values are about `2^-88.0`
and `2^-179.8`.  These rare-event exponents are plainly not losses that can be
subtracted from a decoding work factor.

Esser and Santini, *Not Just Regular Decoding: Asymptotics and Improvements of
Regular Syndrome Decoding Attacks* (ePrint 2023/1568), show for random-code
syndrome decoding that ordinary ISD and regular-ISD have the same leading
asymptotic running time when `t=o(N)`; their difference is confined to
lower-order terms.  Their concrete sparse examples are equally instructive:

| `(n,k,t)` | ordinary ISD | best regular-ISD |
|---:|---:|---:|
| `(245760,245460,15)` | 127 bits | 127 bits |
| `(40960,40660,20)` | 126 bits | 126 bits |
| `(7680,7380,30)` | 127 bits | at least 132 bits |

The paper also describes the exceptions.  Regularity can give substantial
speedups at moderate relative weight, and it creates polynomial-time regimes
when the block equations nearly determine the error (`t >= k`) or when many
solutions exist.  Those are parameter-regime effects, not consequences of the
raw probability that uniform noise happens to be regular.

This is an analogy, not a reduction for our setting.  The cited results concern
random parity-check matrices.  Ring-LPN has a structured multiplication
operator, and Reverse Cuckoo supplies a candidate-dependent soft likelihood
rather than only the fixed block equations of regular noise.  The analogy does
set the right evaluation target: estimate the best attack on regular Ring-LPN
with and without the Reverse-Cuckoo likelihood, instead of treating statistical
distance of the descriptor as a security loss.

### Matching-weighted regular noise

The ideal full channel factors more strongly than the count of 64 unknown
offsets suggests.  Fix the insider's regular offsets `x[a,u]`.  For one unknown
polynomial `y[b]`, define

```
W_b(y[b]) = product over (a,k) of
    Z_{ {x[a,u] + y[b,k-u] : u in [t]} }(lambda[a,b,k]) / E_U[Z].
```

Under the independent ideal-descriptor model,

```
Pr[Y=y | X,Lambda] proportional to product_b W_b(y[b]).
```

The posterior therefore separates into four independent groups.  Each group
has `t=16` unknown offsets and `P*t=64` dense, overlapping matching-count
factors.  This exact factorization does not automatically hold after replacing
the ideal descriptors by the implementation's shared Goldreich feature root.

`rev_cuckoo_regular_noise.py` exactly enumerates one group at reduced sizes.  It
reports full posterior information, the sum of the individual-offset marginal
informations, their difference (posterior total correlation), and a one-shot
regular-Prange diagnostic.  The diagnostic chooses the least likely coordinates
from each block and measures the posterior probability that the resulting
information set is error free.  It is not a repeated-decoding work factor.

For `P=4,t=d=2` and a three-quarter information-set rate:

| block `L` | full information | one-offset information | correlation only | marginal Prange gain |
|---:|---:|---:|---:|---:|
| 8  | 0.2936 | 0.0696 | 0.2240 | 0.724 bits |
| 32 | 0.3177 | 0.0216 | 0.2961 | 0.309 bits |
| 64 | 0.3235 | 0.0115 | 0.3120 | 0.173 bits |

The full information approaches a nonzero value while individual-offset
information and the coordinate-wise information-set advantage decrease.  At
small `L`, an exact jointly coordinated two-block information-set optimizer is
stronger: its three-quarter-rate gain is 1.29 bits at `L=4`, 1.17 at `L=8`,
and 1.04 at `L=12`.  The combinatorial optimizer does not scale far enough to
establish whether that advantage vanishes.  Exact `t=3,4` experiments likewise
move information from individual marginals into pair and higher-order
correlations as `L` grows, but their available block sizes are too small for a
current-parameter extrapolation.

This motivates the following **marginal-quietness conjecture** for the ideal
channel: for fixed `P,t,w,d`,

```
I(Y[b,j] ; Lambda_b | X) -> 0 as L -> infinity,
```

while `I(Y[b] ; Lambda_b | X)` may converge to a positive constant.  The
conjecture says that the transcript asymptotically does not point to individual
coordinates; it can still identify unusually compatible tuples.

For `w=2` the tuple compatibility is far from arbitrary.  The candidate's hash
edges form a bipartite cuckoo graph.  A tree component with `e` edges has `e+1`
injective placements, a unicyclic component has two, and an overfull component
has none.  Relative to isolated edges, shared bins create a repulsive
component-size penalty.  At the current prototype point `t=d=16`, 100,000 uniform and 100,000
planted samples give:

| surrogate for `log2 Z` | uniform feasible `R^2` | planted `R^2` | planted residual stddev |
|---|---:|---:|---:|
| pair-collision count | 0.729 | 0.768 | 0.701 bits |
| occupied-bin deficit | 0.822 | 0.845 | 0.573 bits |

Thus the best current leakage description is **marginally quiet,
collision-correlated matching leakage**.  It resembles a planted repulsive
constraint problem.  This is more structured than the earlier ideal global
rank score and gives a plausible attack path: incorporate collision/component
weights into regular permutations, list construction, or nearest-neighbor
search.

A corresponding, deliberately explicit hardness conjecture is that in the
sparse regular-noise regime the matching-weighted decoder has the same leading
asymptotic exponent as the best regular decoder:

```
log T_matching-weighted / log T_regular -> 1.
```

This would still permit an important additive concrete loss.  It is not implied
by marginal quietness, by the regular-ISD paper, or by these reduced experiments.
The next cryptanalytic task is a collision-aware regular-ISD estimator; the
next channel task remains repeating the analysis with the shared Goldreich
features.

There are 256 product-list descriptors, but they constrain the same 64 unknown
offsets rather than 256 independent secrets.  The ideal rank exponent is about
`166.46/64 = 2.60` bits per hidden offset; it is not 166 bits directly revealed.
The product lists overlap algebraically, and the implementation reuses one
public Goldreich feature root across the batch.  Hence the independent-list
model is neither an upper nor a lower bound on the real joint channel.

These figures assume independent ideal random-function descriptors, disjoint
wrong supports, and an ideal decoder capable of enumerating candidates in
likelihood order.  They are neither an implemented attack nor a measurement of
the shared-Goldreich-feature channel.  The next attack step is to integrate the
matching factors into a list/ISD search; the next channel step is to measure the
joint likelihood with the implementation's shared feature matrices.

## Reproduction

Build `frontend_libOTe`, then run, for example:

```text
frontend_libOTe -invMtx -cuckooLeak -n 16 -w 2 -d 16 -trials 4194304 -seed 10
frontend_libOTe -invMtx -cuckooLeak -real -n 16 -w 2 -d 16 -trials 1048576 -seed 2
frontend_libOTe -invMtx -cuckooLeak -real -n 16 -w 3 -d 109 -batch 256 -trials 131072 -seed 1309
python experiments/rev_cuckoo_w23.py --batch 256 --target 40
python experiments/rev_cuckoo_full_lpn.py --polys 2 --weight 2 --block-length 4 --trials 1024 --seed 1
python experiments/rev_cuckoo_full_lpn_sweep.py
python experiments/rev_cuckoo_full_lpn_field_sweep.py
python experiments/rev_cuckoo_regular_noise.py --polys 4 --weight 2 --block-length 64 --trials 32 --seed 112
python experiments/rev_cuckoo_regular_noise.py --polys 2 --weight 3 --block-length 8 --trials 512 --channel occupancy --selector-k 4 --seed 7
python experiments/rev_cuckoo_regular_noise.py --polys 2 --weight 4 --block-length 6 --trials 256 --channel feasible --seed 11
python experiments/rev_cuckoo_pair_surrogate.py --n 16 --d 16 --trials 100000 --seed 101
```

Current-profile rows, Wilson intervals, and information-density quantiles are
in `rev_cuckoo_results.csv`.  The continuous two/three-choice sweep and its
deterministic seeds are in `rev_cuckoo_w23.csv`; the ordinary-versus-partitioned
sanity sweep is in `rev_cuckoo_normal_control.csv`.  The full-game toy sweeps
are in `rev_cuckoo_full_lpn_results.csv` and
`rev_cuckoo_full_lpn_field_results.csv`; the current-scale rank moments are
in `rev_cuckoo_rank_results.csv`.
