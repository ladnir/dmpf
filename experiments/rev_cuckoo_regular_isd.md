# Leakage-aware regular-Prange diagnostic

`rev_cuckoo_regular_isd.py` measures whether the Reverse-Cuckoo posterior can
increase the success probability of a genuine regular-Prange iteration.  This
goes beyond leakage-only candidate rank: the experiment selects exactly `N`
columns of the structured Ring-LPN matrix and verifies that the selected
square matrix is full rank.

## Experiment

There are `P` hidden regular polynomials, each with one nominal nonzero in
each of `t` blocks of length `L`.  The syndrome matrix is

```text
H = [ I | M_r1 | ... | M_r(P-1) ],
```

where every `M_ri` is the negacyclic multiplication matrix of a fresh uniform
ring element.  A regular-Prange iteration keeps `L/P` positions in every
block, giving exactly `N=tL` columns.  Its uniform support success probability
is

```text
(1/P)^(P*t).
```

For every hidden polynomial, the experiment compares:

- `marginal`: keep the `L/P` most likely positions in every block;
- `pair`: fit the descriptor's log likelihood by its public pair-collision
  count, run damped mean-field inference on the resulting pairwise constraint
  system, and keep the `L/P` largest beliefs in every block;
- `joint`: start from the marginal rectangle and random rectangles, then use
  block-coordinate ascent to maximize full posterior rectangle mass.

The posterior mass of the selected rectangle is the exact conditional support
success probability.  Multiplying it by the selected-matrix rank indicator
gives one-iteration Prange success.  Over `F_65537`, every selected matrix in
the recorded experiments was full rank, so the reported gains below are
support gains rather than small-field rank artifacts.

The `pair` selector is polynomial and uses only public descriptors.  For each
hidden polynomial it groups the pair edges into `Ptw` complete collision
factors (here `w=2`).  Aggregating the per-bin belief mass computes a
mean-field round in `O(Pwt^2L)` time and `O(Pwt^2L)` label storage, rather than
expanding `Pt choose(t,2)w` duplicated edge tables.  Fitting the channel slope
is an offline parameter calculation.  The experiment still enumerates the
exact `L^t` posterior to evaluate the selected rectangle without sampling rare
support hits.  The `marginal` and `joint` selectors themselves require this
posterior, and the joint local optimizer is not guaranteed globally optimal.

## Results

All rows use two-choice Reverse Cuckoo with `d=t`, four joint-selector random
restarts, and the same seeded support and ring-matrix samples across channels.

| `P` | `t` | `L` | channel | marginal gain | joint gain |
|---:|---:|---:|---|---:|---:|
| 2 | 3 | 8 | raw | 1.420 | 1.713 |
| 2 | 3 | 8 | occupancy `K=4` | 0.495 | 0.662 |
| 2 | 3 | 8 | occupancy `K=8` | 0.291 | 0.427 |
| 2 | 3 | 8 | feasibility only | 0.127 | 0.179 |
| 2 | 4 | 6 | raw | 2.565 | 2.884 |
| 2 | 4 | 6 | occupancy `K=4` | 0.930 | 1.170 |
| 2 | 4 | 6 | occupancy `K=8` | 0.563 | 0.799 |
| 2 | 4 | 6 | feasibility only | 0.216 | 0.288 |
| 3 | 3 | 6 | raw | 4.195 | 5.271 |
| 3 | 3 | 6 | occupancy `K=4` | 1.285 | 1.914 |
| 3 | 3 | 6 | occupancy `K=8` | 0.719 | 1.146 |
| 3 | 3 | 6 | feasibility only | 0.304 | 0.430 |
| 4 | 3 | 8 | raw | 7.432 | 9.779 |
| 4 | 3 | 8 | occupancy `K=4` | 2.151 | 3.613 |
| 4 | 3 | 8 | occupancy `K=8` | 1.371 | 2.267 |
| 4 | 3 | 8 | feasibility only | 0.575 | 0.832 |
| 4 | 3 | 16 | occupancy `K=8` | 1.069 | 1.945 |
| 4 | 3 | 16 | feasibility only | 0.484 | 0.804 |
| 4 | 4 | 8 | occupancy `K=4` | 4.130 | 6.537 |
| 4 | 4 | 8 | occupancy `K=8` | 2.526 | 4.634 |
| 4 | 4 | 8 | feasibility only | 1.017 | 1.476 |

The `P=4` rows use the production ratio of four hidden polynomials and a
quarter of every regular block in the solve set.  They are still tiny in `t`
and `L`.

## Interpretation

1. The leakage produces a measurable regular-Prange advantage, not merely an
   exhaustive-ranking signal.  Therefore it is incorrect to say that linear
   attacks cannot use the descriptors.
2. Finite `K` is effective.  At `P=4,t=3,L=8`, `K=8` reduces the measured
   joint gain from 9.78 to 2.27 bits.
3. Hard feasibility remains exploitable: its joint gain is 0.83 bits in that
   row and 1.48 bits at `P=4,t=4,L=8`.
4. Correlation matters.  The joint rectangle consistently improves on the
   independently selected marginal rectangle.
5. Doubling `L` from 8 to 16 at `P=4,t=3,K=8` reduces the marginal gain from
   1.37 to 1.07 bits but the joint gain only from 2.27 to 1.95 bits.  The
   feasibility joint gain is essentially stable (0.83 to 0.80 bits).  This is
   consistent with the marginally-quiet, correlation-persistent hypothesis.

### Scalable pair-collision selector

The polynomial selector gives the following support gains.  The marginal and
joint columns are shown only as exact-posterior reference points from the same
runs; they are not efficient selectors.

| `P` | `t` | `L` | channel | fitted slope | fit `R^2` | marginal | pair | joint |
|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 4 | 3 | 8 | raw | -0.499 | 0.959 | 7.332 | 7.524 | 9.727 |
| 4 | 3 | 8 | occupancy `K=4` | -0.152 | 0.932 | 2.136 | 2.225 | 3.565 |
| 4 | 3 | 8 | occupancy `K=8` | -0.078 | 0.922 | 1.283 | 1.334 | 2.200 |
| 4 | 3 | 8 | feasibility only | 0 | 0 | 0.552 | -0.059 | 0.823 |
| 4 | 3 | 16 | occupancy `K=4` | -0.152 | 0.932 | 1.559 | 1.644 | 3.031 |
| 4 | 3 | 16 | occupancy `K=8` | -0.078 | 0.922 | 1.024 | 1.051 | 1.962 |
| 4 | 4 | 8 | occupancy `K=4` | -0.175 | 0.729 | 4.191 | 4.169 | 6.472 |
| 4 | 4 | 8 | occupancy `K=8` | -0.095 | 0.483 | 2.474 | 2.353 | 4.607 |

The pair selector recovers essentially all of the independently usable
marginal signal and sometimes a little more.  It does not recover most of the
gap to the exponential joint selector.  In particular, feasibility-only
likelihood is constant over feasible graphs, so its fitted pair slope is zero;
the small negative measured gain is sampling noise around the uniform
baseline.  Exploiting that floor requires a higher-order Hall, cycle, or
matching-aware selector.

At these reduced parameters, the demonstrated polynomial attack saves about
`4.2` bits in the most adverse recorded `K=4` row and `2.4` bits for `K=8`.
Those are small additive losses, not evidence of a qualitative break.  They
also exclude the selector's preprocessing cost, which must be charged in any
production work-factor claim.  The decrease from `2.22` to `1.64` bits when
`L` doubles at `P=4,t=3,K=4` is encouraging, but neither that trend nor the
growth with `t` can yet be extrapolated to `t=16,L>=4096`.

### Production-shaped inclusion diagnostic

`rev_cuckoo_regular_isd_mc.py` runs the polynomial selector without enumerating
the posterior.  The robust statistic currently available at large `t` is the
fraction of planted coordinates contained in the selected quarter.  Its
baseline is exactly `0.25`.

| `t` | `L` | channel | trials | coordinate inclusion |
|---:|---:|---|---:|---:|
| 16 | 64 | occupancy `K=4` | 8 | 0.432 (`SE` 0.021 across trials) |
| 16 | 64 | occupancy `K=8` | 8 | 0.346 (`SE` 0.024 across trials) |
| 16 | 256 | occupancy `K=4` | 2 | 0.328 |
| 16 | 256 | occupancy `K=8` | 2 | 0.258 |
| 16 | 1024 | occupancy `K=4` | 4 | 0.285 (`SE` 0.024 across trials) |
| 16 | 1024 | occupancy `K=8` | 4 | 0.266 (`SE` 0.031 across trials) |
| 16 | 4096 | occupancy `K=4` | 4 | 0.277 (`SE` 0.031 across trials) |
| 16 | 4096 | occupancy `K=8` | 4 | 0.297 (`SE` 0.011 across trials) |

The strong small-`L` coordinate signal decays toward the baseline as `L`
grows.  At `L=4096`, neither four-trial result establishes a coordinate
advantage; the trial-level standard errors omit within-trial coordinate
dependence and must not be read as formal confidence intervals.  More trials
are needed, particularly because the apparent ordering of `K=4` and `K=8`
reverses in this small sample.

The same program implements the exact identity
`gain = E[W | rectangle] / E[W]` using defensive importance sampling.  It
matches exact enumeration at `t=3`, but at `t=16,L=64` the observed effective
sample fractions were only about `4e-4`--`6e-4`.  Provisional estimates in the
30--50 bit range are therefore discarded: the proposal misses important
likelihood tails.  This failure motivates the sequential Monte Carlo scorer
below.  Coordinate inclusion alone cannot be exponentiated under an
independence assumption because the central issue is correlation.

### Sequential Monte Carlo scoring

`rev_cuckoo_regular_isd_smc.py` addresses the hard feasibility tail by adding
the `Pt` exact list-likelihood factors sequentially.  It resamples after every
factor and applies exact single-coordinate Metropolis moves against the partial
target.  Every estimate reports replicate uncertainty, minimum and mean stage
ESS, post-resampling uniqueness, and mutation acceptance.

At `P=4,t=3,L=8,K=4`, 512 particles, two moves per stage, and four replicates
estimate the four per-polynomial gains within `0.026` bits of full enumeration
on the identical instance.  This validates the estimator in the regime where
the answer is exactly available.

The first production screen estimated the full-space and rectangle partition
functions separately.  It reported per-polynomial `K=4` gains
`-1.87,-3.36,2.91,3.29` bits and suggested a `2.14`-bit log-mean.  These values
are now **retracted**.  The estimator has a heavy factor-order tail: in a
batched run its dominant component was initially `+5.66` bits, but adding six
replicates changed the identical component to `-3.30` bits despite usable
per-stage ESS.  Subtracting the two global log normalizers understates this
tail uncertainty.

The replacement splitting estimator constructs the exact posterior SMC
population once and imposes the rectangle restrictions one coordinate at a
time.  It estimates sixteen conditional probabilities near `1/4`, resampling
and applying exact-target mutation after each.  On `t=3,L=8`, it recovers exact
rectangle gains within about `0.11` bits.  At `t=16,L=4096`, however, the
independently thresholded mean-field rectangles cause late conditional
survival to fall to `1/256` or zero.  Diversity refreshes and sixteen mutation
moves per restriction do not repair them.  Their precise large negative gains
are unresolved, but they are not useful attacks: an attacker falls back to a
uniform rectangle.

We also added polynomial coordinate ascent that directly optimizes the
fixed-cardinality rectangle under the fitted pair-collision model.  On an
exact toy instance it raises true gain from `0.810` to `0.956` bits.  On the
tested production descriptor it still reaches `1/256` or zero late-stage
survival under the exact likelihood.  Thus pairwise refinement does not yet
convert the surviving higher-order correlation into a production regular-
Prange advantage.

Finally, we tested sampled exact matching-factor messages.  With the other
`t-1` offsets fixed, all `L` values of one offset induce at most `d^2` distinct
endpoint pairs; memoizing those pairs makes an exact component, feasibility,
and placement-multiplicity message practical.  Coordinate updates maximize a
Monte Carlo approximation to the sum of exact factor log-normalizers.  On the
same exact toy, two matching-message rounds raise true gain further from
`0.956` to `1.083` bits.  At `t=16,L=4096,K=4`, the resulting rectangle still
has only `1/256`-scale late conditional survival and no demonstrated positive
gain under the hardened splitting scorer.

The attempted shared-known-support batch consequently provides no evidence
for the earlier speculative `8.6`-bit loss.  This does not prove that every
polynomial-time rectangle attack fails; it shows that the current mean-field
pair-coordinate, and local exact-matching-message selectors fail on the full
higher-order channel.  A more global Hall-aware optimizer remains possible,
but this completes the planned marginal, pairwise, and matching-aware attack
portfolio.

### Shared Goldreich feature channel

`rev_cuckoo_goldreich.py` replaces the independent random-function sampler by
the feature family used in `GoldreichHash.h`.  For each of the two partitions
it samples shared binary maps

```text
a = M0 x
b = M1 x
c = a & b
V(x) = M2 (x || c) + M2 (N || c(N)),
```

so `V(N)=0`, and reuses the same public feature root for all sets in the
attacked batch.  For every set it samples the revealed coefficient solution
uniformly from the affine fiber that programs the secretly assigned endpoint.
The harness matches the implementation's byte-rounded input and intermediate
widths and its effective feature width `c=d+linear_security` (currently
`d+10`).  Finite-`K` proposals are evaluated on their active rows, and only the
selected proposal is expanded over the full domain.

For a candidate active set `A` and a consistent placement `rho`, let `A_0` and
`A_1` be the items assigned to its two partitions.  Relative to uniform public
coefficient matrices, the exact raw likelihood is

```text
  1/(2d)_|A| * sum_rho 2^(log2(d) *
      (rank(V_0(A_0)) + rank(V_1(A_1)))).
```

If every placed feature submatrix has full row rank, the exponent is constant
and this reduces exactly to the ideal placement likelihood.  The reduced
driver enumerates all placements and uses the rank-weighted formula.  The
production driver samples the actual shared-root descriptors but uses the
good-rank simplification after ten items; this limitation is reported as
`approximate_rank_queries`.

At `P=2,t=3,d=4,L=8,K=4`, exact rank-aware runs give:

| $(t,d,L)$ | seed | trials | hash family | marginal support gain | pair support gain | joint support gain |
|---:|---:|---:|---|---:|---:|---:|
| $(3,4,8)$ | 307 | 256 | independent ideal | 0.390 | 0.389 | 0.525 |
| $(3,4,8)$ | 307 | 256 | shared Goldreich | 0.396 | 0.400 | 0.533 |
| $(3,4,8)$ | 911 | 128 | independent ideal | 0.386 | 0.390 | 0.530 |
| $(3,4,8)$ | 911 | 128 | shared Goldreich | 0.384 | 0.389 | 0.526 |
| $(4,4,6)$ | 1201 | 128 | independent ideal | 0.921 | 0.886 | 1.182 |
| $(4,4,6)$ | 1201 | 128 | shared Goldreich | 0.945 | 0.914 | 1.217 |
| $(4,4,6)$ | 2203 | 128 | independent ideal | 0.935 | 0.896 | 1.191 |
| $(4,4,6)$ | 2203 | 128 | shared Goldreich | 0.941 | 0.909 | 1.196 |

At $t=3$, the shared-minus-ideal joint differences are `+0.008` and `-0.004`
bit.  At $t=4$, they are `+0.035` and `+0.005` bit; pooling the two equal-size
runs gives a difference of about `+0.02` bit, below the sampling resolution.
Marginal and pair differences are likewise small and inconsistent across
seeds.  These experiments detect no stable additional leakage from sharing
the Goldreich root as the matching size grows from three to four.  Three of
36,960 attempted $t=3$ systems and two of 32,804 attempted $t=4$ systems were
rank deficient and inconsistent.  This is evidence for the heuristic, not a
joint pseudorandomness proof.

The same two $t=d=4$ seeds give a complete finite-$K$ curve under the shared
Goldreich family.  The pooled gain averages the support probabilities from
the two equal-size runs before taking the logarithm.

| $K$ | relative production OT model | joint gain, seed 1201 | joint gain, seed 2203 | pooled joint gain | gain removed by doubling |
|---:|---:|---:|---:|---:|---:|
| 1 | 1.00x | 2.905 | 2.881 | 2.893 | -- |
| 2 | 1.35x | 1.804 | 1.863 | 1.834 | 1.059 |
| 4 | 1.73x | 1.217 | 1.196 | 1.207 | 0.627 |
| 8 | 2.50x | 0.811 | 0.809 | 0.810 | 0.397 |

Every doubling helps, but its marginal benefit decreases while its marginal
MPC cost increases.  Moving from $K=4$ to $K=8$ adds about 45 percent to the
total modeled OTs and removes only another 0.40 bit in this toy.  Thus $K=4$
is a defensible reduced-attack Pareto knee, not a security threshold: without
an acceptable production leakage budget the table cannot declare either
$K=2$ insufficient or $K=8$ unnecessary.

At `P=4,t=d=16,L=4096,K=4`, the matching-aware splitting run with 256
particles and the implementation width `c=26` estimated `-0.37 +/- 1.00`
bits.  Late rectangle survival was `1/256`, so the estimate is saturated and
the attacker falls back to ordinary Prange.  A separate production-shaped
batch audit found zero rank-deficient and zero inconsistent systems among 512
planted affine solves.  The actual shared family therefore does not rescue the
current local rectangle attack, but the exact toy result means the family
cannot simply be identified with independent random functions in a proof.

We repeated that production screen for raw $K=1$ and occupancy $K=2$.  The
fitted pair slopes were respectively `-0.543`, `-0.440`, and `-0.334` for
$K=1,2,4$, confirming that finite $K$ smooths the local collision signal.  But
all three selected rectangles reached only `1/256` late survival.  Their large
negative splitting estimates are saturation artifacts, not security gains.
The current production rectangle optimizer therefore finds no attack for any
tested $K$ and cannot choose between them.  We did not run the same saturated
screen for $K=8$; the exact finite-$K$ curve is more informative until a
stronger production optimizer is available.

We then replaced the factorwise matching selector by a full-batch conditional
selector.  It first samples the complete posterior with SMC.  For each offset,
posterior particles provide contexts for the other offsets; all $L$ values are
scored by the product of all $Pt$ exact Reverse-Cuckoo likelihood factors and
normalized within each context.  The selector greedily restricts each block
and rejuvenates its training population.  A separate SMC population scores the
result, so the selected quarter cannot merely memorize training particles.

On eight exactly enumerable $P=2,t=3,d=4,L=8,K=4$ instances, its mean paired
exact-gain difference from the previous rectangle was about `-0.001` bit, with
wins and losses in both directions.  Thus the implementation reproduces the
available toy quality without an optimistic training bias.  At the production
$P=4,t=d=16,L=4096$ raw-$K=1$ point, however, both four-context and sixteen-
context versions had zero late survival in independent scoring populations.
This negative result uses complete Hall, matching, and multiplicity likelihoods
across the batch.  It points to the single fixed rectangle as the obstruction,
rather than missing structure in the selector's factor approximation.  A next
attack must use a charged portfolio of mode-specific rectangles or abandon the
rectangular regular-Prange information set.

We tested both continuations.  A charged four-rectangle portfolio was worse
than the ordinary rectangle on exact toys, so merely retrying posterior modes
does not recover the missing gain.  We then abandoned the full rectangle and
restricted only one offset.  The hardened selector chooses either one value or
an exact-width set from a posterior-particle histogram, while an independent
SMC population scores that fixed choice by averaging exact full conditionals.
For a set of size $s$, the charged gain is
`log2(p_set / (s/L))`.  This holdout step is essential: reusing a handful of
training contexts produced spurious 7--10 bit singleton gains by maximizing
over $tL$ noisy estimates.  Those numbers are discarded.

The independent estimator agrees with exact enumeration.  At
`P=2,t=3,d=4,L=8,K=1`, a two-value set gives `0.365` estimated bits versus
`0.387` exact bits.  At `P=4,t=d=16,L=4096` on seed 251, a 256-value set has
posterior mass `0.0569` for occupancy `K=4`, below its uniform baseline
`256/4096 = 0.0625`; the attacker therefore falls back to zero gain.  The same
test at `K=1` gives mass `0.0178`, also no demonstrated gain, but its selector
and evaluator minimum ESS values are only `0.064` and `0.083`.  Thus `K=1` is
inconclusive because independent SMC runs do not mix across the same posterior
modes, whereas the tested `K=4` transcript has usable diagnostics and no
one-offset compression.  This does not prove either parameter secure, and it
does not turn `K=4` into a threshold; it removes an estimator artifact and
narrows the next task to better multimodal sampling or a different decoder.

### Posterior entropy and effective support

The failure of a fixed-set decoder does not mean the leakage is small.  For a
transcript $\Lambda$, let $W_\Lambda(x)$ be its likelihood under hidden offset
vector $x\in[L]^t$.  The SMC normalizing constant and posterior particles give

`D_KL(p(x|Lambda) || U) = E_post[log2 W] - log2 E_U[W]`.

They also give the collision-entropy loss

`D_2(p(x|Lambda) || U) = log2 E_U[W^2] - 2 log2 E_U[W]`.

The second expression is estimated equivalently as
`logmeanexp2_post(log2 W) - log2 E_U[W]`.  On an exactly enumerable
`P=2,t=3,d=4,L=8,K=1` instance, four SMC replicates estimate KL loss `0.300`
and collision loss `0.550` bits, versus exact values `0.294` and `0.547`.

At production `P=4,t=d=16,L=4096`, the prior offset entropy is 192 bits.  The
following table reports residual collision entropy; ranges are over independent
SMC populations for the displayed transcripts.

| channel | transcript | residual collision bits |
|---|---:|---:|
| shared Goldreich, $K=1$ | 251 | 134.2--140.8 (mean 136.2) |
| shared Goldreich, $K=1$ | 252 | 137.3--147.9 (mean 142.6) |
| shared Goldreich, $K=2$ | 251 | 155.1--159.0 (mean 156.7) |
| shared Goldreich, $K=2$ | 252 | 159.83--159.85 |
| shared Goldreich, $K=4$ | 251 | 161.2--163.9 (mean 162.5) |
| shared Goldreich, $K=4$ | 252 | 163.6--165.3 (mean 164.4) |
| shared Goldreich, $K=8$ | 251 | 155.9--165.2 (mean 161.5) |
| ideal hashes, $K=1$ | 251 | 129.8--134.4 (mean 132.2) |
| ideal hashes, $K=2$ | 251 | 149.2--153.8 (mean 151.5) |

Thus the large $K=1$ loss is inherent to Reverse Cuckoo rather than caused by
the shared Goldreich evaluator.  The posterior remains enormous, supporting
the sparsity heuristic, but ideal $K=1$ has as little as a 1.8-bit observed
collision-entropy margin over 128 bits.  By contrast, $K=2$ retains at least
about 149 bits in the current screens and captures most of the improvement of
$K=4$.  The modeled amortized costs are respectively `1.00x`, `1.37x`, and
`1.73x`, making $K=2$ the production support-entropy knee.

This has a precise but limited list interpretation.  For posterior collision
entropy $H_2$, Cauchy--Schwarz gives
$p_\Lambda(S)\leq\sqrt{|S|2^{-H_2}}$.  A leakage-only list containing mass
$\epsilon$ therefore needs at least $\epsilon^2 2^{H_2}$ candidates.  At the
higher-quality shared-$K=1$ range $H_2=134$--141, capturing half the mass needs
at least about $2^{132}$--$2^{139}$ candidates.  This does not lower-bound an
algebraic attack that combines the leakage with the LPN equations.

Increasing regular weight is less predictable.  At fixed $N=65536$, ideal
$K=1,t=17$ leaves roughly 129--147 collision bits across two transcripts, and
$t=18$ leaves roughly 144--159 after correcting the experimental equal-block
domain to $t\log_2(N/t)$.  Although $t=18,K=1$ and $t=16,K=2$ both cost about
`1.35x`, the finite-$K$ change is more stable in the current data.

The preceding collision figures use the posterior-particle moment
$E_{p_\Lambda}[W]$ and can miss rare high-weight modes.  We therefore added a
stronger continuation from target $W$ to target $W^2$.  A second copy of every
descriptor factor is introduced in four quarter-power steps, with resampling
and rejuvenation after every increment.  On the exact toy this tempered
estimate is `0.553` bits versus the exact `0.547`; powered-stage minimum ESS is
above `0.99` there.

At production parameters, two 1024-particle shared-$K=1$ continuations leave
only `138.09` and `131.46` collision bits.  A lower-budget population reaches
`127.73` bits.  Thus the moment-only `134--148` range was optimistic for some
modes, and $K=1,t=16$ has no defensible 128-bit margin.  In contrast, three
quarter-tempered shared-$K=2$ populations across transcripts 251 and 252 leave
`151.64`, `155.89`, and `156.50` bits, with powered minimum ESS between `0.91`
and `0.96`.  Finally, tempering the weaker ideal $t=18,K=1$ transcript leaves
`131.65` equal-block bits, or only `128.59` bits after correcting to fixed
$N=65536$.  Increasing regular weight therefore does not reliably replace the
finite-$K$ smoothing.  The current optimized recommendation is
$t=16,K=2$.

## Reproduction

```text
python experiments/rev_cuckoo_regular_isd.py --polys 4 --weight 3 --block-length 8 --modulus 65537 --trials 64 --channel occupancy --selector-k 8 --seed 61
python experiments/rev_cuckoo_regular_isd.py --polys 4 --weight 3 --block-length 16 --modulus 65537 --trials 32 --channel feasible --seed 81
python experiments/rev_cuckoo_regular_isd.py --polys 4 --weight 4 --block-length 8 --modulus 65537 --trials 32 --channel occupancy --selector-k 8 --seed 71
python experiments/rev_cuckoo_regular_isd.py --polys 4 --weight 3 --block-length 16 --modulus 65537 --trials 32 --channel occupancy --selector-k 4 --seed 101
python experiments/rev_cuckoo_regular_isd_mc.py --polys 4 --weight 16 --block-length 4096 --trials 4 --samples 0 --channel occupancy --selector-k 4 --seed 161
python experiments/rev_cuckoo_regular_isd_smc.py --polys 4 --hidden-polys 1 --weight 16 --block-length 4096 --trials 1 --particles 256 --moves 2 --restriction-moves 16 --refresh-threshold 0.2 --replicates 2 --rectangle-restarts 4 --rectangle-rounds 16 --channel occupancy --selector-k 4 --seed 251 --score-mode splitting
python experiments/rev_cuckoo_regular_isd_smc.py --polys 4 --hidden-polys 1 --weight 16 --block-length 4096 --trials 1 --particles 256 --moves 2 --restriction-moves 16 --refresh-threshold 0.2 --replicates 2 --rectangle-restarts 4 --rectangle-rounds 16 --matching-contexts 4 --matching-rounds 2 --channel occupancy --selector-k 4 --seed 251 --score-mode splitting
python experiments/rev_cuckoo_regular_isd.py --polys 2 --weight 3 --block-length 8 --partition-size 4 --trials 256 --joint-restarts 4 --pair-restarts 4 --pair-rounds 16 --channel occupancy --selector-k 4 --hash-family goldreich --linear-security 10 --exact-rank-limit 3 --seed 307
python experiments/rev_cuckoo_regular_isd.py --polys 2 --weight 4 --block-length 6 --partition-size 4 --modulus 65537 --trials 128 --joint-restarts 4 --pair-restarts 4 --pair-rounds 20 --pair-fit-samples 20000 --channel occupancy --selector-k 4 --hash-family goldreich --linear-security 10 --exact-rank-limit 4 --seed 1201
python experiments/rev_cuckoo_regular_isd_smc.py --polys 4 --hidden-polys 1 --weight 16 --block-length 4096 --partition-size 16 --trials 1 --particles 256 --moves 2 --restriction-moves 16 --refresh-threshold 0.2 --replicates 2 --pair-restarts 0 --pair-rounds 8 --rectangle-restarts 4 --rectangle-rounds 16 --matching-contexts 4 --matching-rounds 2 --pair-fit-samples 50000 --channel occupancy --selector-k 4 --hash-family goldreich --linear-security 10 --exact-rank-limit 10 --seed 251 --score-mode splitting
python experiments/rev_cuckoo_regular_isd_smc.py --polys 4 --hidden-polys 1 --weight 16 --block-length 4096 --partition-size 16 --trials 1 --particles 256 --moves 2 --restriction-moves 16 --refresh-threshold 0.2 --replicates 2 --pair-restarts 0 --pair-rounds 8 --rectangle-restarts 2 --rectangle-rounds 8 --matching-rounds 0 --joint-contexts 16 --joint-particles 256 --joint-moves 2 --joint-restriction-moves 16 --pair-fit-samples 50000 --channel raw --selector-k 1 --hash-family goldreich --linear-security 10 --exact-rank-limit 10 --seed 251 --score-mode splitting
python experiments/rev_cuckoo_regular_isd_portfolio.py --polys 4 --weight 16 --block-length 4096 --partition-size 16 --offset-set --guess-width 256 --guess-contexts 64 --selector-particles 512 --selector-moves 3 --channel occupancy --selector-k 4 --hash-family goldreich --linear-security 10 --exact-rank-limit 10 --seed 251
python experiments/rev_cuckoo_posterior_info.py --polys 4 --weight 16 --block-length 4096 --partition-size 16 --particles 512 --moves 3 --replicates 4 --channel occupancy --selector-k 2 --hash-family goldreich --linear-security 10 --exact-rank-limit 10 --seed 251
```
