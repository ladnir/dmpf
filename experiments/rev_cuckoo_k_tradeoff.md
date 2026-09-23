# Reverse-Cuckoo finite-K cost/leakage frontier

This note joins `rev_cuckoo_k_tradeoff.py` with the implementation-aligned OT
model in `rev_cuckoo_cost.py`.  It is a design diagnostic, not a Ring-LPN
security reduction.  In particular, the reported search gain is the ideal
independent-list likelihood-ranking metric from `rev_cuckoo_leakage.md`; it
must not be subtracted from the nominal ISD exponent.

## Model

For each descriptor, draw `K` ordinary Reverse-Cuckoo proposals.  The exact
selector gives proposal `i` weight `1/Z_i`.  The cheaper occupancy selector
uses the public lookup table

```text
a(D) = 1 / E_U[Z | D],
```

where `D` is the total endpoint occupancy deficit.  A richer moments selector
conditions on `(D,C2)`, with `C2` the number of colliding endpoint pairs.  It
reuses the endpoint-equality matrix and is provisionally assigned the same OT
cost as occupancy; it still needs a separate exact circuit count if retained.

For occupancy, the MPC circuit can use the existing equality and OR-tree
machinery directly.  In each partition the OR tree returns, for every row
after the first, a bit saying that its endpoint appeared earlier.  Summing
these `w(t-1)` bits gives `D` directly.  There is no need to materialize the
complementary first-occurrence flags or subtract the occupied count from
`wt`.

The weighted choice also needs no secret division.  For candidate `i` and
each public deficit value `D`, public randomness supplies an independent
fixed-width exponential-race key with rate `a(D)`.  A secret 5-bit lookup
selects the key belonging to the candidate's `D`, and a tournament retains
the smallest of the four keys.  Conditional on exact continuous keys, this
chooses `i` with probability `a(D_i)/sum_j a(D_j)`.  The implementation should
define the finite-precision race itself as the selector distribution and use
enough key bits that its distance from the continuous diagnostic is below the
desired bound.

The table `a(D)` itself need not be estimated.  `Pr_U[D]` follows from the
standard occupied-bin recurrence.  To obtain `E_U[Z 1_D]`, sum over the
number `a` of planted rows assigned left.  Conditioned on a fixed placement,
the left partition starts with `a` distinct occupied bins and receives
`t-a` independent draws, while the right starts with `t-a` occupied bins and
receives `a` draws.  Therefore

```text
a(D) = Pr_U[D] / E_U[Z 1_D]
```

is computed exactly offline (up to chosen public arithmetic precision).  The
implementation model uses 80-bit race keys.  For `(16,32)` the exact table's
nonzero rates lie between `2^-16` and `1.195`.  Quantizing `[0,2^22)` at
spacing `2^-58` makes the truncation probability at most `4*exp(-64)` and
leaves ample margin below `2^-40` for rounding/ties, provided the public table
generator uses correctly rounded high-precision arithmetic.

The baseline is `t=16`, equal partition sizes `(16,16)`, four noise
polynomials, `N=2^20`, and ten stationary expansions.  Its modeled total is
`3,600,384` OTs.  Costs below are relative to that point.  Non-power-of-two
equal partition sizes are idealized; `(16,32)` is the directly relevant
power-of-two alternative to equal `(24,24)`.

For unequal sizes `(d0,d1)`, the planted injection must not be uniform over
the `d0+d1` slots.  If a placement puts `a` rows on the left, assign that
particular placement probability proportional to

```text
1 / (d0^a d1^(t-a)).
```

Programming its endpoints multiplies the uniform descriptor density by the
inverse factor, so every compatible `(descriptor,placement)` pair again has
the same joint mass and the output density remains proportional to `Z`.  The
induced public distribution of `a` is proportional to

```text
binom(t,a) (d0)_a (d1)_(t-a) / (d0^a d1^(t-a)).
```

The cost model includes this sampler for unequal partitions.  Jointly sample
a secret 64-bit uniform value and invert the fixed 17-point CDF with five comparisons.
Decode the resulting secret count `a` into prefix bits.  Independently apply
hidden random permutations to the 16 rows and to the bins of both partitions.
At permuted rank `j`, use the left bin permutation when `j<a` and the right bin
permutation otherwise.  Conditioned on `a`, every type-`a` injection is
uniform; multiplying by the chosen law of `a` gives each particular injection
the required mass proportional to `1/(d0^a d1^(t-a))`.

For `(16,32)`, one candidate uses 722 OTs for the secret inverse-CDF sample,
512 OTs for the row and two bin permutations, and 32 OTs for the endpoint
muxes, for 1,266 OTs per descriptor.  The 64-bit CDF discretization contributes
at most about `17/2^64` statistical distance per sample; even over four
candidates and 256 descriptors this remains below `2^-40`.

## Preliminary frontier

| partitions | selector | K | modeled search gain | relative OTs |
|---|---|---:|---:|---:|
| `(16,16)` | raw | 1 | 166.5 bits | 1.00x |
| `(16,16)` | occupancy | 2 | 113.8 bits | 1.35x |
| `(16,16)` | occupancy | 4 | 84.4 bits | 1.73x |
| `(16,16)` | occupancy | 8 | 69.0 bits | 2.50x |
| `(20,20)` | moments | 2 | 73.6 bits | 1.67x |
| `(20,20)` | moments | 4 | 48.4 bits | 2.10x |
| `(24,24)` | moments | 2 | 53.4 bits | 2.04x |
| `(24,24)` | moments | 4 | 33.0 bits | 2.49x |
| `(16,32)` | occupancy | 2 | 58.4 bits | 2.22x |
| `(16,32)` | occupancy | 4 | 37.7 bits | 2.76x |
| `(32,32)` | moments | 2 | 34.4 bits | 2.88x |
| `(32,32)` | moments | 4 | 19.4 bits | 3.35x |

The limiting `U | Z>0` search gains are approximately 23.7 bits for
`(16,16)`, 6.5 bits for `(20,20)`, 2.4 bits for `(24,24)`, 3.4 bits for
`(16,32)`, and 0.49 bits for `(32,32)`.  Finite-K selection cannot pass these
limits without also changing the feasibility channel.

## Interpretation

The independent-list diagnostic suggests a fixed-size knee near `K=4`, but it
is not an attack cost and cannot select the security parameter by itself.  The
implementation-shaped exact regular-Prange experiment supplies a more relevant
reduced curve at `P=2,t=d=4,L=6`: the two-seed pooled joint gains for
`K=1,2,4,8` are `2.893, 1.834, 1.207, 0.810` bits.  The gain removed by each
doubling decreases from `1.059` to `0.627` to `0.397` bit, while relative total
OTs grow as `1.00x, 1.35x, 1.73x, 2.50x`.  This makes `K=4` a defensible Pareto
knee in the reduced attack, rather than a proof-derived threshold.  Moving to
`K=8` approximately doubles the incremental debiasing work for the smallest
observed marginal reduction.

A later production one-offset study does not sharpen this choice.  Its
selector chooses an exact-width posterior value set and an independent SMC
population scores the set against the exact `s/L` baseline.  The independent
score is important: maximizing and scoring singleton guesses on the same few
contexts creates spurious 7--10 bit gains.  With the corrected holdout method,
a width-256 set at `P=4,t=d=16,L=4096,K=4` captures mass `0.0569`, below the
uniform `0.0625`, and therefore demonstrates zero gain after fallback.  The
corresponding `K=1` sampler has minimum ESS below 0.09 and does not reproduce
the same posterior modes across independent runs, so it is inconclusive rather
than evidence that `K=1` is either safe or broken.  Consequently the current
case for `K=4` remains the reduced exact curve plus the MPC cost knee.

A production posterior-entropy diagnostic changes the empirical knee.  It
estimates both Shannon loss and Renyi-2 loss from independently replicated SMC
normalizing constants and posterior samples, and matches exact toy values to
within `0.01` bit.  At `P=4,t=d=16,L=4096`, the residual collision-entropy
ranges across two transcripts are about `134--148` bits for shared `K=1`,
`155--160` for shared `K=2`, and `161--165` for shared `K=4`.  Ideal hashes
leave only `130--134` bits for `K=1` and `149--154` for `K=2`, showing that the
large raw loss is not an artifact of the shared Goldreich evaluator.  Since
the modeled amortized costs are `1.00x`, `1.37x`, and `1.73x`, this metric makes
`K=2` the current production cost/support knee.  It does not turn collision
entropy into an LPN attack bound, and the `K=1` posterior remains large, but a
roughly two-bit worst observed ideal margin is too small to absorb modeling
error at a 128-bit target.

A quarter-step $W\to W^2$ continuation supersedes the posterior-moment
collision estimates for the most security-relevant points.  It matches exact
Renyi-2 loss (`0.553` estimated versus `0.547` exact) while maintaining high
incremental ESS.  Two 1024-particle shared-$K=1$ continuations leave only
`131.5` and `138.1` collision bits; a lower-budget run reaches `127.7`.
Three shared-$K=2$ continuations across two transcripts leave `151.6--156.5`
bits.  An ideal $t=18,K=1$ stress transcript leaves only `128.6` bits after
fixed-$N$ correction.  These tail-sensitive results reject both $t=16,K=1$
and the simple $t=18,K=1$ replacement as comfortable 128-bit choices, and
strengthen $t=16,K=2$ as the optimized point.

Exact inverse-`Z` selection is not competitive under the current
`O(t^2 log t)` component-counter estimate: at `(16,16)`, exact `K=2` costs
2.46x and leaves about 108 bits, while occupancy `K=2` costs 1.35x and leaves
about 114 bits.

Increasing the partition capacity is more effective than increasing `K`
alone.  The balanced ideal point `(24,24),K=4` costs 2.49x, while the
power-of-two point `(16,32),K=4` costs 2.76x after accounting for its hidden
side-count sampler.  The latter avoids a non-power-of-two hash range and keeps
the expansion at `3N`.  Equal `(32,32),K=2` is the closest conservative
alternative: it costs 2.88x and has a slightly smaller leakage diagnostic,
but expands to `4N`.

Increasing the regular-noise weight at fixed load `d=t` is not automatically
favorable.  The nominal ISD exponent rises from 128 bits at `t=16` to 144 and
160 bits at `t=18,20`, but the number of leaked product lists also rises from
256 to 288 and 320.  The raw independent-list search diagnostic increases to
roughly 212 and 265 bits.  This does not establish an attack, but it means the
noise knob must be evaluated jointly with `d` and `K`.

## MPC interfaces

The implementation can be factored into three components.  Here `[x]` denotes
an XOR sharing, `[n]={0,...,n-1}`, and `b` is the number of independent
descriptors processed together.

```text
DistinctCount(
    [X] in ([d])^(b*t)
) -> [delta] in ({0,1}^ceil(log2(t)))^b
```

The rows for descriptor `j` are `X[j*t],...,X[(j+1)*t-1]`, and
`delta[j] = t - |{X[j*t],...,X[(j+1)*t-1]}|`.  The implementation performs
all pairwise equality tests in one batch, reduces the equality rows with
batched OR trees, and popcounts the resulting `t-1` duplicate indicators.
The current libOTe prototype processes 1,024 independent `t=16`, 5-bit
instances in one invocation.  Its in-process two-party test completes in
about 125 ms on the development machine, including synthetic base-OT
material setup; this is a smoke benchmark rather than a network measurement.

Let `Bins = ({0} x [d0]) union ({1} x [d1])`.  The asymmetric placement
sampler has interface

```text
AsymmetricPlacement(
    t in N,
    d0 in N,
    d1 in N
) -> [rho] in Bins^t
```

The output is injective.  A particular placement containing `a` left bins has
probability proportional to `1/(d0^a d1^(t-a))`.  Internally the protocol
samples the secret value `a`, a hidden row permutation, and one hidden bin
permutation per side.  None of these values is opened.

Finally, for candidate endpoint arrays
`[L_i] in ([d0])^t`, `[R_i] in ([d1])^t`, and their associated hidden
placements `[rho_i]`, define

```text
OccupancySelect(
    (([L_i], [R_i], [rho_i]))_(i in [K])
) -> ([L_J], [R_J], [rho_J])
```

It invokes `DistinctCount` on every left and right array, adds the two
bounded counts to obtain `D_i`, looks up an 80-bit public-random race key at
secret index `D_i`, and returns the candidate with minimum key.  The lookup
table is public, but `D_i`, `J`, every candidate graph, and every placement
remain secret.  Only the selected Reverse-Cuckoo hash descriptors are later
opened by the ordinary protocol.

The libOTe prototype implements this interface for any power-of-two `K`.
For `t=16`, `(d0,d1)=(16,32)`, and `K=4`, a five-descriptor in-process
two-party test completes in about 16 ms on the development machine, including
synthetic base-OT material setup.  The test exercises all 80 comparison bits,
checks the selected hidden record, and asserts the exact OT count below.  As
with the `DistinctCount` timing, this is a smoke benchmark rather than a
network measurement.

## Implementation status and next step

The occupancy score and fused selector are implemented.  The remaining
preprocessing work is:

1. generate `K=4` secret planted graphs in preprocessing;
2. feed their endpoints and retained placement records to `OccupancySelect`;
3. run the ordinary two bin-solvers once on the selected record;
4. integrate the path into the existing multi-expansion Reverse-Cuckoo code.

The next isolated primitive is `AsymmetricPlacement`.  It should first target
the unequal power-of-two partitions `(16,32)` and emit the endpoint arrays and
the complete retained record in the layout expected by `OccupancySelect`.

For the leading `t=16`, `(16,32)`, `K=4` point, the occupancy and selector OT
count is now:

```text
one proposal score:
    endpoint equalities       120*(7+8) = 1,800 OTs
    two duplicate OR trees    2*105 AND =   420 OTs
    two bounded popcounts        2*22 AND =    88 OTs
    add partition counts            4 AND =     8 OTs
                                             ---------
                                               2,316 OTs

one four-candidate selector (80-bit keys):
    four public-table lookups  4*15 wide products = 120 OTs
    three 80-bit comparisons   3*80 AND          = 480 OTs
    three wide winner muxes    3 wide products   =   6 OTs
                                             ---------
                                                 606 OTs
```

Across the 256 descriptors, proposal scoring costs `2,371,584` OTs and the
weighted selector costs `155,136` OTs.  Including the three extra planted
proposals gives `3,499,008` debiasing OTs and a total modeled cost of
`9,927,168` OTs over ten expansions, or `2.76x` the original `(16,16),K=1`
reference.  The public-table lookup and winner mux use wide OTs; their payload
depends on the race-key and retained-record widths, but not their OT count.
If the retained record is the 160-bit endpoint graph plus planted-side bits,
the selector carries 11,520 correlated-OT payload bits per descriptor, or
360 KiB across all 256 descriptors.  Even retaining a several-kilobit
permutation record keeps this part near one MiB; proposal generation and the
DPF payload remain separate terms.
