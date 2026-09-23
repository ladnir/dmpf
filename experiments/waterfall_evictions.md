# Bounded eviction repair for Waterfall Cuckoo

Cost-accounting update (2026-09-07): setup totals and cost comparisons below
are historical. Use `waterfall_mpc_cost.py` for the corrected DPF correlation
costs, including shared correction selections, sparse activity masks, and
cached sign correlations. The simulation outcomes below are unchanged.

## Simulated router

The simulator starts with the fixed-priority waterfall from
`WaterfallCuckoo.tex`.  It then maintains one current unplaced row and a global
budget of `e` eviction steps.  A step behaves as follows.

1. If any candidate bin of the current row is free, place the row in its
   earliest free partition and continue with the next overflow row.
2. Otherwise, evict one occupant and continue with the evicted row.
3. On the first eviction in a chain, the `biased` policy chooses uniformly from
   partitions `[0,w-1)`.  Later evictions may choose every partition.
4. The bin just traversed is excluded from the next choice, preventing an
   immediate two-step reversal.

The `uniform` control allows the last partition on the first eviction.  The
`oracle` control computes a maximum bipartite matching and therefore gives a
lower bound on the number of unplaced rows achievable by any router.

The global budget is the intended low-cost circuit.  A diagnostic `per-row`
schedule is also implemented, but it schedules `e*t` steps and is much more
expensive.

## Reproduction

```powershell
python experiments\waterfall_evictions.py `
  --d 6 8 10 12 16 `
  --c 0 1 2 `
  --evictions 0 2 4 8 16 32 `
  --schedules global `
  --trials 50000 `
  --batch 256 `
  --seed 789248 `
  --csv experiments\waterfall_evictions_results.csv
```

All policies at a parameter point use the same 50,000 candidate matrices.  The
CSV includes Wilson confidence intervals, the exact unmodified-waterfall
failure probability, and structural MPC-cost columns.

## Representative results

The table reports failures out of 50,000 trials.  These deliberately compact
points put the repaired tail in an observable range; they are not proposed
40-bit parameters.

| stash `c` | partition `d` | columns `3d+c` | global steps `e` | waterfall | biased | uniform | matching oracle |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 12 | 36 | 8 | 3,489 | 8 | 14 | 0 |
| 0 | 16 | 48 | 4 | 801 | 4 | 20 | 0 |
| 1 | 10 | 31 | 8 | 655 | 3 | 4 | 0 |
| 1 | 12 | 37 | 4 | 137 | 3 | 6 | 0 |
| 2 | 8 | 26 | 8 | 295 | 3 | 5 | 0 |
| 2 | 10 | 32 | 4 | 25 | 0 | 1 | 0 |

One eviction cannot reduce the number of unplaced rows: when all candidates
are occupied, it only exchanges the current row for an occupant.  At least two
steps are required to discover and fill a free bin.  The biased first eviction
is most useful at small `e`; its advantage narrows as longer walks forget their
starting tier.

## Cost coordinates

The simulator intentionally does not collapse MPC cost into one number.  At
`t=16,w=3`, the common and variable components are:

- one direct polynomial-hash evaluation, with address powers reused across all
  three partitions;
- `360 = 3*binom(16,2)` initial waterfall collision tests;
- `48d` main-placement one-hot bits;
- `3d+c` output columns;
- `(3+c)N` local expansion work;
- `48e` repair owner-scan units for the global schedule.

One repair scan unit is one candidate-versus-owner test.  A rough binary-index
implementation uses about `2*ceil(log2(d))` bit operations per scan unit, plus
token multiplexing, state updates, and the small random-choice circuit.  These
columns are suitable for fitting weights from an implementation benchmark.

For orientation, replacing the current `(w,c,d)=(3,2,74)` point by the
proof-motivated power-of-two point `(3,2,32)` would reduce the main one-hot
materialization from `3,552` to `1,536` bits and the output width from `224` to
`98`.  A four-step repair adds `192` scan units.  Whether that trade is positive
depends on the relative cost of secret index comparisons, output
materialization/scatter, and the shared direct-hash evaluation.

## Statistical limitation

Ordinary Monte Carlo cannot establish a `2^-40` rejection probability.  With
50,000 trials, even zero failures has a 95% Wilson upper bound near
`7.7e-5`; after a batch of 256 descriptors this is nowhere close to a 40-bit
claim.  The sweep only establishes the shape of the repair curve and identifies
which policies merit rigorous analysis.

This sweep motivated a conditional rare-event experiment to locate the
finite-`e` frontier without extrapolating ordinary Monte Carlo.

## Exact-overflow conditional experiment

The follow-up experiment removes the initial rare-event bottleneck exactly.
The waterfall live count `L_s` is a Markov chain with the occupancy transition
from the paper.  Backward probabilities therefore permit exact sampling from
the uniform candidate-matrix distribution conditioned on `L_w>c`.  For each
conditioned matrix, the experiment enumerates every random eviction choice by
in-place backtracking.  The final estimator is

```text
exact Pr[L_w>c] * estimated Pr[repair fails | L_w>c].
```

This is implemented in `waterfall_evictions_conditional.py`.  At the observable
validation point `(d,c,e)=(8,2,4)`, it estimated `4.792e-4`; the independent
direct sweep observed `24/50000 = 4.8e-4`.

For `(t,w,c)=(16,3,2)`, the compactness/repair frontier is:

| `d` | `e` | columns | failure `log2` per descriptor | batch security bits | repair scan units |
|---:|---:|---:|---:|---:|---:|
| 32 | 12 | 98 | -49.02 | 41.02 | 576 |
| 40 | 8 | 122 | -49.88 | 41.88 | 384 |
| 48 | 4 | 146 | -49.50 | 41.50 | 192 |
| 52 | 2 | 158 | -48.16 | 40.16 | 48 |
| 56 | 2 | 170 | -49.65 | 41.65 | 48 |
| 64 | 2 | 194 | -52.32 | 44.32 | 48 |

For general `e`, the last column uses the conservative `e*w*t` scan count.  At
`e=2`, the circuit needs only `w*t=48` scans: one owner lookup for the initial
overflow row and `w-1` alternative-bin tests for its victim.  The estimates use
between 150,000 and 300,000 conditioned matrices.  They are statistically
precise enough to identify the frontier for `e>2`.  The `e=2` rows in the table
have the exact analysis below and are not estimates.

### Exact two-step failure

Let `(L_0,...,L_w)` be a waterfall live-count path and put
`J_r=L_r-L_{r+1}`, the number of occupied bins created in partition `r`.  If
`L_w>=c+2`, a two-step repair can remove at most one overflow row and therefore
fails.  If `L_w=c+1`, choose the first eviction partition uniformly from
`s in [0,w-1)`.  Its victim lost in every earlier partition, so all of its
earlier candidates are known occupied.  Its later candidates were never used
by the waterfall and remain independent uniform values.  Consequently its
conditional failure probability is exactly

```text
1/(w-1) * sum_{s=0}^{w-2} product_{r=s+1}^{w-1} J_r/d.
```

Multiplying this expression by the exact probability of each live-count path,
summing paths ending at `c+1`, and adding `Pr[L_w>=c+2]` gives the exact
two-step failure probability.  The implementation is
`exact_two_step_biased_failure`; `waterfall_evictions_two_step.py` produces the
parameter table.

For `(t,w,c,d)=(16,3,2,64)`, the exact failure probability is
`1.7769500875e-16 = 2^-52.3214` per descriptor.  A union bound over 256
descriptors is therefore `2^-44.3214`, leaving more than four bits of margin
over the 40-bit target.  The conditional experiment independently estimated
`2^-52.3240`.

The natural implementation point is `(d,e)=(64,2)`.  It is the lowest-cost
power-of-two point in the measured frontier and has roughly four bits of
statistical margin after the 256-descriptor union bound.  The two-step repair
also simplifies: choose partition 0 or 1, replace its occupant by the current
overflow row, and test whether the victim has a free candidate.  If so, place
the victim there; otherwise the repair attempt fails.  A final neutral eviction
is unnecessary.

Relative to the power-of-two unmodified point `d=128`, this reduces the output
width from `386` to `194` columns and the main one-hot materialization from
`6144` to `3072` bits.  The optimized repair contributes `48` owner-scan units,
or roughly `576` binary-index bit operations under the simulator's coarse
proxy.  The resulting proxy is `3648`, about 41% below `6144`, before counting
the additional benefit of the smaller scatter and position map.

The remaining work is circuit specification and implementation benchmarking,
not a probabilistic conjecture for this router: the two-step rejection
probability at `d=64` now has an exact recurrence calculation.

## Minimizing local expansion `w+c`

The literal minimum is `(w,c)=(2,0)`, giving `2N` local expansion, but it is
not a plausible MPC point.  Exact two-step repair first reaches 40 batch bits
at power-of-two `d=65536`, producing 131,072 output columns.  Even the
full-matching Hall-witness proof uses `d=32768` and 65,536 columns.

At `w+c=3`, the relevant choices are:

| `(w,c)` | two-step power-of-two `d` | two-step columns | full-matching proof `d` | matching columns |
|---:|---:|---:|---:|---:|
| `(2,1)` | 2048 | 4097 | 1024 | 2049 |
| `(3,0)` | 1024 | 3072 | 128 | 384 |

Thus `(3,0)` is the only serious three-expansion target, but it needs a much
stronger router than two-step repair.  Exact-overflow conditional sampling at
`(t,w,c,d)=(16,3,0,128)` gives:

| global eviction budget `e` | failure `log2` | batch security bits |
|---:|---:|---:|
| 12 | -42.72 | 34.72 |
| 16 | -46.72 | 38.72 |
| 18 | -48.72 | 40.72 |
| 20 | -50.72 | 42.72 |

The `e=18` point therefore appears sufficient statistically and reduces local
evaluation from `5N` to `3N`, but its MPC routing is far more expensive:
384 output columns, 6,144 main one-hot bits, and 864 generic repair scan units.
Its coarse repair index proxy alone is 12,096 bit operations.  This makes it a
reasonable choice only when repeated local evaluation dominates preprocessing.

There is also a promising middle ground at `w+c=4`: `(w,c)=(4,0)` has a
Hall-witness proof at `d=16`, only 64 columns, but reaching the matching floor
requires a substantially stronger routing circuit than the current biased
walk.  It adds routing work for a fourth partition but no additional field
powers, trades one full-domain DPF away, and may
be the best next configuration to investigate if `3N` is not mandatory.

### The `(w,c,d)=(4,0,16)` routing experiment

The biased random walk converges too slowly on this dense graph:

| global eviction budget `e` | failure `log2` | batch security bits |
|---:|---:|---:|
| 4 | -25.23 | 17.23 |
| 8 | -31.22 | 23.22 |
| 14 | -37.50 | 29.50 |
| 20 | -42.28 | 34.28 |

At larger `e`, a small number of cyclic configurations dominate the estimator,
so extending the same walk is not attractive.

A deterministic breadth-first augmenting-path control is dramatically better:
among 300,000 matrices sampled exactly conditional on waterfall overflow,
depth-2 BFS found a path in every sample.  This observation is not itself a
tail bound.  For a live-count path ending in exactly one overflow row, let
`J_r=L_r-L_{r+1}`.  Depth-2 BFS fails only if every later candidate of every
possible first victim is occupied, which has exact path-conditional
probability

```text
product_{r=1}^{w-1} (J_r/d)^r.
```

Summing this over paths for `(4,0,16)` gives `2^-36.9156` per descriptor, so
depth two is not sufficient despite producing no failures in the sample.

A rigorous fallback is to support at most three overflow rows and run complete
bounded BFS for each.  The exact probability of at least four waterfall
overflows is `2^-53.2361`.  Combining this with the full-matching Hall-witness
bound `2^-51.8263` gives total failure below `2^-51.3654`, or 43.37 bits after
the 256-descriptor union bound.

The fallback is structurally compact but MPC-heavy.  A direct fixed circuit
with three roots, at most 16 row levels, and 64 row-candidate edges per level
has roughly 3,072 edge-scan units, in addition to the direct hash circuit.  This is far
more preprocessing than `(3,2,64,e=2)`, even though it reduces the position map
to 64 columns and local evaluation to `4N`.  A shallow-BFS proof, ideally depth
three or four, is required for `(4,0,16)` to become the overall winner.

#### Targeted shallow-BFS experiment

The ordinary conditioned experiment still almost never encounters a depth-2
BFS obstruction.  A second sampler therefore conditions exactly on the
dominant one-overflow depth-2 obstruction.  It samples live-count paths
proportional to their exact obstruction mass and then programs the unused
victim candidates as uniform occupied bins.  The total probability of the
target event is the exact `2^-36.9156` value above.

Among 300,000 such targeted matrices:

| BFS depth | failures | implied one-overflow contribution |
|---:|---:|---:|
| 3 | 18 | `2^-50.94` |
| 4 | 12 | `2^-51.53` |
| 8 | 12 | `2^-51.53` |
| 16 | 12 | `2^-51.53` |

Depth 16 is maximum matching for these instances.  Thus the 12 persistent
failures are genuine no-matching configurations, while the other six have a
shortest augmenting path of depth four.  No targeted matrix required a path
longer than four while still admitting a matching.

The other overflow strata agree with this picture:

- conditioned on exactly two overflow rows, depth 2 failed on 33 of 500,000
  matrices and depth 3 repaired all 33;
- conditioned on exactly three overflow rows, depth 2 failed on 153 of 200,000
  matrices and depth 3 repaired all 153;
- four or more overflow rows already have exact probability `2^-53.2361`.

These experiments suggested that three overflow-row slots and BFS depth four
might achieve the full-matching bound.  A fixed depth-4 circuit costs
approximately `3*4*(w*t)=768` edge-scan units, rather than 3,072 for three full
depth-16 searches.  Together with 1,024 main one-hot bits and the direct hash
circuit, it
is still more MPC setup than `(3,2,64,e=2)`, but it has only 64 output columns
and `4N` local expansion.

#### A depth-five counterexample

Depth four is not universally sufficient.  Exhaustive row-subset minimization
of a random counterexample gives the following 11-row matrix at `w=4,d=4`:

```text
(0,2,0,3)  (1,0,0,1)  (0,2,3,0)  (3,3,0,1)
(3,2,3,1)  (0,3,3,1)  (3,3,0,1)  (3,2,0,1)
(3,3,3,1)  (0,2,3,0)  (1,3,0,3)
```

Waterfall leaves row 8 unplaced, yet the candidate graph has a perfect
matching.  Its shortest augmenting path is

```text
row 8  -> partition 1, bin 3   (occupied by row 5)
row 5  -> partition 0, bin 0   (occupied by row 0)
row 0  -> partition 3, bin 3   (occupied by row 10)
row 10 -> partition 0, bin 1   (occupied by row 1)
row 1  -> partition 1, bin 0   (free)
```

Thus depth four fails and depth five succeeds.  The matrix embeds into the
actual `t=16,w=4,d=16` parameters by retaining these labels and adding five
isolated rows on fresh labels 4 through 8.  The same overflow row and shortest
path remain.  No proper row subset of this particular 11-row matrix is still a
depth-five witness; this is subset minimality, not a proof of global minimality.

This counterexample disproves only a structural impossibility claim.  It is
extremely collision-heavy: among its 11 rows, the four partitions use only
`(3,3,2,3)` distinct labels.  It therefore does not show that depth-five
witnesses have noticeable probability at `d=16`, and it does not contradict
the targeted experiment, which simply did not sample one.

The remaining mathematical task is to upper-bound the probability of a
*matchable* Waterfall instance with at most three overflow rows but no
depth-four augmenting repair.  Protocol failure is then the sum of this new
deep-witness event, the no-perfect-matching event, and `L_w>=4`.  The targeted
sampler is implemented in `waterfall_bfs_targeted.py`; reproducible witness
search, minimization, path extraction, and embedding are implemented in
`waterfall_bfs_witness.py`.

#### Collision-minimized deep witnesses

To distinguish essential structure from accidental collisions in the small
`d=4` search space, embed each witness at `t=d=16` and mutate individual
candidate labels while preserving all of the following properties: Waterfall
has between one and three overflow rows, the graph has a perfect matching,
and BFS at the rejected depth still fails.  The optimization score is

```text
collision excess = sum_s (16 - number of distinct labels in partition s).
```

This score counts repetitions in the complete concrete matrix.  Entries that
are irrelevant to the routing can mutate to fresh labels, so the local search
removes most accidental equalities.  It is a structural diagnostic, not by
itself a probability or a probability bound.

Deterministic searches and 750,000--1,000,000 mutation steps per witness give:

| minimum successful BFS depth | rejected depth | collision excess | log2 mass of complete equality pattern |
|---:|---:|---:|---:|
| 5 | 4 | 22 | -111.80 |
| 6 | 5 | 25 | -120.00 |
| 7 | 6 | 28 | -128.39 |
| 8 | 7 | 31 | -137.37 |

Thus neither depth four nor depths five through seven are universal cutoffs.
More importantly, the best witnesses found obey the exact progression
`3*depth+7`.  This has a natural BFS interpretation for `w=4`: after following
the currently occupied placement edge, extending a blocked chain by one row
requires its other three candidate edges to remain inside the occupied
frontier.  The experiment therefore suggests that each extra forbidden BFS
level costs three additional collision constraints.

The multiplicity of possible occupied targets means that three collision
constraints do **not** automatically contribute exactly `3 log2(d)=12`
security bits.  Even the last column is only the probability of one complete
equality pattern and overconstrains routing-irrelevant entries; a proof must
sum over the possible patterns, paths, and occupied targets.  The regular `+3`
pattern nevertheless gives the right prospective lemma: bound
the mass of depth-`ell` blocked alternating paths, conditioned on the exact
Waterfall live-count path, and union-bound over at most three overflow rows.
The `--optimize-collisions` mode in `waterfall_bfs_witness.py` reproduces this
search.

### A proved shallow-router point

There is a simple rigorous alternative to proving the full depth-four tail.
Give the router stash capacity `c`, accept immediately when `L_w<=c`, and,
when `L_w=c+1`, run depth-two BFS for one fixed overflow row.  Conservatively
reject every proposal with `L_w>=c+2`.  Only one augmenting path is attempted,
so there is no interaction between multiple repairs.

The circuit is particularly small.  Let `x` be the first overflow row.  Its
candidate in partition `p` is occupied by a distinct row `v_p`.  In fixed
lexicographic order, inspect the `w` candidates of `v_0`, then those of
`v_1`, and so on, and select the first free candidate.  If `v_p` has a free
candidate `b`, move `v_p` to `b` and place `x` in the bin vacated by `v_p`.
Otherwise set the rejection bit.  This is one layer of `w` owner lookups and
at most `w^2` occupied/free tests, followed by two one-hot placement updates;
it is not a general BFS or matching circuit.  After a successful update,
exactly `c` rows remain and stable-compaction places them in the stash.

Fix a live-count path `(L_0,...,L_w)` and write
`J_s=L_s-L_{s+1}` for the number of rows placed in partition `s`.  The chosen
overflow row has one occupied candidate in every partition.  Its victim in
partition `p` was live in every earlier partition, so its candidates through
partition `p` are already known to be occupied.  Its candidates in later
partitions were never inspected by Waterfall and remain mutually independent
uniform values in `[d]`.  Consequently depth-two BFS fails with the exact
path-conditional probability

```text
product_{q=1}^{w-1} (J_q/d)^q.
```

Indeed, partition `q` receives one still-unseen candidate from each of the
victims placed in partitions `0,...,q-1`, giving exactly `q` independent
tests against its `J_q` occupied bins.  The victims are distinct because one
placed row owns only one physical bin.

Combining this observation with the exact Stirling transition recurrence gives
the following rejection bound for the proved router:

```text
Pr[reject]
 <= Pr[L_w >= c+2]
    + sum over paths ending at L_w=c+1 of
        Pr[path] * product_{q=1}^{w-1} (J_q/d)^q.
```

For `t=16,w=4,d=16,c=2`, exact rational arithmetic gives

```text
Pr[L_4 >= 4]                         = 2^-53.236062418
Pr[L_4 = 3 and depth-two BFS fails] = 2^-60.665137797
Pr[reject]                           = 2^-53.227715176
256 * Pr[reject]                     = 2^-45.227715176.
```

Thus `(w,d,c)=(4,16,2)` has a complete, assumption-free routing proof with
only one depth-two repair circuit.  It has `m=66` position-map columns and
local expansion `(w+c)N=6N`.  In the ideal arbitrary-range model, `d=14`
already gives 40.679 batch bits and `m=58`; the paper's binary hash output
rounds this to the proved `d=16` point.

For comparison, retaining `w+c=5` with `c=1` and using the same one-repair
proof requires `d=22` in the arbitrary-range model, hence `d=32` with the
paper's binary hash output.  It uses less local expansion but substantially
more main-bin state.  The exact proof and parameter sweep are implemented in
`waterfall_shallow_proof.py`.  Its result was independently checked by the
exact conditional importance sampler: 50,000 samples in the `L_4=3` stratum
estimated the obstruction at `2^-60.6765`, consistent with the exact
`2^-60.6651` value.

Some relevant power-of-two points are:

| `w` | `c` | `d` | columns `wd+c` | local factor `w+c` | 256-descriptor bits |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 64 | 194 | 5 | 51.8300 |
| 3 | 3 | 32 | 99 | 6 | 47.8553 |
| 4 | 0 | 64 | 256 | 4 | 49.4276 |
| 4 | 1 | 32 | 129 | 5 | 51.2110 |
| 4 | 2 | 16 | 66 | 6 | 45.2277 |
| 5 | 0 | 16 | 80 | 5 | 40.0183 |
| 5 | 1 | 16 | 81 | 6 | 57.2351 |

The `(4,2,16)` point is the strongest MPC candidate in this proved family:
relative to `(5,0,16)`, it removes one complete main partition and 14 position
columns at the cost of two cheap stash DPFs and one extra unit of local
expansion.  It also leaves more than five bits of margin beyond the 40-bit
batch target, whereas `(5,0,16)` has essentially no budget for other errors.
