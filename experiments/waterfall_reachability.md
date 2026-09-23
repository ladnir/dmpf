# Complete Waterfall repair

Cost-accounting update (2026-09-07): historical setup totals in this note
are superseded by `waterfall_mpc_cost.py` and the paper's current tables.
The model now charges shared correction selections, sparse activity masks,
and cached sign correlations explicitly. The main 4N / 3N totals are
39,052 / 55,358 OTs and 942,888 / 3,623,410 payload bits for G=Z_(2^64).
The algorithms and correctness bounds below are unchanged.

This note records the complete-repair algorithm implemented by
`WaterfallReachability.h` and analyzed in `WaterfallCuckooAppendix.tex`.
The main profiles are

```text
t=16, d=(16,16,16,16), c=0, r=3
t=16, d=(128,128,64),   c=0, r=2
```

They use the same deterministic repair circuit. The exact occupancy and Hall
bounds give 43.37 and 41.13 batch-correctness bits, respectively, for a batch
of 256 independently routed states.

## Matching representation

The generator stores the current main-bin matching as row-to-partition bits

```text
M[i,s] = 1 iff row i occupies its candidate in partition s.
```

Each row of `M` has weight zero or one. An all-zero row is unmatched. The
candidate decoders `D[i,s]` are cached from the initial waterfall scan. At the
start of each repair round, the circuit materializes an occupancy bit and an
owner record for every bin from `M` and `D`. A secret-index mux can therefore
read the status and owner of any candidate.

## Multi-source reachability search

Every unmatched row is an initial root. Initialize `visited` with all roots
and `expanded` with zero. For exactly `t` public steps:

1. Select the first visited but unexpanded row in public row order.
2. Mark it expanded and inspect its candidates in public partition order.
3. If a candidate is free and no target has been selected, record it as the
   terminal candidate.
4. For every occupied candidate whose owner is unseen, mark the owner visited
   and record the current row and partition as its parent.

After a terminal candidate is selected, all remaining forward steps are
gated no-ops. If no terminal is selected after `t` steps, every row reachable
from any unmatched root has been expanded. Otherwise the parent records form
a forest rooted at the unmatched rows, and the selected terminal has a unique
parent chain to one root. The circuit follows that chain for `t` fixed steps
and flips the corresponding assignments. Unused backtracking steps are again
gated no-ops.

This traversal is deliberately called a *multi-source reachability search*.
The public row-order priority does not process distance layers and therefore
does not promise a shortest path.

## Correctness reduction

Let `L` be the number of rows left unmatched by the initial waterfall stage.
If `L <= c+r` and the candidate graph has a size-`t` matching, then every
round in which more than `c` rows remain unmatched has an augmenting path.
The complete multi-source search finds one and increases the main matching by
one. Therefore

```text
Pr[reject] <= Pr[L > c+r] + Pr[no size-t matching].
```

The first term is computed by the exact occupancy recurrence. The second is
bounded by the Hall-witness sum. These are the bounds reported in the paper.

## Reference implementations

`waterfall_reachability.py` contains two references:

- `augment_first_reachable_root` is a simple exact-distance oracle that finds
  a canonical shortest augmenting path.
- `augment_from_all_roots_dfs` emulates the optimized multi-source traversal
  used by the C++ circuit. The historical function name is retained for test
  compatibility; the traversal is not presented as a depth-first search in
  the paper.

The exhaustive tests compare their success behavior on small candidate
graphs and verify that bounded complete repair attains the expected matching
criterion. `waterfall_mpc_cost.py` contains the OT and payload model used for
the concrete table.
