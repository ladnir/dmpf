#!/usr/bin/env python3
"""Importance-sample multi-overflow shallow-BFS obstructions exactly.

Matrices are sampled from the exact distribution conditioned on ``L_w=k``.
For a selected prefix of the initial overflow rows, determine the later,
unused candidate cells of their first-level victims.  Each selected row's
depth-two search in the initial Waterfall matching fails exactly when all its
distinct cells hit occupied bins.  We sample those cells uniformly from the
occupied bins and attach the exact likelihood ratio
``product_s (J_s/d)^a_s``, where ``a_s`` is the number of distinct conditioned
cells in partition ``s``.

For multiple selected rows this is their simultaneous obstruction in the
*initial* matching.  It is a strict subset of sequential-router failure,
because an earlier successful repair can change the matching seen by a later
row.  The weighted estimate for the selected event is unbiased for the
original uniform hash distribution.  It is statistical evidence, not a
deterministic upper bound.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from dataclasses import dataclass
from pathlib import Path

from waterfall_evictions import (
    Configuration,
    bounded_bfs_repair,
    maximum_matching_unplaced,
    waterfall,
)
from waterfall_evictions_conditional import (
    backward_event_probabilities,
    occupancy_transitions,
    sample_conditioned_matrix,
)


@dataclass
class WeightedMoments:
    count: int = 0
    total: float = 0.0
    total_square: float = 0.0
    hits: int = 0

    def add(self, weight: float, hit: bool) -> None:
        value = weight if hit else 0.0
        self.count += 1
        self.total += value
        self.total_square += value * value
        self.hits += hit

    @property
    def mean(self) -> float:
        return self.total / self.count

    @property
    def standard_error(self) -> float:
        if self.count < 2:
            return math.nan
        centered = self.total_square - self.total * self.total / self.count
        return math.sqrt(max(0.0, centered / (self.count - 1) / self.count))


def force_initial_depth2_obstructions(
    candidates: list[int],
    configuration: Configuration,
    rng: random.Random,
    forced_roots: int,
) -> tuple[float, list[int], list[int], list[int], tuple[int, ...]]:
    """Condition a prefix of roots to fail in the initial matching."""
    t, w, d = configuration.t, configuration.w, configuration.d
    owner, placement, overflow = waterfall(candidates, configuration)
    occupied_labels: list[list[int]] = [[] for _ in range(w)]
    occupied_counts = [0] * w
    for partition in range(w):
        offset = partition * d
        labels = [
            label
            for label in range(d)
            if owner[offset + label] != -1
        ]
        occupied_labels[partition] = labels
        occupied_counts[partition] = len(labels)

    forced_cells: set[int] = set()
    if not 1 <= forced_roots <= len(overflow):
        raise ValueError("forced_roots must select a nonempty overflow prefix")
    for root in overflow[:forced_roots]:
        base = root * w
        for partition in range(w):
            slot = partition * d + candidates[base + partition]
            victim = owner[slot]
            if victim == -1:
                raise AssertionError("Waterfall overflow candidate must be occupied")
            victim_partition = placement[victim] // d
            for later_partition in range(victim_partition + 1, w):
                forced_cells.add(victim * w + later_partition)

    exponents = [0] * w
    for cell in forced_cells:
        partition = cell % w
        labels = occupied_labels[partition]
        if not labels:
            raise AssertionError("forced partition has no occupied bins")
        candidates[cell] = labels[rng.randrange(len(labels))]
        exponents[partition] += 1

    weight = 1.0
    for partition, exponent in enumerate(exponents):
        weight *= (occupied_counts[partition] / d) ** exponent

    conditioned_owner, conditioned_placement, conditioned_overflow = waterfall(
        candidates, configuration
    )
    if conditioned_overflow != overflow:
        raise AssertionError("conditioning changed the Waterfall placement")
    for root in conditioned_overflow[:forced_roots]:
        if bounded_bfs_repair(
            candidates,
            configuration,
            conditioned_owner,
            conditioned_placement,
            [root],
            2,
        ) == 0:
            raise AssertionError("selected root does not defeat depth-two BFS")
    return (
        weight,
        conditioned_owner,
        conditioned_placement,
        conditioned_overflow,
        tuple(exponents),
    )


def log2_or_negative_infinity(value: float) -> float:
    return -math.inf if value == 0.0 else math.log2(value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument("--w", type=int, default=4)
    parser.add_argument("--d", type=int, default=16)
    parser.add_argument("--stash", type=int, default=0)
    parser.add_argument("--overflow", type=int, default=2)
    parser.add_argument(
        "--forced-roots",
        type=int,
        default=1,
        help="initial overflow-row prefix to obstruct simultaneously",
    )
    parser.add_argument("--depth", type=int, nargs="+", default=[3, 4, 5, 16])
    parser.add_argument("--samples", type=int, default=500_000)
    parser.add_argument("--seed", type=int, default=0x2B5F)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    if args.overflow < 1 or args.samples < 2 or not 0 <= args.stash < args.overflow:
        raise SystemExit("require overflow>=1, samples>=2, and 0<=stash<overflow")

    sampling_configuration = Configuration(args.t, args.w, args.d, 0)
    evaluation_configuration = Configuration(args.t, args.w, args.d, args.stash)
    transitions = occupancy_transitions(args.t, args.d)
    backward = backward_event_probabilities(
        sampling_configuration, transitions, terminal_live=args.overflow
    )
    exact_stratum = backward[0][args.t]
    if exact_stratum == 0.0:
        raise SystemExit("requested overflow stratum has zero probability")

    depths = sorted(set(args.depth))
    rng = random.Random(args.seed)
    moments = {depth: WeightedMoments() for depth in depths}
    no_matching = WeightedMoments()
    obstruction = WeightedMoments()
    exponent_histogram: dict[tuple[int, ...], int] = {}

    for _ in range(args.samples):
        candidates, _ = sample_conditioned_matrix(
            sampling_configuration,
            transitions,
            backward,
            rng,
            terminal_live=args.overflow,
        )
        weight, owner, placement, overflow, exponents = (
            force_initial_depth2_obstructions(
                candidates, sampling_configuration, rng, args.forced_roots
            )
        )
        exponent_histogram[exponents] = exponent_histogram.get(exponents, 0) + 1
        obstruction.add(weight, True)
        matching_failure = (
            maximum_matching_unplaced(candidates, evaluation_configuration)
            > args.stash
        )
        no_matching.add(weight, matching_failure)
        for depth in depths:
            failure = bounded_bfs_repair(
                candidates,
                evaluation_configuration,
                owner,
                placement,
                overflow,
                depth,
            ) > args.stash
            moments[depth].add(weight, failure)

    rows: list[dict[str, object]] = []
    obstruction_label = f"initial_{args.forced_roots}_root_obstruction"
    for label, accumulator in [
        (obstruction_label, obstruction),
        ("no_matching", no_matching),
    ]:
        estimate = exact_stratum * accumulator.mean
        standard_error = exact_stratum * accumulator.standard_error
        rows.append(
            {
                "event": label,
                "hits": accumulator.hits,
                "samples": accumulator.count,
                "estimate": estimate,
                "estimate_log2": log2_or_negative_infinity(estimate),
                "standard_error": standard_error,
            }
        )
    for depth in depths:
        accumulator = moments[depth]
        estimate = exact_stratum * accumulator.mean
        standard_error = exact_stratum * accumulator.standard_error
        rows.append(
            {
                "event": f"bfs_depth_{depth}",
                "hits": accumulator.hits,
                "samples": accumulator.count,
                "estimate": estimate,
                "estimate_log2": log2_or_negative_infinity(estimate),
                "standard_error": standard_error,
            }
        )

    print(
        f"Pr[L_w={args.overflow}]={exact_stratum:.12e} "
        f"log2={math.log2(exact_stratum):.6f}"
    )
    for row in rows:
        print(
            f"{row['event']} hits={row['hits']}/{row['samples']} "
            f"estimate={row['estimate']:.12e} "
            f"log2={row['estimate_log2']:.6f} se={row['standard_error']:.3e}"
        )
    print("forced-cell exponents:")
    for exponents, count in sorted(exponent_histogram.items()):
        print(f"  {exponents}: {count}")

    if args.csv:
        with args.csv.open("w", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
