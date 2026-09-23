#!/usr/bin/env python3
"""Target the single-overflow depth-2 BFS obstruction exactly.

For a waterfall live-count path ending at L_w=1, let J_r be the number of
occupied bins created in partition r.  A depth-2 augmenting-path search fails
when every later candidate of every possible first victim is occupied.  Its
path-conditional probability is product_{r=1}^{w-1}(J_r/d)^r.

This program samples live-count paths proportional to their exact probability
times that obstruction probability, samples a uniform matrix with the chosen
path, and conditions the unused victim candidates to hit uniform occupied
bins.  It therefore samples exactly from the one-overflow depth-2 failure
event, allowing deeper deterministic BFS to be measured without first waiting
for a roughly 2^-37 event.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import random
from dataclasses import dataclass
from pathlib import Path

from waterfall_evictions import Configuration, bounded_bfs_repair, waterfall
from waterfall_evictions_conditional import (
    occupancy_transitions,
    uniform_surjection_labels,
)


@dataclass
class PathEntry:
    path: tuple[int, ...]
    probability: float
    obstruction_probability: float
    cumulative_weight: float = 0.0


def enumerate_target_paths(
    configuration: Configuration,
    transitions: list[list[tuple[int, float]]],
) -> tuple[list[PathEntry], float]:
    t, w, d = configuration.t, configuration.w, configuration.d
    path = [t]
    entries: list[PathEntry] = []

    def recurse(partition: int, probability: float) -> None:
        if partition == w:
            if path[-1] != 1:
                return
            occupied = [path[index] - path[index + 1] for index in range(w)]
            obstruction = 1.0
            for later_partition in range(1, w):
                obstruction *= (occupied[later_partition] / d) ** later_partition
            if obstruction:
                entries.append(PathEntry(tuple(path), probability, obstruction))
            return
        for next_live, transition_probability in transitions[path[-1]]:
            path.append(next_live)
            recurse(partition + 1, probability * transition_probability)
            path.pop()

    recurse(0, 1.0)
    cumulative = 0.0
    for entry in entries:
        cumulative += entry.probability * entry.obstruction_probability
        entry.cumulative_weight = cumulative
    return entries, cumulative


def sample_matrix_for_path(
    configuration: Configuration,
    path: tuple[int, ...],
    rng: random.Random,
) -> list[int]:
    t, w, d = configuration.t, configuration.w, configuration.d
    candidates = [0] * (t * w)
    live_rows = list(range(t))
    for partition in range(w):
        next_live_count = path[partition + 1]
        occupied = len(live_rows) - next_live_count
        for row in range(t):
            candidates[row * w + partition] = rng.randrange(d)
        bins = rng.sample(range(d), occupied)
        labels = uniform_surjection_labels(len(live_rows), occupied, rng)
        seen: set[int] = set()
        next_live_rows: list[int] = []
        for row, label in zip(live_rows, labels):
            value = bins[label]
            candidates[row * w + partition] = value
            if value in seen:
                next_live_rows.append(row)
            else:
                seen.add(value)
        if len(next_live_rows) != next_live_count:
            raise AssertionError("sampled path has the wrong live count")
        live_rows = next_live_rows
    return candidates


def force_depth2_failure(
    candidates: list[int],
    configuration: Configuration,
    rng: random.Random,
) -> tuple[list[int], list[int], list[int]]:
    t, w, d = configuration.t, configuration.w, configuration.d
    owner, placement, unplaced = waterfall(candidates, configuration)
    if len(unplaced) != 1:
        raise AssertionError("targeted sampler requires exactly one overflow row")
    starting_row = unplaced[0]
    base = starting_row * w
    occupied_bins = [
        [slot - partition * d for slot in range(partition * d, (partition + 1) * d) if owner[slot] != -1]
        for partition in range(w)
    ]

    victims: set[int] = set()
    for first_partition in range(w - 1):
        slot = first_partition * d + candidates[base + first_partition]
        victim = owner[slot]
        if victim == -1 or victim in victims:
            raise AssertionError("first-level victims must be distinct and placed")
        victims.add(victim)
        for later_partition in range(first_partition + 1, w):
            choices = occupied_bins[later_partition]
            if not choices:
                raise AssertionError("later waterfall partition has no occupied bin")
            candidates[victim * w + later_partition] = choices[
                rng.randrange(len(choices))
            ]

    conditioned_owner, conditioned_placement, conditioned_unplaced = waterfall(
        candidates, configuration
    )
    if conditioned_unplaced != unplaced:
        raise AssertionError("conditioning changed the waterfall path")
    if bounded_bfs_repair(
        candidates,
        configuration,
        conditioned_owner,
        conditioned_placement,
        conditioned_unplaced,
        2,
    ) == 0:
        raise AssertionError("conditioned matrix does not defeat depth-2 BFS")
    return conditioned_owner, conditioned_placement, conditioned_unplaced


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument("--w", type=int, default=4)
    parser.add_argument("--d", type=int, default=16)
    parser.add_argument("--depth", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--samples", type=int, default=500_000)
    parser.add_argument("--seed", type=int, default=0xBF52)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    if args.t < 1 or args.w < 2 or args.d < 1 or args.samples < 2:
        raise SystemExit("require t>=1, w>=2, d>=1, and samples>=2")
    if any(depth < 3 for depth in args.depth):
        raise SystemExit("every --depth must be at least three")

    configuration = Configuration(args.t, args.w, args.d, 0)
    transitions = occupancy_transitions(args.t, args.d)
    paths, exact_depth2_single = enumerate_target_paths(configuration, transitions)
    cumulative = [entry.cumulative_weight for entry in paths]
    rng = random.Random(args.seed)
    depths = sorted(set(args.depth))
    failures = {depth: 0 for depth in depths}

    for _ in range(args.samples):
        draw = rng.random() * exact_depth2_single
        entry = paths[bisect.bisect_left(cumulative, draw)]
        candidates = sample_matrix_for_path(configuration, entry.path, rng)
        owner, placement, unplaced = force_depth2_failure(
            candidates, configuration, rng
        )
        for depth in depths:
            failures[depth] += bounded_bfs_repair(
                candidates,
                configuration,
                owner,
                placement,
                unplaced,
                depth,
            ) > 0

    rows = []
    for depth in depths:
        conditional = failures[depth] / args.samples
        estimate = exact_depth2_single * conditional
        standard_error = exact_depth2_single * math.sqrt(
            conditional * (1.0 - conditional) / args.samples
        )
        row = {
            "t": args.t,
            "w": args.w,
            "d": args.d,
            "terminal_live": 1,
            "target_event": "depth2_bfs_failure",
            "depth": depth,
            "samples": args.samples,
            "target_probability": exact_depth2_single,
            "target_probability_log2": math.log2(exact_depth2_single),
            "conditional_failures": failures[depth],
            "conditional_failure_rate": conditional,
            "failure_estimate": estimate,
            "failure_estimate_log2": (
                -math.inf if estimate == 0.0 else math.log2(estimate)
            ),
            "standard_error": standard_error,
        }
        rows.append(row)
        print(
            f"depth={depth} failures={failures[depth]}/{args.samples} "
            f"conditional={conditional:.8e} estimate={estimate:.8e} "
            f"log2={row['failure_estimate_log2']}"
        )
    print(
        f"exact one-overflow depth-2 failure={exact_depth2_single:.12e} "
        f"log2={math.log2(exact_depth2_single):.6f} paths={len(paths)}"
    )

    if args.csv is not None:
        with args.csv.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=tuple(rows[0]),
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
