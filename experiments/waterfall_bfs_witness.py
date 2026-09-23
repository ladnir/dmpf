#!/usr/bin/env python3
"""Search for matchable Waterfall instances needing a deep BFS repair.

The search first runs the ordinary fixed-priority waterfall placement.  It
keeps matrices with a bounded number of overflow rows, a perfect matching,
and no repair using the requested BFS depth.  A found witness can then be
minimized over subsets of its rows while preserving their order.

The built-in witness is a small, reproducible counterexample to the conjecture
that depth-four BFS always reaches every matchable Waterfall instance with at
most three overflow rows.  It has one overflow row and needs depth five.
"""

from __future__ import annotations

import argparse
import itertools
import json
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


KNOWN_WITNESS = [
    0, 2, 0, 3,
    1, 0, 0, 1,
    0, 2, 3, 0,
    3, 3, 0, 1,
    3, 2, 3, 1,
    0, 3, 3, 1,
    3, 3, 0, 1,
    3, 2, 0, 1,
    3, 3, 3, 1,
    0, 2, 3, 0,
    1, 3, 0, 3,
]


@dataclass(frozen=True)
class Evaluation:
    overflow: tuple[int, ...]
    placement: tuple[int, ...]
    matching_deficit: int
    minimum_depth: int | None


def evaluate(
    candidates: list[int],
    configuration: Configuration,
    rejected_depth: int,
    max_overflow: int,
) -> Evaluation | None:
    """Return diagnostics iff this is a qualifying deep, matchable witness."""
    owner, placement, overflow = waterfall(candidates, configuration)
    if not 1 <= len(overflow) <= max_overflow:
        return None
    if bounded_bfs_repair(
        candidates,
        configuration,
        owner,
        placement,
        overflow,
        rejected_depth,
    ) == 0:
        return None
    matching_deficit = maximum_matching_unplaced(candidates, configuration)
    if matching_deficit:
        return None

    minimum_depth = None
    for depth in range(rejected_depth + 1, configuration.t + 1):
        if bounded_bfs_repair(
            candidates,
            configuration,
            owner,
            placement,
            overflow,
            depth,
        ) == 0:
            minimum_depth = depth
            break
    if minimum_depth is None:
        raise AssertionError("perfect matching was not reached by full-depth BFS")
    return Evaluation(
        tuple(overflow),
        tuple(placement),
        matching_deficit,
        minimum_depth,
    )


def shortest_path_for_row(
    candidates: list[int],
    configuration: Configuration,
    owner: list[int],
    starting_row: int,
    max_depth: int | None = None,
) -> list[tuple[int, int]] | None:
    """Return one shortest augmenting path as (row, target-slot) pairs."""
    t, w, d = configuration.t, configuration.w, configuration.d
    queue = [starting_row]
    queue_index = 0
    seen_row = [False] * t
    seen_row[starting_row] = True
    depth = [-1] * t
    depth[starting_row] = 0
    seen_slot = [False] * (w * d)
    parent_row = [-1] * t
    parent_slot = [-1] * t

    terminal_row = -1
    terminal_slot = -1
    while queue_index < len(queue) and terminal_slot == -1:
        row = queue[queue_index]
        queue_index += 1
        if max_depth is not None and depth[row] >= max_depth:
            continue
        base = row * w
        for partition in range(w):
            slot = partition * d + candidates[base + partition]
            if seen_slot[slot]:
                continue
            seen_slot[slot] = True
            victim = owner[slot]
            if victim == -1:
                terminal_row = row
                terminal_slot = slot
                break
            if not seen_row[victim]:
                seen_row[victim] = True
                parent_row[victim] = row
                parent_slot[victim] = slot
                depth[victim] = depth[row] + 1
                queue.append(victim)

    if terminal_slot == -1:
        return None

    reverse_path = [(terminal_row, terminal_slot)]
    row = terminal_row
    while row != starting_row:
        row = parent_row[row]
        if row == -1:
            raise AssertionError("broken BFS parent chain")
        reverse_path.append((row, parent_slot[reverse_path[-1][0]]))
    reverse_path.reverse()
    return reverse_path


def sequential_repair_paths(
    candidates: list[int],
    configuration: Configuration,
    max_depth: int,
) -> tuple[list[list[tuple[int, int]]], list[int]]:
    """Run the sequential router and return its applied augmenting paths."""
    owner, placement, overflow = waterfall(candidates, configuration)
    paths: list[list[tuple[int, int]]] = []
    unrepaired: list[int] = []
    for starting_row in overflow:
        path = shortest_path_for_row(
            candidates,
            configuration,
            owner,
            starting_row,
            max_depth,
        )
        if path is None:
            unrepaired.append(starting_row)
            continue

        for row, slot in reversed(path):
            if owner[slot] != -1:
                raise AssertionError("repair path did not terminate at a free bin")
            old_slot = placement[row]
            owner[slot] = row
            placement[row] = slot
            if old_slot != -1:
                owner[old_slot] = -1
        paths.append(path)
    return paths, unrepaired


def select_rows(candidates: list[int], w: int, rows: tuple[int, ...]) -> list[int]:
    selected: list[int] = []
    for row in rows:
        base = row * w
        selected.extend(candidates[base : base + w])
    return selected


def minimize_rows(
    candidates: list[int],
    configuration: Configuration,
    rejected_depth: int,
    max_overflow: int,
    exhaustive_limit: int,
) -> tuple[list[int], Configuration, tuple[int, ...]]:
    """Find a minimum-row subwitness when the subset space is affordable."""
    t, w = configuration.t, configuration.w
    if 1 << t > exhaustive_limit:
        rows = list(range(t))
        changed = True
        while changed:
            changed = False
            for row in rows.copy():
                trial_rows = tuple(candidate for candidate in rows if candidate != row)
                trial = select_rows(candidates, w, trial_rows)
                trial_configuration = Configuration(len(trial_rows), w, configuration.d, 0)
                if evaluate(trial, trial_configuration, rejected_depth, max_overflow):
                    rows = list(trial_rows)
                    changed = True
                    break
        chosen = tuple(rows)
        return (
            select_rows(candidates, w, chosen),
            Configuration(len(chosen), w, configuration.d, 0),
            chosen,
        )

    for size in range(1, t + 1):
        for rows in itertools.combinations(range(t), size):
            trial = select_rows(candidates, w, rows)
            trial_configuration = Configuration(size, w, configuration.d, 0)
            if evaluate(trial, trial_configuration, rejected_depth, max_overflow):
                return trial, trial_configuration, rows
    raise AssertionError("the input witness ceased to be a witness")


def embed_witness(
    candidates: list[int],
    configuration: Configuration,
    target_t: int,
    target_d: int,
) -> tuple[list[int], Configuration]:
    """Embed a witness and add isolated rows that do not change its routing."""
    if target_t < configuration.t or target_d < configuration.d:
        raise ValueError("target dimensions must contain the source witness")
    if target_d - configuration.d < target_t - configuration.t:
        raise ValueError("need one fresh bin label for every appended row")

    embedded = candidates.copy()
    for row in range(configuration.t, target_t):
        fresh_bin = configuration.d + row - configuration.t
        embedded.extend([fresh_bin] * configuration.w)
    return embedded, Configuration(target_t, configuration.w, target_d, 0)


def random_search(
    configuration: Configuration,
    rejected_depth: int,
    max_overflow: int,
    trials: int,
    rng: random.Random,
) -> tuple[list[int], int] | None:
    """Return the first uniform random witness and its one-based trial index."""
    entries = configuration.t * configuration.w
    d = configuration.d
    for trial in range(1, trials + 1):
        candidates = [rng.randrange(d) for _ in range(entries)]
        if evaluate(candidates, configuration, rejected_depth, max_overflow):
            return candidates, trial
    return None


def matrix_rows(candidates: list[int], configuration: Configuration) -> list[list[int]]:
    w = configuration.w
    return [candidates[row * w : (row + 1) * w] for row in range(configuration.t)]


def equality_pattern_statistics(
    candidates: list[int], configuration: Configuration
) -> tuple[list[int], int, float]:
    """Return distinct counts, collision excess, and exact pattern log-mass.

    Collision excess is ``sum(t-k_s)``, where partition ``s`` uses ``k_s``
    distinct labels.  The log-mass is for the complete unlabeled equality
    pattern, including entries irrelevant to routing; it is consequently a
    diagnostic score rather than a probability bound for the witness event.
    """
    t, w, d = configuration.t, configuration.w, configuration.d
    distinct = [
        len({candidates[row * w + partition] for row in range(t)})
        for partition in range(w)
    ]
    collision_excess = sum(t - count for count in distinct)
    log2_mass = 0.0
    for count in distinct:
        log2_mass += sum(math.log2(d - index) for index in range(count))
        log2_mass -= t * math.log2(d)
    return distinct, collision_excess, log2_mass


def optimize_collision_excess(
    candidates: list[int],
    configuration: Configuration,
    rejected_depth: int,
    max_overflow: int,
    iterations: int,
    rng: random.Random,
) -> tuple[list[int], int]:
    """Locally remove accidental equalities while preserving the witness.

    The hot loop mutates one matrix entry in place and evaluates only matrices
    that change the selected label.  Equal-score moves are accepted to cross
    plateaus; worse moves are rejected.  Periodic restarts return to the best
    matrix found so a long neutral walk cannot lose progress.
    """
    if evaluate(candidates, configuration, rejected_depth, max_overflow) is None:
        raise ValueError("collision optimization requires a qualifying witness")
    current = candidates.copy()
    _, current_score, _ = equality_pattern_statistics(current, configuration)
    best = current.copy()
    best_score = current_score
    entries = len(current)
    d = configuration.d

    for iteration in range(iterations):
        if iteration and iteration % (entries * 64) == 0:
            current = best.copy()
            current_score = best_score
        index = rng.randrange(entries)
        old_label = current[index]
        new_label = rng.randrange(d - 1)
        if new_label >= old_label:
            new_label += 1
        current[index] = new_label
        statistics = equality_pattern_statistics(current, configuration)
        score = statistics[1]
        if score <= current_score and evaluate(
            current, configuration, rejected_depth, max_overflow
        ) is not None:
            current_score = score
            if score < best_score:
                best = current.copy()
                best_score = score
        else:
            current[index] = old_label
    return best, best_score


def describe(
    candidates: list[int],
    configuration: Configuration,
    rejected_depth: int,
    max_overflow: int,
) -> dict[str, object]:
    evaluation = evaluate(candidates, configuration, rejected_depth, max_overflow)
    if evaluation is None:
        raise ValueError("matrix is not a qualifying witness")
    if evaluation.minimum_depth is None:
        raise AssertionError("qualifying witness has no successful repair depth")
    paths, unrepaired = sequential_repair_paths(
        candidates, configuration, evaluation.minimum_depth
    )
    if unrepaired:
        raise AssertionError("reported minimum depth did not repair every row")
    distinct, collision_excess, equality_pattern_log2_mass = (
        equality_pattern_statistics(candidates, configuration)
    )
    return {
        "t": configuration.t,
        "w": configuration.w,
        "d": configuration.d,
        "rejected_depth": rejected_depth,
        "overflow_rows": list(evaluation.overflow),
        "waterfall_placement_slots": list(evaluation.placement),
        "minimum_bfs_depth": evaluation.minimum_depth,
        "augmenting_paths": [
            [
                {
                    "row": row,
                    "partition": slot // configuration.d,
                    "bin": slot % configuration.d,
                }
                for row, slot in path
            ]
            for path in paths
        ],
        "distinct_labels_by_partition": distinct,
        "collision_excess": collision_excess,
        "complete_equality_pattern_log2_mass": equality_pattern_log2_mass,
        "matrix": matrix_rows(candidates, configuration),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", action="store_true", help="search instead of checking the built-in witness")
    parser.add_argument("--t", type=int, default=12)
    parser.add_argument("--w", type=int, default=4)
    parser.add_argument("--d", type=int, default=4)
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--max-overflow", type=int, default=3)
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=12_345)
    parser.add_argument("--exhaustive-limit", type=int, default=1 << 20)
    parser.add_argument("--embed-t", type=int, default=16)
    parser.add_argument("--embed-d", type=int, default=16)
    parser.add_argument("--json", type=Path)
    parser.add_argument(
        "--optimize-collisions",
        type=int,
        default=0,
        metavar="ITERATIONS",
        help="run a neutral local search minimizing collision excess",
    )
    args = parser.parse_args()

    if args.search:
        configuration = Configuration(args.t, args.w, args.d, 0)
        result = random_search(
            configuration,
            args.depth,
            args.max_overflow,
            args.trials,
            random.Random(args.seed),
        )
        if result is None:
            raise SystemExit(f"no witness found in {args.trials:,} trials")
        candidates, trial = result
    else:
        configuration = Configuration(11, 4, 4, 0)
        candidates = KNOWN_WITNESS.copy()
        trial = None

    minimized, minimized_configuration, source_rows = minimize_rows(
        candidates,
        configuration,
        args.depth,
        args.max_overflow,
        args.exhaustive_limit,
    )
    embedded, embedded_configuration = embed_witness(
        minimized,
        minimized_configuration,
        args.embed_t,
        args.embed_d,
    )
    preoptimization = None
    if args.optimize_collisions:
        preoptimization = describe(
            embedded,
            embedded_configuration,
            args.depth,
            args.max_overflow,
        )
        embedded, _ = optimize_collision_excess(
            embedded,
            embedded_configuration,
            args.depth,
            args.max_overflow,
            args.optimize_collisions,
            random.Random(args.seed ^ 0xC0111510),
        )
    output: dict[str, object] = {
        "search_trial": trial,
        "source_rows_retained": list(source_rows),
        "minimized": describe(
            minimized,
            minimized_configuration,
            args.depth,
            args.max_overflow,
        ),
        "embedded": describe(
            embedded,
            embedded_configuration,
            args.depth,
            args.max_overflow,
        ),
    }
    if preoptimization is not None:
        output["embedded_before_collision_optimization"] = preoptimization
    rendered = json.dumps(output, indent=2)
    print(rendered)
    if args.json:
        args.json.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
