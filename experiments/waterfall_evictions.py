#!/usr/bin/env python3
"""Simulate bounded eviction repair for Waterfall Cuckoo.

The initial placement is the fixed-priority waterfall from WaterfallCuckoo.tex.
Repair maintains a current unplaced row and performs bounded eviction chains:

1. place the current row in its first free candidate, if one exists;
2. otherwise evict an occupant and continue with that occupant;
3. never immediately traverse the previously used bin in reverse.

The default ``global`` schedule shares one eviction-step budget across every
overflow row; the diagnostic ``per-row`` schedule gives that entire budget to
each row.  The ``biased`` policy avoids the final partition on the first
eviction in a chain, but permits every partition on later evictions.
``uniform`` permits every partition from the first eviction.  All policies use
the same candidate matrices in a run.  The ``oracle`` row reports maximum
bipartite matching and is an unattainable lower bound on the number of
unplaced rows.

Monte Carlo cannot certify a 2^-40 tail.  This program is intended to compare
router shapes in their observable regime and to identify configurations for a
subsequent rigorous analysis.  It reports structural MPC-cost components
rather than combining them using unjustified implementation-dependent weights.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(frozen=True)
class Configuration:
    t: int
    w: int
    d: int
    c: int


@dataclass(frozen=True)
class Experiment:
    policy: str
    evictions: int
    schedule: str


@dataclass
class Accumulator:
    trials: int = 0
    failures: int = 0
    unplaced_sum: int = 0
    histogram: Counter[int] | None = None

    def __post_init__(self) -> None:
        if self.histogram is None:
            self.histogram = Counter()

    def add(self, unplaced: int, stash: int) -> None:
        self.trials += 1
        self.failures += unplaced > stash
        self.unplaced_sum += unplaced
        assert self.histogram is not None
        self.histogram[unplaced] += 1


def waterfall(
    candidates: list[int], configuration: Configuration
) -> tuple[list[int], list[int], list[int]]:
    """Return (owner by main bin, placement by row, unplaced rows)."""
    t, w, d = configuration.t, configuration.w, configuration.d
    owner = [-1] * (w * d)
    placement = [-1] * t
    live = list(range(t))

    for partition in range(w):
        offset = partition * d
        candidate_column = partition
        next_live: list[int] = []
        for row in live:
            slot = offset + candidates[row * w + candidate_column]
            if owner[slot] == -1:
                owner[slot] = row
                placement[row] = slot
            else:
                next_live.append(row)
        live = next_live

    return owner, placement, live


def repair_one_chain(
    candidates: list[int],
    configuration: Configuration,
    owner: list[int],
    placement: list[int],
    starting_row: int,
    policy: str,
    evictions: int,
    rng: random.Random,
) -> tuple[int, int]:
    """Run one chain and return (current token or -1, steps consumed)."""
    _, w, d = configuration.t, configuration.w, configuration.d
    token = starting_row
    previous_slot = -1

    for step in range(evictions):
        # Any free candidate completes the augmenting path.  Earliest-tier
        # placement preserves the waterfall preference for later chains.
        free_slot = -1
        base = token * w
        for partition in range(w):
            slot = partition * d + candidates[base + partition]
            if owner[slot] == -1:
                free_slot = slot
                break

        if free_slot != -1:
            owner[free_slot] = token
            placement[token] = free_slot
            return -1, step + 1

        if policy == "biased" and step == 0 and w > 1:
            partitions: Iterable[int] = range(w - 1)
        elif policy in ("biased", "uniform"):
            partitions = range(w)
        else:
            raise ValueError(f"unknown repair policy: {policy}")

        eviction_slots = [
            partition * d + candidates[base + partition]
            for partition in partitions
            if partition * d + candidates[base + partition] != previous_slot
        ]
        # This is reachable only for a degenerate one-choice configuration.
        if not eviction_slots:
            eviction_slots = [
                partition * d + candidates[base + partition]
                for partition in range(w)
            ]

        slot = eviction_slots[rng.randrange(len(eviction_slots))]
        victim = owner[slot]
        if victim == -1:
            raise AssertionError("free candidates must be handled first")

        owner[slot] = token
        placement[token] = slot
        placement[victim] = -1
        token = victim
        previous_slot = slot

    return token, evictions


def repair(
    candidates: list[int],
    configuration: Configuration,
    initial_owner: list[int],
    initial_placement: list[int],
    initial_unplaced: list[int],
    experiment: Experiment,
    rng: random.Random,
) -> int:
    """Return the number of rows unplaced after bounded eviction repair."""
    if experiment.evictions == 0:
        return len(initial_unplaced)
    if len(initial_unplaced) <= configuration.c:
        return len(initial_unplaced)

    owner = initial_owner.copy()
    placement = initial_placement.copy()
    unplaced_count = len(initial_unplaced)

    if experiment.schedule == "per-row":
        for starting_row in initial_unplaced:
            if unplaced_count <= configuration.c:
                break
            if placement[starting_row] != -1:
                raise AssertionError("an initially unplaced row became placed early")
            token, _ = repair_one_chain(
                candidates,
                configuration,
                owner,
                placement,
                starting_row,
                experiment.policy,
                experiment.evictions,
                rng,
            )
            if token == -1:
                unplaced_count -= 1
    elif experiment.schedule == "global":
        remaining_budget = experiment.evictions
        pending_index = 0
        while pending_index < len(initial_unplaced) and remaining_budget > 0:
            starting_row = initial_unplaced[pending_index]
            pending_index += 1
            token, consumed = repair_one_chain(
                candidates,
                configuration,
                owner,
                placement,
                starting_row,
                experiment.policy,
                remaining_budget,
                rng,
            )
            remaining_budget -= consumed
            if token == -1:
                unplaced_count -= 1
                if unplaced_count <= configuration.c:
                    break
            # A capped chain leaves exactly one current row unplaced.  There is
            # no budget with which to continue it or start another chain.
            if token != -1:
                break
    else:
        raise ValueError(f"unknown eviction schedule: {experiment.schedule}")

    return sum(slot == -1 for slot in placement)


def exact_repair_failure_probability(
    candidates: list[int],
    configuration: Configuration,
    initial_owner: list[int],
    initial_placement: list[int],
    initial_unplaced: list[int],
    experiment: Experiment,
) -> float:
    """Integrate all router coins for a global bounded-eviction schedule.

    This enumerates at most ``w*(w-1)^(e-1)`` eviction-choice paths and uses
    in-place backtracking to avoid allocating placement state at every node.
    It is intended for rare-event analysis at small ``e``, not the main Monte
    Carlo sweep.
    """
    if experiment.schedule != "global":
        raise ValueError("exact router enumeration requires the global schedule")
    if experiment.policy not in ("biased", "uniform"):
        raise ValueError(f"unknown repair policy: {experiment.policy}")
    if len(initial_unplaced) <= configuration.c:
        return 0.0
    if experiment.evictions == 0:
        return 1.0

    _, w, d = configuration.t, configuration.w, configuration.d
    owner = initial_owner.copy()
    placement = initial_placement.copy()

    def recurse(
        token: int,
        pending_index: int,
        unplaced_count: int,
        remaining_budget: int,
        first_eviction: bool,
        previous_slot: int,
    ) -> float:
        if unplaced_count <= configuration.c:
            return 0.0
        if remaining_budget == 0:
            return 1.0

        base = token * w
        for partition in range(w):
            slot = partition * d + candidates[base + partition]
            if owner[slot] != -1:
                continue
            owner[slot] = token
            placement[token] = slot
            next_unplaced = unplaced_count - 1
            if next_unplaced <= configuration.c:
                probability = 0.0
            else:
                if pending_index >= len(initial_unplaced):
                    raise AssertionError("missing pending overflow row")
                next_token = initial_unplaced[pending_index]
                if placement[next_token] != -1:
                    raise AssertionError("pending overflow row became placed")
                probability = recurse(
                    next_token,
                    pending_index + 1,
                    next_unplaced,
                    remaining_budget - 1,
                    True,
                    -1,
                )
            owner[slot] = -1
            placement[token] = -1
            return probability

        if experiment.policy == "biased" and first_eviction and w > 1:
            partitions: Iterable[int] = range(w - 1)
        else:
            partitions = range(w)
        eviction_slots = [
            partition * d + candidates[base + partition]
            for partition in partitions
            if partition * d + candidates[base + partition] != previous_slot
        ]
        if not eviction_slots:
            eviction_slots = [
                partition * d + candidates[base + partition]
                for partition in range(w)
            ]

        probability = 0.0
        branch_weight = 1.0 / len(eviction_slots)
        for slot in eviction_slots:
            victim = owner[slot]
            if victim == -1:
                raise AssertionError("free candidates must be handled first")
            owner[slot] = token
            placement[token] = slot
            placement[victim] = -1
            probability += branch_weight * recurse(
                victim,
                pending_index,
                unplaced_count,
                remaining_budget - 1,
                False,
                slot,
            )
            owner[slot] = victim
            placement[token] = -1
            placement[victim] = slot
        return probability

    return recurse(
        initial_unplaced[0],
        1,
        len(initial_unplaced),
        experiment.evictions,
        True,
        -1,
    )


def maximum_matching_unplaced(
    candidates: list[int], configuration: Configuration
) -> int:
    """Return t minus the maximum matching size in the partitioned graph."""
    t, w, d = configuration.t, configuration.w, configuration.d
    owner = [-1] * (w * d)

    def augment(row: int, seen: list[bool]) -> bool:
        base = row * w
        for partition in range(w):
            slot = partition * d + candidates[base + partition]
            if seen[slot]:
                continue
            seen[slot] = True
            victim = owner[slot]
            if victim == -1 or augment(victim, seen):
                owner[slot] = row
                return True
        return False

    matched = 0
    for row in range(t):
        matched += augment(row, [False] * (w * d))
    return t - matched


def bounded_bfs_repair(
    candidates: list[int],
    configuration: Configuration,
    initial_owner: list[int],
    initial_placement: list[int],
    initial_unplaced: list[int],
    max_steps: int,
) -> int:
    """Repair with shortest augmenting paths of at most ``max_steps`` bins.

    The search alternates from a row to its candidate bins and then, for an
    occupied bin, to its owner.  On finding a free bin, parent pointers apply
    the augmenting path backwards in place.  This is a deterministic control
    for evaluating stronger routing circuits; it is not used by the biased
    random-walk protocol.
    """
    if max_steps < 0:
        raise ValueError("max_steps must be nonnegative")
    t, w, d, c = (
        configuration.t,
        configuration.w,
        configuration.d,
        configuration.c,
    )
    owner = initial_owner.copy()
    placement = initial_placement.copy()
    unplaced_count = len(initial_unplaced)

    for starting_row in initial_unplaced:
        if unplaced_count <= c:
            break
        if placement[starting_row] != -1:
            raise AssertionError("an initially unplaced row became placed early")

        queue = [starting_row]
        queue_index = 0
        depth = [-1] * t
        depth[starting_row] = 0
        parent_row = [-1] * t
        parent_slot = [-1] * t
        seen_slot = [False] * (w * d)
        terminal_row = -1
        terminal_slot = -1

        while queue_index < len(queue) and terminal_slot == -1:
            row = queue[queue_index]
            queue_index += 1
            if depth[row] >= max_steps:
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
                if depth[victim] == -1:
                    depth[victim] = depth[row] + 1
                    parent_row[victim] = row
                    parent_slot[victim] = slot
                    queue.append(victim)

        if terminal_slot == -1:
            continue

        row = terminal_row
        slot = terminal_slot
        while True:
            if owner[slot] != -1:
                raise AssertionError("augmenting path did not free its next bin")
            old_slot = placement[row]
            owner[slot] = row
            placement[row] = slot
            if row == starting_row:
                if old_slot != -1:
                    raise AssertionError("starting row was unexpectedly placed")
                break
            if old_slot == -1:
                raise AssertionError("internal augmenting-path row was unplaced")
            owner[old_slot] = -1
            next_row = parent_row[row]
            next_slot = parent_slot[row]
            if next_slot != old_slot:
                raise AssertionError("broken augmenting-path parent pointer")
            row = next_row
            slot = next_slot
        unplaced_count -= 1

    return unplaced_count


def wilson_interval(successes: int, trials: int, z: float = 1.95996398454) -> tuple[float, float]:
    """Two-sided Wilson interval for a binomial probability."""
    if trials <= 0:
        raise ValueError("trials must be positive")
    estimate = successes / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials
    center = (estimate + z2 / (2.0 * trials)) / denominator
    radius = z * math.sqrt(
        estimate * (1.0 - estimate) / trials + z2 / (4.0 * trials * trials)
    ) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)


def log2_or_negative_infinity(value: float) -> float:
    return -math.inf if value == 0.0 else math.log2(value)


def exact_waterfall_failure(configuration: Configuration) -> float:
    """Evaluate the paper's exact occupancy recurrence in binary64."""
    t, w, d, c = (
        configuration.t,
        configuration.w,
        configuration.d,
        configuration.c,
    )
    stirling = [[0] * (t + 1) for _ in range(t + 1)]
    stirling[0][0] = 1
    for rows in range(1, t + 1):
        for occupied in range(1, rows + 1):
            stirling[rows][occupied] = (
                stirling[rows - 1][occupied - 1]
                + occupied * stirling[rows - 1][occupied]
            )

    live_distribution = [0.0] * (t + 1)
    live_distribution[t] = 1.0
    for _ in range(w):
        next_distribution = [0.0] * (t + 1)
        for live, live_probability in enumerate(live_distribution):
            if live_probability == 0.0:
                continue
            if live == 0:
                next_distribution[0] += live_probability
                continue
            falling = 1
            denominator = d**live
            for occupied in range(1, min(live, d) + 1):
                falling *= d - occupied + 1
                probability = falling * stirling[live][occupied] / denominator
                next_distribution[live - occupied] += live_probability * probability
        live_distribution = next_distribution
    return sum(live_distribution[c + 1 :])


def exact_two_step_biased_failure(configuration: Configuration) -> float:
    """Return the exact failure probability of global biased repair with e=2.

    If L_w>=c+2, one successful augmentation still leaves an overflow.  When
    L_w=c+1 and the first eviction uses partition s<w-1, its victim previously
    lost in every earlier partition, so those candidates are occupied.  The
    victim was placed in s, hence its later candidates were never inspected
    and remain independent uniform values.  For a live-count path, failure in
    branch s is therefore the product of the later occupied-bin counts divided
    by the corresponding powers of d.  Averaging s and summing path
    probabilities gives the exact expression below.
    """
    t, w, d, c = (
        configuration.t,
        configuration.w,
        configuration.d,
        configuration.c,
    )
    if w < 2:
        raise ValueError("biased two-step repair requires at least two partitions")

    stirling = [[0] * (t + 1) for _ in range(t + 1)]
    stirling[0][0] = 1
    for rows in range(1, t + 1):
        for occupied in range(1, rows + 1):
            stirling[rows][occupied] = (
                stirling[rows - 1][occupied - 1]
                + occupied * stirling[rows - 1][occupied]
            )

    transitions: list[list[tuple[int, float]]] = [[] for _ in range(t + 1)]
    transitions[0].append((0, 1.0))
    for live in range(1, t + 1):
        falling = 1
        denominator = d**live
        for occupied in range(1, min(live, d) + 1):
            falling *= d - occupied + 1
            probability = falling * stirling[live][occupied] / denominator
            transitions[live].append((live - occupied, probability))

    total_failure = 0.0
    path = [t]

    def enumerate_paths(partition: int, probability: float) -> None:
        nonlocal total_failure
        live = path[-1]
        if partition == w:
            if live >= c + 2:
                total_failure += probability
                return
            if live != c + 1:
                return

            occupied = [path[index] - path[index + 1] for index in range(w)]
            conditional_failure = 0.0
            for first_partition in range(w - 1):
                branch_failure = 1.0
                for later_partition in range(first_partition + 1, w):
                    branch_failure *= occupied[later_partition] / d
                conditional_failure += branch_failure
            conditional_failure /= w - 1
            total_failure += probability * conditional_failure
            return

        for next_live, transition_probability in transitions[live]:
            path.append(next_live)
            enumerate_paths(
                partition + 1, probability * transition_probability
            )
            path.pop()

    enumerate_paths(0, 1.0)
    return total_failure


def percentile_from_histogram(histogram: Counter[int], trials: int, probability: float) -> int:
    threshold = math.ceil(probability * trials)
    cumulative = 0
    for value in sorted(histogram):
        cumulative += histogram[value]
        if cumulative >= threshold:
            return value
    raise AssertionError("incomplete histogram")


def structural_cost(
    configuration: Configuration, experiment: Experiment
) -> dict[str, int]:
    """Return implementation-independent circuit-size components.

    Each step scans ``w`` candidates against up to ``t`` owners.  A concrete
    circuit will add token multiplexing, state updates, and selection gates, so
    ``repair_scan_units`` is a comparison unit rather than a gate count.  The
    global schedule uses exactly ``e`` circuit steps; the per-row diagnostic
    schedules ``e*t`` steps and masks rows that were not initially unplaced.
    """
    t, w, d, c = (
        configuration.t,
        configuration.w,
        configuration.d,
        configuration.c,
    )
    index_bits = max(1, (d - 1).bit_length())
    waterfall_pair_tests = w * math.comb(t, 2)
    if experiment.schedule == "global":
        repair_steps = experiment.evictions
    elif experiment.schedule == "per-row":
        repair_steps = experiment.evictions * t
    else:
        repair_steps = 0
    repair_scan_units = repair_steps * w * t
    optimized_repair_scan_units = repair_scan_units
    if experiment.schedule == "global" and experiment.evictions == 2:
        # The initial overflow row has no free candidate, so step one only
        # looks up the selected victim.  On step two, the traversed bin is
        # occupied by the first row, so only the victim's w-1 alternatives
        # need occupancy tests.  No final neutral eviction is required.
        optimized_repair_scan_units = w * t
    return {
        "output_columns": w * d + c,
        "local_expansion_N": w + c,
        "bin_solver_calls": w,
        "index_bits": index_bits,
        "waterfall_pair_tests": waterfall_pair_tests,
        "main_one_hot_bits": w * t * d,
        "repair_steps": repair_steps,
        "repair_scan_units": repair_scan_units,
        "rough_repair_index_bit_ops": 2 * index_bits * repair_scan_units,
        "optimized_repair_scan_units": optimized_repair_scan_units,
        "rough_optimized_repair_index_bit_ops": (
            2 * index_bits * optimized_repair_scan_units
        ),
    }


FIELDNAMES = (
    "t",
    "w",
    "d",
    "c",
    "policy",
    "schedule",
    "evictions",
    "trials",
    "failures",
    "failure_rate",
    "failure_log2",
    "failure_ci95_low",
    "failure_ci95_high",
    "exact_waterfall_failure",
    "exact_waterfall_failure_log2",
    "batch_failure_log2",
    "batch_ci95_high_log2",
    "mean_unplaced",
    "p95_unplaced",
    "output_columns",
    "local_expansion_N",
    "bin_solver_calls",
    "index_bits",
    "waterfall_pair_tests",
    "main_one_hot_bits",
    "repair_steps",
    "repair_scan_units",
    "rough_repair_index_bit_ops",
    "optimized_repair_scan_units",
    "rough_optimized_repair_index_bit_ops",
)


def result_row(
    configuration: Configuration,
    experiment: Experiment,
    accumulator: Accumulator,
    batch: int,
) -> dict[str, int | float | str]:
    low, high = wilson_interval(accumulator.failures, accumulator.trials)
    rate = accumulator.failures / accumulator.trials
    exact_failure = exact_waterfall_failure(configuration)
    assert accumulator.histogram is not None
    cost = structural_cost(configuration, experiment)
    return {
        "t": configuration.t,
        "w": configuration.w,
        "d": configuration.d,
        "c": configuration.c,
        "policy": experiment.policy,
        "schedule": experiment.schedule,
        "evictions": experiment.evictions,
        "trials": accumulator.trials,
        "failures": accumulator.failures,
        "failure_rate": rate,
        "failure_log2": log2_or_negative_infinity(rate),
        "failure_ci95_low": low,
        "failure_ci95_high": high,
        "exact_waterfall_failure": exact_failure,
        "exact_waterfall_failure_log2": log2_or_negative_infinity(exact_failure),
        "batch_failure_log2": min(0.0, log2_or_negative_infinity(batch * rate)),
        "batch_ci95_high_log2": min(0.0, math.log2(batch * high)),
        "mean_unplaced": accumulator.unplaced_sum / accumulator.trials,
        "p95_unplaced": percentile_from_histogram(
            accumulator.histogram, accumulator.trials, 0.95
        ),
        **cost,
    }


def run_configuration(
    configuration: Configuration,
    experiments: tuple[Experiment, ...],
    trials: int,
    batch: int,
    seed: int,
) -> list[dict[str, int | float | str]]:
    descriptor_rng = random.Random(seed)
    policy_rngs = {
        experiment: random.Random(
            seed
            ^ (0x9E3779B97F4A7C15 * (index + 1))
            ^ (experiment.evictions << 32)
        )
        for index, experiment in enumerate(experiments)
    }
    accumulators = {experiment: Accumulator() for experiment in experiments}
    oracle = Accumulator()
    t, w, d = configuration.t, configuration.w, configuration.d

    for _ in range(trials):
        candidates = [descriptor_rng.randrange(d) for _ in range(t * w)]
        owner, placement, unplaced = waterfall(candidates, configuration)
        for experiment in experiments:
            remaining = repair(
                candidates,
                configuration,
                owner,
                placement,
                unplaced,
                experiment,
                policy_rngs[experiment],
            )
            accumulators[experiment].add(remaining, configuration.c)
        oracle.add(
            maximum_matching_unplaced(candidates, configuration), configuration.c
        )

    rows = [
        result_row(configuration, experiment, accumulators[experiment], batch)
        for experiment in experiments
    ]
    rows.append(
        result_row(
            configuration,
            Experiment("oracle", 0, "none"),
            oracle,
            batch,
        )
    )
    return rows


def validate_args(args: argparse.Namespace) -> None:
    if args.t < 1:
        raise SystemExit("--t must be positive")
    if args.w < 2:
        raise SystemExit("--w must be at least two")
    if args.trials < 1:
        raise SystemExit("--trials must be positive")
    if args.batch < 1:
        raise SystemExit("--batch must be positive")
    if any(d < 1 for d in args.d):
        raise SystemExit("every --d must be positive")
    if any(c < 0 for c in args.c):
        raise SystemExit("every --c must be nonnegative")
    if any(value < 0 for value in args.evictions):
        raise SystemExit("every --evictions value must be nonnegative")


def output_stream(path: Path | None) -> tuple[TextIO, bool]:
    if path is None:
        return sys.stdout, False
    return path.open("w", newline=""), True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument("--w", type=int, default=3)
    parser.add_argument("--d", type=int, nargs="+", default=[8, 10, 12, 16, 24, 32])
    parser.add_argument("--c", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--evictions", type=int, nargs="+", default=[0, 1, 2, 4, 8])
    parser.add_argument("--policies", choices=("biased", "uniform"), nargs="+", default=["biased", "uniform"])
    parser.add_argument(
        "--schedules",
        choices=("global", "per-row"),
        nargs="+",
        default=["global"],
        help="global eviction budget (cheap) or that budget for every row",
    )
    parser.add_argument("--trials", type=int, default=100_000)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--seed", type=int, default=0xC0C00)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    validate_args(args)

    # e=0 is policy-independent, so emit a single waterfall baseline.
    experiments = [Experiment("waterfall", 0, "none")]
    experiments.extend(
        Experiment(policy, evictions, schedule)
        for evictions in sorted(set(args.evictions))
        if evictions > 0
        for schedule in args.schedules
        for policy in args.policies
    )
    experiments_tuple = tuple(experiments)

    stream, should_close = output_stream(args.csv)
    try:
        writer = csv.DictWriter(
            stream,
            fieldnames=FIELDNAMES,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for d in args.d:
            for c in args.c:
                configuration = Configuration(args.t, args.w, d, c)
                configuration_seed = args.seed ^ (d << 16) ^ c
                for row in run_configuration(
                    configuration,
                    experiments_tuple,
                    args.trials,
                    args.batch,
                    configuration_seed,
                ):
                    writer.writerow(row)
    finally:
        if should_close:
            stream.close()


if __name__ == "__main__":
    main()
