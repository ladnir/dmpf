#!/usr/bin/env python3
"""Exact-overflow conditional sampling for bounded Waterfall repair.

The ordinary waterfall has a Markov state L_s: the number of rows entering
partition s.  Its transition probabilities are exact occupancy probabilities.
This program uses the corresponding backward probabilities to sample a
uniform candidate matrix C conditioned on L_w>c.  Given a sampled transition,
it samples a uniform mapping of the live labeled rows onto the required number
of occupied bins, so the resulting matrix has exactly the target conditional
distribution rather than a heuristic collision tilt.

For every conditioned matrix, all random choices of the bounded eviction
router are enumerated exactly.  If q_e is the resulting conditional failure
estimate, the unconditional failure estimate is

    Pr[L_w>c] * q_e,

where the first factor is obtained from the exact recurrence.  This is still a
statistical estimate over conditioned matrices, not a security proof, but it
removes both the rare initial-overflow event and router-coin noise.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from waterfall_evictions import (
    bounded_bfs_repair,
    Configuration,
    Experiment,
    exact_repair_failure_probability,
    exact_waterfall_failure,
    waterfall,
)


NEG_INF = -math.inf


@dataclass
class Moments:
    count: int = 0
    total: float = 0.0
    total_square: float = 0.0

    def add(self, value: float) -> None:
        self.count += 1
        self.total += value
        self.total_square += value * value

    @property
    def mean(self) -> float:
        return self.total / self.count

    @property
    def variance_of_mean(self) -> float:
        if self.count < 2:
            return math.nan
        centered = self.total_square - self.total * self.total / self.count
        return max(0.0, centered / (self.count - 1) / self.count)

    @property
    def effective_sample_size(self) -> float:
        if self.total_square == 0.0:
            return 0.0
        return self.total * self.total / self.total_square


def occupancy_transitions(t: int, d: int) -> list[list[tuple[int, float]]]:
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
    return transitions


def backward_event_probabilities(
    configuration: Configuration,
    transitions: list[list[tuple[int, float]]],
    terminal_live: int | None = None,
) -> list[list[float]]:
    t, w, c = configuration.t, configuration.w, configuration.c
    backward = [[0.0] * (t + 1) for _ in range(w + 1)]
    for live in range(t + 1):
        if terminal_live is None:
            backward[w][live] = float(live > c)
        else:
            backward[w][live] = float(live == terminal_live)
    for partition in range(w - 1, -1, -1):
        for live in range(t + 1):
            backward[partition][live] = sum(
                probability * backward[partition + 1][next_live]
                for next_live, probability in transitions[live]
            )
    return backward


def uniform_surjection_labels(rows: int, occupied: int, rng: random.Random) -> list[int]:
    """Sample a uniform onto map [rows] -> [occupied] exactly."""
    ways = [[0] * (occupied + 1) for _ in range(rows + 1)]
    ways[0][occupied] = 1
    for remaining in range(1, rows + 1):
        for used in range(occupied, -1, -1):
            existing = used * ways[remaining - 1][used]
            new = 0
            if used < occupied:
                new = (occupied - used) * ways[remaining - 1][used + 1]
            ways[remaining][used] = existing + new

    labels: list[int] = []
    used_labels: list[int] = []
    unused_labels = list(range(occupied))
    for position in range(rows):
        remaining = rows - position
        used = len(used_labels)
        existing_ways = ways[remaining - 1][used]
        new_ways = ways[remaining - 1][used + 1] if used < occupied else 0
        existing_total = used * existing_ways
        total = ways[remaining][used]
        draw = rng.randrange(total)
        if draw < existing_total:
            labels.append(used_labels[draw // existing_ways])
        else:
            new_index = (draw - existing_total) // new_ways
            label = unused_labels.pop(new_index)
            used_labels.append(label)
            labels.append(label)
    if len(used_labels) != occupied:
        raise AssertionError("surjection sampler omitted an occupied label")
    return labels


def choose_weighted_transition(
    options: list[tuple[int, float]],
    next_backward: list[float],
    total_probability: float,
    rng: random.Random,
) -> int:
    draw = rng.random() * total_probability
    cumulative = 0.0
    fallback = -1
    for next_live, probability in options:
        weight = probability * next_backward[next_live]
        if weight == 0.0:
            continue
        fallback = next_live
        cumulative += weight
        if draw <= cumulative:
            return next_live
    if fallback == -1:
        raise AssertionError("conditioned transition has no support")
    return fallback


def sample_conditioned_matrix(
    configuration: Configuration,
    transitions: list[list[tuple[int, float]]],
    backward: list[list[float]],
    rng: random.Random,
    terminal_live: int | None = None,
) -> tuple[list[int], tuple[int, ...]]:
    t, w, d = configuration.t, configuration.w, configuration.d
    candidates = [0] * (t * w)
    live_rows = list(range(t))
    path = [t]

    for partition in range(w):
        live = len(live_rows)
        next_live_count = choose_weighted_transition(
            transitions[live],
            backward[partition + 1],
            backward[partition][live],
            rng,
        )
        occupied = live - next_live_count

        # Rows that already left the waterfall have unused, independent target
        # entries in this column.  Fill every row uniformly, then override the
        # live rows with the exactly conditioned occupancy map.
        for row in range(t):
            candidates[row * w + partition] = rng.randrange(d)
        bins = rng.sample(range(d), occupied)
        labels = uniform_surjection_labels(live, occupied, rng)

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
            raise AssertionError("sampled transition has the wrong live count")
        live_rows = next_live_rows
        path.append(next_live_count)

    if terminal_live is None:
        if len(live_rows) <= configuration.c:
            raise AssertionError("conditioned sampler produced no overflow")
    elif len(live_rows) != terminal_live:
        raise AssertionError("conditioned sampler produced wrong terminal state")
    return candidates, tuple(path)


def matrix_collision_signature(
    candidates: list[int], configuration: Configuration
) -> tuple[int, ...]:
    t, w = configuration.t, configuration.w
    return tuple(
        t - len({candidates[row * w + partition] for row in range(t)})
        for partition in range(w)
    )


def log2_or_negative_infinity(value: float) -> float:
    return NEG_INF if value == 0.0 else math.log2(value)


def validate_args(args: argparse.Namespace) -> None:
    if args.t < 1 or args.w < 2 or args.d < 1 or args.c < 0:
        raise SystemExit("require t>=1, w>=2, d>=1, and c>=0")
    if args.samples < 2:
        raise SystemExit("--samples must be at least two")
    if any(value < 1 for value in args.evictions):
        raise SystemExit("every --evictions value must be positive")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument("--w", type=int, default=3)
    parser.add_argument("--d", type=int, default=32)
    parser.add_argument("--c", type=int, default=2)
    parser.add_argument("--evictions", type=int, nargs="+", default=[4, 8])
    parser.add_argument(
        "--policy",
        choices=("biased", "uniform", "bfs"),
        default="biased",
    )
    parser.add_argument("--samples", type=int, default=200_000)
    parser.add_argument("--seed", type=int, default=0xC01D17)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--top", type=int, default=12)
    parser.add_argument(
        "--terminal-live",
        type=int,
        help="condition on exactly this many final live rows instead of L_w>c",
    )
    args = parser.parse_args()
    validate_args(args)

    configuration = Configuration(args.t, args.w, args.d, args.c)
    transitions = occupancy_transitions(args.t, args.d)
    if args.terminal_live is not None and not 0 <= args.terminal_live <= args.t:
        raise SystemExit("--terminal-live must lie in [0,t]")
    backward = backward_event_probabilities(
        configuration, transitions, args.terminal_live
    )
    exact_event = backward[0][args.t]
    if args.terminal_live is None:
        exact_overflow = exact_waterfall_failure(configuration)
        if not math.isclose(exact_event, exact_overflow, rel_tol=1e-12):
            raise AssertionError("backward recurrence disagrees with exact waterfall tail")
        conditioning_event = f"L_w>{args.c}"
    else:
        conditioning_event = f"L_w={args.terminal_live}"

    experiments = tuple(
        Experiment(args.policy, value, "global")
        for value in sorted(set(args.evictions))
    )
    rng = random.Random(args.seed)
    moments = {experiment: Moments() for experiment in experiments}
    nonzero_matrices = Counter({experiment: 0 for experiment in experiments})
    structural_matrices = Counter({experiment: 0 for experiment in experiments})
    sampled_paths: Counter[tuple[int, ...]] = Counter()
    path_contribution = {
        experiment: defaultdict(float) for experiment in experiments
    }
    signature_contribution = {
        experiment: defaultdict(float) for experiment in experiments
    }

    for _ in range(args.samples):
        candidates, path = sample_conditioned_matrix(
            configuration,
            transitions,
            backward,
            rng,
            args.terminal_live,
        )
        sampled_paths[path] += 1
        owner, placement, unplaced = waterfall(candidates, configuration)
        if len(unplaced) != path[-1]:
            raise AssertionError("router disagrees with conditioned live path")
        signature: tuple[int, ...] | None = None

        for experiment in experiments:
            if experiment.policy == "bfs":
                conditional_failure = float(
                    bounded_bfs_repair(
                        candidates,
                        configuration,
                        owner,
                        placement,
                        unplaced,
                        experiment.evictions,
                    )
                    > configuration.c
                )
            else:
                conditional_failure = exact_repair_failure_probability(
                    candidates,
                    configuration,
                    owner,
                    placement,
                    unplaced,
                    experiment,
                )
            moments[experiment].add(conditional_failure)
            if conditional_failure > 0.0:
                nonzero_matrices[experiment] += 1
                path_contribution[experiment][path] += conditional_failure
                if signature is None:
                    signature = matrix_collision_signature(candidates, configuration)
                signature_contribution[experiment][signature] += conditional_failure
            if conditional_failure == 1.0:
                structural_matrices[experiment] += 1

    rows: list[dict[str, int | float | str]] = []
    for experiment in experiments:
        conditional = moments[experiment].mean
        conditional_se = math.sqrt(moments[experiment].variance_of_mean)
        probability = exact_event * conditional
        standard_error = exact_event * conditional_se
        relative_error = math.inf if probability == 0.0 else standard_error / probability
        row: dict[str, int | float | str] = {
            "t": args.t,
            "w": args.w,
            "d": args.d,
            "c": args.c,
            "policy": args.policy,
            "evictions": experiment.evictions,
            "conditioned_samples": args.samples,
            "conditioning_event": conditioning_event,
            "exact_initial_overflow": exact_event,
            "exact_initial_overflow_log2": log2_or_negative_infinity(exact_event),
            "conditional_failure": conditional,
            "conditional_failure_log2": log2_or_negative_infinity(conditional),
            "failure_estimate": probability,
            "failure_estimate_log2": log2_or_negative_infinity(probability),
            "standard_error": standard_error,
            "relative_standard_error": relative_error,
            "log2_standard_error": (
                math.inf if probability == 0.0 else relative_error / math.log(2.0)
            ),
            "failure_effective_sample_size": moments[experiment].effective_sample_size,
            "nonzero_failure_matrices": nonzero_matrices[experiment],
            "structural_failure_matrices": structural_matrices[experiment],
        }
        rows.append(row)

        print(
            f"e={experiment.evictions} conditional={conditional:.8e} "
            f"failure={probability:.8e} "
            f"log2={log2_or_negative_infinity(probability):.4f} "
            f"rel_se={relative_error:.3%} "
            f"failure_ess={moments[experiment].effective_sample_size:.1f} "
            f"nonzero_matrices={nonzero_matrices[experiment]} "
            f"structural={structural_matrices[experiment]}"
        )
        total_contribution = moments[experiment].total
        print("  top live-count paths by failure contribution:")
        for path, contribution in sorted(
            path_contribution[experiment].items(),
            key=lambda item: item[1],
            reverse=True,
        )[: args.top]:
            print(f"    {path}: {contribution / total_contribution:.3%}")
        print("  top collision signatures by failure contribution:")
        for signature, contribution in sorted(
            signature_contribution[experiment].items(),
            key=lambda item: item[1],
            reverse=True,
        )[: args.top]:
            print(f"    {signature}: {contribution / total_contribution:.3%}")

    print(
        f"exact conditioning event {conditioning_event}={exact_event:.8e} "
        f"log2={log2_or_negative_infinity(exact_event):.4f}"
    )
    print("sampled conditioned terminal live counts:")
    terminal_counts: Counter[int] = Counter()
    for path, count in sampled_paths.items():
        terminal_counts[path[-1]] += count
    for live, count in sorted(terminal_counts.items()):
        print(f"  L_w={live}: {count / args.samples:.6%}")

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
