"""Measure pairwise surrogates for the exact w=2 matching-count likelihood."""

from __future__ import annotations

import argparse
import math
import random
import statistics
from collections.abc import Sequence

from rev_cuckoo_full_lpn import plant_descriptor_w2, placement_count_w2


def graph_statistics(
    active: Sequence[int],
    descriptor: tuple[Sequence[int], Sequence[int]],
    partition_size: int,
) -> tuple[int, int]:
    left, right = descriptor
    left_counts = [0] * partition_size
    right_counts = [0] * partition_size
    for item in active:
        left_counts[left[item]] += 1
        right_counts[right[item]] += 1
    collision_pairs = sum(count * (count - 1) // 2 for count in left_counts)
    collision_pairs += sum(count * (count - 1) // 2 for count in right_counts)
    occupancy_deficit = 2 * len(active) - sum(count != 0 for count in left_counts) - sum(
        count != 0 for count in right_counts
    )
    return collision_pairs, occupancy_deficit


def linear_fit(feature: Sequence[float], target: Sequence[float]) -> tuple[float, float, float, float]:
    feature_mean = statistics.mean(feature)
    target_mean = statistics.mean(target)
    variance = sum((value - feature_mean) ** 2 for value in feature)
    covariance = sum(
        (value - feature_mean) * (result - target_mean)
        for value, result in zip(feature, target)
    )
    slope = covariance / variance
    intercept = target_mean - slope * feature_mean
    residuals = [
        result - (intercept + slope * value)
        for value, result in zip(feature, target)
    ]
    total = sum((result - target_mean) ** 2 for result in target)
    residual = sum(value * value for value in residuals)
    r_squared = 1.0 - residual / total
    return intercept, slope, r_squared, math.sqrt(residual / len(residuals))


def sample(
    active_count: int,
    partition_size: int,
    trials: int,
    planted: bool,
    rng: random.Random,
) -> None:
    active = tuple(range(active_count))
    collision_pairs = []
    occupancy_deficits = []
    log_placements = []
    failures = 0
    for _ in range(trials):
        if planted:
            descriptor = plant_descriptor_w2(
                active, active_count, partition_size, rng
            )
        else:
            descriptor = (
                tuple(rng.randrange(partition_size) for _ in active),
                tuple(rng.randrange(partition_size) for _ in active),
            )
        placements = placement_count_w2(active, descriptor, partition_size)
        if placements == 0:
            failures += 1
            continue
        pairs, deficit = graph_statistics(active, descriptor, partition_size)
        collision_pairs.append(pairs)
        occupancy_deficits.append(deficit)
        log_placements.append(math.log2(placements))

    print(
        f"distribution={'planted' if planted else 'uniform'} "
        f"n={active_count} d={partition_size} trials={trials} "
        f"failures={failures} feasible={len(log_placements)}"
    )
    for name, feature in (
        ("collision_pairs", collision_pairs),
        ("occupancy_deficit", occupancy_deficits),
    ):
        intercept, slope, r_squared, residual_stddev = linear_fit(
            feature, log_placements
        )
        print(
            f"{name}: intercept={intercept:.6f} slope={slope:.6f} "
            f"r2={r_squared:.6f} residual_stddev={residual_stddev:.6f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=16)
    parser.add_argument("--d", type=int, default=16)
    parser.add_argument("--trials", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if min(args.n, args.d, args.trials) < 1:
        parser.error("n, d, and trials must be positive")
    rng = random.Random(args.seed)
    sample(args.n, args.d, args.trials, False, rng)
    sample(args.n, args.d, args.trials, True, rng)


if __name__ == "__main__":
    main()
