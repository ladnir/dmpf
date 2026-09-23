"""Sequential reduced sweep of coefficient-aware basic Ring-LPN Prange."""

from __future__ import annotations

import csv
import itertools
import math
import random
import time
from dataclasses import dataclass
from pathlib import Path

from rev_cuckoo_basic_ring_lpn_isd import (
    best_mixed_expected_work_gain,
    best_weighted_expected_work_gain,
    best_weighted_half_success_gain,
    fit_pair_slope,
    mean,
    one_trial,
    standard_error,
    states_by_support,
)
from rev_cuckoo_full_lpn_field import make_model
from rev_cuckoo_regular_noise import LeakageChannel


@dataclass(frozen=True)
class Config:
    name: str
    weight: int
    block_length: int
    trials: int


CONFIGS = (
    Config("t2-L4", 2, 4, 256),
    Config("t3-L4", 3, 4, 128),
    Config("t3-L8", 3, 8, 64),
    Config("t4-L4", 4, 4, 64),
)
SELECTOR_COUNTS = (1, 2, 4, 8)
SELECTORS = (
    ("marginal", "marginal_support", "marginal_rank"),
    ("pair", "pair_support", "pair_rank"),
    ("joint", "joint_support", "joint_rank"),
)


def main() -> None:
    output = Path(__file__).with_name(
        "rev_cuckoo_basic_ring_lpn_isd_results.csv"
    )
    field_modulus = 3
    rank_modulus = 65537
    seed = 20260810
    rows = []
    for config_index, config in enumerate(CONFIGS):
        model = make_model(
            2,
            config.weight,
            config.block_length,
            field_modulus,
            enumerate_full_states=False,
        )
        candidates = tuple(
            itertools.product(
                range(config.block_length), repeat=config.weight
            )
        )
        grouped_states = states_by_support(model)
        for selector_k in SELECTOR_COUNTS:
            channel = LeakageChannel.build(
                "occupancy",
                selector_k,
                config.weight,
                config.weight,
            )
            pair_slope, pair_r_squared = fit_pair_slope(
                channel,
                config.weight,
                5000,
                random.Random(seed + config_index * 100 + selector_k),
            )
            start = time.perf_counter()
            results = [
                one_trial(
                    model,
                    candidates,
                    grouped_states,
                    channel,
                    joint_restarts=2,
                    pair_slope=pair_slope,
                    pair_restarts=2,
                    pair_rounds=16,
                    pair_damping=0.5,
                    rng=random.Random(
                        seed + config_index * 1_000_000 + trial_index
                    ),
                    rank_modulus=rank_modulus,
                )
                for trial_index in range(config.trials)
            ]
            seconds = time.perf_counter() - start
            base_support = 2.0 ** (-2 * config.weight)
            base_rank = mean([result.base_rank for result in results])
            base_success = base_support * base_rank
            for selector, support_field, rank_field in SELECTORS:
                support = [
                    getattr(result, support_field) for result in results
                ]
                ranks = [getattr(result, rank_field) for result in results]
                success = [
                    probability * rank_ok
                    for probability, rank_ok in zip(support, ranks)
                ]
                support_mean = mean(support)
                success_mean = mean(success)
                mixed_gain, mixed_alpha = best_mixed_expected_work_gain(
                    [
                        getattr(result, f"{selector}_hits")
                        for result in results
                    ],
                    2 * config.weight,
                    tuple(index / 100.0 for index in range(0, 100, 5)),
                )
                if selector in ("marginal", "pair"):
                    weighted_successes = [
                        getattr(
                            result,
                            f"{selector}_weighted_success",
                        )
                        for result in results
                    ]
                    weighted_gain, weighted_beta = (
                        best_weighted_expected_work_gain(
                            weighted_successes,
                            2 * config.weight,
                        )
                    )
                    weighted_half_gain, weighted_half_beta = (
                        best_weighted_half_success_gain(
                            weighted_successes,
                            2 * config.weight,
                        )
                    )
                else:
                    weighted_gain, weighted_beta = math.nan, math.nan
                    weighted_half_gain, weighted_half_beta = math.nan, math.nan
                rows.append(
                    {
                        "name": config.name,
                        "selector_k": selector_k,
                        "selector": selector,
                        "field_modulus": field_modulus,
                        "rank_modulus": rank_modulus,
                        "weight": config.weight,
                        "block_length": config.block_length,
                        "partition_size": config.weight,
                        "trials": config.trials,
                        "pair_slope": f"{pair_slope:.9f}",
                        "pair_r_squared": f"{pair_r_squared:.9f}",
                        "base_support": f"{base_support:.9f}",
                        "base_rank": f"{base_rank:.9f}",
                        "support": f"{support_mean:.9f}",
                        "support_se": f"{standard_error(support):.9f}",
                        "support_gain_bits": (
                            f"{math.log2(support_mean / base_support):.9f}"
                        ),
                        "rank": f"{mean(ranks):.9f}",
                        "success": f"{success_mean:.9f}",
                        "success_se": f"{standard_error(success):.9f}",
                        "gain_bits": (
                            f"{math.log2(success_mean / base_success):.9f}"
                        ),
                        "mixed_expected_work_gain_bits": f"{mixed_gain:.9f}",
                        "mixed_alpha": f"{mixed_alpha:.2f}",
                        "weighted_expected_work_gain_bits": (
                            f"{weighted_gain:.9f}"
                        ),
                        "weighted_beta": f"{weighted_beta:.2f}",
                        "weighted_half_success_gain_bits": (
                            f"{weighted_half_gain:.9f}"
                        ),
                        "weighted_half_beta": f"{weighted_half_beta:.2f}",
                        "seconds": f"{seconds:.6f}",
                        "seed": seed,
                    }
                )
            print(
                f"{config.name} K={selector_k} trials={config.trials} "
                + " ".join(
                    f"{selector}={rows[-3 + index]['gain_bits']}"
                    for index, (selector, _, _) in enumerate(SELECTORS)
                )
                + f" seconds={seconds:.3f}",
                flush=True,
            )

    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
