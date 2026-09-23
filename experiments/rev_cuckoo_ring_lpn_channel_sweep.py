"""Compare finite-K Reverse-Cuckoo channels in exact basic Ring-LPN toys.

Every row conditions on the complete coefficient-aware syndrome
``f + r*e`` over ``F_q[X]/(X^N+1)``.  The experiment enumerates all regular
supports and all nonzero coefficients, then reports how much the descriptor
improves an optimal ranking of syndrome-compatible ``(f,e)`` candidates.

The same Ring-LPN instance seeds are reused for every K.  Descriptor randomness
is kept separate so that consuming K proposals does not change later public
Ring-LPN instances.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import time
from dataclasses import dataclass
from pathlib import Path

from rev_cuckoo_full_lpn import mean, quantile, standard_error
from rev_cuckoo_full_lpn_field import field_trial, make_model, self_test
from rev_cuckoo_regular_noise import LeakageChannel


@dataclass(frozen=True)
class Config:
    name: str
    weight: int
    block_length: int
    partition_size: int
    modulus: int
    trials: int


CONFIGS = (
    Config("base", 2, 2, 2, 3, 4096),
    Config("longer-block", 2, 3, 2, 3, 1024),
    Config("larger-field", 2, 2, 2, 5, 512),
    Config("larger-weight", 3, 2, 3, 3, 512),
)

SELECTOR_COUNTS = (1, 2, 4, 8)

FIELDNAMES = (
    "name",
    "leakage_input",
    "channel",
    "selector_k",
    "modulus",
    "polys",
    "weight",
    "block_length",
    "ring_dimension",
    "partition_size",
    "lists",
    "position_candidates",
    "full_candidates",
    "trials",
    "seed",
    "seconds",
    "conditional_information_mean",
    "conditional_information_se",
    "conditional_information_p99",
    "posterior_entropy_loss_mean",
    "posterior_min_entropy_loss_mean",
    "mean_log_rank_gain",
    "expected_rank_gain",
    "support_mean_log_rank_gain",
    "support_expected_rank_gain",
    "base_entropy_mean",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
    )
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--leakage-input",
        choices=("support", "product-values"),
        default="product-values",
    )
    args = parser.parse_args()
    self_test()
    output = args.output or Path(__file__).with_name(
        "rev_cuckoo_ring_lpn_product_results.csv"
        if args.leakage_input == "product-values"
        else "rev_cuckoo_ring_lpn_channel_results.csv"
    )

    rows = []
    for config_index, config in enumerate(CONFIGS):
        trials = min(config.trials, 16) if args.quick else config.trials
        model = make_model(2, config.weight, config.block_length, config.modulus)
        for selector_k in SELECTOR_COUNTS:
            channel = LeakageChannel.build(
                "occupancy",
                selector_k,
                config.weight,
                config.partition_size,
            )
            start = time.perf_counter()
            results = []
            for trial_index in range(trials):
                instance_seed = args.seed + config_index * 1_000_000 + trial_index
                descriptor_seed = instance_seed + selector_k * 100_000_000
                results.append(
                    field_trial(
                        model,
                        channel,
                        random.Random(instance_seed),
                        random.Random(descriptor_seed),
                        product_values=args.leakage_input == "product-values",
                    )
                )
            seconds = time.perf_counter() - start

            conditional_information = [
                result.decoding_information for result in results
            ]
            entropy_loss = [
                result.base_entropy - result.posterior_entropy
                for result in results
            ]
            min_entropy_loss = [
                result.base_entropy - result.posterior_min_entropy
                for result in results
            ]
            mean_log_rank_gain = mean(
                [
                    math.log2(result.base_rank / result.posterior_rank)
                    for result in results
                ]
            )
            expected_rank_gain = -math.log2(
                mean(
                    [
                        result.posterior_rank / result.base_rank
                        for result in results
                    ]
                )
            )
            support_mean_log_rank_gain = mean(
                [
                    math.log2(
                        result.search_base_rank / result.search_rank
                    )
                    for result in results
                ]
            )
            support_expected_rank_gain = -math.log2(
                mean(
                    [
                        result.search_rank / result.search_base_rank
                        for result in results
                    ]
                )
            )
            row = {
                "name": config.name,
                "leakage_input": args.leakage_input,
                "channel": "occupancy",
                "selector_k": selector_k,
                "modulus": config.modulus,
                "polys": 2,
                "weight": config.weight,
                "block_length": config.block_length,
                "ring_dimension": config.weight * config.block_length,
                "partition_size": config.partition_size,
                "lists": 4 * config.weight,
                "position_candidates": len(model.position_candidates),
                "full_candidates": len(model.full_states),
                "trials": trials,
                "seed": args.seed,
                "seconds": f"{seconds:.6f}",
                "conditional_information_mean": (
                    f"{mean(conditional_information):.9f}"
                ),
                "conditional_information_se": (
                    f"{standard_error(conditional_information):.9f}"
                ),
                "conditional_information_p99": (
                    f"{quantile(conditional_information, 0.99):.9f}"
                ),
                "posterior_entropy_loss_mean": f"{mean(entropy_loss):.9f}",
                "posterior_min_entropy_loss_mean": (
                    f"{mean(min_entropy_loss):.9f}"
                ),
                "mean_log_rank_gain": f"{mean_log_rank_gain:.9f}",
                "expected_rank_gain": f"{expected_rank_gain:.9f}",
                "support_mean_log_rank_gain": (
                    f"{support_mean_log_rank_gain:.9f}"
                ),
                "support_expected_rank_gain": (
                    f"{support_expected_rank_gain:.9f}"
                ),
                "base_entropy_mean": (
                    f"{mean([result.base_entropy for result in results]):.9f}"
                ),
            }
            rows.append(row)
            print(
                f"{config.name} K={selector_k} trials={trials} "
                f"input={args.leakage_input} "
                f"I_dec={mean(conditional_information):.6f} "
                f"conditional-rank={expected_rank_gain:.6f} "
                f"support-rank={support_expected_rank_gain:.6f} "
                f"seconds={seconds:.3f}",
                flush=True,
            )

    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
