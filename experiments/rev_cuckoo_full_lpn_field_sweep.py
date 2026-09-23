"""Sequential coefficient-aware sweep for the exact toy Ring-LPN game."""

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
    polys: int
    weight: int
    block_length: int
    partition_size: int
    modulus: int
    trials: int
    channel: str = "raw"
    selector_k: int = 1


CONFIGS = (
    Config(2, 2, 2, 2, 3, 2048),
    Config(2, 2, 3, 2, 3, 512),
    Config(2, 2, 4, 2, 3, 128),
    Config(2, 2, 3, 3, 3, 512),
    Config(2, 2, 3, 4, 3, 512),
    Config(2, 2, 3, 8, 3, 512),
    Config(2, 2, 2, 2, 5, 128),
    Config(2, 3, 2, 3, 3, 128),
    Config(3, 2, 2, 2, 3, 128),
    Config(4, 2, 2, 2, 3, 8),
)


FIELDNAMES = (
    "mode",
    "channel",
    "selector_k",
    "modulus",
    "polys",
    "weight",
    "block_length",
    "partition_size",
    "lists",
    "position_candidates",
    "full_candidates",
    "trials",
    "seed",
    "seconds",
    "total_mean",
    "total_se",
    "total_p99",
    "dist_mean",
    "dist_se",
    "dec_mean",
    "dec_se",
    "entropy_loss_mean",
    "min_entropy_loss_mean",
    "posterior_mean_log_rank_gain",
    "decoding_expected_work_gain",
    "search_mean_log_gain",
    "search_expected_work_gain",
    "base_entropy_mean",
    "decomposition_error",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("rev_cuckoo_full_lpn_field_results.csv"),
    )
    parser.add_argument("--seed", type=int, default=20260831)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    self_test()

    rows = []
    for index, original in enumerate(CONFIGS):
        config = Config(
            original.polys,
            original.weight,
            original.block_length,
            original.partition_size,
            original.modulus,
            min(original.trials, 16) if args.quick else original.trials,
            original.channel,
            original.selector_k,
        )
        seed = args.seed + index
        model = make_model(
            config.polys, config.weight, config.block_length, config.modulus
        )
        rng = random.Random(seed)
        channel = LeakageChannel.build(
            config.channel,
            config.selector_k,
            config.weight,
            config.partition_size,
        )
        start = time.perf_counter()
        results = [
            field_trial(model, channel, rng)
            for _ in range(config.trials)
        ]
        seconds = time.perf_counter() - start
        total = [result.total_information for result in results]
        dist = [result.distinguishing_information for result in results]
        dec = [result.decoding_information for result in results]
        entropy_loss = [
            result.base_entropy - result.posterior_entropy for result in results
        ]
        min_entropy_loss = [
            result.base_entropy - result.posterior_min_entropy
            for result in results
        ]
        posterior_gain = [
            math.log2(result.base_rank / result.posterior_rank)
            for result in results
        ]
        search_gain = [
            math.log2(result.search_base_rank / result.search_rank)
            for result in results
        ]
        decoding_expected_work_gain = -math.log2(
            mean(
                [
                    result.posterior_rank / result.base_rank
                    for result in results
                ]
            )
        )
        expected_work_gain = -math.log2(
            mean(
                [
                    result.search_rank / result.search_base_rank
                    for result in results
                ]
            )
        )
        row = {
            "mode": "prime-field",
            "channel": config.channel,
            "selector_k": config.selector_k,
            "modulus": config.modulus,
            "polys": config.polys,
            "weight": config.weight,
            "block_length": config.block_length,
            "partition_size": config.partition_size,
            "lists": config.polys * config.polys * config.weight,
            "position_candidates": len(model.position_candidates),
            "full_candidates": len(model.full_states),
            "trials": config.trials,
            "seed": seed,
            "seconds": f"{seconds:.6f}",
            "total_mean": f"{mean(total):.9f}",
            "total_se": f"{standard_error(total):.9f}",
            "total_p99": f"{quantile(total, 0.99):.9f}",
            "dist_mean": f"{mean(dist):.9f}",
            "dist_se": f"{standard_error(dist):.9f}",
            "dec_mean": f"{mean(dec):.9f}",
            "dec_se": f"{standard_error(dec):.9f}",
            "entropy_loss_mean": f"{mean(entropy_loss):.9f}",
            "min_entropy_loss_mean": f"{mean(min_entropy_loss):.9f}",
            "posterior_mean_log_rank_gain": f"{mean(posterior_gain):.9f}",
            "decoding_expected_work_gain": f"{decoding_expected_work_gain:.9f}",
            "search_mean_log_gain": f"{mean(search_gain):.9f}",
            "search_expected_work_gain": f"{expected_work_gain:.9f}",
            "base_entropy_mean": f"{mean([r.base_entropy for r in results]):.9f}",
            "decomposition_error": f"{mean(total) - mean(dist) - mean(dec):.3e}",
        }
        rows.append(row)
        print(
            f"row={index + 1}/{len(CONFIGS)} q={config.modulus} "
            f"P={config.polys} t={config.weight} L={config.block_length} "
            f"d={config.partition_size} states={len(model.full_states)} "
            f"channel={config.channel} K={config.selector_k} "
            f"trials={config.trials} decode={decoding_expected_work_gain:.6f} "
            f"search={expected_work_gain:.6f} "
            f"seconds={seconds:.6f}",
            flush=True,
        )

    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
