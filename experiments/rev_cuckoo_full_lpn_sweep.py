"""Reproducible sequential sweep for the exact toy full-LPN leakage game."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import random
import time
from dataclasses import dataclass
from pathlib import Path

from rev_cuckoo_full_lpn import (
    TrialResult,
    mean,
    one_trial,
    quantile,
    self_test,
    standard_error,
)


@dataclass(frozen=True)
class Config:
    polys: int
    weight: int
    block_length: int
    partition_size: int
    trials: int


CONFIGS = (
    # Block-length scaling at the current d=t shape.
    Config(2, 2, 2, 2, 4096),
    Config(2, 2, 3, 2, 2048),
    Config(2, 2, 4, 2, 1024),
    Config(2, 2, 5, 2, 512),
    Config(2, 2, 6, 2, 256),
    Config(2, 2, 8, 2, 64),
    # Partition-size scaling for the shared reference instance.
    Config(2, 2, 4, 3, 1024),
    Config(2, 2, 4, 4, 1024),
    Config(2, 2, 4, 8, 1024),
    # Weight and number-of-polynomials scaling at enumerable sizes.
    Config(2, 3, 2, 3, 2048),
    Config(2, 3, 3, 3, 256),
    Config(2, 4, 2, 4, 512),
    Config(3, 2, 2, 2, 1024),
    Config(4, 2, 2, 2, 256),
)


FIELDNAMES = (
    "mode",
    "polys",
    "weight",
    "block_length",
    "partition_size",
    "lists",
    "candidates",
    "trials",
    "seed",
    "seconds",
    "total_mean",
    "total_se",
    "total_p99",
    "dist_mean",
    "dist_se",
    "dist_p99",
    "dec_mean",
    "dec_se",
    "dec_p99",
    "entropy_loss_mean",
    "entropy_loss_se",
    "min_entropy_loss_mean",
    "min_entropy_loss_se",
    "rank_gain_mean",
    "rank_gain_se",
    "search_mean_log_gain",
    "search_mean_log_gain_se",
    "search_expected_work_gain",
    "base_entropy_mean",
    "posterior_entropy_mean",
    "decomposition_error",
)


def summarize(
    config: Config, seed: int, seconds: float, candidates: int,
    results: list[TrialResult]
) -> dict[str, str | int]:
    total = [result.total_information for result in results]
    dist = [result.distinguishing_information for result in results]
    dec = [result.decoding_information for result in results]
    entropy_loss = [
        result.base_entropy - result.posterior_entropy for result in results
    ]
    min_entropy_loss = [
        result.base_entropy - result.posterior_min_entropy for result in results
    ]
    rank_gain = [
        math.log2(result.base_rank / result.posterior_rank) for result in results
    ]
    search_rank_gain = [
        math.log2(result.search_base_rank / result.search_rank)
        for result in results
    ]
    return {
        "mode": "binary-unit",
        "polys": config.polys,
        "weight": config.weight,
        "block_length": config.block_length,
        "partition_size": config.partition_size,
        "lists": config.polys * config.polys * config.weight,
        "candidates": candidates,
        "trials": config.trials,
        "seed": seed,
        "seconds": f"{seconds:.6f}",
        "total_mean": f"{mean(total):.9f}",
        "total_se": f"{standard_error(total):.9f}",
        "total_p99": f"{quantile(total, 0.99):.9f}",
        "dist_mean": f"{mean(dist):.9f}",
        "dist_se": f"{standard_error(dist):.9f}",
        "dist_p99": f"{quantile(dist, 0.99):.9f}",
        "dec_mean": f"{mean(dec):.9f}",
        "dec_se": f"{standard_error(dec):.9f}",
        "dec_p99": f"{quantile(dec, 0.99):.9f}",
        "entropy_loss_mean": f"{mean(entropy_loss):.9f}",
        "entropy_loss_se": f"{standard_error(entropy_loss):.9f}",
        "min_entropy_loss_mean": f"{mean(min_entropy_loss):.9f}",
        "min_entropy_loss_se": f"{standard_error(min_entropy_loss):.9f}",
        "rank_gain_mean": f"{mean(rank_gain):.9f}",
        "rank_gain_se": f"{standard_error(rank_gain):.9f}",
        "search_mean_log_gain": f"{mean(search_rank_gain):.9f}",
        "search_mean_log_gain_se": f"{standard_error(search_rank_gain):.9f}",
        "search_expected_work_gain": f"{-math.log2(mean([r.search_rank / r.search_base_rank for r in results])):.9f}",
        "base_entropy_mean": f"{mean([r.base_entropy for r in results]):.9f}",
        "posterior_entropy_mean": f"{mean([r.posterior_entropy for r in results]):.9f}",
        "decomposition_error": f"{mean(total) - mean(dist) - mean(dec):.3e}",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("rev_cuckoo_full_lpn_results.csv"),
    )
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument(
        "--quick", action="store_true", help="cap every row at 32 trials"
    )
    args = parser.parse_args()

    self_test()
    rows = []
    for index, original in enumerate(CONFIGS):
        config = Config(
            original.polys,
            original.weight,
            original.block_length,
            original.partition_size,
            min(original.trials, 32) if args.quick else original.trials,
        )
        seed = args.seed + index
        candidates = tuple(
            itertools.product(
                range(config.block_length),
                repeat=config.polys * config.weight,
            )
        )
        rng = random.Random(seed)
        start = time.perf_counter()
        results = [
            one_trial(
                candidates,
                config.polys,
                config.weight,
                config.block_length,
                config.partition_size,
                rng,
            )
            for _ in range(config.trials)
        ]
        seconds = time.perf_counter() - start
        row = summarize(config, seed, seconds, len(candidates), results)
        rows.append(row)
        print(
            f"row={index + 1}/{len(CONFIGS)} P={config.polys} "
            f"t={config.weight} L={config.block_length} "
            f"d={config.partition_size} trials={config.trials} "
            f"total={row['total_mean']} dec={row['dec_mean']} "
            f"seconds={row['seconds']}",
            flush=True,
        )

    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
