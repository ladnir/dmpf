"""Scalable pair-weighted Prange diagnostic for regular Ring-LPN.

This script does not enumerate a posterior.  It samples the complete ideal
coefficient-deduplicated product leakage for ``P`` known and ``P`` hidden regular
polynomials, runs polynomial collision-aware selectors, and evaluates the
planted support under repeatable fixed-cardinality information-set samplers.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import time
from pathlib import Path

from rev_cuckoo_basic_ring_lpn_isd import (
    BALANCED_BETAS,
    balanced_stratified_coordinate_probabilities,
    best_weighted_expected_work_gain,
    best_weighted_half_success_gain,
    iterations_for_average_success,
)
from rev_cuckoo_full_lpn_field import FieldModel, PolyState, product_active_list
from rev_cuckoo_regular_isd import (
    collision_factors,
    fit_pair_slope,
    pair_center_coordinate_ascent,
    pair_mean_field_distributions,
)
from rev_cuckoo_regular_noise import LeakageChannel


MIXTURE_ALPHAS = (0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0)


def best_uniform_mixture_gains(
    successes: tuple[float, ...] | list[float], baseline_bits: float
) -> tuple[float, float, float, float]:
    """Optimize a fixed mixture of a selector and uniform regular Prange."""

    baseline = 2.0 ** (-baseline_bits)
    baseline_iterations = math.log(0.5) / math.log1p(-baseline)
    best_expected = (-math.inf, 0.0)
    best_half = (-math.inf, 0.0)
    for alpha in MIXTURE_ALPHAS:
        mixed = [
            (1.0 - alpha) * baseline + alpha * success
            for success in successes
        ]
        expected_work = sum(1.0 / success for success in mixed) / len(mixed)
        expected_gain = math.log2((1.0 / baseline) / expected_work)
        if expected_gain > best_expected[0]:
            best_expected = (expected_gain, alpha)
        iterations = iterations_for_average_success(mixed)
        half_gain = math.log2(baseline_iterations / iterations)
        if half_gain > best_half[0]:
            best_half = (half_gain, alpha)
    return best_expected[0], best_expected[1], best_half[0], best_half[1]


def sample_poly(
    weight: int,
    block_length: int,
    modulus: int,
    rng: random.Random,
) -> PolyState:
    return PolyState(
        tuple(rng.randrange(block_length) for _ in range(weight)),
        tuple(rng.randrange(1, modulus) for _ in range(weight)),
        (),
    )


def active_product_lists(
    model: FieldModel,
    known: tuple[PolyState, ...],
    hidden: tuple[PolyState, ...],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        product_active_list(model, known_poly, hidden_poly, output_block)
        for known_poly in known
        for hidden_poly in hidden
        for output_block in range(model.weight)
    )


def factor_descriptors(
    descriptors, hidden_index: int, polys: int, weight: int
):
    return tuple(
        descriptors[
            (known_index * polys + hidden_index) * weight + output_block
        ]
        for known_index in range(polys)
        for output_block in range(weight)
    )


def one_trial(
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    pair_rounds: int,
    center_rounds: int,
    rng: random.Random,
) -> tuple[tuple[float, ...], int, float]:
    known = tuple(
        sample_poly(
            model.weight, model.block_length, model.modulus, rng
        )
        for _ in range(model.polys)
    )
    hidden = tuple(
        sample_poly(
            model.weight, model.block_length, model.modulus, rng
        )
        for _ in range(model.polys)
    )
    active_lists = active_product_lists(model, known, hidden)
    descriptors = tuple(
        channel.sample(active, 2 * model.block_length, rng)
        for active in active_lists
    )
    known_positions = tuple(
        position for state in known for position in state.positions
    )
    kept = model.block_length // model.polys
    successes = [1.0] * len(BALANCED_BETAS)
    center_hits = 0
    for hidden_index, hidden_poly in enumerate(hidden):
        hidden_descriptors = factor_descriptors(
            descriptors, hidden_index, model.polys, model.weight
        )
        factors = collision_factors(
            known_positions,
            hidden_descriptors,
            polys=model.polys,
            weight=model.weight,
            block_length=model.block_length,
        )
        distributions = pair_mean_field_distributions(
            known_positions,
            hidden_descriptors,
            polys=model.polys,
            weight=model.weight,
            block_length=model.block_length,
            kept=kept,
            slope=pair_slope,
            restarts=0,
            rounds=pair_rounds,
            damping=0.5,
            rng=rng,
            factors=factors,
        )
        for block, position in enumerate(hidden_poly.positions):
            inclusions = balanced_stratified_coordinate_probabilities(
                distributions[block], position, kept, BALANCED_BETAS
            )
            for beta_index, inclusion in enumerate(inclusions):
                successes[beta_index] *= inclusion
        initial_center = tuple(
            max(range(model.block_length), key=distribution.__getitem__)
            for distribution in distributions
        )
        center = pair_center_coordinate_ascent(
            factors, initial_center, pair_slope, center_rounds
        )
        center_hits += sum(
            guess == position
            for guess, position in zip(center, hidden_poly.positions)
        )
    filler_probability = (
        (kept - 1) / (model.block_length - 1)
    )
    center_success = filler_probability ** (
        model.polys * model.weight - center_hits
    )
    return tuple(successes), center_hits, center_success


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=2)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=4096)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--modulus", type=int, default=65537)
    parser.add_argument("--selector-k", type=int, default=1)
    parser.add_argument("--trials", type=int, default=8)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--center-rounds", type=int, default=2)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trial-output", type=Path)
    args = parser.parse_args()
    if args.polys < 2:
        parser.error("the number of polynomials must be at least two")
    if args.block_length % args.polys:
        parser.error("block length must be divisible by the number of polynomials")
    partition_size = args.partition_size or args.weight
    model = FieldModel(
        polys=args.polys,
        weight=args.weight,
        block_length=args.block_length,
        modulus=args.modulus,
        poly_states=(),
        full_states=(),
        position_candidates=(),
    )
    channel = LeakageChannel.build(
        "occupancy",
        args.selector_k,
        args.weight,
        partition_size,
    )
    pair_slope, pair_r_squared = fit_pair_slope(
        channel,
        args.weight,
        args.pair_fit_samples,
        random.Random(args.seed ^ 0xA5A5A5A5),
    )
    start = time.perf_counter()
    trials = [
        one_trial(
            model,
            channel,
            pair_slope,
            args.pair_rounds,
            args.center_rounds,
            random.Random(args.seed + trial),
        )
        for trial in range(args.trials)
    ]
    successes = [trial[0] for trial in trials]
    center_hits = [trial[1] for trial in trials]
    seconds = time.perf_counter() - start
    support_work_bits = (
        args.polys * args.weight * math.log2(args.polys)
    )
    expected_gain, expected_beta = best_weighted_expected_work_gain(
        successes,
        support_work_bits,
        BALANCED_BETAS,
    )
    half_gain, half_beta = best_weighted_half_success_gain(
        successes,
        support_work_bits,
        BALANCED_BETAS,
    )
    (
        center_expected_gain,
        center_expected_alpha,
        center_half_gain,
        center_half_alpha,
    ) = best_uniform_mixture_gains(
        [trial[2] for trial in trials], support_work_bits
    )
    row = {
        "polys": args.polys,
        "weight": args.weight,
        "block_length": args.block_length,
        "ring_dimension": args.weight * args.block_length,
        "kept_per_block": args.block_length // args.polys,
        "support_work_bits": support_work_bits,
        "partition_size": partition_size,
        "modulus": args.modulus,
        "selector_k": args.selector_k,
        "trials": args.trials,
        "pair_rounds": args.pair_rounds,
        "center_rounds": args.center_rounds,
        "pair_slope": pair_slope,
        "pair_r_squared": pair_r_squared,
        "expected_work_gain_bits": expected_gain,
        "expected_beta": expected_beta,
        "half_success_gain_bits": half_gain,
        "half_beta": half_beta,
        "center_mean_hits": sum(center_hits) / len(center_hits),
        "center_expected_work_gain_bits": center_expected_gain,
        "center_expected_alpha": center_expected_alpha,
        "center_half_success_gain_bits": center_half_gain,
        "center_half_alpha": center_half_alpha,
        "seconds": seconds,
        "seed": args.seed,
    }
    print(" ".join(f"{key}={value}" for key, value in row.items()))
    for beta_index, beta in enumerate(BALANCED_BETAS):
        beta_successes = [
            (trial[0][beta_index],) for trial in trials
        ]
        beta_expected, _ = best_weighted_expected_work_gain(
            beta_successes, support_work_bits, (beta,)
        )
        beta_half, _ = best_weighted_half_success_gain(
            beta_successes, support_work_bits, (beta,)
        )
        print(
            f"balanced_beta={beta} "
            f"expected_work_gain_bits={beta_expected} "
            f"half_success_gain_bits={beta_half}"
        )
    if args.output:
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=tuple(row))
            writer.writeheader()
            writer.writerow(row)
        print(f"wrote {args.output}")
    if args.trial_output:
        fields = (
            "trial",
            *(f"balanced_beta_{beta:g}" for beta in BALANCED_BETAS),
            "center_hits",
            "center_success",
        )
        with args.trial_output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for trial_index, trial in enumerate(trials):
                trial_row = {
                    "trial": trial_index,
                    **{
                        f"balanced_beta_{beta:g}": trial[0][beta_index]
                        for beta_index, beta in enumerate(BALANCED_BETAS)
                    },
                    "center_hits": trial[1],
                    "center_success": trial[2],
                }
                writer.writerow(trial_row)
        print(f"wrote {args.trial_output}")


if __name__ == "__main__":
    main()
