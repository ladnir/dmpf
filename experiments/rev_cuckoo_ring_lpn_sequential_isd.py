"""Sequential collision-conditioned regular-Prange attack.

The attack samples one fixed-cardinality block subset at a time.  Later block
beliefs are corrected using the collision-label mass of the subsets already
chosen.  A bootstrap particle filter evaluates the planted-support success by
sampling the attack conditioned on success so far, avoiding a 2^-128 rare-hit
experiment.  The planted support is used only by the evaluator, never by the
attack's transition rule.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from rev_cuckoo_basic_ring_lpn_isd import (
    balanced_stratified_coordinate_probabilities,
    best_weighted_expected_work_gain,
    best_weighted_half_success_gain,
    sample_balanced_stratified_subset,
    sample_balanced_stratified_subset_conditioned,
)
from rev_cuckoo_basic_ring_lpn_isd_production import (
    active_product_lists,
    factor_descriptors,
    sample_poly,
)
from rev_cuckoo_full_lpn_field import FieldModel
from rev_cuckoo_regular_isd import (
    CollisionFactor,
    collision_factors,
    fit_pair_slope,
    pair_mean_field_distributions,
)
from rev_cuckoo_regular_noise import LeakageChannel


@dataclass(frozen=True)
class TrialModel:
    true_positions: tuple[tuple[int, ...], ...]
    distributions: tuple[tuple[tuple[float, ...], ...], ...]
    factor_sets: tuple[tuple[CollisionFactor, ...], ...]
    known_positions: tuple[int, ...]
    descriptors: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]


def normalized_log2_weights(log_weights: np.ndarray) -> np.ndarray:
    maximum = float(np.max(log_weights))
    weights = np.exp2(log_weights - maximum)
    return weights / np.sum(weights)


def entropy(probabilities: np.ndarray) -> float:
    positive = probabilities[probabilities > 0.0]
    return float(-np.sum(positive * np.log2(positive)))


def systematic_resample(
    weights: list[float], count: int, rng: random.Random
) -> tuple[int, ...]:
    total = sum(weights)
    if total <= 0.0:
        raise ValueError("particle weights must have positive mass")
    step = total / count
    point = rng.random() * step
    result = []
    prefix = weights[0]
    index = 0
    for _ in range(count):
        while point > prefix and index + 1 < len(weights):
            index += 1
            prefix += weights[index]
        result.append(index)
        point += step
    return tuple(result)


def label_arrays(
    factors: tuple[CollisionFactor, ...],
) -> tuple[np.ndarray, int]:
    labels = np.asarray(
        [factor.labels for factor in factors], dtype=np.int32
    )
    label_count = max(factor.label_count for factor in factors)
    return labels, label_count


def add_distribution_masses(
    masses: np.ndarray,
    labels: np.ndarray,
    variable: int,
    probabilities: np.ndarray,
    label_count: int,
) -> None:
    for factor in range(labels.shape[0]):
        masses[factor] += np.bincount(
            labels[factor, variable],
            weights=probabilities,
            minlength=label_count,
        )


def add_subset_masses(
    masses: np.ndarray,
    labels: np.ndarray,
    variable: int,
    subset: frozenset[int],
    kept: int,
    label_count: int,
) -> None:
    coordinates = np.fromiter(subset, dtype=np.int64, count=kept)
    for factor in range(labels.shape[0]):
        masses[factor] += np.bincount(
            labels[factor, variable, coordinates],
            minlength=label_count,
        ) / kept


def corrected_distribution(
    base: np.ndarray,
    labels: np.ndarray,
    variable: int,
    particle_masses: np.ndarray,
    baseline_masses: np.ndarray,
    slope: float,
) -> np.ndarray:
    log_weights = np.log2(base)
    for factor in range(labels.shape[0]):
        delta = particle_masses[factor] - baseline_masses[factor]
        log_weights += slope * delta[labels[factor, variable]]
    return normalized_log2_weights(log_weights)


def sequential_polynomial_success(
    true_positions: tuple[int, ...],
    distributions: tuple[tuple[float, ...], ...],
    factors: tuple[CollisionFactor, ...],
    kept: int,
    beta: float,
    slope: float,
    particles: int,
    rng: random.Random,
) -> float:
    """Estimate one hidden polynomial's conditional support success."""

    base = np.asarray(distributions, dtype=np.float64)
    labels, label_count = label_arrays(factors)
    order = tuple(
        sorted(range(len(distributions)), key=lambda index: entropy(base[index]))
    )
    baseline_masses = np.zeros(
        (labels.shape[0], label_count), dtype=np.float64
    )
    particle_masses = [np.zeros_like(baseline_masses) for _ in range(particles)]
    log_success = 0.0
    for variable in order:
        children = []
        inclusion_probabilities = []
        for masses in particle_masses:
            conditioned = corrected_distribution(
                base[variable],
                labels,
                variable,
                masses,
                baseline_masses,
                slope=slope,
            )
            subset, inclusion = sample_balanced_stratified_subset_conditioned(
                conditioned,
                kept,
                beta,
                true_positions[variable],
                rng,
            )
            child = masses.copy()
            add_subset_masses(
                child,
                labels,
                variable,
                subset,
                kept,
                label_count,
            )
            children.append(child)
            inclusion_probabilities.append(inclusion)
        mean_inclusion = sum(inclusion_probabilities) / particles
        log_success += math.log2(mean_inclusion)
        indices = systematic_resample(
            inclusion_probabilities, particles, rng
        )
        particle_masses = [children[index].copy() for index in indices]
        add_distribution_masses(
            baseline_masses,
            labels,
            variable,
            base[variable],
            label_count,
        )
    return 2.0 ** log_success


def sample_sequential_polynomial_selection(
    distributions: tuple[tuple[float, ...], ...],
    factors: tuple[CollisionFactor, ...],
    kept: int,
    beta: float,
    slope: float,
    rng: random.Random,
) -> tuple[frozenset[int], ...]:
    """Sample the attack's actual correlated block subsets."""

    base = np.asarray(distributions, dtype=np.float64)
    labels, label_count = label_arrays(factors)
    order = tuple(
        sorted(range(len(distributions)), key=lambda index: entropy(base[index]))
    )
    baseline_masses = np.zeros(
        (labels.shape[0], label_count), dtype=np.float64
    )
    selected_masses = np.zeros_like(baseline_masses)
    selections: list[frozenset[int] | None] = [None] * len(distributions)
    for variable in order:
        conditioned = corrected_distribution(
            base[variable],
            labels,
            variable,
            selected_masses,
            baseline_masses,
            slope,
        )
        subset = sample_balanced_stratified_subset(
            conditioned, kept, beta, rng
        )
        selections[variable] = subset
        add_subset_masses(
            selected_masses,
            labels,
            variable,
            subset,
            kept,
            label_count,
        )
        add_distribution_masses(
            baseline_masses,
            labels,
            variable,
            base[variable],
            label_count,
        )
    if any(selection is None for selection in selections):
        raise AssertionError("sequential selector omitted a block")
    return tuple(
        selection for selection in selections if selection is not None
    )


def independent_polynomial_success(
    true_positions: tuple[int, ...],
    distributions: tuple[tuple[float, ...], ...],
    kept: int,
    beta: float,
) -> float:
    result = 1.0
    for position, distribution in zip(true_positions, distributions):
        result *= balanced_stratified_coordinate_probabilities(
            distribution, position, kept, (beta,)
        )[0]
    return result


def sample_trial_instance(
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    pair_rounds: int,
    rng: random.Random,
) -> TrialModel:
    known = tuple(
        sample_poly(model.weight, model.block_length, model.modulus, rng)
        for _ in range(model.polys)
    )
    hidden = tuple(
        sample_poly(model.weight, model.block_length, model.modulus, rng)
        for _ in range(model.polys)
    )
    descriptors = tuple(
        channel.sample(active, 2 * model.block_length, rng)
        for active in active_product_lists(model, known, hidden)
    )
    known_positions = tuple(
        position for state in known for position in state.positions
    )
    kept = model.block_length // model.polys
    distributions = []
    factor_sets = []
    for hidden_index in range(model.polys):
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
        factor_sets.append(factors)
        distributions.append(
            pair_mean_field_distributions(
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
        )
    return TrialModel(
        true_positions=tuple(state.positions for state in hidden),
        distributions=tuple(distributions),
        factor_sets=tuple(factor_sets),
        known_positions=known_positions,
        descriptors=descriptors,
    )


def sample_trial_model(
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    pair_rounds: int,
    rng: random.Random,
) -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[tuple[float, ...], ...], ...],
    tuple[tuple[CollisionFactor, ...], ...],
]:
    """Backward-compatible view used by the pair-only experiment."""

    instance = sample_trial_instance(
        model, channel, pair_slope, pair_rounds, rng
    )
    return (
        instance.true_positions,
        instance.distributions,
        instance.factor_sets,
    )


def one_trial(
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    pair_rounds: int,
    beta: float,
    particles: int,
    particle_repeats: int,
    seed: int,
) -> tuple[float, float]:
    true_positions, distributions, factor_sets = sample_trial_model(
        model,
        channel,
        pair_slope,
        pair_rounds,
        random.Random(seed),
    )
    kept = model.block_length // model.polys
    independent = 1.0
    for positions, marginal in zip(true_positions, distributions):
        independent *= independent_polynomial_success(
            positions, marginal, kept, beta
        )
    sequential_repeats = []
    for repeat in range(particle_repeats):
        success = 1.0
        for hidden_index in range(model.polys):
            success *= sequential_polynomial_success(
                true_positions[hidden_index],
                distributions[hidden_index],
                factor_sets[hidden_index],
                kept,
                beta,
                pair_slope,
                particles,
                random.Random(seed ^ (repeat << 32) ^ hidden_index ^ 0x51D),
            )
        sequential_repeats.append(success)
    return independent, sum(sequential_repeats) / particle_repeats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=4096)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--modulus", type=int, default=65537)
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--particles", type=int, default=16)
    parser.add_argument("--particle-repeats", type=int, default=2)
    parser.add_argument("--trials", type=int, default=4)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys < 2 or args.block_length % args.polys:
        parser.error("P must be at least two and divide the block length")
    if min(args.particles, args.particle_repeats, args.trials) < 1:
        parser.error("particles, repeats, and trials must be positive")
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
        "occupancy", args.selector_k, args.weight, partition_size
    )
    pair_slope, pair_r_squared = fit_pair_slope(
        channel,
        args.weight,
        args.pair_fit_samples,
        random.Random(args.seed ^ 0xA5A5A5A5),
    )
    start = time.perf_counter()
    results = [
        one_trial(
            model,
            channel,
            pair_slope,
            args.pair_rounds,
            args.beta,
            args.particles,
            args.particle_repeats,
            args.seed + trial,
        )
        for trial in range(args.trials)
    ]
    seconds = time.perf_counter() - start
    baseline_bits = args.polys * args.weight * math.log2(args.polys)
    independent_rows = [(result[0],) for result in results]
    sequential_rows = [(result[1],) for result in results]
    independent_gain, _ = best_weighted_expected_work_gain(
        independent_rows, baseline_bits, (args.beta,)
    )
    sequential_gain, _ = best_weighted_expected_work_gain(
        sequential_rows, baseline_bits, (args.beta,)
    )
    independent_half, _ = best_weighted_half_success_gain(
        independent_rows, baseline_bits, (args.beta,)
    )
    sequential_half, _ = best_weighted_half_success_gain(
        sequential_rows, baseline_bits, (args.beta,)
    )
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} K={args.selector_k} "
        f"beta={args.beta} particles={args.particles} "
        f"particle_repeats={args.particle_repeats} trials={args.trials} "
        f"pair_slope={pair_slope} pair_r_squared={pair_r_squared} "
        f"independent_expected_work_gain_bits={independent_gain} "
        f"sequential_expected_work_gain_bits={sequential_gain} "
        f"independent_half_success_gain_bits={independent_half} "
        f"sequential_half_success_gain_bits={sequential_half} "
        f"seconds={seconds} seed={args.seed}"
    )
    for trial, (independent, sequential) in enumerate(results):
        print(
            f"trial={trial} "
            f"independent_log_gain={math.log2(independent) + baseline_bits} "
            f"sequential_log_gain={math.log2(sequential) + baseline_bits}"
        )
    if args.output:
        fields = (
            "trial",
            "independent_success",
            "sequential_success",
            "independent_log_gain",
            "sequential_log_gain",
        )
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for trial, (independent, sequential) in enumerate(results):
                writer.writerow(
                    {
                        "trial": trial,
                        "independent_success": independent,
                        "sequential_success": sequential,
                        "independent_log_gain": math.log2(independent)
                        + baseline_bits,
                        "sequential_log_gain": math.log2(sequential)
                        + baseline_bits,
                    }
                )
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
