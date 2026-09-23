"""Leakage-aware support selection after sparse-factor Ring-LPN reduction.

The Ring-LPN attack of Boyle et al. first reduces the regular error modulo a
sparse factor ``X^n + c`` of ``X^N + 1``.  With ``P`` error polynomials, the
reduced syndrome is a length-``P*n`` code with redundancy ``n``.  A minimum-
weight dual check may prescribe any ``n-1`` zero coordinates: the resulting
homogeneous system has ``n-1`` equations in the ``n``-dimensional dual space.
The main diagnostic therefore chooses a global ``n-1``-subset and asks for
the probability that it contains the complete reduced error support.  We also
retain the earlier, more restrictive balanced selector as a control.

Reverse-Cuckoo beliefs live on the original regular blocks.  We project each
block belief modulo ``n``, combine the blocks into reduced-coordinate
occupancy scores, and sample a fixed-size subset with probability proportional
to the product of its scores.  The probability of containing the planted
support is evaluated exactly with elementary symmetric polynomials; no rare
support-hit Monte Carlo is required.

The global selector is a concrete statistical-decoding strategy under the
paper's favorable assumption that arbitrary minimum-weight dual checks can be
constructed or precomputed.  Its exact hypergeometric baseline is stricter
than the independent-coordinate lower bound used in the paper's parameter
table; paired leakage gains remain directly comparable.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from rev_cuckoo_basic_ring_lpn_isd import BALANCED_BETAS
from rev_cuckoo_basic_ring_lpn_isd_production import (
    active_product_lists,
    factor_descriptors,
    sample_poly,
)
from rev_cuckoo_full_lpn_field import FieldModel, PolyState
from rev_cuckoo_regular_isd import (
    collision_factors,
    fit_pair_slope,
    pair_mean_field_distributions,
)
from rev_cuckoo_regular_noise import LeakageChannel


@dataclass(frozen=True)
class FactorTrial:
    trial: int
    reduced_weights: tuple[int, ...]
    balanced_success_probabilities: tuple[float, ...]
    dual_success_probabilities: tuple[float, ...]


def elementary_symmetric(weights: Sequence[float], degree: int) -> float:
    """Return the elementary symmetric polynomial ``e_degree(weights)``."""

    if not 0 <= degree <= len(weights):
        return 0.0
    coefficients = [0.0] * (degree + 1)
    coefficients[0] = 1.0
    for weight in weights:
        for index in range(degree, 0, -1):
            coefficients[index] += weight * coefficients[index - 1]
    return coefficients[degree]


def product_weighted_support_probability(
    scores: Sequence[float],
    planted_support: frozenset[int],
    kept: int,
    beta: float,
) -> float:
    """Probability a product-weighted ``kept``-subset contains the support.

    A subset ``C`` has probability proportional to
    ``prod(scores[index] ** beta for index in C)``.  At ``beta == 0`` this is
    the uniform fixed-cardinality distribution, including when a score is
    zero.
    """

    width = len(scores)
    if not 0 <= kept <= width:
        raise ValueError("bad subset size")
    if beta < 0.0:
        raise ValueError("beta must be nonnegative")
    if any(not 0 <= index < width for index in planted_support):
        raise ValueError("support coordinate outside the score domain")
    support_size = len(planted_support)
    if support_size > kept:
        return 0.0
    if beta == 0.0:
        return math.comb(width - support_size, kept - support_size) / math.comb(
            width, kept
        )

    maximum = max(scores, default=0.0)
    if maximum <= 0.0:
        raise ValueError("scores have zero total mass")
    weights = tuple(
        (score / maximum) ** beta if score > 0.0 else 0.0
        for score in scores
    )
    planted_weight = math.prod(weights[index] for index in planted_support)
    if planted_weight == 0.0:
        return 0.0
    outside = tuple(
        weight
        for index, weight in enumerate(weights)
        if index not in planted_support
    )
    numerator = planted_weight * elementary_symmetric(
        outside, kept - support_size
    )
    denominator = elementary_symmetric(weights, kept)
    if denominator <= 0.0:
        raise ValueError("fewer than kept coordinates have positive weight")
    return numerator / denominator


def project_block_distribution(
    probabilities: Sequence[float],
    block: int,
    block_length: int,
    factor_degree: int,
) -> tuple[float, ...]:
    """Project one original block distribution to exponents modulo ``n``."""

    if len(probabilities) != block_length:
        raise ValueError("block distribution has the wrong length")
    projected = [0.0] * factor_degree
    base = block * block_length
    for offset, probability in enumerate(probabilities):
        projected[(base + offset) % factor_degree] += probability
    return tuple(projected)


def projected_occupancy_scores(
    distributions: Sequence[Sequence[float]],
    block_length: int,
    factor_degree: int,
) -> tuple[float, ...]:
    """Return marginal probabilities that reduced coordinates are occupied."""

    empty = [1.0] * factor_degree
    for block, probabilities in enumerate(distributions):
        projected = project_block_distribution(
            probabilities, block, block_length, factor_degree
        )
        for residue, probability in enumerate(projected):
            empty[residue] *= 1.0 - probability
    return tuple(1.0 - probability for probability in empty)


def reduced_support(
    state: PolyState,
    block_length: int,
    factor_degree: int,
    modulus: int,
    factor_constant: int = 1,
) -> frozenset[int]:
    """Reduce a sparse polynomial modulo ``X^n + factor_constant``.

    The constant affects only the nonzero scalar applied to each folded
    coefficient.  Choosing one is sufficient for support experiments; random
    nonzero coefficients still model the negligible cancellation event.
    """

    coefficients = [0] * factor_degree
    negative_constant = (-factor_constant) % modulus
    for block, (offset, coefficient) in enumerate(
        zip(state.positions, state.coefficients)
    ):
        quotient, residue = divmod(
            block * block_length + offset, factor_degree
        )
        multiplier = pow(negative_constant, quotient, modulus)
        coefficients[residue] = (
            coefficients[residue] + coefficient * multiplier
        ) % modulus
    return frozenset(
        residue for residue, coefficient in enumerate(coefficients) if coefficient
    )


def expected_reduced_weight(polys: int, weight: int, factor_degree: int) -> float:
    """Expected total weight after throwing each block point into ``n`` bins."""

    return polys * factor_degree * (
        1.0 - (1.0 - 1.0 / factor_degree) ** weight
    )


def statistical_decoding_bits(
    polys: int, reduced_weight: float, factor_degree: int
) -> float:
    """Sparse-factor statistical-decoding estimate from Boyle et al."""

    return (
        reduced_weight
        * math.log2(polys * factor_degree / (factor_degree - 1))
        + math.log2(factor_degree)
    )


def exact_dual_check_bits(
    polys: int, reduced_weight: int, factor_degree: int
) -> float:
    """Exact uniform zero-set containment cost, including ``n`` checks."""

    width = polys * factor_degree
    zeros = factor_degree - 1
    if not 0 <= reduced_weight <= zeros:
        return math.inf
    probability = math.comb(
        width - reduced_weight, zeros - reduced_weight
    ) / math.comb(width, zeros)
    return -math.log2(probability) + math.log2(factor_degree)


def reduced_weight_distribution(
    polys: int, weight: int, factor_degree: int
) -> tuple[float, ...]:
    """Exact occupancy distribution, ignoring negligible field cancellation."""

    one = [0.0] * (weight + 1)
    one[0] = 1.0
    for _ in range(weight):
        following = [0.0] * (weight + 1)
        for occupied, probability in enumerate(one):
            if probability == 0.0:
                continue
            following[occupied] += (
                probability * occupied / factor_degree
            )
            if occupied < weight:
                following[occupied + 1] += (
                    probability
                    * (factor_degree - occupied)
                    / factor_degree
                )
        one = following

    total = [1.0]
    for _ in range(polys):
        following = [0.0] * (len(total) + weight)
        for left_weight, left_probability in enumerate(total):
            for right_weight, right_probability in enumerate(one):
                following[left_weight + right_weight] += (
                    left_probability * right_probability
                )
        total = following
    return tuple(total)


def mixed_statistical_decoding_bits(
    polys: int,
    factor_degree: int,
    distribution: Sequence[float],
    minimum_weight: int = 0,
) -> tuple[float, float]:
    """Average the statistical test over a (possibly truncated) weight law."""

    accepted_mass = sum(distribution[minimum_weight:])
    if accepted_mass <= 0.0:
        raise ValueError("minimum weight rejects the complete distribution")
    zero_probability = (factor_degree - 1) / (polys * factor_degree)
    success = sum(
        probability * zero_probability**weight
        for weight, probability in enumerate(distribution)
        if weight >= minimum_weight
    ) / (accepted_mass * factor_degree)
    return -math.log2(success), accepted_mass


def mixed_exact_dual_bits(
    polys: int,
    factor_degree: int,
    distribution: Sequence[float],
    minimum_weight: int = 0,
) -> tuple[float, float]:
    """Average the exact uniform dual-check signal over the weight law."""

    accepted_mass = sum(distribution[minimum_weight:])
    if accepted_mass <= 0.0:
        raise ValueError("minimum weight rejects the complete distribution")
    width = polys * factor_degree
    zeros = factor_degree - 1
    denominator = math.comb(width, zeros)
    success = sum(
        probability
        * (
            math.comb(width - weight, zeros - weight) / denominator
            if weight <= zeros
            else 0.0
        )
        for weight, probability in enumerate(distribution)
        if weight >= minimum_weight
    ) / (accepted_mass * factor_degree)
    return -math.log2(success), accepted_mass


def one_trial(
    trial: int,
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    factor_degree: int,
    betas: Sequence[float],
    pair_rounds: int,
    rng: random.Random,
) -> FactorTrial:
    known = tuple(
        sample_poly(model.weight, model.block_length, model.modulus, rng)
        for _ in range(model.polys)
    )
    hidden = tuple(
        sample_poly(model.weight, model.block_length, model.modulus, rng)
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
    balanced_kept = factor_degree // model.polys
    balanced_successes = [1.0] * len(betas)
    global_scores = []
    global_support: set[int] = set()
    supports = []
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
            kept=model.block_length // model.polys,
            slope=pair_slope,
            restarts=0,
            rounds=pair_rounds,
            damping=0.5,
            rng=rng,
            factors=factors,
        )
        scores = projected_occupancy_scores(
            distributions, model.block_length, factor_degree
        )
        support = reduced_support(
            hidden_poly,
            model.block_length,
            factor_degree,
            model.modulus,
        )
        supports.append(len(support))
        base = hidden_index * factor_degree
        global_scores.extend(scores)
        global_support.update(base + residue for residue in support)
        for beta_index, beta in enumerate(betas):
            balanced_successes[beta_index] *= (
                product_weighted_support_probability(
                    scores, support, balanced_kept, beta
                )
            )
    dual_successes = tuple(
        product_weighted_support_probability(
            global_scores,
            frozenset(global_support),
            factor_degree - 1,
            beta,
        )
        for beta in betas
    )
    return FactorTrial(
        trial,
        tuple(supports),
        tuple(balanced_successes),
        dual_successes,
    )


def log2_mean_inverse(probabilities: Sequence[float]) -> float:
    """Return ``log2(mean(1 / p))`` without overflowing."""

    inverse_logs = [-math.log2(probability) for probability in probabilities]
    maximum = max(inverse_logs)
    return maximum + math.log2(
        sum(2.0 ** (value - maximum) for value in inverse_logs)
        / len(inverse_logs)
    )


def log2_inverse_mean(probabilities: Sequence[float]) -> float:
    """Return ``-log2(mean(p))`` for average distinguishing signal."""

    return -math.log2(sum(probabilities) / len(probabilities))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=4096)
    parser.add_argument("--factor-degree", type=int, default=128)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--modulus", type=int, default=65537)
    parser.add_argument("--selector-k", type=int, default=1)
    parser.add_argument("--trials", type=int, default=4)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trial-output", type=Path)
    args = parser.parse_args()
    if args.polys < 2:
        parser.error("the number of polynomials must be at least two")
    if args.factor_degree % args.polys:
        parser.error("factor degree must be divisible by the polynomial count")
    if args.block_length <= 0 or args.factor_degree <= 1:
        parser.error("block length and factor degree must be positive")
    if args.factor_degree - 1 < args.polys * args.weight:
        parser.error(
            "factor degree minus one must contain the unreduced support"
        )
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
    trials = tuple(
        one_trial(
            trial,
            model,
            channel,
            pair_slope,
            args.factor_degree,
            BALANCED_BETAS,
            args.pair_rounds,
            random.Random(args.seed + trial),
        )
        for trial in range(args.trials)
    )
    seconds = time.perf_counter() - start
    balanced_work_bits = tuple(
        log2_mean_inverse(
            [
                trial.balanced_success_probabilities[index]
                for trial in trials
            ]
        )
        for index in range(len(BALANCED_BETAS))
    )
    balanced_baseline_bits = balanced_work_bits[0]
    balanced_gains = tuple(
        balanced_baseline_bits - bits for bits in balanced_work_bits
    )
    balanced_best_index = max(
        range(len(balanced_gains)), key=balanced_gains.__getitem__
    )
    dual_work_bits = tuple(
        log2_mean_inverse(
            [trial.dual_success_probabilities[index] for trial in trials]
        )
        + math.log2(args.factor_degree)
        for index in range(len(BALANCED_BETAS))
    )
    dual_baseline_bits = dual_work_bits[0]
    dual_gains = tuple(
        dual_baseline_bits - bits for bits in dual_work_bits
    )
    dual_best_index = max(
        range(len(dual_gains)), key=dual_gains.__getitem__
    )
    dual_signal_bits = tuple(
        log2_inverse_mean(
            [trial.dual_success_probabilities[index] for trial in trials]
        )
        + math.log2(args.factor_degree)
        for index in range(len(BALANCED_BETAS))
    )
    dual_signal_baseline_bits = dual_signal_bits[0]
    dual_signal_gains = tuple(
        dual_signal_baseline_bits - bits for bits in dual_signal_bits
    )
    dual_signal_best_index = max(
        range(len(dual_signal_gains)), key=dual_signal_gains.__getitem__
    )
    expected_weight = expected_reduced_weight(
        args.polys, args.weight, args.factor_degree
    )
    weight_distribution = reduced_weight_distribution(
        args.polys, args.weight, args.factor_degree
    )
    weak_weight_threshold = math.ceil(expected_weight)
    mixed_stat_bits, _ = mixed_statistical_decoding_bits(
        args.polys, args.factor_degree, weight_distribution
    )
    conditioned_stat_bits, rejection_acceptance = (
        mixed_statistical_decoding_bits(
            args.polys,
            args.factor_degree,
            weight_distribution,
            weak_weight_threshold,
        )
    )
    exact_mixed_dual_bits, _ = mixed_exact_dual_bits(
        args.polys, args.factor_degree, weight_distribution
    )
    exact_conditioned_dual_bits, _ = mixed_exact_dual_bits(
        args.polys,
        args.factor_degree,
        weight_distribution,
        weak_weight_threshold,
    )
    row = {
        "polys": args.polys,
        "weight": args.weight,
        "block_length": args.block_length,
        "ring_dimension": args.weight * args.block_length,
        "factor_degree": args.factor_degree,
        "kept_per_polynomial": args.factor_degree // args.polys,
        "expected_reduced_weight": expected_weight,
        "mean_sampled_reduced_weight": statistics.mean(
            sum(trial.reduced_weights) for trial in trials
        ),
        "statistical_decoding_bits": statistical_decoding_bits(
            args.polys, expected_weight, args.factor_degree
        ),
        "unconditioned_statistical_decoding_bits": mixed_stat_bits,
        "weak_weight_threshold": weak_weight_threshold,
        "weak_weight_probability": 1.0 - rejection_acceptance,
        "conditioned_statistical_decoding_bits": conditioned_stat_bits,
        "unconditioned_exact_dual_bits": exact_mixed_dual_bits,
        "conditioned_exact_dual_bits": exact_conditioned_dual_bits,
        "rejection_expected_samples": 1.0 / rejection_acceptance,
        "partition_size": partition_size,
        "selector_k": args.selector_k,
        "trials": args.trials,
        "pair_rounds": args.pair_rounds,
        "pair_slope": pair_slope,
        "pair_r_squared": pair_r_squared,
        "balanced_uniform_selector_bits": balanced_baseline_bits,
        "balanced_best_beta": BALANCED_BETAS[balanced_best_index],
        "balanced_best_selector_bits": balanced_work_bits[
            balanced_best_index
        ],
        "balanced_best_gain_bits": balanced_gains[balanced_best_index],
        "dual_uniform_selector_bits": dual_baseline_bits,
        "dual_best_beta": BALANCED_BETAS[dual_best_index],
        "dual_best_selector_bits": dual_work_bits[dual_best_index],
        "dual_best_gain_bits": dual_gains[dual_best_index],
        "dual_uniform_signal_bits": dual_signal_baseline_bits,
        "dual_best_signal_beta": BALANCED_BETAS[
            dual_signal_best_index
        ],
        "dual_best_signal_bits": dual_signal_bits[
            dual_signal_best_index
        ],
        "dual_best_signal_gain_bits": dual_signal_gains[
            dual_signal_best_index
        ],
        "seconds": seconds,
        "seed": args.seed,
    }
    print(" ".join(f"{key}={value}" for key, value in row.items()))
    for (
        beta,
        balanced_bits,
        balanced_gain,
        dual_bits,
        dual_gain,
        signal_bits,
        signal_gain,
    ) in zip(
        BALANCED_BETAS,
        balanced_work_bits,
        balanced_gains,
        dual_work_bits,
        dual_gains,
        dual_signal_bits,
        dual_signal_gains,
    ):
        print(
            f"beta={beta} "
            f"balanced_bits={balanced_bits} "
            f"balanced_gain_bits={balanced_gain} "
            f"dual_bits={dual_bits} "
            f"dual_gain_bits={dual_gain} "
            f"dual_signal_bits={signal_bits} "
            f"dual_signal_gain_bits={signal_gain}"
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
            "reduced_weights",
            "total_reduced_weight",
            *(f"balanced_success_beta_{beta:g}" for beta in BALANCED_BETAS),
            *(f"dual_success_beta_{beta:g}" for beta in BALANCED_BETAS),
        )
        with args.trial_output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for trial in trials:
                writer.writerow(
                    {
                        "trial": trial.trial,
                        "reduced_weights": " ".join(
                            str(weight) for weight in trial.reduced_weights
                        ),
                        "total_reduced_weight": sum(trial.reduced_weights),
                        **{
                            f"balanced_success_beta_{beta:g}": (
                                trial.balanced_success_probabilities[index]
                            )
                            for index, beta in enumerate(BALANCED_BETAS)
                        },
                        **{
                            f"dual_success_beta_{beta:g}": (
                                trial.dual_success_probabilities[index]
                            )
                            for index, beta in enumerate(BALANCED_BETAS)
                        },
                    }
                )
        print(f"wrote {args.trial_output}")


if __name__ == "__main__":
    main()
