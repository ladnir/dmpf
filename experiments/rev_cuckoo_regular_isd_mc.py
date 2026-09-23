"""Monte Carlo scorer for the scalable Reverse-Cuckoo regular-ISD selector.

Unlike ``rev_cuckoo_regular_isd.py``, this program never enumerates ``L^t``
support tuples.  For a rectangle R chosen by the polynomial pair selector it
uses

    Pr[Y in R | descriptor]
      = P^(-t) E[W(Y) | Y uniform in R] / E[W(Y) | Y uniform],

where W is the exact finite-K descriptor likelihood.  The ratio of Monte
Carlo means is therefore an estimate of the support-success gain over regular
Prange.  Selected-matrix rank is omitted; the large-field reduced experiments
found rank one in every recorded trial.
"""

from __future__ import annotations

import argparse
import math
import random
import statistics
from collections.abc import Sequence

from rev_cuckoo_regular_isd import (
    fit_pair_slope,
    marginal_selection,
    pair_mean_field_distributions,
)
from rev_cuckoo_regular_noise import (
    LeakageChannel,
    support_lists_one_unknown,
)


def logmeanexp2(values: Sequence[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return -math.inf
    maximum = max(finite)
    return maximum + math.log2(
        sum(2.0 ** (value - maximum) for value in finite) / len(values)
    )


def sample_uniform_candidate(
    weight: int, block_length: int, rng: random.Random
) -> tuple[int, ...]:
    return tuple(rng.randrange(block_length) for _ in range(weight))


def sample_rectangle_candidate(
    selection: Sequence[frozenset[int]], rng: random.Random
) -> tuple[int, ...]:
    choices = tuple(tuple(block) for block in selection)
    return tuple(rng.choice(block) for block in choices)


def sample_product_candidate(
    distributions: Sequence[Sequence[float]], rng: random.Random
) -> tuple[int, ...]:
    return tuple(
        rng.choices(range(len(distribution)), weights=distribution, k=1)[0]
        for distribution in distributions
    )


def log2_product_probability(
    candidate: Sequence[int], distributions: Sequence[Sequence[float]]
) -> float:
    return sum(
        math.log2(distribution[value])
        for value, distribution in zip(candidate, distributions)
    )


def logaddexp2(left: float, right: float) -> float:
    maximum = max(left, right)
    return maximum + math.log2(
        2.0 ** (left - maximum) + 2.0 ** (right - maximum)
    )


def conditional_distributions(
    distributions: Sequence[Sequence[float]],
    selection: Sequence[frozenset[int]],
) -> tuple[tuple[float, ...], ...]:
    result = []
    for distribution, block in zip(distributions, selection):
        total = sum(distribution[value] for value in block)
        result.append(
            tuple(
                distribution[value] / total if value in block else 0.0
                for value in range(len(distribution))
            )
        )
    return tuple(result)


def candidate_log_weight(
    known: Sequence[int],
    candidate: Sequence[int],
    descriptors,
    polys: int,
    weight: int,
    channel: LeakageChannel,
) -> float:
    lists = support_lists_one_unknown(known, candidate, polys, weight)
    return channel.candidate_log_weight(lists, descriptors)


def score_rectangle(
    known: Sequence[int],
    selection: Sequence[frozenset[int]],
    distributions: Sequence[Sequence[float]],
    descriptors,
    polys: int,
    weight: int,
    block_length: int,
    channel: LeakageChannel,
    samples: int,
    proposal_mix: float,
    rng: random.Random,
) -> tuple[float, float, float, float, float]:
    """Score with a defensive product-belief/uniform importance mixture."""

    full_weights = []
    rectangle_weights = []
    rectangle_distributions = conditional_distributions(
        distributions, selection
    )
    uniform_full_log_probability = -weight * math.log2(block_length)
    kept = len(selection[0])
    uniform_rectangle_log_probability = -weight * math.log2(kept)
    log_mix = math.log2(proposal_mix)
    log_uniform_mix = math.log2(1.0 - proposal_mix)
    for _ in range(samples):
        candidate = (
            sample_product_candidate(distributions, rng)
            if rng.random() < proposal_mix
            else sample_uniform_candidate(weight, block_length, rng)
        )
        candidate_probability = log2_product_probability(
            candidate, distributions
        )
        proposal_probability = logaddexp2(
            log_mix + candidate_probability,
            log_uniform_mix + uniform_full_log_probability,
        )
        likelihood = candidate_log_weight(
            known, candidate, descriptors, polys, weight, channel
        )
        full_weights.append(
            likelihood + uniform_full_log_probability - proposal_probability
        )

        candidate = (
            sample_product_candidate(rectangle_distributions, rng)
            if rng.random() < proposal_mix
            else sample_rectangle_candidate(selection, rng)
        )
        candidate_probability = log2_product_probability(
            candidate, rectangle_distributions
        )
        proposal_probability = logaddexp2(
            log_mix + candidate_probability,
            log_uniform_mix + uniform_rectangle_log_probability,
        )
        likelihood = candidate_log_weight(
            known, candidate, descriptors, polys, weight, channel
        )
        rectangle_weights.append(
            likelihood
            + uniform_rectangle_log_probability
            - proposal_probability
        )
    full_log_mean = logmeanexp2(full_weights)
    rectangle_log_mean = logmeanexp2(rectangle_weights)
    return (
        rectangle_log_mean - full_log_mean,
        full_log_mean,
        rectangle_log_mean,
        effective_sample_fraction(full_weights),
        effective_sample_fraction(rectangle_weights),
    )


def effective_sample_fraction(log_weights: Sequence[float]) -> float:
    finite = [value for value in log_weights if math.isfinite(value)]
    if not finite:
        return 0.0
    maximum = max(finite)
    weights = [2.0 ** (value - maximum) for value in finite]
    total = sum(weights)
    return total * total / (sum(value * value for value in weights) * len(log_weights))


def logmean_bits(values: Sequence[float]) -> float:
    return logmeanexp2(values)


def standard_error(values: Sequence[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values) / math.sqrt(len(values))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=64)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--trials", type=int, default=4)
    parser.add_argument("--samples", type=int, default=4096)
    parser.add_argument("--pair-restarts", type=int, default=0)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--pair-damping", type=float, default=0.5)
    parser.add_argument("--pair-slope", type=float)
    parser.add_argument("--pair-fit-samples", type=int, default=50000)
    parser.add_argument("--proposal-mix", type=float, default=0.9)
    parser.add_argument(
        "--channel", choices=("raw", "occupancy"), default="occupancy"
    )
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if min(
        args.polys,
        args.weight,
        args.block_length,
        args.trials,
        args.pair_rounds,
    ) < 1:
        parser.error("all size, trial, and round parameters must be positive")
    if args.samples < 0:
        parser.error("sample count must be nonnegative")
    if args.block_length % args.polys:
        parser.error("block length must be divisible by the polynomial count")
    if args.pair_restarts < 0:
        parser.error("pair restarts must be nonnegative")
    if not 0.0 < args.proposal_mix < 1.0:
        parser.error("proposal mix must lie strictly between zero and one")

    partition_size = args.partition_size or args.weight
    channel = LeakageChannel.build(
        args.channel,
        args.selector_k,
        args.weight,
        partition_size,
    )
    fit_rng = random.Random(args.seed ^ 0x4B1D)
    if args.pair_slope is None:
        pair_slope, pair_r2 = fit_pair_slope(
            channel, args.weight, args.pair_fit_samples, fit_rng
        )
    else:
        pair_slope = args.pair_slope
        pair_r2 = math.nan

    rng = random.Random(args.seed)
    kept = args.block_length // args.polys
    trial_gains = []
    polynomial_gains = []
    inclusion_rates = []
    full_log_means = []
    rectangle_log_means = []
    full_ess = []
    rectangle_ess = []
    for trial in range(args.trials):
        known = tuple(
            rng.randrange(args.block_length)
            for _ in range(args.polys * args.weight)
        )
        trial_gain = 0.0
        included = 0
        coordinates = 0
        for _ in range(args.polys):
            true_candidate = sample_uniform_candidate(
                args.weight, args.block_length, rng
            )
            true_lists = support_lists_one_unknown(
                known, true_candidate, args.polys, args.weight
            )
            descriptors = tuple(
                channel.sample(active, 2 * args.block_length, rng)
                for active in true_lists
            )
            distributions = pair_mean_field_distributions(
                known,
                descriptors,
                args.polys,
                args.weight,
                args.block_length,
                kept,
                pair_slope,
                args.pair_restarts,
                args.pair_rounds,
                args.pair_damping,
                rng,
            )
            selection = marginal_selection(distributions, kept)
            if args.samples:
                (
                    gain,
                    full_log_mean,
                    rectangle_log_mean,
                    full_ess_fraction,
                    rectangle_ess_fraction,
                ) = score_rectangle(
                    known,
                    selection,
                    distributions,
                    descriptors,
                    args.polys,
                    args.weight,
                    args.block_length,
                    channel,
                    args.samples,
                    args.proposal_mix,
                    rng,
                )
                polynomial_gains.append(gain)
                full_log_means.append(full_log_mean)
                rectangle_log_means.append(rectangle_log_mean)
                full_ess.append(full_ess_fraction)
                rectangle_ess.append(rectangle_ess_fraction)
                trial_gain += gain
            included += sum(
                value in block
                for value, block in zip(true_candidate, selection)
            )
            coordinates += args.weight
        if args.samples:
            trial_gains.append(trial_gain)
        inclusion_rates.append(included / coordinates)
        gain_text = f"gain_bits={trial_gain:.6f} " if args.samples else ""
        print(
            f"trial={trial} {gain_text}"
            f"coordinate_inclusion={inclusion_rates[-1]:.6f}"
        )

    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} d={partition_size} "
        f"channel={args.channel} selector_k={args.selector_k} "
        f"pair_slope={pair_slope:.9f} pair_r2={pair_r2:.9f} "
        f"trials={args.trials} samples={args.samples} seed={args.seed}"
    )
    if args.samples:
        print(
            f"gain_bits_logmean={logmean_bits(trial_gains):.6f} "
            f"gain_bits_meanlog={statistics.mean(trial_gains):.6f} "
            f"gain_bits_meanlog_se={standard_error(trial_gains):.6f} "
            f"per_poly_meanlog={statistics.mean(polynomial_gains):.6f} "
            f"coordinate_inclusion={statistics.mean(inclusion_rates):.6f} "
            f"full_log_mean={statistics.mean(full_log_means):.6f} "
            f"rectangle_log_mean={statistics.mean(rectangle_log_means):.6f} "
            f"full_ess_fraction={statistics.mean(full_ess):.6f} "
            f"rectangle_ess_fraction={statistics.mean(rectangle_ess):.6f}"
        )
    else:
        print(
            f"coordinate_inclusion={statistics.mean(inclusion_rates):.6f} "
            f"coordinate_inclusion_trial_se={standard_error(inclusion_rates):.6f}"
        )


if __name__ == "__main__":
    main()
