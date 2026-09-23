"""Leakage-weighted regular-Stern estimator for one Ring-LPN instance.

The model is the q-ary analogue of enumeration-based regular ISD.  A regular
permutation places ``v`` coordinates from each error block in the information
set J.  The attack allows ``p/2`` error blocks in each half of J, enumerates
their positions and nonzero field coefficients, and collides the two lists on
``ell`` field equations.  Reverse-Cuckoo leakage changes only the probability
that the planted support has this block distribution; all list, coefficient,
and collision costs are charged independently of the leakage.

This is a single-syndrome attack.  It receives no benefit from stationary
support reuse or from additional Ring-LPN samples.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
import time
from dataclasses import dataclass
from pathlib import Path

from rev_cuckoo_basic_ring_lpn_isd import (
    BALANCED_BETAS,
    balanced_stratified_coordinate_probabilities,
)
from rev_cuckoo_basic_ring_lpn_isd_production import (
    active_product_lists,
    factor_descriptors,
    sample_poly,
)
from rev_cuckoo_full_lpn_field import FieldModel
from rev_cuckoo_regular_isd import (
    collision_factors,
    fit_pair_slope,
    pair_mean_field_distributions,
)
from rev_cuckoo_regular_noise import LeakageChannel


@dataclass(frozen=True)
class SternTrial:
    # beta-major, followed by the P*t regular blocks in public fixed order.
    complement_inclusions: tuple[tuple[float, ...], ...]


@dataclass(frozen=True)
class SternEstimate:
    selector: str
    beta: float
    p: int
    ell: int
    list_log2: float
    collision_list_log2: float
    iteration_log2: float
    mean_inverse_success_log2: float
    expected_work_log2: float
    gain_vs_uniform_prange_bits: float


def log2_binomial(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
    ) / math.log(2.0)


def poisson_binomial_probability(
    probabilities: tuple[float, ...], target: int
) -> float:
    """Probability of exactly target successes for independent Bernoullis."""

    if not 0 <= target <= len(probabilities):
        return 0.0
    coefficients = [0.0] * (target + 1)
    coefficients[0] = 1.0
    for probability in probabilities:
        for degree in range(target, 0, -1):
            coefficients[degree] = (
                coefficients[degree] * (1.0 - probability)
                + coefficients[degree - 1] * probability
            )
        coefficients[0] *= 1.0 - probability
    return coefficients[target]


def regular_stern_success(
    information_inclusions: tuple[float, ...], p: int
) -> float:
    if p % 2:
        raise ValueError("regular Stern requires even p")
    if len(information_inclusions) % 2:
        raise ValueError("the regular blocks must split into equal halves")
    midpoint = len(information_inclusions) // 2
    target = p // 2
    return poisson_binomial_probability(
        information_inclusions[:midpoint], target
    ) * poisson_binomial_probability(
        information_inclusions[midpoint:], target
    )


def uniform_regular_stern_success(
    blocks: int, information_fraction: float, p: int
) -> float:
    return 2.0 ** (
        2 * log2_binomial(blocks // 2, p // 2)
        + p * math.log2(information_fraction)
        + (blocks - p) * math.log2(1.0 - information_fraction)
    )


def log2_add(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    maximum = max(left, right)
    return maximum + math.log2(
        2.0 ** (left - maximum) + 2.0 ** (right - maximum)
    )


def log2_elementary_symmetric(
    log2_weights: tuple[float, ...], degree: int
) -> float:
    """Log2 of the requested elementary symmetric polynomial."""

    if not 0 <= degree <= len(log2_weights):
        return -math.inf
    coefficients = [-math.inf] * (degree + 1)
    coefficients[0] = 0.0
    for weight in log2_weights:
        for index in range(degree, 0, -1):
            coefficients[index] = log2_add(
                coefficients[index], coefficients[index - 1] + weight
            )
    return coefficients[degree]


def regular_information_widths(
    blocks: int, block_length: int, polys: int, ell: int
) -> tuple[int, ...]:
    """Balanced integral block widths whose sum is k + ell."""

    if blocks % 2:
        raise ValueError("the regular blocks must split into equal halves")
    base = (polys - 1) * block_length // polys
    quotient, remainder = divmod(ell, blocks)
    if base + quotient + (remainder != 0) > block_length:
        raise ValueError("ell makes the information set exceed the code")
    widths = [base + quotient] * blocks
    half = blocks // 2
    left_extra = (remainder + 1) // 2
    right_extra = remainder // 2
    for index in range(left_extra):
        widths[index] += 1
    for index in range(right_extra):
        widths[half + index] += 1
    assert sum(widths) == blocks * base + ell
    return tuple(widths)


def stern_iteration_cost(
    information_widths: tuple[int, ...],
    modulus: int,
    p: int,
    ell: int,
) -> tuple[float, float, float]:
    """Return largest leaf, collision-list, and iteration exponents."""

    if p % 2:
        raise ValueError("regular Stern requires even p")
    if len(information_widths) % 2:
        raise ValueError("the regular blocks must split into equal halves")
    half_weight = p // 2
    midpoint = len(information_widths) // 2
    coefficient_log2 = math.log2(modulus - 1)
    left_log2 = log2_elementary_symmetric(
        tuple(
            math.log2(width) + coefficient_log2
            for width in information_widths[:midpoint]
        ),
        half_weight,
    )
    right_log2 = log2_elementary_symmetric(
        tuple(
            math.log2(width) + coefficient_log2
            for width in information_widths[midpoint:]
        ),
        half_weight,
    )
    list_log2 = max(left_log2, right_log2)
    collision_log2 = left_log2 + right_log2 - ell * math.log2(modulus)
    return list_log2, collision_log2, max(0.0, list_log2, collision_log2)


def log2_mean_inverse(probabilities: list[float]) -> float:
    if any(probability <= 0.0 for probability in probabilities):
        return math.inf
    exponents = [-math.log2(probability) for probability in probabilities]
    maximum = max(exponents)
    return maximum + math.log2(
        statistics.fmean(2.0 ** (value - maximum) for value in exponents)
    )


def one_trial(
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    pair_rounds: int,
    rng: random.Random,
) -> SternTrial:
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
        position for polynomial in known for position in polynomial.positions
    )
    complement_size = model.block_length // model.polys
    inclusions = [[] for _ in BALANCED_BETAS]
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
            kept=complement_size,
            slope=pair_slope,
            restarts=0,
            rounds=pair_rounds,
            damping=0.5,
            rng=rng,
            factors=factors,
        )
        for block, true_position in enumerate(hidden_poly.positions):
            complement_probabilities = (
                balanced_stratified_coordinate_probabilities(
                    distributions[block],
                    true_position,
                    complement_size,
                    BALANCED_BETAS,
                )
            )
            for beta_index, probability in enumerate(
                complement_probabilities
            ):
                inclusions[beta_index].append(probability)
    return SternTrial(tuple(tuple(values) for values in inclusions))


def estimate_attack(
    trials: list[SternTrial],
    polys: int,
    weight: int,
    block_length: int,
    modulus: int,
    p: int,
    ell: int,
) -> tuple[SternEstimate, ...]:
    blocks = polys * weight
    information_widths = regular_information_widths(
        blocks, block_length, polys, ell
    )
    list_log2, collision_log2, iteration_log2 = stern_iteration_cost(
        information_widths,
        modulus,
        p,
        ell,
    )
    estimates = []
    for beta_index, beta in enumerate(BALANCED_BETAS):
        successes = []
        for trial in trials:
            # The production selector is evaluated at complement size L/P.
            # Scale each block by its exact integral complement width.  For
            # ell << N this correction is tiny, but it avoids pretending that
            # ell/(P*t) is a fractional coordinate in every regular block.
            complement = tuple(
                min(
                    1.0,
                    max(
                        0.0,
                        probability
                        * (block_length - width)
                        / (block_length / polys),
                    ),
                )
                for probability, width in zip(
                    trial.complement_inclusions[beta_index],
                    information_widths,
                )
            )
            information = tuple(1.0 - value for value in complement)
            successes.append(regular_stern_success(information, p))
        inverse_success = log2_mean_inverse(successes)
        expected_work = iteration_log2 + inverse_success
        baseline = blocks * math.log2(polys)
        estimates.append(
            SternEstimate(
                selector="leakage",
                beta=beta,
                p=p,
                ell=ell,
                list_log2=list_log2,
                collision_list_log2=collision_log2,
                iteration_log2=iteration_log2,
                mean_inverse_success_log2=inverse_success,
                expected_work_log2=expected_work,
                gain_vs_uniform_prange_bits=baseline - expected_work,
            )
        )

    uniform_success = regular_stern_success(
        tuple(width / block_length for width in information_widths), p
    )
    uniform_inverse = -math.log2(uniform_success)
    estimates.append(
        SternEstimate(
            selector="uniform",
            beta=0.0,
            p=p,
            ell=ell,
            list_log2=list_log2,
            collision_list_log2=collision_log2,
            iteration_log2=iteration_log2,
            mean_inverse_success_log2=uniform_inverse,
            expected_work_log2=iteration_log2 + uniform_inverse,
            gain_vs_uniform_prange_bits=(
                blocks * math.log2(polys)
                - iteration_log2
                - uniform_inverse
            ),
        )
    )
    return tuple(estimates)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=4096)
    parser.add_argument("--modulus", type=int, default=65537)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument("--trials", type=int, default=6)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--p-max", type=int, default=8)
    parser.add_argument("--ell-max", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys < 2 or args.block_length % args.polys:
        parser.error("P must divide L and be at least two")
    if args.p_max < 0 or args.ell_max < 0 or args.trials < 1:
        parser.error("search limits must be nonnegative and trials positive")
    blocks = args.polys * args.weight
    if blocks % 2 or args.p_max > blocks:
        parser.error("P*t must be even and p_max cannot exceed it")
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
    trials = [
        one_trial(
            model,
            channel,
            pair_slope,
            args.pair_rounds,
            random.Random(args.seed + trial),
        )
        for trial in range(args.trials)
    ]
    rows = [
        estimate
        for p in range(0, args.p_max + 1, 2)
        for ell in range(args.ell_max + 1)
        for estimate in estimate_attack(
            trials,
            args.polys,
            args.weight,
            args.block_length,
            args.modulus,
            p,
            ell,
        )
    ]
    seconds = time.perf_counter() - start
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} modulus={args.modulus} "
        f"K={args.selector_k} trials={args.trials} pair_slope={pair_slope} "
        f"pair_r_squared={pair_r_squared} seconds={seconds} seed={args.seed}"
    )
    selectors = ["uniform", *(f"beta={beta}" for beta in BALANCED_BETAS)]
    for selector_name in selectors:
        if selector_name == "uniform":
            selected = [row for row in rows if row.selector == "uniform"]
        else:
            beta = float(selector_name.split("=", 1)[1])
            selected = [
                row
                for row in rows
                if row.selector == "leakage" and row.beta == beta
            ]
        best = min(selected, key=lambda row: row.expected_work_log2)
        print(
            f"selector={selector_name} best_p={best.p} best_ell={best.ell} "
            f"list_log2={best.list_log2} "
            f"collision_list_log2={best.collision_list_log2} "
            f"iteration_log2={best.iteration_log2} "
            f"inverse_success_log2={best.mean_inverse_success_log2} "
            f"expected_work_log2={best.expected_work_log2} "
            f"gain_vs_uniform_prange_bits={best.gain_vs_uniform_prange_bits}"
        )
    if args.output:
        fields = tuple(SternEstimate.__dataclass_fields__)
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {field: getattr(row, field) for field in fields}
                )
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
