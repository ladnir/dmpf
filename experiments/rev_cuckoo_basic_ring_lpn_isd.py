"""Coefficient-aware leakage-assisted regular Prange for basic Ring-LPN.

The hidden error is ``(f,e)`` in the regular syndrome equation

    [I | M_r] (f,e)^T = y.

Each polynomial has one nonzero in every one of ``t`` blocks of length ``L``.
A regular-Prange iteration retains ``L/2`` coordinates per block, producing
exactly ``N=tL`` columns.  Three leakage-aware selectors are compared:

* marginal: retain the largest exact one-coordinate posterior marginals;
* pair: a polynomial collision-feature mean-field heuristic;
* joint: coordinate-ascent on the exact posterior rectangle mass.

The exact posterior marginalizes the hidden nonzero coefficients and models
all four coefficient-deduplicated product descriptor families.  Marginal and
joint are reduced-size oracles; pair is the scalable attack heuristic.  Every
selector is charged the actual full-rank event of its selected ``[I | M_r]``
matrix.  Hash descriptors use only the ideal independent channel.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
from dataclasses import dataclass
from typing import Sequence

from rev_cuckoo_full_lpn import logsumexp2
from rev_cuckoo_full_lpn_field import (
    FieldModel,
    PolyState,
    make_model,
    product_active_lists,
    product_active_lists_one_hidden,
)
from rev_cuckoo_regular_isd import (
    Selection,
    fit_pair_slope,
    joint_selection,
    marginal_selection,
    mean,
    pair_mean_field_distributions,
    random_selection,
    rectangle_mass,
    selected_matrix_full_rank,
    standard_error,
)
from rev_cuckoo_regular_noise import LeakageChannel, marginal_tables, normalize


Candidate = tuple[int, ...]


@dataclass(frozen=True)
class BasicTrialResult:
    base_rank: bool
    marginal_support: float
    marginal_rank: bool
    pair_support: float
    pair_rank: bool
    joint_support: float
    joint_rank: bool
    marginal_hits: int
    pair_hits: int
    joint_hits: int
    marginal_weighted_success: tuple[float, ...]
    pair_weighted_success: tuple[float, ...]
    marginal_balanced_success: tuple[float, ...]
    pair_balanced_success: tuple[float, ...]
    posterior_center_success: float


WEIGHT_BETAS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)
BALANCED_BETAS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)


def elementary_symmetric(weights: Sequence[float], degree: int) -> float:
    coefficients = [0.0] * (degree + 1)
    coefficients[0] = 1.0
    for weight in weights:
        for index in range(degree, 0, -1):
            coefficients[index] += weight * coefficients[index - 1]
    return coefficients[degree]


def weighted_subset_inclusion_probabilities(
    probabilities: Sequence[float], kept: int, beta: float
) -> tuple[float, ...]:
    """Inclusion probabilities for a product-weighted fixed-size subset."""

    if not 0 <= kept <= len(probabilities):
        raise ValueError("bad subset size")
    if beta < 0.0:
        raise ValueError("beta must be nonnegative")
    if kept == 0:
        return (0.0,) * len(probabilities)
    maximum = max(probabilities)
    if maximum <= 0.0:
        raise ValueError("probabilities have zero mass")
    weights = tuple(
        (probability / maximum) ** beta
        if probability > 0.0
        else 0.0
        for probability in probabilities
    )
    normalizer = elementary_symmetric(weights, kept)
    if normalizer <= 0.0:
        raise ValueError("fewer than kept coordinates have positive weight")
    return tuple(
        weight
        * elementary_symmetric(
            weights[:index] + weights[index + 1 :], kept - 1
        )
        / normalizer
        for index, weight in enumerate(weights)
    )


def balanced_half_inclusion_probabilities(
    probabilities: Sequence[float], beta: float
) -> tuple[float, ...]:
    """Select one coordinate from each upper/lower belief pair.

    The coordinates are sorted by belief, rank ``i`` is paired with rank
    ``i+L/2``, and exactly one member of every pair is sampled.  At beta zero
    every coordinate has inclusion probability one half.  Positive beta
    smoothly favors the higher-belief member without excluding either side.
    """

    if len(probabilities) % 2:
        raise ValueError("balanced-half selection requires an even domain")
    return balanced_stratified_inclusion_probabilities(
        probabilities, len(probabilities) // 2, beta
    )


def balanced_stratified_inclusion_probabilities(
    probabilities: Sequence[float], kept: int, beta: float
) -> tuple[float, ...]:
    """Select one coordinate from each cross-rank stratum.

    For ``P = L / kept``, rank ``i`` is grouped with ranks
    ``i + kept, ..., i + (P-1) kept``.  Sampling one member of every group
    gives exactly ``kept`` coordinates.  At beta zero every coordinate has
    inclusion probability ``1/P``.
    """

    if kept <= 0 or len(probabilities) % kept:
        raise ValueError("kept must be a positive divisor of the domain")
    if beta < 0.0:
        raise ValueError("beta must be nonnegative")
    order = sorted(
        range(len(probabilities)),
        key=lambda index: probabilities[index],
        reverse=True,
    )
    strata = len(order) // kept
    inclusion = [0.0] * len(order)
    for group_index in range(kept):
        group = tuple(
            order[group_index + stratum * kept]
            for stratum in range(strata)
        )
        if beta == 0.0:
            weights = [1.0] * strata
        else:
            positive = [probabilities[index] for index in group]
            maximum = max(positive)
            if maximum <= 0.0:
                weights = [1.0] * strata
            else:
                weights = [
                    (probability / maximum) ** beta
                    for probability in positive
                ]
        normalizer = sum(weights)
        for index, weight in zip(group, weights):
            inclusion[index] = weight / normalizer
    return tuple(inclusion)


def balanced_half_coordinate_probabilities(
    probabilities: Sequence[float],
    coordinate: int,
    betas: Sequence[float],
) -> tuple[float, ...]:
    """Return one coordinate's inclusion probability for several betas.

    This production evaluator sorts the belief vector once, identifies the
    coordinate's paired alternative, and evaluates every requested
    temperature.
    """

    if len(probabilities) % 2:
        raise ValueError("balanced-half selection requires an even domain")
    return balanced_stratified_coordinate_probabilities(
        probabilities,
        coordinate,
        len(probabilities) // 2,
        betas,
    )


def balanced_stratified_coordinate_probabilities(
    probabilities: Sequence[float],
    coordinate: int,
    kept: int,
    betas: Sequence[float],
) -> tuple[float, ...]:
    """Evaluate one coordinate in the stratified sampler for many betas."""

    if kept <= 0 or len(probabilities) % kept:
        raise ValueError("kept must be a positive divisor of the domain")
    if not 0 <= coordinate < len(probabilities):
        raise ValueError("coordinate is outside the domain")
    if any(beta < 0.0 for beta in betas):
        raise ValueError("beta must be nonnegative")
    order = sorted(
        range(len(probabilities)),
        key=lambda index: probabilities[index],
        reverse=True,
    )
    rank = order.index(coordinate)
    group_index = rank % kept
    strata = len(order) // kept
    group = tuple(
        order[group_index + stratum * kept]
        for stratum in range(strata)
    )
    positive = [probabilities[index] for index in group]
    coordinate_index = group.index(coordinate)
    maximum = max(positive)
    result = []
    for beta in betas:
        if beta == 0.0 or maximum <= 0.0:
            result.append(1.0 / strata)
            continue
        weights = [
            (probability / maximum) ** beta
            for probability in positive
        ]
        result.append(weights[coordinate_index] / sum(weights))
    return tuple(result)


def sample_balanced_stratified_subset(
    probabilities: Sequence[float],
    kept: int,
    beta: float,
    rng: random.Random,
) -> frozenset[int]:
    """Sample one fixed-cardinality subset from the stratified distribution."""

    if kept <= 0 or len(probabilities) % kept:
        raise ValueError("kept must be a positive divisor of the domain")
    if beta < 0.0:
        raise ValueError("beta must be nonnegative")
    order = sorted(
        range(len(probabilities)),
        key=lambda index: probabilities[index],
        reverse=True,
    )
    strata = len(order) // kept
    selected = []
    for group_index in range(kept):
        group = tuple(
            order[group_index + stratum * kept]
            for stratum in range(strata)
        )
        if beta == 0.0:
            selected.append(group[rng.randrange(strata)])
            continue
        positive = [probabilities[index] for index in group]
        maximum = max(positive)
        if maximum <= 0.0:
            selected.append(group[rng.randrange(strata)])
            continue
        weights = [
            (probability / maximum) ** beta
            for probability in positive
        ]
        threshold = rng.random() * sum(weights)
        prefix = 0.0
        for index, weight in zip(group, weights):
            prefix += weight
            if threshold <= prefix:
                selected.append(index)
                break
        else:
            selected.append(group[-1])
    return frozenset(selected)


def sample_balanced_stratified_subset_conditioned(
    probabilities: Sequence[float],
    kept: int,
    beta: float,
    required: int,
    rng: random.Random,
) -> tuple[frozenset[int], float]:
    """Sample a stratified subset conditioned on containing ``required``.

    The returned probability is the unconditional inclusion probability of
    ``required``.  For all other rank groups the conditional distribution is
    identical to the original sampler.
    """

    if kept <= 0 or len(probabilities) % kept:
        raise ValueError("kept must be a positive divisor of the domain")
    if not 0 <= required < len(probabilities):
        raise ValueError("required coordinate is outside the domain")
    if beta < 0.0:
        raise ValueError("beta must be nonnegative")
    order = sorted(
        range(len(probabilities)),
        key=lambda index: probabilities[index],
        reverse=True,
    )
    rank = order.index(required)
    required_group = rank % kept
    strata = len(order) // kept
    selected = []
    required_probability = None
    for group_index in range(kept):
        group = tuple(
            order[group_index + stratum * kept]
            for stratum in range(strata)
        )
        if beta == 0.0:
            weights = [1.0] * strata
        else:
            positive = [probabilities[index] for index in group]
            maximum = max(positive)
            weights = (
                [1.0] * strata
                if maximum <= 0.0
                else [
                    (probability / maximum) ** beta
                    for probability in positive
                ]
            )
        normalizer = sum(weights)
        if group_index == required_group:
            required_index = group.index(required)
            required_probability = weights[required_index] / normalizer
            selected.append(required)
            continue
        threshold = rng.random() * normalizer
        prefix = 0.0
        for index, weight in zip(group, weights):
            prefix += weight
            if threshold <= prefix:
                selected.append(index)
                break
        else:
            selected.append(group[-1])
    if required_probability is None:
        raise AssertionError("required rank group was not sampled")
    return frozenset(selected), required_probability


def states_by_support(
    model: FieldModel,
) -> dict[Candidate, tuple[PolyState, ...]]:
    grouped: dict[Candidate, list[PolyState]] = {}
    for state in model.poly_states:
        grouped.setdefault(state.positions, []).append(state)
    return {
        support: tuple(states)
        for support, states in grouped.items()
    }


def coefficient_marginalized_posterior(
    model: FieldModel,
    candidates: Sequence[Candidate],
    grouped_states: dict[Candidate, tuple[PolyState, ...]],
    known_states: Sequence[PolyState],
    descriptors,
    channel: LeakageChannel,
) -> tuple[list[float], list[list[float]]]:
    log_weights = []
    for support in candidates:
        states = grouped_states[support]
        state_weights = [
            channel.candidate_log_weight(
                product_active_lists_one_hidden(
                    model, known_states, hidden_state
                ),
                descriptors,
            )
            for hidden_state in states
        ]
        log_weights.append(
            logsumexp2(state_weights) - math.log2(len(state_weights))
        )
    probabilities = normalize(log_weights)
    marginals = marginal_tables(
        candidates,
        probabilities,
        model.weight,
        model.block_length,
    )
    return probabilities, marginals


def posterior_center_success_probability(
    candidates: Sequence[Candidate],
    probabilities: Sequence[float],
    kept: int,
    block_length: int,
) -> float:
    """Success when a posterior-drawn center is forced into every block.

    The remaining ``kept-1`` coordinates are uniform.  This exact quadratic
    diagnostic measures whether a correlated center sampler could convert the
    joint posterior into a repeatable Prange iteration; it is not itself a
    scalable sampler.
    """

    filler = (kept - 1) / (block_length - 1)
    return sum(
        left_probability
        * right_probability
        * filler
        ** sum(left != right for left, right in zip(left_candidate, right_candidate))
        for left_candidate, left_probability in zip(candidates, probabilities)
        for right_candidate, right_probability in zip(candidates, probabilities)
    )


def factor_descriptors(
    descriptors,
    hidden_poly_index: int,
    polys: int,
    weight: int,
):
    return tuple(
        descriptors[
            (known_poly_index * polys + hidden_poly_index) * weight
            + output_block
        ]
        for known_poly_index in range(polys)
        for output_block in range(weight)
    )


def one_trial(
    model: FieldModel,
    candidates: Sequence[Candidate],
    grouped_states: dict[Candidate, tuple[PolyState, ...]],
    channel: LeakageChannel,
    joint_restarts: int,
    pair_slope: float,
    pair_restarts: int,
    pair_rounds: int,
    pair_damping: float,
    rng: random.Random,
    rank_modulus: int | None = None,
) -> BasicTrialResult:
    support_rng = random.Random(rng.getrandbits(128))
    descriptor_rng = random.Random(rng.getrandbits(128))
    selector_rng = random.Random(rng.getrandbits(128))
    rank_rng = random.Random(rng.getrandbits(128))
    known_states = tuple(
        model.poly_states[support_rng.randrange(len(model.poly_states))]
        for _ in range(model.polys)
    )
    hidden_state = tuple(
        support_rng.randrange(len(model.poly_states))
        for _ in range(model.polys)
    )
    active_lists = product_active_lists(model, known_states, hidden_state)
    descriptors = tuple(
        channel.sample(active, 2 * model.block_length, descriptor_rng)
        for active in active_lists
    )
    known_positions = tuple(
        position
        for state in known_states
        for position in state.positions
    )

    kept = model.block_length // model.polys
    matrix_modulus = rank_modulus or model.modulus
    marginal_selections: list[Selection] = []
    pair_selections: list[Selection] = []
    joint_selections: list[Selection] = []
    marginal_support = 1.0
    pair_support = 1.0
    joint_support = 1.0
    marginal_hits = 0
    pair_hits = 0
    joint_hits = 0
    marginal_weighted_success = [1.0] * len(WEIGHT_BETAS)
    pair_weighted_success = [1.0] * len(WEIGHT_BETAS)
    marginal_balanced_success = [1.0] * len(BALANCED_BETAS)
    pair_balanced_success = [1.0] * len(BALANCED_BETAS)
    posterior_center_success = 1.0
    for hidden_poly_index in range(model.polys):
        descriptors_for_poly = factor_descriptors(
            descriptors,
            hidden_poly_index,
            model.polys,
            model.weight,
        )
        probabilities, marginals = coefficient_marginalized_posterior(
            model,
            candidates,
            grouped_states,
            known_states,
            descriptors_for_poly,
            channel,
        )
        posterior_center_success *= posterior_center_success_probability(
            candidates,
            probabilities,
            kept,
            model.block_length,
        )
        marginal = marginal_selection(marginals, kept)
        pair_distributions = pair_mean_field_distributions(
            known_positions,
            descriptors_for_poly,
            model.polys,
            model.weight,
            model.block_length,
            kept,
            pair_slope,
            pair_restarts,
            pair_rounds,
            pair_damping,
            selector_rng,
        )
        pair = marginal_selection(pair_distributions, kept)
        joint, joint_mass = joint_selection(
            candidates,
            probabilities,
            marginals,
            kept,
            joint_restarts,
            selector_rng,
        )
        marginal_selections.append(marginal)
        pair_selections.append(pair)
        joint_selections.append(joint)
        true_positions = model.poly_states[
            hidden_state[hidden_poly_index]
        ].positions
        marginal_hits += sum(
            position in marginal[block]
            for block, position in enumerate(true_positions)
        )
        pair_hits += sum(
            position in pair[block]
            for block, position in enumerate(true_positions)
        )
        joint_hits += sum(
            position in joint[block]
            for block, position in enumerate(true_positions)
        )
        for beta_index, beta in enumerate(WEIGHT_BETAS):
            for block, position in enumerate(true_positions):
                marginal_inclusion = weighted_subset_inclusion_probabilities(
                    marginals[block], kept, beta
                )
                pair_inclusion = weighted_subset_inclusion_probabilities(
                    pair_distributions[block], kept, beta
                )
                marginal_weighted_success[beta_index] *= (
                    marginal_inclusion[position]
                )
                pair_weighted_success[beta_index] *= pair_inclusion[position]
        for beta_index, beta in enumerate(BALANCED_BETAS):
            for block, position in enumerate(true_positions):
                marginal_inclusion = balanced_half_inclusion_probabilities(
                    marginals[block], beta
                )
                pair_inclusion = balanced_half_inclusion_probabilities(
                    pair_distributions[block], beta
                )
                marginal_balanced_success[beta_index] *= (
                    marginal_inclusion[position]
                )
                pair_balanced_success[beta_index] *= pair_inclusion[position]
        marginal_support *= rectangle_mass(
            candidates, probabilities, marginal
        )
        pair_support *= rectangle_mass(candidates, probabilities, pair)
        joint_support *= joint_mass

    uniform_selections = tuple(
        random_selection(
            model.weight,
            model.block_length,
            kept,
            rank_rng,
        )
        for _ in range(model.polys)
    )
    width = model.weight * model.block_length
    multipliers = tuple(
        tuple(rank_rng.randrange(matrix_modulus) for _ in range(width))
        for _ in range(model.polys - 1)
    )
    return BasicTrialResult(
        base_rank=selected_matrix_full_rank(
            uniform_selections,
            multipliers,
            model.block_length,
            matrix_modulus,
        ),
        marginal_support=marginal_support,
        marginal_rank=selected_matrix_full_rank(
            marginal_selections,
            multipliers,
            model.block_length,
            matrix_modulus,
        ),
        pair_support=pair_support,
        pair_rank=selected_matrix_full_rank(
            pair_selections,
            multipliers,
            model.block_length,
            matrix_modulus,
        ),
        joint_support=joint_support,
        joint_rank=selected_matrix_full_rank(
            joint_selections,
            multipliers,
            model.block_length,
            matrix_modulus,
        ),
        marginal_hits=marginal_hits,
        pair_hits=pair_hits,
        joint_hits=joint_hits,
        marginal_weighted_success=tuple(marginal_weighted_success),
        pair_weighted_success=tuple(pair_weighted_success),
        marginal_balanced_success=tuple(marginal_balanced_success),
        pair_balanced_success=tuple(pair_balanced_success),
        posterior_center_success=posterior_center_success,
    )


def mixed_rectangle_success(hits: int, variables: int, alpha: float) -> float:
    """Per-iteration support success of exploit-or-uniform block sampling."""

    if not 0.0 <= alpha < 1.0:
        raise ValueError("alpha must be in [0,1)")
    return (
        ((1.0 + alpha) / 2.0) ** hits
        * ((1.0 - alpha) / 2.0) ** (variables - hits)
    )


def best_mixed_expected_work_gain(
    hit_counts: Sequence[int],
    variables: int,
    alphas: Sequence[float],
) -> tuple[float, float]:
    """Return ``(best gain, alpha)`` under average per-instance work."""

    baseline_work = 2.0 ** variables
    best_gain = -math.inf
    best_alpha = 0.0
    for alpha in alphas:
        expected_work = mean(
            [
                1.0 / mixed_rectangle_success(hits, variables, alpha)
                for hits in hit_counts
            ]
        )
        gain = math.log2(baseline_work / expected_work)
        if gain > best_gain:
            best_gain = gain
            best_alpha = alpha
    return best_gain, best_alpha


def best_weighted_expected_work_gain(
    successes: Sequence[Sequence[float]],
    variables: int,
    betas: Sequence[float] = WEIGHT_BETAS,
) -> tuple[float, float]:
    """Return expected-work gain of the best product-weight temperature."""

    baseline_work = 2.0 ** variables
    best_gain = -math.inf
    best_beta = 0.0
    for beta_index, beta in enumerate(betas):
        expected_work = mean(
            [1.0 / trial[beta_index] for trial in successes]
        )
        gain = math.log2(baseline_work / expected_work)
        if gain > best_gain:
            best_gain = gain
            best_beta = beta
    return best_gain, best_beta


def iterations_for_average_success(
    per_iteration_success: Sequence[float], target: float = 0.5
) -> float:
    """Iterations reaching ``target`` success over the instance distribution."""

    if not per_iteration_success:
        raise ValueError("success list must be nonempty")
    if not 0.0 < target < 1.0:
        raise ValueError("target must be in (0,1)")

    def success(iterations: float) -> float:
        return mean(
            [
                1.0 - math.exp(iterations * math.log1p(-probability))
                if probability < 1.0
                else 1.0
                for probability in per_iteration_success
            ]
        )

    lower = 0.0
    upper = 1.0
    while success(upper) < target:
        upper *= 2.0
        if not math.isfinite(upper):
            return math.inf
    for _ in range(80):
        midpoint = (lower + upper) / 2.0
        if success(midpoint) < target:
            lower = midpoint
        else:
            upper = midpoint
    return upper


def best_weighted_half_success_gain(
    successes: Sequence[Sequence[float]],
    variables: int,
    betas: Sequence[float] = WEIGHT_BETAS,
) -> tuple[float, float]:
    """Return the best 50%-success work gain and its temperature."""

    baseline_probability = 2.0 ** (-variables)
    baseline_iterations = math.log(0.5) / math.log1p(-baseline_probability)
    best_gain = -math.inf
    best_beta = 0.0
    for beta_index, beta in enumerate(betas):
        iterations = iterations_for_average_success(
            [trial[beta_index] for trial in successes]
        )
        gain = math.log2(baseline_iterations / iterations)
        if gain > best_gain:
            best_gain = gain
            best_beta = beta
    return best_gain, best_beta


def report_selector(
    name: str,
    support: Sequence[float],
    rank: Sequence[bool],
    base_support: float,
    base_rank: float,
    base_success: float,
) -> None:
    successes = [
        probability * rank_ok
        for probability, rank_ok in zip(support, rank)
    ]
    success = mean(successes)
    print(
        f"{name}: support={mean(support):.9g} "
        f"support_se={standard_error(support):.3g} "
        f"support_gain_bits={math.log2(mean(support) / base_support):.6f} "
        f"rank={mean(rank):.6f} "
        f"total={success:.9g} total_se={standard_error(successes):.3g} "
        f"work_bits={-math.log2(success):.6f} "
        f"gain_bits={math.log2(success / base_success):.6f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", type=int, default=3)
    parser.add_argument("--block-length", type=int, default=8)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--modulus", type=int, default=3)
    parser.add_argument(
        "--rank-modulus",
        type=int,
        help="prime field for the Prange matrix rank audit (defaults to modulus)",
    )
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--joint-restarts", type=int, default=4)
    parser.add_argument("--pair-restarts", type=int, default=4)
    parser.add_argument("--pair-rounds", type=int, default=24)
    parser.add_argument("--pair-damping", type=float, default=0.5)
    parser.add_argument("--pair-fit-samples", type=int, default=20000)
    parser.add_argument(
        "--channel",
        choices=("raw", "occupancy", "feasible"),
        default="raw",
    )
    parser.add_argument("--selector-k", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if args.block_length % 2:
        parser.error("basic P=2 Ring-LPN requires an even block length")
    if min(
        args.weight,
        args.block_length,
        args.modulus - 1,
        (args.rank_modulus or args.modulus) - 1,
        args.trials,
        args.selector_k,
        args.pair_rounds,
        args.pair_fit_samples,
    ) < 1:
        parser.error("all dimensions, trials, and selector K must be positive")

    partition_size = args.partition_size or args.weight
    model = make_model(
        2,
        args.weight,
        args.block_length,
        args.modulus,
        enumerate_full_states=False,
    )
    candidates = tuple(
        itertools.product(range(args.block_length), repeat=args.weight)
    )
    grouped_states = states_by_support(model)
    channel = LeakageChannel.build(
        args.channel,
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
    rng = random.Random(args.seed)
    results = [
        one_trial(
            model,
            candidates,
            grouped_states,
            channel,
            args.joint_restarts,
            pair_slope,
            args.pair_restarts,
            args.pair_rounds,
            args.pair_damping,
            rng,
            args.rank_modulus,
        )
        for _ in range(args.trials)
    ]

    kept = args.block_length // 2
    base_support = (kept / args.block_length) ** (2 * args.weight)
    base_rank = mean([result.base_rank for result in results])
    base_success = base_support * base_rank
    print(
        f"basic-ring-lpn field={args.modulus} weight={args.weight} "
        f"rank_field={args.rank_modulus or args.modulus} "
        f"block_length={args.block_length} d={partition_size} kept={kept} "
        f"channel={args.channel} K={args.selector_k} "
        f"pair_slope={pair_slope:.6f} pair_r2={pair_r_squared:.6f} "
        f"candidates_per_poly={len(candidates)} trials={args.trials} "
        f"seed={args.seed}"
    )
    print(
        f"uniform: support={base_support:.9g} rank={base_rank:.6f} "
        f"total={base_success:.9g} work_bits={-math.log2(base_success):.6f} "
        "gain_bits=0.000000"
    )
    report_selector(
        "marginal",
        [result.marginal_support for result in results],
        [result.marginal_rank for result in results],
        base_support,
        base_rank,
        base_success,
    )
    report_selector(
        "pair",
        [result.pair_support for result in results],
        [result.pair_rank for result in results],
        base_support,
        base_rank,
        base_success,
    )
    report_selector(
        "joint",
        [result.joint_support for result in results],
        [result.joint_rank for result in results],
        base_support,
        base_rank,
        base_success,
    )
    alpha_grid = tuple(index / 100.0 for index in range(0, 100, 5))
    for selector in ("marginal", "pair", "joint"):
        gain, alpha = best_mixed_expected_work_gain(
            [getattr(result, f"{selector}_hits") for result in results],
            2 * args.weight,
            alpha_grid,
        )
        print(
            f"{selector}-mixed: expected_work_gain_bits={gain:.6f} "
            f"alpha={alpha:.2f}"
        )
    for selector in ("marginal", "pair"):
        weighted_successes = [
            getattr(result, f"{selector}_weighted_success")
            for result in results
        ]
        gain, beta = best_weighted_expected_work_gain(
            weighted_successes,
            2 * args.weight,
        )
        half_gain, half_beta = best_weighted_half_success_gain(
            weighted_successes,
            2 * args.weight,
        )
        print(
            f"{selector}-weighted: expected_work_gain_bits={gain:.6f} "
            f"beta={beta:.2f} half_success_gain_bits={half_gain:.6f} "
            f"half_beta={half_beta:.2f}"
        )
        balanced_successes = [
            getattr(result, f"{selector}_balanced_success")
            for result in results
        ]
        balanced_expected_gain, balanced_expected_beta = (
            best_weighted_expected_work_gain(
                balanced_successes,
                2 * args.weight,
                BALANCED_BETAS,
            )
        )
        balanced_half_gain, balanced_half_beta = (
            best_weighted_half_success_gain(
                balanced_successes,
                2 * args.weight,
                BALANCED_BETAS,
            )
        )
        print(
            f"{selector}-balanced: "
            f"expected_work_gain_bits={balanced_expected_gain:.6f} "
            f"beta={balanced_expected_beta:.2f} "
            f"half_success_gain_bits={balanced_half_gain:.6f} "
            f"half_beta={balanced_half_beta:.2f}"
        )
    posterior_center_successes = [
        (result.posterior_center_success,) for result in results
    ]
    posterior_center_expected, _ = best_weighted_expected_work_gain(
        posterior_center_successes,
        2 * args.weight,
        (1.0,),
    )
    posterior_center_half, _ = best_weighted_half_success_gain(
        posterior_center_successes,
        2 * args.weight,
        (1.0,),
    )
    print(
        "posterior-center: "
        f"expected_work_gain_bits={posterior_center_expected:.6f} "
        f"half_success_gain_bits={posterior_center_half:.6f}"
    )


if __name__ == "__main__":
    main()
