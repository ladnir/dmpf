"""Toy leakage-aware regular-Prange estimator for Reverse Cuckoo.

The experiment models the honest party's ``P`` regular sparse polynomials as
the error vector in

    syndrome = E[0] + R[1] E[1] + ... + R[P-1] E[P-1]

over F_q[X]/(X^N+1).  Every polynomial has one nonzero in each of ``t``
blocks of length ``L``.  A regular Prange iteration chooses ``L/P`` columns
from every block of every polynomial, for exactly ``N`` selected columns.

Four selectors are compared on the same leakage instance:

* uniform: the ordinary regular-Prange baseline;
* marginal: choose the most likely coordinates independently in every block;
* pair: fit a pair-collision likelihood and optimize its mean-field marginals;
* joint: coordinate-ascent rectangles scored with the complete correlated
  Reverse-Cuckoo posterior.

For every selector the experiment checks the rank of the actual selected
negacyclic multiplication matrix.  Posterior rectangle mass is used as a
Rao-Blackwellized support-success probability, so rare support hits need not
be observed directly.  The pair selector is polynomial in ``t`` and ``L``;
the exact ``L^t`` posterior is used only to score it in this reduced diagnostic.
The marginal and joint selectors themselves require that exact posterior.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
import statistics
from collections.abc import Sequence
from dataclasses import dataclass

from rev_cuckoo_goldreich import SharedGoldreichChannel
from rev_cuckoo_regular_noise import (
    Descriptor,
    LeakageChannel,
    marginal_tables,
    normalize,
    support_lists_one_unknown,
)


Candidate = tuple[int, ...]
Selection = tuple[frozenset[int], ...]


@dataclass(frozen=True)
class PairConstraint:
    left_variable: int
    right_variable: int
    left_labels: tuple[int, ...]
    right_labels: tuple[int, ...]


@dataclass(frozen=True)
class CollisionFactor:
    labels: tuple[tuple[int, ...], ...]
    label_count: int


def rectangle_mass(
    candidates: Sequence[Candidate],
    probabilities: Sequence[float],
    selection: Selection,
) -> float:
    return sum(
        probability
        for candidate, probability in zip(candidates, probabilities)
        if all(value in selection[block] for block, value in enumerate(candidate))
    )


def marginal_selection(
    marginals: Sequence[Sequence[float]], kept: int
) -> Selection:
    return tuple(
        frozenset(
            sorted(
                range(len(marginal)),
                key=lambda value: marginal[value],
                reverse=True,
            )[:kept]
        )
        for marginal in marginals
    )


def random_selection(
    variable_count: int,
    block_length: int,
    kept: int,
    rng: random.Random,
) -> Selection:
    return tuple(
        frozenset(rng.sample(range(block_length), kept))
        for _ in range(variable_count)
    )


def coordinate_ascent_selection(
    candidates: Sequence[Candidate],
    probabilities: Sequence[float],
    initial: Selection,
    kept: int,
    max_rounds: int = 32,
) -> tuple[Selection, float]:
    """Optimize a fixed-cardinality posterior rectangle one block at a time."""

    selection = [set(values) for values in initial]
    block_length = max(max(candidate) for candidate in candidates) + 1
    for _ in range(max_rounds):
        changed = False
        for block in range(len(selection)):
            scores = [0.0] * block_length
            for candidate, probability in zip(candidates, probabilities):
                if all(
                    other == block or candidate[other] in selection[other]
                    for other in range(len(selection))
                ):
                    scores[candidate[block]] += probability
            following = set(
                sorted(
                    range(block_length),
                    key=lambda value: scores[value],
                    reverse=True,
                )[:kept]
            )
            if following != selection[block]:
                selection[block] = following
                changed = True
        if not changed:
            break
    result = tuple(frozenset(values) for values in selection)
    return result, rectangle_mass(candidates, probabilities, result)


def joint_selection(
    candidates: Sequence[Candidate],
    probabilities: Sequence[float],
    marginals: Sequence[Sequence[float]],
    kept: int,
    restarts: int,
    rng: random.Random,
) -> tuple[Selection, float]:
    starts = [marginal_selection(marginals, kept)]
    starts.extend(
        random_selection(
            len(marginals), len(marginals[0]), kept, rng
        )
        for _ in range(restarts)
    )
    best_selection = starts[0]
    best_mass = -1.0
    for initial in starts:
        selection, mass = coordinate_ascent_selection(
            candidates, probabilities, initial, kept
        )
        if mass > best_mass:
            best_selection = selection
            best_mass = mass
    return best_selection, best_mass


def collision_pair_count(
    active: Sequence[int], descriptor: Descriptor
) -> int:
    result = 0
    for partition in descriptor:
        loads: dict[int, int] = {}
        for item in active:
            label = partition[item]
            result += loads.get(label, 0)
            loads[label] = loads.get(label, 0) + 1
    return result


def fit_pair_slope(
    channel: LeakageChannel,
    items: int,
    samples: int,
    rng: random.Random,
) -> tuple[float, float]:
    """Fit log2 channel density by the public pair-collision count."""

    active = tuple(range(items))
    features = []
    targets = []
    for _ in range(samples):
        descriptor = (
            tuple(
                rng.randrange(channel.partition_size) for _ in active
            ),
            tuple(
                rng.randrange(channel.partition_size) for _ in active
            ),
        )
        log_density = channel.log_likelihood(active, descriptor)
        if math.isfinite(log_density):
            features.append(collision_pair_count(active, descriptor))
            targets.append(log_density)
    feature_mean = mean(features)
    target_mean = mean(targets)
    variance = sum((value - feature_mean) ** 2 for value in features)
    if variance == 0.0:
        return 0.0, 0.0
    covariance = sum(
        (feature - feature_mean) * (target - target_mean)
        for feature, target in zip(features, targets)
    )
    slope = covariance / variance
    total = sum((target - target_mean) ** 2 for target in targets)
    residual = sum(
        (
            target
            - target_mean
            - slope * (feature - feature_mean)
        )
        ** 2
        for feature, target in zip(features, targets)
    )
    r_squared = 0.0 if total == 0.0 else 1.0 - residual / total
    return slope, r_squared


def pair_constraints(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    polys: int,
    weight: int,
    block_length: int,
) -> tuple[PairConstraint, ...]:
    constraints = []
    for known_poly in range(polys):
        for output_block in range(weight):
            descriptor = descriptors[known_poly * weight + output_block]
            for left_item in range(weight):
                left_variable = (output_block - left_item) % weight
                left_base = known[known_poly * weight + left_item]
                for right_item in range(left_item + 1, weight):
                    right_variable = (output_block - right_item) % weight
                    right_base = known[known_poly * weight + right_item]
                    for partition in descriptor:
                        constraints.append(
                            PairConstraint(
                                left_variable,
                                right_variable,
                                tuple(
                                    partition[left_base + value]
                                    for value in range(block_length)
                                ),
                                tuple(
                                    partition[right_base + value]
                                    for value in range(block_length)
                                ),
                            )
                        )
    return tuple(constraints)


def collision_factors(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    polys: int,
    weight: int,
    block_length: int,
) -> tuple[CollisionFactor, ...]:
    """Build complete collision factors without duplicating pair label vectors."""

    factors = []
    for known_poly in range(polys):
        for output_block in range(weight):
            descriptor = descriptors[known_poly * weight + output_block]
            for partition in descriptor:
                labels_by_variable: list[tuple[int, ...] | None] = [
                    None
                ] * weight
                label_count = 0
                for item in range(weight):
                    variable = (output_block - item) % weight
                    base = known[known_poly * weight + item]
                    labels = tuple(
                        partition[base + value]
                        for value in range(block_length)
                    )
                    labels_by_variable[variable] = labels
                    label_count = max(label_count, max(labels) + 1)
                if any(labels is None for labels in labels_by_variable):
                    raise AssertionError("collision factor omitted a variable")
                factors.append(
                    CollisionFactor(
                        tuple(labels for labels in labels_by_variable if labels is not None),
                        label_count,
                    )
                )
    return tuple(factors)


def normalized_exponents(log_weights: Sequence[float]) -> list[float]:
    maximum = max(log_weights)
    weights = [2.0 ** (value - maximum) for value in log_weights]
    total = sum(weights)
    return [value / total for value in weights]


def pair_variational_objective(
    distributions: Sequence[Sequence[float]],
    factors: Sequence[CollisionFactor],
    slope: float,
) -> float:
    expected_collisions = 0.0
    for factor in factors:
        total_mass = [0.0] * factor.label_count
        square_mass = [0.0] * factor.label_count
        for variable, labels in enumerate(factor.labels):
            variable_mass = [0.0] * factor.label_count
            for value, probability in enumerate(distributions[variable]):
                variable_mass[labels[value]] += probability
            for label, mass in enumerate(variable_mass):
                total_mass[label] += mass
                square_mass[label] += mass * mass
        expected_collisions += 0.5 * sum(
            total * total - squares
            for total, squares in zip(total_mass, square_mass)
        )
    entropy = -sum(
        probability * math.log2(probability)
        for distribution in distributions
        for probability in distribution
        if probability
    )
    return slope * expected_collisions + entropy


def pair_mean_field_distributions(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    polys: int,
    weight: int,
    block_length: int,
    kept: int,
    slope: float,
    restarts: int,
    rounds: int,
    damping: float,
    rng: random.Random,
    factors: Sequence[CollisionFactor] | None = None,
) -> tuple[tuple[float, ...], ...]:
    """Infer product beliefs using only a pair-collision likelihood surrogate."""

    if factors is None:
        factors = collision_factors(
            known, descriptors, polys, weight, block_length
        )
    best_objective = -math.inf
    best_distributions = None
    for restart in range(restarts + 1):
        if restart == 0:
            distributions = [
                [1.0 / block_length] * block_length
                for _ in range(weight)
            ]
        else:
            distributions = []
            for _ in range(weight):
                values = [rng.expovariate(1.0) for _ in range(block_length)]
                total = sum(values)
                distributions.append([value / total for value in values])

        for _ in range(rounds):
            log_updates = [[0.0] * block_length for _ in range(weight)]
            for factor in factors:
                total_mass = [0.0] * factor.label_count
                variable_masses = []
                for variable, labels in enumerate(factor.labels):
                    variable_mass = [0.0] * factor.label_count
                    for value in range(block_length):
                        variable_mass[labels[value]] += distributions[variable][value]
                    variable_masses.append(variable_mass)
                    for label, mass in enumerate(variable_mass):
                        total_mass[label] += mass
                for variable, labels in enumerate(factor.labels):
                    variable_mass = variable_masses[variable]
                    updates = log_updates[variable]
                    for value in range(block_length):
                        label = labels[value]
                        updates[value] += slope * (
                            total_mass[label] - variable_mass[label]
                        )

            maximum_change = 0.0
            for variable in range(weight):
                updated = normalized_exponents(log_updates[variable])
                for value in range(block_length):
                    following = (
                        (1.0 - damping) * distributions[variable][value]
                        + damping * updated[value]
                    )
                    maximum_change = max(
                        maximum_change,
                        abs(following - distributions[variable][value]),
                    )
                    distributions[variable][value] = following
            if maximum_change < 1e-9:
                break

        objective = pair_variational_objective(
            distributions, factors, slope
        )
        if objective > best_objective:
            best_objective = objective
            best_distributions = distributions

    if best_distributions is None:
        raise AssertionError("pair selector produced no variational solution")
    return tuple(tuple(distribution) for distribution in best_distributions)


def pair_center_coordinate_ascent(
    factors: Sequence[CollisionFactor],
    initial: Sequence[int],
    slope: float,
    rounds: int,
) -> Candidate:
    """Find a collision-model center without enumerating the joint domain.

    Each coordinate update maximizes the complete pair-collision surrogate
    with all other coordinates fixed.  One round has the same asymptotic cost
    as one mean-field round, ``O(len(factors) * t * L)``.
    """

    if not factors:
        return tuple(initial)
    state = list(initial)
    variable_count = len(factors[0].labels)
    block_length = len(factors[0].labels[0])
    if len(state) != variable_count:
        raise ValueError("center dimension does not match collision factors")
    for _ in range(rounds):
        changed = False
        for variable in range(variable_count):
            scores = [0] * block_length
            for factor in factors:
                other_counts = [0] * factor.label_count
                for other, labels in enumerate(factor.labels):
                    if other != variable:
                        other_counts[labels[state[other]]] += 1
                labels = factor.labels[variable]
                for value, label in enumerate(labels):
                    scores[value] += other_counts[label]
            if slope < 0.0:
                following = min(range(block_length), key=scores.__getitem__)
            else:
                following = max(range(block_length), key=scores.__getitem__)
            if following != state[variable]:
                state[variable] = following
                changed = True
        if not changed:
            break
    return tuple(state)


def pair_mean_field_selection(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    polys: int,
    weight: int,
    block_length: int,
    kept: int,
    slope: float,
    restarts: int,
    rounds: int,
    damping: float,
    rng: random.Random,
) -> Selection:
    """Select a rectangle using only a pair-collision likelihood surrogate."""

    distributions = pair_mean_field_distributions(
        known,
        descriptors,
        polys,
        weight,
        block_length,
        kept,
        slope,
        restarts,
        rounds,
        damping,
        rng,
    )
    return marginal_selection(distributions, kept)


def pair_rectangle_objective(
    selection: Selection,
    factors: Sequence[CollisionFactor],
    slope: float,
) -> float:
    kept = len(selection[0])
    block_length = len(factors[0].labels[0])
    distributions = [
        [
            1.0 / kept if value in block else 0.0
            for value in range(block_length)
        ]
        for block in selection
    ]
    return pair_variational_objective(distributions, factors, slope)


def pair_rectangle_coordinate_ascent(
    initial: Selection,
    factors: Sequence[CollisionFactor],
    kept: int,
    slope: float,
    rounds: int,
) -> Selection:
    """Optimize a fixed-cardinality rectangle under the pair surrogate."""

    selection = [set(block) for block in initial]
    block_length = len(factors[0].labels[0])
    for _ in range(rounds):
        changed = False
        for variable in range(len(selection)):
            scores = [0.0] * block_length
            for factor in factors:
                other_mass = [0.0] * factor.label_count
                for other, labels in enumerate(factor.labels):
                    if other == variable:
                        continue
                    mass = 1.0 / kept
                    for value in selection[other]:
                        other_mass[labels[value]] += mass
                labels = factor.labels[variable]
                for value in range(block_length):
                    scores[value] += slope * other_mass[labels[value]]
            following = set(
                sorted(
                    range(block_length),
                    key=lambda value: scores[value],
                    reverse=True,
                )[:kept]
            )
            if following != selection[variable]:
                selection[variable] = following
                changed = True
        if not changed:
            break
    return tuple(frozenset(block) for block in selection)


def pair_refined_rectangle_selection(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    polys: int,
    weight: int,
    block_length: int,
    kept: int,
    slope: float,
    distributions: Sequence[Sequence[float]],
    restarts: int,
    rounds: int,
    rng: random.Random,
) -> Selection:
    factors = collision_factors(
        known, descriptors, polys, weight, block_length
    )
    starts = [marginal_selection(distributions, kept)]
    starts.extend(
        random_selection(weight, block_length, kept, rng)
        for _ in range(restarts)
    )
    best = starts[0]
    best_objective = -math.inf
    for initial in starts:
        following = pair_rectangle_coordinate_ascent(
            initial, factors, kept, slope, rounds
        )
        objective = pair_rectangle_objective(following, factors, slope)
        if objective > best_objective:
            best = following
            best_objective = objective
    return best


def negacyclic_column(
    multiplier: Sequence[int], column: int, modulus: int
) -> list[int]:
    width = len(multiplier)
    result = [0] * width
    for row, value in enumerate(multiplier):
        output = row + column
        if output >= width:
            output -= width
            value = -value
        result[output] = value % modulus
    return result


def matrix_rank(matrix: list[list[int]], modulus: int) -> int:
    if not matrix:
        return 0
    rows = len(matrix)
    columns = len(matrix[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(rank, rows)
                if matrix[row][column] % modulus
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column] % modulus, -1, modulus)
        matrix[rank] = [
            value * inverse % modulus for value in matrix[rank]
        ]
        for row in range(rows):
            if row == rank:
                continue
            coefficient = matrix[row][column] % modulus
            if coefficient:
                matrix[row] = [
                    (value - coefficient * pivot_value) % modulus
                    for value, pivot_value in zip(matrix[row], matrix[rank])
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def selected_matrix_full_rank(
    selections: Sequence[Selection],
    multipliers: Sequence[Sequence[int]],
    block_length: int,
    modulus: int,
) -> bool:
    polys = len(selections)
    weight = len(selections[0])
    width = weight * block_length
    columns: list[list[int]] = []
    for poly in range(polys):
        for block, offsets in enumerate(selections[poly]):
            for offset in sorted(offsets):
                column = block * block_length + offset
                if poly == 0:
                    vector = [0] * width
                    vector[column] = 1
                else:
                    vector = negacyclic_column(
                        multipliers[poly - 1], column, modulus
                    )
                columns.append(vector)
    if len(columns) != width:
        raise ValueError(
            f"regular solve set has {len(columns)} columns, expected {width}"
        )
    matrix = [
        [columns[column][row] for column in range(width)]
        for row in range(width)
    ]
    return matrix_rank(matrix, modulus) == width


def posterior_for_polynomial(
    candidates: Sequence[Candidate],
    candidate_lists: Sequence[Sequence[Sequence[int]]],
    descriptors,
    weight: int,
    channel: LeakageChannel,
) -> tuple[list[float], list[list[float]]]:
    log_weights = []
    for lists in candidate_lists:
        log_weights.append(channel.candidate_log_weight(lists, descriptors))
    probabilities = normalize(log_weights)
    block_length = max(max(candidate) for candidate in candidates) + 1
    marginals = marginal_tables(
        candidates, probabilities, weight, block_length
    )
    return probabilities, marginals


@dataclass(frozen=True)
class TrialResult:
    base_rank: bool
    marginal_support: float
    marginal_rank: bool
    joint_support: float
    joint_rank: bool
    pair_support: float
    pair_rank: bool
    joint_over_marginal_support_gain: float


def one_trial(
    candidates: Sequence[Candidate],
    polys: int,
    weight: int,
    block_length: int,
    kept: int,
    modulus: int,
    channel: LeakageChannel,
    joint_restarts: int,
    pair_slope: float,
    pair_restarts: int,
    pair_rounds: int,
    pair_damping: float,
    rng: random.Random,
) -> TrialResult:
    support_rng = random.Random(rng.getrandbits(128))
    descriptor_rng = random.Random(rng.getrandbits(128))
    selector_rng = random.Random(rng.getrandbits(128))
    rank_rng = random.Random(rng.getrandbits(128))
    known = tuple(
        support_rng.randrange(block_length) for _ in range(polys * weight)
    )
    candidate_lists = tuple(
        support_lists_one_unknown(known, candidate, polys, weight)
        for candidate in candidates
    )
    marginal_selections = []
    joint_selections = []
    pair_selections = []
    marginal_support = 1.0
    joint_support = 1.0
    pair_support = 1.0

    for _ in range(polys):
        true_index = support_rng.randrange(len(candidates))
        true_lists = candidate_lists[true_index]
        if isinstance(channel, SharedGoldreichChannel):
            descriptors = channel.sample_batch(
                true_lists, 2 * block_length, descriptor_rng
            )
        else:
            descriptors = tuple(
                channel.sample(active, 2 * block_length, descriptor_rng)
                for active in true_lists
            )
        probabilities, marginals = posterior_for_polynomial(
            candidates,
            candidate_lists,
            descriptors,
            weight,
            channel,
        )
        marginal = marginal_selection(marginals, kept)
        marginal_mass = rectangle_mass(
            candidates, probabilities, marginal
        )
        joint, joint_mass = joint_selection(
            candidates,
            probabilities,
            marginals,
            kept,
            joint_restarts,
            selector_rng,
        )
        pair = pair_mean_field_selection(
            known,
            descriptors,
            polys,
            weight,
            block_length,
            kept,
            pair_slope,
            pair_restarts,
            pair_rounds,
            pair_damping,
            selector_rng,
        )
        pair_mass = rectangle_mass(candidates, probabilities, pair)
        marginal_selections.append(marginal)
        joint_selections.append(joint)
        pair_selections.append(pair)
        marginal_support *= marginal_mass
        joint_support *= joint_mass
        pair_support *= pair_mass

    uniform_selections = [
        random_selection(weight, block_length, kept, rank_rng)
        for _ in range(polys)
    ]
    width = weight * block_length
    multipliers = [
        tuple(rank_rng.randrange(modulus) for _ in range(width))
        for _ in range(polys - 1)
    ]
    base_rank = selected_matrix_full_rank(
        uniform_selections, multipliers, block_length, modulus
    )
    marginal_rank = selected_matrix_full_rank(
        marginal_selections, multipliers, block_length, modulus
    )
    joint_rank = selected_matrix_full_rank(
        joint_selections, multipliers, block_length, modulus
    )
    pair_rank = selected_matrix_full_rank(
        pair_selections, multipliers, block_length, modulus
    )
    return TrialResult(
        base_rank=base_rank,
        marginal_support=marginal_support,
        marginal_rank=marginal_rank,
        joint_support=joint_support,
        joint_rank=joint_rank,
        pair_support=pair_support,
        pair_rank=pair_rank,
        joint_over_marginal_support_gain=math.log2(
            joint_support / marginal_support
        ),
    )


def mean(values: Sequence[float]) -> float:
    return statistics.mean(values)


def standard_error(values: Sequence[float]) -> float:
    return (
        statistics.stdev(values) / math.sqrt(len(values))
        if len(values) > 1
        else 0.0
    )


def report_selector(
    name: str,
    support: Sequence[float],
    rank: Sequence[bool],
    base_support: float,
    base_rank: float,
    base_success: float,
) -> None:
    total = [probability * rank_ok for probability, rank_ok in zip(support, rank)]
    success = mean(total)
    support_mean = mean(support)
    rank_mean = mean(rank)
    print(
        f"{name}: support={support_mean:.9g} "
        f"support_se={standard_error(support):.3g} "
        f"support_gain_bits={math.log2(support_mean / base_support):.6f} "
        f"rank={rank_mean:.6f} rank_se={standard_error(rank):.3g} "
        f"rank_gain_bits={math.log2(rank_mean / base_rank):.6f} "
        f"total={success:.9g} total_se={standard_error(total):.3g} "
        f"work_bits={-math.log2(success):.6f} "
        f"gain_bits={math.log2(success / base_success):.6f}"
    )


def self_test() -> None:
    candidates = tuple(itertools.product(range(2), repeat=2))
    probabilities = [0.7, 0.1, 0.1, 0.1]
    marginals = marginal_tables(candidates, probabilities, 2, 2)
    marginal = marginal_selection(marginals, 1)
    joint, joint_mass = joint_selection(
        candidates,
        probabilities,
        marginals,
        1,
        2,
        random.Random(0),
    )
    if joint_mass + 1e-12 < rectangle_mass(
        candidates, probabilities, marginal
    ):
        raise AssertionError("joint selection is worse than its marginal start")
    if rectangle_mass(candidates, probabilities, joint) != joint_mass:
        raise AssertionError("reported joint mass is inconsistent")

    if negacyclic_column((1, 2, 3, 4), 3, 5) != [3, 2, 1, 1]:
        raise AssertionError("negacyclic multiplication column is wrong")

    descriptor = ((0, 1, 0, 1), (1, 0, 1, 0))
    constraints = pair_constraints((0, 1), (descriptor, descriptor), 1, 2, 2)
    if len(constraints) != 4:
        raise AssertionError("wrong number of pair constraints")

    identity_selection = (
        (frozenset({0, 1}), frozenset({0, 1})),
    )
    if not selected_matrix_full_rank(identity_selection, (), 2, 3):
        raise AssertionError("identity solve matrix is singular")
    zero_multiplier_selection = (
        (frozenset({0}), frozenset({0})),
        (frozenset({0}), frozenset({0})),
    )
    if selected_matrix_full_rank(
        zero_multiplier_selection, ((0, 0, 0, 0),), 2, 3
    ):
        raise AssertionError("zero multiplication columns were accepted")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=2)
    parser.add_argument("--weight", type=int, default=3)
    parser.add_argument("--block-length", type=int, default=8)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--modulus", type=int, default=3)
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--joint-restarts", type=int, default=4)
    parser.add_argument("--pair-restarts", type=int, default=4)
    parser.add_argument("--pair-rounds", type=int, default=24)
    parser.add_argument("--pair-damping", type=float, default=0.5)
    parser.add_argument("--pair-slope", type=float)
    parser.add_argument("--pair-fit-samples", type=int, default=20000)
    parser.add_argument(
        "--channel",
        choices=("raw", "occupancy", "feasible"),
        default="raw",
    )
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument(
        "--hash-family",
        choices=("ideal", "goldreich"),
        default="ideal",
    )
    parser.add_argument("--linear-security", type=int, default=10)
    parser.add_argument("--exact-rank-limit", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if min(
        args.polys,
        args.weight,
        args.block_length,
        args.modulus - 1,
        args.trials,
        args.selector_k,
        args.pair_rounds,
        args.pair_fit_samples,
    ) < 1:
        parser.error("dimensions, trials, and selector K must be positive")
    if min(args.joint_restarts, args.pair_restarts) < 0:
        parser.error("selector restarts must be nonnegative")
    if not 0.0 < args.pair_damping <= 1.0:
        parser.error("pair damping must be in (0,1]")
    if args.block_length % args.polys:
        parser.error("block length must be divisible by the polynomial count")
    if args.linear_security < 0:
        parser.error("linear-security must be nonnegative")
    if args.exact_rank_limit < 0:
        parser.error("exact-rank-limit must be nonnegative")

    self_test()
    partition_size = args.partition_size or args.weight
    kept = args.block_length // args.polys
    candidates = tuple(
        itertools.product(range(args.block_length), repeat=args.weight)
    )
    base_channel = LeakageChannel.build(
        args.channel,
        args.selector_k,
        args.weight,
        partition_size,
    )
    channel = (
        SharedGoldreichChannel(
            base_channel,
            linear_security=args.linear_security,
            exact_rank_limit=args.exact_rank_limit,
        )
        if args.hash_family == "goldreich"
        else base_channel
    )
    if args.pair_slope is None:
        pair_slope, pair_r_squared = fit_pair_slope(
            channel,
            args.weight,
            args.pair_fit_samples,
            random.Random(args.seed ^ 0xA5A5A5A5),
        )
    else:
        pair_slope = args.pair_slope
        pair_r_squared = math.nan
    rng = random.Random(args.seed)
    results = [
        one_trial(
            candidates,
            args.polys,
            args.weight,
            args.block_length,
            kept,
            args.modulus,
            channel,
            args.joint_restarts,
            pair_slope,
            args.pair_restarts,
            args.pair_rounds,
            args.pair_damping,
            rng,
        )
        for _ in range(args.trials)
    ]

    base_support = (kept / args.block_length) ** (
        args.polys * args.weight
    )
    base_rank = mean([result.base_rank for result in results])
    base_success = base_support * base_rank
    print(
        f"field={args.modulus} polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} d={partition_size} "
        f"kept={kept} candidates_per_poly={len(candidates)} "
        f"channel={args.channel} "
        f"hash_family={args.hash_family} "
        f"linear_security={args.linear_security} "
        f"exact_rank_limit={args.exact_rank_limit} "
        f"selector_k={args.selector_k if args.channel == 'occupancy' else '-'} "
        f"joint_restarts={args.joint_restarts} "
        f"pair_slope={pair_slope:.6f} pair_r2={pair_r_squared:.6f} "
        f"pair_restarts={args.pair_restarts} pair_rounds={args.pair_rounds} "
        f"trials={args.trials} "
        f"seed={args.seed}"
    )
    if isinstance(channel, SharedGoldreichChannel):
        print(
            f"goldreich_feature_bits={channel.feature_bits} "
            f"sampled_systems={channel.sampled_systems} "
            f"rank_deficient_systems={channel.rank_deficient_systems} "
            f"rank_failures={channel.rank_failures} "
            f"approximate_rank_queries={channel.approximate_rank_queries}"
        )
    print(
        f"uniform: support={base_support:.9g} rank={base_rank:.6f} "
        f"total={base_success:.9g} "
        f"work_bits={-math.log2(base_success):.6f} gain_bits=0.000000"
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
    joint_gain = [
        result.joint_over_marginal_support_gain for result in results
    ]
    print(
        "joint_over_marginal_support_gain: "
        f"mean={mean(joint_gain):.6f} se={standard_error(joint_gain):.6f}"
    )


if __name__ == "__main__":
    main()
