import itertools
import math
import random

from rev_cuckoo_regular_isd import (
    collision_factors,
    collision_pair_count,
    coordinate_ascent_selection,
    marginal_selection,
    marginal_tables,
    matrix_rank,
    negacyclic_column,
    pair_constraints,
    pair_mean_field_selection,
    pair_rectangle_coordinate_ascent,
    pair_rectangle_objective,
    pair_variational_objective,
    rectangle_mass,
    selected_matrix_full_rank,
)


def test_coordinate_ascent_does_not_decrease_marginal_rectangle() -> None:
    candidates = tuple(itertools.product(range(3), repeat=3))
    raw = [1 + (17 * index) % 23 for index in range(len(candidates))]
    total = sum(raw)
    probabilities = [value / total for value in raw]
    marginals = marginal_tables(candidates, probabilities, 3, 3)
    initial = marginal_selection(marginals, 1)
    following, mass = coordinate_ascent_selection(
        candidates, probabilities, initial, 1
    )
    assert mass >= rectangle_mass(candidates, probabilities, initial)
    assert mass == rectangle_mass(candidates, probabilities, following)


def test_negacyclic_column_wraps_with_a_sign() -> None:
    assert negacyclic_column((1, 2, 3, 4), 3, 5) == [3, 2, 1, 1]


def test_modular_rank() -> None:
    assert matrix_rank([[1, 2], [2, 1]], 5) == 2
    assert matrix_rank([[1, 2], [2, 4]], 5) == 1


def test_selected_identity_matrix_is_full_rank() -> None:
    selection = ((frozenset({0, 1}), frozenset({0, 1})),)
    assert selected_matrix_full_rank(selection, (), 2, 65537)


def test_zero_multiplier_forces_rank_failure() -> None:
    selection = (
        (frozenset({0}), frozenset({0})),
        (frozenset({0}), frozenset({0})),
    )
    assert not selected_matrix_full_rank(
        selection, ((0, 0, 0, 0),), 2, 65537
    )


def test_pair_collision_count_sums_both_partitions() -> None:
    descriptor = ((0, 0, 1), (2, 2, 2))
    assert collision_pair_count((0, 1, 2), descriptor) == 4


def test_pair_constraints_have_the_expected_shape() -> None:
    descriptor = ((0, 1, 0, 1), (1, 0, 1, 0))
    constraints = pair_constraints((0, 1), (descriptor, descriptor), 1, 2, 2)
    assert len(constraints) == 4
    assert all(len(constraint.left_labels) == 2 for constraint in constraints)
    assert all(len(constraint.right_labels) == 2 for constraint in constraints)


def test_collision_factors_compress_complete_pair_graphs() -> None:
    descriptor = ((0, 1, 0, 1), (1, 0, 1, 0))
    factors = collision_factors((0, 1), (descriptor, descriptor), 1, 2, 2)
    assert len(factors) == 4
    assert all(len(factor.labels) == 2 for factor in factors)


def test_collision_factor_objective_matches_explicit_pairs() -> None:
    descriptor = ((0, 1, 0, 1), (1, 0, 1, 0))
    distributions = ((0.25, 0.75), (0.6, 0.4))
    constraints = pair_constraints((0, 1), (descriptor, descriptor), 1, 2, 2)
    expected_collisions = 0.0
    for constraint in constraints:
        for left, left_probability in enumerate(
            distributions[constraint.left_variable]
        ):
            for right, right_probability in enumerate(
                distributions[constraint.right_variable]
            ):
                if constraint.left_labels[left] == constraint.right_labels[right]:
                    expected_collisions += left_probability * right_probability
    entropy = -sum(
        probability * math.log2(probability)
        for distribution in distributions
        for probability in distribution
    )
    factors = collision_factors((0, 1), (descriptor, descriptor), 1, 2, 2)
    assert math.isclose(
        pair_variational_objective(distributions, factors, -0.25),
        entropy - 0.25 * expected_collisions,
    )


def test_pair_rectangle_coordinate_ascent_is_monotone() -> None:
    descriptor = (
        (0, 1, 0, 1, 2, 2),
        (0, 0, 1, 1, 2, 2),
    )
    factors = collision_factors((0, 3), (descriptor, descriptor), 1, 2, 3)
    initial = (frozenset({0}), frozenset({0}))
    following = pair_rectangle_coordinate_ascent(
        initial, factors, kept=1, slope=-0.5, rounds=8
    )
    assert pair_rectangle_objective(
        following, factors, -0.5
    ) >= pair_rectangle_objective(initial, factors, -0.5)


def test_pair_selector_returns_fixed_cardinality_rectangle() -> None:
    descriptor = (
        (0, 1, 0, 1, 2, 2),
        (0, 0, 1, 1, 2, 2),
    )
    selection = pair_mean_field_selection(
        (0, 3),
        (descriptor, descriptor),
        polys=1,
        weight=2,
        block_length=3,
        kept=1,
        slope=-0.5,
        restarts=1,
        rounds=8,
        damping=0.5,
        rng=random.Random(7),
    )
    assert len(selection) == 2
    assert all(len(block) == 1 for block in selection)
