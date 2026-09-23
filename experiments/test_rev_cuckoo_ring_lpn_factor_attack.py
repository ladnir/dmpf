import math

from rev_cuckoo_full_lpn_field import PolyState
from rev_cuckoo_ring_lpn_factor_attack import (
    elementary_symmetric,
    exact_dual_check_bits,
    expected_reduced_weight,
    product_weighted_support_probability,
    project_block_distribution,
    projected_occupancy_scores,
    reduced_support,
    reduced_weight_distribution,
    mixed_statistical_decoding_bits,
    mixed_exact_dual_bits,
    statistical_decoding_bits,
)


def test_elementary_symmetric_is_binomial_for_unit_weights() -> None:
    for width in range(1, 9):
        for degree in range(width + 1):
            assert elementary_symmetric([1.0] * width, degree) == math.comb(
                width, degree
            )


def test_uniform_joint_inclusion_is_hypergeometric() -> None:
    scores = (0.05, 0.1, 0.2, 0.25, 0.4, 0.0)
    support = frozenset((1, 4))
    actual = product_weighted_support_probability(scores, support, 4, 0.0)
    expected = math.comb(4, 2) / math.comb(6, 4)
    assert math.isclose(actual, expected, rel_tol=1e-15)


def test_weighted_joint_inclusion_matches_enumeration() -> None:
    scores = (0.1, 0.2, 0.3, 0.4, 0.5)
    support = frozenset((1, 3))
    kept = 3
    beta = 1.5
    numerator = 0.0
    denominator = 0.0
    for mask in range(1 << len(scores)):
        subset = frozenset(
            index for index in range(len(scores)) if mask & (1 << index)
        )
        if len(subset) != kept:
            continue
        weight = math.prod(scores[index] ** beta for index in subset)
        denominator += weight
        if support <= subset:
            numerator += weight
    actual = product_weighted_support_probability(
        scores, support, kept, beta
    )
    assert math.isclose(actual, numerator / denominator, rel_tol=1e-14)


def test_projection_uses_absolute_ring_coordinate() -> None:
    probabilities = (0.1, 0.2, 0.3, 0.4)
    assert project_block_distribution(probabilities, 1, 4, 3) == (
        0.3,
        0.1 + 0.4,
        0.2,
    )


def test_projected_occupancy_combines_independent_blocks() -> None:
    distributions = ((1.0, 0.0), (0.5, 0.5))
    actual = projected_occupancy_scores(distributions, 2, 2)
    assert actual == (1.0, 0.5)


def test_reduced_support_folds_and_cancels_coefficients() -> None:
    state = PolyState(
        positions=(0, 0),
        coefficients=(3, 3),
        polynomial=(),
    )
    assert reduced_support(state, 2, 2, 17) == frozenset()
    state = PolyState(
        positions=(0, 1),
        coefficients=(3, 4),
        polynomial=(),
    )
    assert reduced_support(state, 2, 2, 17) == frozenset((0, 1))


def test_paper_parameter_calibration() -> None:
    reduced_weight = expected_reduced_weight(4, 16, 128)
    assert math.isclose(reduced_weight, 60.3833115418203, rel_tol=1e-14)
    assert math.isclose(
        statistical_decoding_bits(4, reduced_weight, 128),
        128.44987916747,
        rel_tol=1e-14,
    )
    assert math.isclose(
        exact_dual_check_bits(4, 60, 128),
        146.77355886479714,
        rel_tol=1e-14,
    )


def test_reduced_weight_distribution_and_weak_tail() -> None:
    distribution = reduced_weight_distribution(4, 16, 128)
    assert math.isclose(sum(distribution), 1.0, rel_tol=1e-14)
    assert math.isclose(
        sum(weight * probability for weight, probability in enumerate(distribution)),
        expected_reduced_weight(4, 16, 128),
        rel_tol=1e-14,
    )
    assert math.isclose(sum(distribution[:61]), 0.49808554174791797)
    unconditioned, mass = mixed_statistical_decoding_bits(
        4, 128, distribution
    )
    assert math.isclose(mass, 1.0, rel_tol=1e-14)
    assert math.isclose(unconditioned, 122.4817691212011)
    conditioned, mass = mixed_statistical_decoding_bits(
        4, 128, distribution, 61
    )
    assert math.isclose(mass, 0.5019144582520833)
    assert math.isclose(conditioned, 130.5717350459001)
    exact_unconditioned, mass = mixed_exact_dual_bits(
        4, 128, distribution
    )
    assert math.isclose(mass, 1.0, rel_tol=1e-14)
    assert math.isclose(exact_unconditioned, 136.43055398206982)
    exact_conditioned, mass = mixed_exact_dual_bits(
        4, 128, distribution, 61
    )
    assert math.isclose(mass, 0.5019144582520833)
    assert math.isclose(exact_conditioned, 150.52635859933227)
