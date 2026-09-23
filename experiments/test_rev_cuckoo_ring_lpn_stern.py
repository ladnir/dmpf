import math

from rev_cuckoo_basic_ring_lpn_isd import BALANCED_BETAS
from rev_cuckoo_ring_lpn_stern import (
    SternTrial,
    estimate_attack,
    poisson_binomial_probability,
    regular_information_widths,
    regular_stern_success,
    stern_iteration_cost,
    uniform_regular_stern_success,
)


def test_poisson_binomial_matches_binomial_distribution() -> None:
    probabilities = (0.3,) * 8
    expected = math.comb(8, 3) * 0.3**3 * 0.7**5
    assert math.isclose(
        poisson_binomial_probability(probabilities, 3),
        expected,
        rel_tol=1e-12,
    )


def test_integral_information_widths_are_balanced() -> None:
    widths = regular_information_widths(64, 4096, 4, 17)
    assert sum(widths) == 64 * 3072 + 17
    assert max(widths) - min(widths) == 1
    assert sum(widths[:32]) - sum(widths[32:]) == 1


def test_uniform_success_matches_closed_form_without_extension() -> None:
    blocks = 64
    p = 4
    information_fraction = 0.75
    direct = regular_stern_success(
        (information_fraction,) * blocks, p
    )
    closed = uniform_regular_stern_success(
        blocks, information_fraction, p
    )
    assert math.isclose(direct, closed, rel_tol=1e-12)


def test_beta_zero_estimate_is_exactly_uniform() -> None:
    polys = 4
    weight = 2
    block_length = 16
    blocks = polys * weight
    trial = SternTrial(
        tuple(
            (1.0 / polys,) * blocks
            for _ in BALANCED_BETAS
        )
    )
    rows = estimate_attack(
        [trial], polys, weight, block_length, 17, p=2, ell=3
    )
    leakage = next(
        row
        for row in rows
        if row.selector == "leakage" and row.beta == 0.0
    )
    uniform = next(row for row in rows if row.selector == "uniform")
    assert math.isclose(
        leakage.mean_inverse_success_log2,
        uniform.mean_inverse_success_log2,
        abs_tol=1e-12,
    )
    assert math.isclose(
        leakage.expected_work_log2,
        uniform.expected_work_log2,
        abs_tol=1e-12,
    )


def test_production_stern_split_is_costlier_than_prange() -> None:
    widths = regular_information_widths(64, 4096, 4, 0)
    _, _, iteration = stern_iteration_cost(
        widths, modulus=65537, p=2, ell=0
    )
    success = regular_stern_success((0.75,) * 64, p=2)
    assert iteration - math.log2(success) > 140.0


if __name__ == "__main__":
    test_poisson_binomial_matches_binomial_distribution()
    test_integral_information_widths_are_balanced()
    test_uniform_success_matches_closed_form_without_extension()
    test_beta_zero_estimate_is_exactly_uniform()
    test_production_stern_split_is_costlier_than_prange()
