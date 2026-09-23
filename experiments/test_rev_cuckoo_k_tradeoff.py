import math
import random

from rev_cuckoo_k_tradeoff import (
    conditional_inverse_mean_weights,
    exact_deficit_masses,
    exact_occupancy_weights,
    leakage_metrics,
    occupancy_weights,
    proposal_distribution,
    selected_distribution,
    expected_placements,
    finite_race_factors,
    exact_occupancy_race_factors,
)


def toy_uniform():
    return {(0, 2, 3): 0.1, (2, 1, 1): 0.3, (4, 0, 0): 0.6}


def test_reverse_distribution_is_size_biased() -> None:
    proposal = proposal_distribution(toy_uniform())
    assert proposal.get((0, 2, 3), 0) == 0
    assert math.isclose(proposal[(2, 1, 1)], 0.2)
    assert math.isclose(proposal[(4, 0, 0)], 0.8)


def test_occupancy_lookup_inverts_conditional_mean() -> None:
    weights = occupancy_weights(toy_uniform())
    assert math.isclose(weights[1], 0.5)
    assert math.isclose(weights[0], 0.25)


def test_conditional_lookup_accepts_richer_scores() -> None:
    weights = conditional_inverse_mean_weights(
        toy_uniform(), lambda feature: (feature[1], feature[2])
    )
    assert math.isclose(weights[(1, 1)], 0.5)
    assert math.isclose(weights[(0, 0)], 0.25)


def test_k_one_is_exactly_raw() -> None:
    proposal = proposal_distribution(toy_uniform())
    selected = selected_distribution(
        proposal, 1, 10, lambda feature: 1 / feature[0], random.Random(0)
    )
    assert selected == proposal


def test_feasible_limit_has_expected_kl() -> None:
    uniform = toy_uniform()
    feasible = {(2, 1, 1): 1 / 3, (4, 0, 0): 2 / 3}
    metrics = leakage_metrics(uniform, feasible, 8)
    assert math.isclose(metrics.kl_bits, -math.log2(0.9))
    assert math.isclose(metrics.search_gain_bits, 8 * metrics.pair_chernoff_bits)


def test_equal_partition_mean_matches_closed_form() -> None:
    assert math.isclose(
        expected_placements(4, (3, 3)),
        math.prod((6 - index) / 3 for index in range(4)),
    )


def test_exact_occupancy_table_matches_two_item_enumeration() -> None:
    bins = (2, 3)
    uniform = {}
    graph_count = (bins[0] * bins[1]) ** 2
    for left0 in range(bins[0]):
        for right0 in range(bins[1]):
            for left1 in range(bins[0]):
                for right1 in range(bins[1]):
                    placements = (
                        4
                        - int(left0 == left1)
                        - int(right0 == right1)
                    )
                    deficit = int(left0 == left1) + int(right0 == right1)
                    feature = (placements, deficit, deficit)
                    uniform[feature] = uniform.get(feature, 0.0) + 1 / graph_count
    sampled_formula = occupancy_weights(uniform)
    exact = exact_occupancy_weights(2, bins)
    assert exact.keys() == sampled_formula.keys()
    for deficit in exact:
        assert math.isclose(exact[deficit], sampled_formula[deficit])


def test_exact_weighted_deficit_mass_has_correct_mean() -> None:
    _, weighted = exact_deficit_masses(4, (3, 5))
    assert math.isclose(sum(weighted.values()), expected_placements(4, (3, 5)))


def test_finite_race_factor_matches_two_point_closed_form() -> None:
    proposal = {"low": 0.25, "high": 0.75}
    rates = {"low": 1.0, "high": 3.0}
    factors = finite_race_factors(proposal, 2, rates)
    expected_low = 2.0 * (0.25 * 0.5 + 0.75 * 0.25)
    expected_high = 2.0 * (0.25 * 0.75 + 0.75 * 0.5)
    assert math.isclose(factors["low"], expected_low, rel_tol=1e-9)
    assert math.isclose(factors["high"], expected_high, rel_tol=1e-9)
    assert math.isclose(
        sum(proposal[score] * factors[score] for score in proposal),
        1.0,
        rel_tol=1e-9,
    )


def test_occupancy_race_k_one_is_raw() -> None:
    factors = exact_occupancy_race_factors(4, (3, 3), 1)
    assert all(factor == 1.0 for factor in factors.values())


def test_occupancy_race_factors_normalize() -> None:
    deficit_mass, placement_mass = exact_deficit_masses(4, (4, 4))
    del deficit_mass
    total = sum(placement_mass.values())
    proposal = {
        deficit: mass / total
        for deficit, mass in placement_mass.items()
        if mass > 0.0
    }
    factors = exact_occupancy_race_factors(4, (4, 4), 8)
    assert math.isclose(
        sum(proposal[deficit] * factors[deficit] for deficit in proposal),
        1.0,
        rel_tol=1e-9,
    )
