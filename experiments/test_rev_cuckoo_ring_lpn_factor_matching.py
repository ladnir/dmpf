import math

from rev_cuckoo_ring_lpn_factor_matching import (
    completion_probability,
    particle_occupancy_scores,
    planted_completion_probability,
    projected_candidate_support,
    uniform_zero_set_probability,
)


def test_projected_candidate_support_folds_collisions() -> None:
    assert projected_candidate_support((0, 1, 2, 3), 4, 4) == frozenset(
        (0, 1, 2, 3)
    )
    assert projected_candidate_support((0, 0, 0, 0), 4, 4) == frozenset((0,))


def test_particle_occupancy_scores_count_each_residue_once() -> None:
    scores = particle_occupancy_scores(((0, 0), (0, 1)), 2, 2)
    assert scores == (1.0, 0.5)


def test_completion_probability_matches_direct_count() -> None:
    # Width six, plant two coordinates including one of the two true ones,
    # then choose one of the four remaining coordinates.
    assert completion_probability(6, 3, 2, 2, 1) == 0.25


def test_empty_plant_is_the_uniform_selector() -> None:
    assert math.isclose(
        completion_probability(8, 3, 0, 2, 0),
        uniform_zero_set_probability(8, 3, 2),
    )


def test_planted_completion_averages_the_empirical_population() -> None:
    populations = (((0,), (1,)),)
    truth = ((0,),)
    # Half the particles plant the true coordinate and succeed certainly.
    # The other half plant a false coordinate and hit truth with probability
    # 1/3 while filling one of the three remaining coordinates.
    assert math.isclose(
        planted_completion_probability(populations, truth, 4, 4, 2),
        2.0 / 3.0,
    )
