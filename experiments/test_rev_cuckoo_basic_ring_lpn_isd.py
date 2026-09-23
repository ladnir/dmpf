import itertools
import math
import random

import numpy as np

from rev_cuckoo_basic_ring_lpn_isd import (
    best_mixed_expected_work_gain,
    best_weighted_half_success_gain,
    balanced_half_coordinate_probabilities,
    balanced_half_inclusion_probabilities,
    balanced_stratified_coordinate_probabilities,
    balanced_stratified_inclusion_probabilities,
    coefficient_marginalized_posterior,
    one_trial,
    posterior_center_success_probability,
    sample_balanced_stratified_subset,
    sample_balanced_stratified_subset_conditioned,
    states_by_support,
    weighted_subset_inclusion_probabilities,
)
from rev_cuckoo_full_lpn_field import make_model, product_active_lists
from rev_cuckoo_regular_noise import LeakageChannel
from rev_cuckoo_regular_isd import (
    CollisionFactor,
    matrix_rank,
    pair_center_coordinate_ascent,
    selected_matrix_full_rank,
)
from rev_cuckoo_ring_lpn_rank_audit import (
    matrix_rank_mod,
    structured_selected_matrix,
)
from rev_cuckoo_ring_lpn_matching_isd import (
    condition_contexts,
    matching_distribution,
    matching_sequential_polynomial_success,
    pair_proposal_target_particles,
)
from rev_cuckoo_regular_isd_smc import factor_active_set
from rev_cuckoo_ring_lpn_sequential_isd import (
    sample_sequential_polynomial_selection,
    sequential_polynomial_success,
)


def test_coefficient_marginalized_posterior_normalizes() -> None:
    model = make_model(polys=2, weight=2, block_length=2, modulus=3)
    candidates = tuple(itertools.product(range(2), repeat=2))
    grouped = states_by_support(model)
    known = (model.poly_states[2], model.poly_states[9])
    hidden = (5, 11)
    channel = LeakageChannel.build("occupancy", 4, model.weight, 2)
    descriptors = tuple(
        channel.sample(active, 4, random.Random(101 + index))
        for index, active in enumerate(
            product_active_lists(model, known, hidden)
        )
    )
    factor = tuple(descriptors[index] for index in (0, 1, 4, 5))
    probabilities, marginals = coefficient_marginalized_posterior(
        model,
        candidates,
        grouped,
        known,
        factor,
        channel,
    )

    assert math.isclose(sum(probabilities), 1.0, abs_tol=1e-12)
    assert all(math.isclose(sum(row), 1.0, abs_tol=1e-12) for row in marginals)


def test_basic_ring_lpn_isd_trial_has_valid_probabilities() -> None:
    model = make_model(polys=2, weight=2, block_length=2, modulus=3)
    candidates = tuple(itertools.product(range(2), repeat=2))
    channel = LeakageChannel.build("raw", 1, model.weight, 2)
    result = one_trial(
        model,
        candidates,
        states_by_support(model),
        channel,
        joint_restarts=1,
        pair_slope=-0.5,
        pair_restarts=1,
        pair_rounds=8,
        pair_damping=0.5,
        rng=random.Random(13),
    )

    assert 0.0 <= result.marginal_support <= 1.0
    assert 0.0 <= result.pair_support <= 1.0
    assert result.joint_support + 1e-12 >= result.marginal_support
    assert 0 <= result.pair_hits <= 2 * model.weight


def test_mixed_policy_prefers_uniform_without_signal() -> None:
    hits = [0, 1, 2, 1] * 16
    gain, alpha = best_mixed_expected_work_gain(
        hits,
        variables=2,
        alphas=(0.0, 0.25, 0.5, 0.75),
    )
    assert math.isclose(gain, 0.0, abs_tol=1e-12)
    assert alpha == 0.0


def test_product_weighted_subset_inclusion_probabilities() -> None:
    uniform = weighted_subset_inclusion_probabilities(
        (0.1, 0.2, 0.3, 0.4), kept=2, beta=0.0
    )
    assert all(math.isclose(value, 0.5, abs_tol=1e-12) for value in uniform)

    biased = weighted_subset_inclusion_probabilities(
        (0.1, 0.2, 0.3, 0.4), kept=2, beta=1.0
    )
    assert math.isclose(sum(biased), 2.0, abs_tol=1e-12)
    assert tuple(sorted(biased)) == biased


def test_uniform_weighted_policy_has_zero_half_success_gain() -> None:
    successes = [(0.25, 0.25)] * 8
    gain, beta = best_weighted_half_success_gain(
        successes,
        variables=2,
        betas=(0.0, 1.0),
    )
    assert math.isclose(gain, 0.0, abs_tol=1e-12)
    assert beta == 0.0


def test_balanced_half_selector_has_exact_cardinality() -> None:
    probabilities = (0.4, 0.3, 0.2, 0.1)
    uniform = balanced_half_inclusion_probabilities(probabilities, beta=0.0)
    assert uniform == (0.5, 0.5, 0.5, 0.5)

    biased = balanced_half_inclusion_probabilities(probabilities, beta=1.0)
    assert math.isclose(sum(biased), 2.0, abs_tol=1e-12)
    assert biased[0] > biased[2]
    assert biased[1] > biased[3]
    for coordinate in range(len(probabilities)):
        direct = balanced_half_coordinate_probabilities(
            probabilities, coordinate, (0.0, 1.0, 2.0)
        )
        for beta_index, beta in enumerate((0.0, 1.0, 2.0)):
            full = balanced_half_inclusion_probabilities(probabilities, beta)
            assert abs(direct[beta_index] - full[coordinate]) < 1e-15


def test_balanced_stratified_selector_for_four_polynomials() -> None:
    probabilities = (0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
    uniform = balanced_stratified_inclusion_probabilities(
        probabilities, kept=2, beta=0.0
    )
    assert uniform == (0.25,) * 8

    biased = balanced_stratified_inclusion_probabilities(
        probabilities, kept=2, beta=1.0
    )
    assert math.isclose(sum(biased), 2.0, abs_tol=1e-12)
    assert biased[0] > biased[2] > biased[4] > biased[6]
    assert biased[1] > biased[3] > biased[5] > biased[7]
    for coordinate in range(len(probabilities)):
        direct = balanced_stratified_coordinate_probabilities(
            probabilities, coordinate, kept=2, betas=(0.0, 1.0, 2.0)
        )
        for beta_index, beta in enumerate((0.0, 1.0, 2.0)):
            full = balanced_stratified_inclusion_probabilities(
                probabilities, kept=2, beta=beta
            )
            assert abs(direct[beta_index] - full[coordinate]) < 1e-15
    rng = random.Random(101)
    for beta in (0.0, 0.5, 2.0):
        for _ in range(20):
            selected = sample_balanced_stratified_subset(
                probabilities, kept=2, beta=beta, rng=rng
            )
            assert len(selected) == 2
        for required in range(len(probabilities)):
            selected, inclusion = (
                sample_balanced_stratified_subset_conditioned(
                    probabilities,
                    kept=2,
                    beta=beta,
                    required=required,
                    rng=rng,
                )
            )
            direct = balanced_stratified_coordinate_probabilities(
                probabilities, required, kept=2, betas=(beta,)
            )[0]
            assert required in selected
            assert len(selected) == 2
            assert math.isclose(inclusion, direct, abs_tol=1e-15)


def test_pair_center_coordinate_ascent_avoids_collision() -> None:
    factor = CollisionFactor(
        labels=((0, 1, 2), (0, 1, 2)),
        label_count=3,
    )
    center = pair_center_coordinate_ascent(
        (factor,), initial=(0, 0), slope=-1.0, rounds=2
    )
    assert center[0] != center[1]


def test_uniform_posterior_center_is_uniform_prange() -> None:
    candidates = tuple(itertools.product(range(4), repeat=2))
    probabilities = [1.0 / len(candidates)] * len(candidates)
    success = posterior_center_success_probability(
        candidates, probabilities, kept=2, block_length=4
    )
    assert math.isclose(success, 0.25, abs_tol=1e-12)


def test_vectorized_structured_rank_matches_reference() -> None:
    selections = (
        (frozenset((0,)), frozenset((1,))),
        (frozenset((2,)), frozenset((3,))),
        (frozenset((1,)), frozenset((0,))),
        (frozenset((3,)), frozenset((2,))),
    )
    multipliers = (
        (1, 3, 5, 7, 9, 11, 13, 15),
        (2, 4, 6, 8, 10, 12, 14, 16),
        (4, 1, 9, 2, 13, 5, 7, 11),
    )
    modulus = 17
    reference = selected_matrix_full_rank(
        selections, multipliers, block_length=4, modulus=modulus
    )
    matrix = structured_selected_matrix(
        selections,
        tuple(np.asarray(values, dtype=np.uint64) for values in multipliers),
        block_length=4,
        modulus=modulus,
    )
    assert (matrix_rank_mod(matrix, modulus) == 8) == reference
    rng = random.Random(707)
    for _ in range(20):
        random_matrix = np.asarray(
            [
                [rng.randrange(modulus) for _ in range(12)]
                for _ in range(12)
            ],
            dtype=np.uint64,
        )
        expected = matrix_rank(random_matrix.tolist(), modulus)
        assert matrix_rank_mod(random_matrix, modulus) == expected


def test_sequential_selector_beta_zero_is_exact_baseline() -> None:
    distributions = tuple(
        tuple([1.0 / 8.0] * 8) for _ in range(3)
    )
    factor = CollisionFactor(
        labels=tuple(tuple(value % 2 for value in range(8)) for _ in range(3)),
        label_count=2,
    )
    success = sequential_polynomial_success(
        true_positions=(1, 3, 5),
        distributions=distributions,
        factors=(factor,),
        kept=2,
        beta=0.0,
        slope=-0.5,
        particles=8,
        rng=random.Random(909),
    )
    assert math.isclose(success, 4.0 ** -3, abs_tol=1e-15)
    selection = sample_sequential_polynomial_selection(
        distributions=distributions,
        factors=(factor,),
        kept=2,
        beta=0.5,
        slope=-0.5,
        rng=random.Random(910),
    )
    assert len(selection) == 3
    assert all(len(block) == 2 for block in selection)


def test_matching_selector_preserves_full_support_and_beta_zero() -> None:
    pair = np.asarray([0.05, 0.15, 0.3, 0.5], dtype=np.float64)
    contexts = ((0, 1), (0, 2), (3, 1), (3, 2))
    mixed = matching_distribution(pair, contexts, 0, 4, 0.75)
    assert math.isclose(float(np.sum(mixed)), 1.0, abs_tol=1e-15)
    assert np.all(mixed > 0.0)

    conditioned = condition_contexts(
        contexts, 0, frozenset((0, 2)), random.Random(911)
    )
    assert len(conditioned) == len(contexts)
    assert all(context[0] == 0 for context in conditioned)

    distributions = tuple(tuple([1.0 / 8.0] * 8) for _ in range(3))
    factor = CollisionFactor(
        labels=tuple(tuple(value % 2 for value in range(8)) for _ in range(3)),
        label_count=2,
    )
    success = matching_sequential_polynomial_success(
        true_positions=(1, 3, 5),
        distributions=distributions,
        factors=(factor,),
        posterior_contexts=((0, 0, 0), (7, 7, 7)),
        kept=2,
        beta=0.0,
        slope=-0.5,
        posterior_mix=0.75,
        particles=8,
        rng=random.Random(912),
    )
    assert math.isclose(success, 4.0 ** -3, abs_tol=1e-15)


def test_pair_proposal_matching_population_is_well_formed() -> None:
    rng = random.Random(913)
    weight = 2
    known = (0, 1, 2, 3)
    planted = (1, 2)
    channel = LeakageChannel.build("raw", 1, weight, 2)
    descriptors = tuple(
        channel.sample(
            factor_active_set(
                known,
                planted,
                known_poly,
                output_block,
                weight,
            ),
            8,
            rng,
        )
        for known_poly in range(2)
        for output_block in range(weight)
    )
    diagnostic, particles = pair_proposal_target_particles(
        known=known,
        descriptors=descriptors,
        distributions=((0.1, 0.2, 0.3, 0.4),) * weight,
        weight=weight,
        channel=channel,
        particle_count=64,
        moves_per_stage=2,
        rng=random.Random(914),
    )
    assert len(particles) == 64
    assert math.isfinite(diagnostic.log_mean_weight)
    assert 0.0 < diagnostic.minimum_ess_fraction <= 1.0
    assert all(
        len(particle.candidate) == weight
        and all(0 <= value < 4 for value in particle.candidate)
        for particle in particles
    )
