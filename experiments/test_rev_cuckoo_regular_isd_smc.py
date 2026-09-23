import math
import random

from rev_cuckoo_regular_isd_smc import (
    Particle,
    SmcResult,
    anchored_rectangle_selection,
    combine_replicates,
    conditional_value_logs,
    diverse_posterior_anchors,
    factor_active_set,
    factor_log_likelihood,
    normalized_log_distribution,
    normalized_weights,
    posterior_single_offset_set,
    smc_continue_to_squared_target,
    systematic_indices,
)
from rev_cuckoo_regular_noise import LeakageChannel, support_lists_one_unknown


def test_factor_active_set_matches_regular_product_geometry() -> None:
    polys = 2
    weight = 3
    known = (0, 3, 6, 1, 4, 7)
    candidate = (2, 4, 7)
    lists = support_lists_one_unknown(known, candidate, polys, weight)
    for factor, expected in enumerate(lists):
        known_poly, output_block = divmod(factor, weight)
        assert set(
            factor_active_set(
                known,
                candidate,
                known_poly,
                output_block,
                weight,
            )
        ) == set(expected)


def test_normalized_weights_preserve_hard_zeros() -> None:
    weights = normalized_weights((0.0, -math.inf, -1.0))
    assert weights == [2.0 / 3.0, 0.0, 1.0 / 3.0]
    assert normalized_weights((-math.inf, -math.inf)) == []


def test_systematic_resampling_is_identity_for_uniform_weights() -> None:
    assert systematic_indices([0.25] * 4, random.Random(3)) == [0, 1, 2, 3]


def test_failed_smc_replicate_contributes_zero_not_missing_mass() -> None:
    results = (
        SmcResult(1.0, 1.0, 1.0, 1.0, 1.0),
        SmcResult(-math.inf, 0.0, 0.0, 0.0, 0.0),
    )
    estimate, _ = combine_replicates(results)
    assert estimate == 0.0


def test_normalized_log_distribution_preserves_hard_zeros() -> None:
    distribution = normalized_log_distribution((0.0, -math.inf, -1.0))
    assert distribution == [2.0 / 3.0, 0.0, 1.0 / 3.0]
    assert normalized_log_distribution((-math.inf, -math.inf)) == [0.0, 0.0]


def test_conditional_value_logs_match_direct_full_likelihood() -> None:
    polys = 2
    weight = 3
    block_length = 8
    known = (0, 3, 6, 1, 4, 7)
    context = (2, 4, 7)
    channel = LeakageChannel.build("occupancy", 4, weight, 4)
    rng = random.Random(17)
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng)
        for active in support_lists_one_unknown(
            known, context, polys, weight
        )
    )
    variable = 1
    actual = conditional_value_logs(
        known,
        context,
        variable,
        descriptors,
        weight,
        block_length,
        channel,
    )
    expected = []
    for value in range(block_length):
        candidate = list(context)
        candidate[variable] = value
        factor_logs = [
            factor_log_likelihood(
                known,
                candidate,
                factor,
                descriptors,
                weight,
                channel,
            )
            for factor in range(len(descriptors))
        ]
        expected.append(
            -math.inf
            if any(not math.isfinite(item) for item in factor_logs)
            else sum(factor_logs)
        )
    for observed, reference in zip(actual, expected):
        if math.isfinite(reference):
            assert math.isclose(observed, reference, abs_tol=1e-12)
        else:
            assert observed == reference


def test_diverse_anchors_prefer_frequency_then_hamming_distance() -> None:
    particles = [
        Particle((0, 0, 0), [], 0.0),
        Particle((0, 0, 0), [], 0.0),
        Particle((0, 0, 1), [], 0.0),
        Particle((1, 1, 1), [], 0.0),
    ]
    anchors = diverse_posterior_anchors(particles, 2)
    assert anchors[0] == (0, 0, 0)
    assert anchors[1] == (1, 1, 1)


def test_anchored_rectangle_has_fixed_block_cardinality() -> None:
    polys = 2
    weight = 3
    block_length = 8
    known = (0, 3, 6, 1, 4, 7)
    anchor = (2, 4, 7)
    channel = LeakageChannel.build("occupancy", 4, weight, 4)
    rng = random.Random(23)
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng)
        for active in support_lists_one_unknown(
            known, anchor, polys, weight
        )
    )
    rectangle = anchored_rectangle_selection(
        known,
        descriptors,
        anchor,
        weight,
        block_length,
        block_length // polys,
        channel,
    )
    assert len(rectangle) == weight
    assert all(len(block) == block_length // polys for block in rectangle)


def test_offset_set_selector_has_exact_width_and_holdout_score() -> None:
    polys = 2
    weight = 3
    block_length = 8
    known = (0, 3, 6, 1, 4, 7)
    true_candidate = (2, 4, 7)
    channel = LeakageChannel.build("occupancy", 1, weight, 4)
    rng = random.Random(29)
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng)
        for active in support_lists_one_unknown(
            known, true_candidate, polys, weight
        )
    )
    selection, guess, probability, selector, evaluator = (
        posterior_single_offset_set(
            known,
            descriptors,
            weight,
            block_length,
            channel,
            128,
            2,
            16,
            5,
            random.Random(31),
        )
    )
    variable, values = guess
    assert len(selection) == weight
    assert selection[variable] == frozenset(values)
    assert len(selection[variable]) == 5
    assert 0.0 <= probability <= 1.0
    assert selector.mean_ess_fraction > 0.0
    assert evaluator.mean_ess_fraction > 0.0


def test_squared_target_continuation_adds_one_likelihood_copy() -> None:
    polys = 2
    weight = 3
    block_length = 8
    known = (0, 3, 6, 1, 4, 7)
    candidate = (2, 4, 7)
    channel = LeakageChannel.build("occupancy", 1, weight, 4)
    rng = random.Random(41)
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng)
        for active in support_lists_one_unknown(
            known, candidate, polys, weight
        )
    )
    logs = [
        factor_log_likelihood(
            known, candidate, factor, descriptors, weight, channel
        )
        for factor in range(len(descriptors))
    ]
    total = sum(logs)
    diagnostic, particles = smc_continue_to_squared_target(
        known,
        descriptors,
        tuple((value,) for value in candidate),
        weight,
        channel,
        (Particle(candidate, logs, total),),
        0,
        4,
        random.Random(43),
    )
    assert math.isclose(diagnostic.log_mean_weight, total, abs_tol=1e-12)
    assert len(particles) == 1
    assert math.isclose(particles[0].total_log, 2.0 * total, abs_tol=1e-12)
