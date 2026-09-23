import math
import random

from rev_cuckoo_full_lpn import support_lists
from rev_cuckoo_full_lpn_field import (
    PolyState,
    field_trial,
    make_model,
    product_active_lists,
    product_active_lists_one_hidden,
    regular_poly,
)
from rev_cuckoo_regular_noise import LeakageChannel


def test_occupancy_k1_is_raw_ring_lpn_channel() -> None:
    model = make_model(polys=2, weight=2, block_length=2, modulus=3)
    raw = LeakageChannel.build("raw", 1, model.weight, 2)
    occupancy = LeakageChannel.build("occupancy", 1, model.weight, 2)
    known = (0, 1, 1, 0)
    candidate = (1, 0, 0, 1)
    lists = support_lists(known, candidate, model.polys, model.weight)
    rng = random.Random(29)
    descriptors = tuple(
        raw.sample(active, 2 * model.block_length, rng) for active in lists
    )

    assert raw.candidate_log_weight(lists, descriptors) == (
        occupancy.candidate_log_weight(lists, descriptors)
    )


def test_ring_lpn_information_decomposition() -> None:
    model = make_model(polys=2, weight=2, block_length=2, modulus=3)
    channel = LeakageChannel.build("occupancy", 4, model.weight, 2)
    result = field_trial(model, channel, random.Random(71))

    assert math.isfinite(result.total_information)
    assert math.isfinite(result.decoding_information)
    assert math.isfinite(result.distinguishing_information)
    assert math.isclose(
        result.total_information,
        result.decoding_information + result.distinguishing_information,
        abs_tol=1e-10,
    )
    assert result.posterior_rank >= 1.0
    assert result.base_rank >= 1.0

    product_result = field_trial(
        model,
        channel,
        random.Random(71),
        product_values=True,
    )
    assert math.isclose(
        product_result.total_information,
        product_result.decoding_information
        + product_result.distinguishing_information,
        abs_tol=1e-10,
    )


def test_product_rows_apply_deduplication_and_zero_cancellation() -> None:
    model = make_model(polys=2, weight=2, block_length=2, modulus=3)
    polynomial = regular_poly((0, 0), (1, 1), 2, 3)
    known_poly = PolyState((0, 0), (1, 1), polynomial)
    hidden_index = model.poly_states.index(known_poly)

    lists = product_active_lists(
        model,
        (known_poly, known_poly),
        (hidden_index, hidden_index),
    )
    assert lists == ((), (0,), (), (0,), (), (0,), (), (0,))

    channel = LeakageChannel.build("occupancy", 4, model.weight, 2)
    descriptor = channel.sample((), 4, random.Random(97))
    assert channel.log_likelihood((), descriptor) == 0.0


def test_product_rows_match_support_geometry_without_cancellation() -> None:
    model = make_model(polys=2, weight=2, block_length=4, modulus=3)
    known_poly = PolyState(
        (0, 0),
        (1, 1),
        regular_poly((0, 0), (1, 1), 4, 3),
    )
    hidden_poly = PolyState(
        (0, 2),
        (1, 1),
        regular_poly((0, 2), (1, 1), 4, 3),
    )
    hidden_index = model.poly_states.index(hidden_poly)
    known_positions = known_poly.positions * 2
    hidden_positions = hidden_poly.positions * 2

    assert product_active_lists(
        model,
        (known_poly, known_poly),
        (hidden_index, hidden_index),
    ) == support_lists(
        known_positions,
        hidden_positions,
        model.polys,
        model.weight,
    )


def test_one_hidden_factorization_preserves_descriptor_order() -> None:
    model = make_model(polys=2, weight=2, block_length=3, modulus=3)
    known = (model.poly_states[3], model.poly_states[11])
    hidden_state = (7, 19)
    full_lists = product_active_lists(model, known, hidden_state)

    for hidden_poly_index, state_index in enumerate(hidden_state):
        expected = tuple(
            full_lists[
                (known_poly_index * model.polys + hidden_poly_index)
                * model.weight
                + output_block
            ]
            for known_poly_index in range(model.polys)
            for output_block in range(model.weight)
        )
        assert product_active_lists_one_hidden(
            model,
            known,
            model.poly_states[state_index],
        ) == expected
