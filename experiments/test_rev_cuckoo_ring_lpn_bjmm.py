import itertools
import math

from rev_cuckoo_ring_lpn_bjmm import (
    estimate_candidate,
    estimate_mmt_from_representation_count,
    log2_binomial,
    representation_log2,
)


def test_binary_count_reduces_to_regular_isd_formula() -> None:
    blocks = 12
    target_weight = 6
    epsilon = 2
    addend_weight = target_weight // 2 + epsilon
    width = 17
    actual = representation_log2(
        target_weight, addend_weight, blocks, width, modulus=2
    )
    expected = (
        log2_binomial(target_weight, target_weight // 2)
        + log2_binomial(blocks - target_weight, epsilon)
        + epsilon * math.log2(width)
    )
    assert math.isclose(actual, expected, abs_tol=1e-12)


def block_states(width: int, modulus: int):
    yield (0,) * width
    for position in range(width):
        for value in range(1, modulus):
            block = [0] * width
            block[position] = value
            yield tuple(block)


def test_qary_count_matches_exhaustive_regular_representations() -> None:
    blocks = 3
    width = 2
    modulus = 3
    target = (
        (1, 0),
        (0, 2),
        (0, 0),
    )
    states = tuple(block_states(width, modulus))
    for addend_weight in range(blocks + 1):
        count = 0
        for left in itertools.product(states, repeat=blocks):
            if sum(any(block) for block in left) != addend_weight:
                continue
            right = tuple(
                tuple(
                    (target_value - left_value) % modulus
                    for target_value, left_value in zip(
                        target_block, left_block
                    )
                )
                for target_block, left_block in zip(target, left)
            )
            if any(sum(value != 0 for value in block) > 1 for block in right):
                continue
            if sum(any(block) for block in right) != addend_weight:
                continue
            count += 1
        modeled = representation_log2(
            target_weight=2,
            addend_weight=addend_weight,
            blocks=blocks,
            block_width=width,
            modulus=modulus,
        )
        if count == 0:
            assert modeled == -math.inf
        else:
            assert math.isclose(2.0**modeled, count, abs_tol=1e-10)


def test_zero_level_is_regular_prange_at_production_parameters() -> None:
    row = estimate_candidate(
        polys=4,
        weight=16,
        block_length=4096,
        modulus=65537,
        p=0,
        p_x=0,
        p_y=0,
        ell=0,
    )
    assert row is not None
    assert row.iteration_log2 == 0.0
    assert row.inverse_success_log2 == 128.0
    assert row.expected_work_log2 == 128.0


def test_production_mmt_screen_matches_recorded_optimum() -> None:
    blocks = 64
    ell = 5
    width = (3 * 16 * 4096 + ell) / blocks
    log2_rx = representation_log2(
        target_weight=2,
        addend_weight=2,
        blocks=blocks,
        block_width=width,
        modulus=65537,
    )
    row = estimate_mmt_from_representation_count(
        blocks=blocks,
        block_length=4096,
        modulus=65537,
        information_width=width,
        p=2,
        p_x=2,
        ell=ell,
        log2_rx=log2_rx,
    )
    assert row is not None
    assert math.isclose(row.expected_work_log2, 146.444545, abs_tol=1e-6)


if __name__ == "__main__":
    test_binary_count_reduces_to_regular_isd_formula()
    test_qary_count_matches_exhaustive_regular_representations()
    test_zero_level_is_regular_prange_at_production_parameters()
    test_production_mmt_screen_matches_recorded_optimum()
