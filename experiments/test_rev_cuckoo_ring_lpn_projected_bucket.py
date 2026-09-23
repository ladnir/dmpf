import itertools
import random

import numpy as np

from rev_cuckoo_regular_isd import fit_pair_slope
from rev_cuckoo_regular_noise import LeakageChannel
from rev_cuckoo_ring_lpn_projected_bucket import (
    add_projection,
    cyclic_mass_convolution,
    one_trial,
    subtract_projection,
)


def test_projection_arithmetic_round_trip() -> None:
    left = (0, 16, 8)
    right = (3, 2, 10)
    total = add_projection(left, right, 17)
    assert subtract_projection(total, right, 17) == left


def test_fft_cyclic_convolution_matches_brute_force() -> None:
    rng = random.Random(2001)
    modulus = 3
    shape = (modulus, modulus)
    left = np.asarray(
        [rng.random() for _ in range(modulus**2)]
    ).reshape(shape)
    right = np.asarray(
        [rng.random() for _ in range(modulus**2)]
    ).reshape(shape)
    actual = cyclic_mass_convolution(left, right)
    expected = np.zeros(shape)
    for left_key in itertools.product(range(modulus), repeat=2):
        for right_key in itertools.product(range(modulus), repeat=2):
            target = add_projection(left_key, right_key, modulus)
            expected[target] += left[left_key] * right[right_key]
    assert np.allclose(actual, expected, atol=1e-12)


def test_projected_and_flat_decoders_recover_planted_solution() -> None:
    polys = 4
    weight = 2
    block_length = 3
    channel = LeakageChannel.build("occupancy", 4, weight, weight)
    pair_slope, _ = fit_pair_slope(
        channel, weight, 500, random.Random(2002)
    )
    candidates = tuple(
        itertools.product(range(block_length), repeat=weight)
    )
    rows = one_trial(
        trial=0,
        candidates=candidates,
        polys=polys,
        weight=weight,
        block_length=block_length,
        modulus=17,
        projection_width=1,
        channel=channel,
        pair_slope=pair_slope,
        seed=20260811,
    )
    assert {row.selector for row in rows} == {
        "uniform_projected",
        "exact_projected",
        "uniform_flat",
        "exact_flat",
    }
    assert all(row.recovered_planted for row in rows)
    assert all(row.charged_work > 0 for row in rows)


if __name__ == "__main__":
    test_projection_arithmetic_round_trip()
    test_fft_cyclic_convolution_matches_brute_force()
    test_projected_and_flat_decoders_recover_planted_solution()
