"""Checks of the reusable-DPF output-sharing boundary, not a privacy attack.

The current Expand computes a local vector u from saved seeds, then returns
u + (-1)^p gamma t for one opened group element gamma and a local bit vector t.
These tests check the resulting same-tag invariant. They do not simulate the
DPF, prove privacy, or inspect hidden state of an ideal subprotocol.
"""

from itertools import combinations, product

import pytest

from test_sparse_dpf_schedule import setup


def equal_tag_pair(tags):
    return next((i, j) for i, j in combinations(range(len(tags)), 2)
                if tags[i] == tags[j])


def corrected_output(provisional, tags, gamma, party, modulus):
    return tuple((value + (1 - 2 * party) * gamma * tag) % modulus
                 for value, tag in zip(provisional, tags, strict=True))


@pytest.mark.parametrize("modulus", (2, 3, 4, 5, 8, 17))
def test_any_scalar_correction_preserves_same_tag_difference(modulus):
    for party in (0, 1):
        for tags in product((0, 1), repeat=3):
            i, j = equal_tag_pair(tags)
            for provisional in ((0, 0, 0), (2, 4, 7), (6, 1, 3)):
                predicted = (provisional[i] - provisional[j]) % modulus
                for gamma in range(modulus):
                    output = corrected_output(provisional, tags, gamma, party, modulus)
                    assert (output[i] - output[j]) % modulus == predicted


@pytest.mark.parametrize("modulus", (2, 3, 4, 5, 8, 17))
def test_fresh_uniform_map_matches_fixed_prediction_with_probability_one_over_order(modulus):
    # The pair and predicted difference are fixed before the fresh ideal output.
    # Allow the simulator any choice of the subsequent scalar gamma: it cancels.
    for tags in product((0, 1), repeat=3):
        i, j = equal_tag_pair(tags)
        for predicted in range(modulus):
            accepted = sum((output[i] - output[j]) % modulus == predicted
                           for output in product(range(modulus), repeat=3))
            assert accepted == modulus ** 2


@pytest.mark.parametrize("modulus", (2, 3, 4, 5, 8, 17))
def test_invariant_on_actual_small_sparse_setup_states(modulus):
    support = (0, 1, 4)
    for alpha in support:
        for salt in (0, 1, 17):
            leaves, _, _, _ = setup(support, alpha, salt=salt)
            for party in (0, 1):
                seeds = tuple(leaves[x][party] for x in support)
                tags = tuple(seed & 1 for seed in seeds)
                i, j = equal_tag_pair(tags)
                for round_number in range(3):
                    # Match the existing scheduling model's toy leaf conversion.
                    provisional = tuple((1 - 2 * party) * (
                        ((seed >> 1) * 11 + round_number * 7) ^ (seed >> 3)
                    ) % modulus for seed in seeds)
                    predicted = (provisional[i] - provisional[j]) % modulus
                    for gamma in range(modulus):
                        output = corrected_output(provisional, tags, gamma, party, modulus)
                        assert (output[i] - output[j]) % modulus == predicted
