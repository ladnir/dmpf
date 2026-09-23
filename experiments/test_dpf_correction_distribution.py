"""Exact finite checks for ReusableDpfAppendix's random-value hybrids.

Fresh child blocks are independent by experiment definition. These tests do
not justify replacing a real PRG or AES evaluation by those blocks.
"""

from collections import Counter
from itertools import product

import pytest

from test_reusable_dpf_sharing_contract import real_expansion, simulated_expansion


def binary_correction(children, digit):
    prefix = children[1 - digit] >> 1
    return ((prefix << 1) | ((children[0] & 1) ^ digit ^ 1),
            (prefix << 1) | ((children[1] & 1) ^ digit))


@pytest.mark.parametrize("bits", (2, 3, 4))
def test_binary_correction_law_and_hidden_active_prefix(bits):
    blocks, prefixes = 1 << bits, 1 << (bits - 1)
    expected_corrections = Counter({
        ((a << 1) | t0, (a << 1) | t1): prefixes
        for a, t0, t1 in product(range(prefixes), range(2), range(2))})
    for digit, corrupt_tag in product(range(2), repeat=2):
        corrupt_children = (3 % blocks, 6 % blocks)
        law, joint = Counter(), Counter()
        for sums in product(range(blocks), repeat=2):
            correction = binary_correction(sums, digit)
            law[correction] += 1
            honest_raw = sums[digit] ^ corrupt_children[digit]
            honest_corrected = honest_raw ^ ((1 - corrupt_tag) * correction[digit])
            joint[correction, honest_corrected >> 1] += 1
        assert law == expected_corrections
        assert len(joint) == len(expected_corrections) * prefixes
        assert set(joint.values()) == {1}


@pytest.mark.parametrize("width,bits", ((3, 2), (3, 3), (4, 2)))
def test_wide_correction_is_uniform_independent_of_active_mask(width, bits):
    blocks, prefixes = 1 << bits, 1 << (bits - 1)
    expected = Counter({values: 1 for values in product(range(blocks), repeat=width)})
    for digit, mask in product(range(width), range(prefixes)):
        law = Counter()
        for sums in product(range(blocks), repeat=width):
            correction = tuple(value ^ (((mask << 1) | 1) if j == digit else 0)
                               for j, value in enumerate(sums))
            law[correction] += 1
        # Identical for every fixed mask: conditioning on the public tuple
        # leaves the originally uniform active mask uniform.
        assert law == expected


@pytest.mark.parametrize("bits", (2, 3, 4))
def test_inactive_sparse_level_has_the_same_correction_law(bits):
    # The inactive reconstructed sums are zero; XOR a fresh two-block mask.
    law = [Counter(binary_correction(mask, digit)
                   for mask in product(range(1 << bits), repeat=2))
           for digit in range(2)]
    assert law[0] == law[1]
    assert len(law[0]) == 1 << (bits + 1)
    assert set(law[0].values()) == {1 << (bits - 1)}


@pytest.mark.parametrize("modulus", (2, 3, 4, 5, 8))
def test_opening_messages_and_sign_multiplication_shares(modulus):
    local, tags = (2, 4, 7), (0, 1, 0)
    corrupt_beta = 2 % modulus
    for party, alpha, beta in product(range(2), range(3), range(modulus)):
        real, ideal = Counter(), Counter()
        # In characteristic two there is no sign-multiplication output.
        local_openings = ((corrupt_beta - sum(local)) % modulus,) if modulus == 2 else range(modulus)
        for random_value, local_opening in product(range(modulus), local_openings):
            gamma, corrupt, honest = real_expansion(
                local, tags, alpha, beta, random_value, party, modulus)
            honest_message = (gamma - local_opening) % modulus
            real[local_opening, honest_message, corrupt, honest] += 1

            # The simulator samples gamma, never the hidden active value.
            sim_gamma, sim_corrupt = simulated_expansion(
                local, tags, random_value, party, modulus)
            sim_message = (sim_gamma - local_opening) % modulus
            completed = tuple(((beta if x == alpha else 0) - value) % modulus
                              for x, value in enumerate(sim_corrupt))
            ideal[local_opening, sim_message, sim_corrupt, completed] += 1
        assert real == ideal
