"""Exact transcript checks for the inspected C++ sparse correction encoding.

Before the approved fix, SparseDpf.h opened the full selected block and two
separate tag bits. These tests retain that disclosure's equations, not AES
or a network execution. The fixed encoder has a separate C++ socket regression.
Uniform reconstructed child sums occur at the directly sampled root; the same
experiment also models independent-child hybrids and masked inactive levels.
The packed encoding is the paper's existing correction format.
"""

from collections import Counter
from fractions import Fraction
from itertools import product

import pytest


def raw_correction(z0, z1, digit):
    sigma = (z0, z1)[1 - digit]
    tau0 = (z0 & 1) ^ digit ^ 1
    tau1 = (z1 & 1) ^ digit
    return sigma, tau0, tau1


def infer_digit(opening):
    sigma, tau0, tau1 = opening
    if tau0 == tau1:
        return None
    return int((sigma & 1) == tau0)


def packed_correction(z0, z1, digit):
    sigma, tau0, tau1 = raw_correction(z0, z1, digit)
    return (sigma & ~1) | tau0, (sigma & ~1) | tau1


@pytest.mark.parametrize("bits", (1, 2, 3, 4))
def test_full_sigma_gives_exact_digit_on_half_of_root_openings(bits):
    count = 1 << (2 * bits)
    laws = []
    for digit in (0, 1):
        law = Counter()
        recovered = 0
        for z0, z1 in product(range(1 << bits), repeat=2):
            opening = raw_correction(z0, z1, digit)
            law[opening] += 1
            inferred = infer_digit(opening)
            if inferred is not None:
                assert inferred == digit
                recovered += 1
        assert recovered == count // 2
        laws.append(law)

    distance = sum(abs(laws[0][x] - laws[1][x])
                   for x in laws[0].keys() | laws[1].keys())
    assert Fraction(distance, 2 * count) == Fraction(1, 2)


@pytest.mark.parametrize("bits", (1, 2, 3, 4))
def test_paper_packed_encoding_removes_the_extra_disclosure(bits):
    laws = [Counter(packed_correction(z0, z1, digit)
                    for z0, z1 in product(range(1 << bits), repeat=2))
            for digit in (0, 1)]
    assert laws[0] == laws[1]
    assert len(laws[0]) == 1 << (bits + 1)


def test_separate_seed_tag_representation_requires_repacking_generator_input():
    # Removing sigma's low bit alone is not a representation conversion.
    # An off-path child must become equal in both seed and logical tag.
    for digit, parent0, child0, child1 in product(
            range(2), range(2), range(4), range(4)):
        parent1 = parent0 ^ 1
        branch = 1 - digit
        # Choose the other reconstructed child sum arbitrarily.
        sums = [2, 2]
        sums[branch] = child0 ^ child1
        packed = packed_correction(*sums, digit)[branch]
        corrected0 = child0 ^ (parent0 * packed)
        corrected1 = child1 ^ (parent1 * packed)
        assert corrected0 == corrected1

        prefix_correction = packed & ~1
        seed0 = (child0 ^ (parent0 * prefix_correction)) & ~1
        seed1 = (child1 ^ (parent1 * prefix_correction)) & ~1
        tag0 = (child0 & 1) ^ (parent0 * (packed & 1))
        tag1 = (child1 & 1) ^ (parent1 * (packed & 1))
        assert seed0 == seed1
        assert tag0 == tag1
        assert seed0 | tag0 == corrected0
        assert seed1 | tag1 == corrected1
