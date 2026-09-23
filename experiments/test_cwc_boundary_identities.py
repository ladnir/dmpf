"""Algebraic regression checks for the CWC corrections, not security tests."""

from collections import Counter
from itertools import product
from random import Random

import pytest

from test_sparse_dpf_schedule import check_invariant, correction, expand


def standard_leaves(n, branching, alpha, salt):
    """Model the standard box with local seeds, tags, and public corrections."""
    rng = Random(salt)
    depth, capacity = 0, 1
    while capacity < n:
        capacity *= branching
        depth += 1
    assert depth > 0
    nodes = [(rng.randrange(256), rng.randrange(256), 0, 1)
             for _ in range(branching)]
    for level in range(depth - 1, -1, -1):
        sums = [0] * branching
        for i, (left, right, _, _) in enumerate(nodes):
            sums[i % branching] ^= left ^ right
        digit = (alpha // branching**level) % branching
        if branching == 2:
            words = correction(sums, digit)
        else:
            active = rng.randrange(128) * 2 + 1
            words = [s ^ (active if j == digit else 0)
                     for j, s in enumerate(sums)]
        corrected = [(a ^ (words[i % branching] if ta else 0),
                      b ^ (words[i % branching] if tb else 0))
                     for i, (a, b, ta, tb) in enumerate(nodes)]
        if not level:
            return dict(enumerate(corrected[:n]))
        nodes = []
        for left, right in corrected:
            for j in range(branching):
                def prg(seed):
                    return ((seed * (13 + 2 * j)) ^ (seed >> 1) ^ (73 * j)) & 255
                nodes.append((prg(left), prg(right), left & 1, right & 1))


@pytest.mark.parametrize("branching", (2, 3, 4, 5))
def test_standard_non_power_domains(branching):
    for n in range(2, 34):
        for alpha in range(n):
            leaves = standard_leaves(n, branching, alpha, n + alpha)
            check_invariant(leaves, alpha)
            for round_number, beta in enumerate((0, 1, 16)):
                assert expand(leaves, beta, round_number, 17) == {
                    i: beta if i == alpha else 0 for i in range(n)}


def test_singleton_identity_resharing():
    for modulus in (2, 3, 17):
        for beta, left in product(range(modulus), repeat=2):
            right = (beta - left) % modulus
            assert (left + right) % modulus == beta


def test_group_by_counts_and_joint_opening():
    for n, t, capacity in ((2, 4, 4), (3, 4, 2), (2, 8, 8), (3, 5, 2)):
        count_modulus = 1 << t.bit_length()
        for addresses in product(range(n), repeat=t):
            counts = Counter(addresses)
            opened = Counter((a, 1) for a in addresses)
            for key in range(n):
                count = counts[key] % count_modulus
                assert count == counts[key]
                selected = capacity - min(count, capacity)
                opened[key, 1] += selected
                opened[0, 0] += capacity - selected
                assert opened[key, 1] >= capacity
            if max(counts.values()) <= capacity:
                expected = Counter({(key, 1): capacity for key in range(n)})
                expected[0, 0] = t
                assert opened == expected


def negacyclic_product(left, right, modulus):
    n = len(left)
    out = [0] * n
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[(i + j) % n] += a * b * (-1 if i + j >= n else 1)
    return [x % modulus for x in out]


def test_regular_product_list_wraparound():
    rng = Random(830)
    for blocks, length in ((2, 2), (3, 3), (4, 3)):
        n, modulus = blocks * length, 17
        for _ in range(100):
            offsets = [[rng.randrange(length) for _ in range(blocks)] for _ in range(2)]
            coeffs = [[rng.randrange(1, modulus) for _ in range(blocks)] for _ in range(2)]
            polys = [[0] * n for _ in range(2)]
            for party in range(2):
                for j in range(blocks):
                    polys[party][j * length + offsets[party][j]] = coeffs[party][j]
            reconstructed = [0] * n
            for k, j in product(range(blocks), repeat=2):
                other = (k - j) % blocks
                pos = k * length + offsets[0][j] + offsets[1][other]
                sign = (-1)**((j + other) // blocks + pos // n)
                reconstructed[pos % n] += sign * coeffs[0][j] * coeffs[1][other]
            assert [x % modulus for x in reconstructed] == negacyclic_product(*polys, modulus)


def test_ring_mask_weights():
    rng = Random(93)
    n, modulus = 8, 17
    for _ in range(100):
        e00, e01, e10, e11, mask = [[rng.randrange(modulus) for _ in range(n)] for _ in range(5)]
        mul = lambda a, b: negacyclic_product(a, b, modulus)
        add = lambda *xs: [sum(v) % modulus for v in zip(*xs)]
        actual = mul(add(e00, mul(mask, e01)), add(e10, mul(mask, e11)))
        weighted = add(mul(e00, e10), mul(mask, add(mul(e00, e11), mul(e01, e10))),
                       mul(mul(mask, mask), mul(e01, e11)))
        assert actual == weighted


def test_four_polynomial_ring_mask_weights():
    rng = Random(940)
    n, modulus, count = 8, 17, 4
    mul = lambda a, b: negacyclic_product(a, b, modulus)
    add = lambda xs: [sum(v) % modulus for v in zip(*xs)]
    for _ in range(30):
        noises = [[[rng.randrange(modulus) for _ in range(n)]
                   for _ in range(count)] for _ in range(2)]
        masks = [[1] + [0] * (n - 1)] + [
            [rng.randrange(modulus) for _ in range(n)] for _ in range(count - 1)]
        inputs = [add([mul(masks[i], noises[p][i]) for i in range(count)])
                  for p in range(2)]
        weighted = add([mul(mul(masks[i], masks[j]), mul(noises[0][i], noises[1][j]))
                        for i, j in product(range(count), repeat=2)])
        assert mul(*inputs) == weighted
