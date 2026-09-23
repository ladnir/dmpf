"""Executable small-tree model of the paper's sparse-DPF scheduling.

Both parties' seed shares are evaluated locally. Shared circuit operations
are reconstructed and freshly reshared. The toy PRG is deterministic, not
cryptographic: the tests establish algebraic/scheduling identities only.
"""

from dataclasses import dataclass
from itertools import combinations, product
from random import Random

import pytest


def partition(support, interval):
    """Use the paper's one-based inclusive interval convention."""
    first, last = interval
    assert 1 <= first <= last <= len(support)
    if first == last:
        return 0, (interval, None)
    depth = (support[first - 1] ^ support[last - 1]).bit_length()
    count = sum(not ((x >> (depth - 1)) & 1)
                for x in support[first - 1:last])
    assert 0 < count < last - first + 1
    return depth, ((first, first + count - 1), (first + count, last))


def xor_pair(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def correction(sums, bit):
    """One-based child 2-bit is off-path; upper correction bits coincide."""
    upper = sums[1 - bit] & ~1
    return upper | ((sums[0] & 1) ^ bit ^ 1), upper | ((sums[1] & 1) ^ bit)


def toy_prg(seed, width):
    mask = (1 << width) - 1
    return ((seed * 13) ^ (seed >> 1) ^ 0xA5) & mask, (
        (seed * 29) ^ (seed >> 2) ^ 0x53) & mask


@dataclass
class Node:
    child: int
    parent_depth: int
    intervals: tuple
    seeds: tuple
    tags: tuple


def setup(support, alpha, *, depth=None, width=8, salt=0,
          mask_inactive=True):
    support = tuple(support)
    assert support == tuple(sorted(set(support))) and len(support) > 1
    assert alpha in support
    depth = depth if depth is not None else max(1, support[-1].bit_length())
    rng = Random(salt)
    mask = (1 << width) - 1

    def share(value, bits=width):
        left = rng.getrandbits(bits)
        return left, left ^ value

    queues = [[] for _ in range(depth + 1)]
    # sums[d][child] is a pair of local XOR shares.
    sums = [[(0, 0), (0, 0)] for _ in range(depth + 1)]
    active = [0] * (depth + 1)
    words, events, leaves = {}, [], {}
    root, intervals = partition(support, (1, len(support)))
    active[root] = 1
    for child in (0, 1):
        target, children = partition(support, intervals[child])
        seeds = rng.getrandbits(width), rng.getrandbits(width)
        sums[root][child] = seeds
        queues[target].append(Node(child, root, children, seeds, share(1, 1)))

    def apply(node, target):
        assert node.parent_depth in words
        assert node.parent_depth > target
        events.append(("consume", node.parent_depth, target))
        word = words[node.parent_depth][node.child]
        return tuple(s ^ (word if t else 0) for s, t in zip(node.seeds, node.tags))

    calls = 0
    for d in range(depth, 0, -1):
        if d != root and not queues[d]:
            continue
        for node in queues[d]:
            corrected = apply(node, d)
            expanded = tuple(toy_prg(s, width) for s in corrected)
            calls += 1
            tags = tuple(s & 1 for s in corrected)
            active[d] ^= tags[0] ^ tags[1]
            for child in (0, 1):
                seeds = expanded[0][child], expanded[1][child]
                sums[d][child] = xor_pair(sums[d][child], seeds)
                target, children = partition(support, node.intervals[child])
                assert target < d
                queues[target].append(Node(child, d, children, seeds, tags))
        reconstructed = tuple(a ^ b for a, b in sums[d])
        if mask_inactive and not active[d]:
            reconstructed = tuple(s ^ rng.getrandbits(width) for s in reconstructed)
        words[d] = correction(reconstructed, (alpha >> (d - 1)) & 1)
        events.append(("produce", d))

    for node in queues[0]:
        address = support[node.intervals[0][0] - 1]
        assert address not in leaves
        leaves[address] = apply(node, 0)
    assert set(leaves) == set(support)
    assert all(0 <= s <= mask for pair in leaves.values() for s in pair)
    return leaves, words, events, calls


def check_invariant(leaves, alpha):
    for x, (left, right) in leaves.items():
        assert (left ^ right) & 1 == (x == alpha)
        if x != alpha:
            assert left == right


def expand(leaves, beta, round_number, modulus):
    """Prime and composite cyclic groups; signs matter outside char. two."""
    tags = [[pair[p] & 1 for pair in leaves.values()] for p in (0, 1)]
    values = [[((1 - 2 * p) * (
        ((pair[p] >> 1) * 11 + round_number * 7) ^ (pair[p] >> 3)
    )) % modulus for pair in leaves.values()] for p in (0, 1)]
    counts = [sum(row) for row in tags]
    d = ((counts[0] // 2) % 2) ^ ((counts[1] // 2 + counts[1]) % 2)
    delta = (beta - sum(map(sum, values))) % modulus
    opened = (1 - 2 * d) * delta % modulus
    outputs = [[(value + (1 - 2 * p) * opened * tag) % modulus
                for value, tag in zip(values[p], tags[p])] for p in (0, 1)]
    return {x: (outputs[0][i] + outputs[1][i]) % modulus
            for i, x in enumerate(leaves)}


def test_partition_boundaries_and_singletons():
    support = (1, 2, 4, 5, 8, 9, 10, 13, 15)
    assert partition(support, (2, 7)) == (4, ((2, 4), (5, 7)))
    for first in range(1, len(support) + 1):
        for last in range(first, len(support) + 1):
            d, ranges = partition(support, (first, last))
            if first == last:
                assert (d, ranges) == (0, ((first, last), None))
                continue
            left, right = ranges
            assert left[0] == first and right[1] == last
            assert left[1] + 1 == right[0]
            assert all(not ((x >> (d - 1)) & 1)
                       for x in support[left[0] - 1:left[1]])
            assert all((x >> (d - 1)) & 1
                       for x in support[right[0] - 1:right[1]])


def test_two_leaf_root_correction_with_empty_expansion_queues():
    for alpha in (0, 1):
        leaves, words, events, calls = setup((0, 1), alpha)
        check_invariant(leaves, alpha)
        assert set(words) == {1} and calls == 0
        assert events == [("produce", 1), ("consume", 1, 0), ("consume", 1, 0)]


def test_deferred_parent_correction_and_leaf_address():
    leaves, words, events, calls = setup((0, 1, 4), 1, depth=6)
    check_invariant(leaves, 1)
    assert set(words) == {1, 3} and calls == 1
    assert events[:3] == [("produce", 3), ("consume", 3, 1), ("produce", 1)]
    assert ("consume", 3, 0) in events


def test_publicly_present_inactive_split_requires_secret_activity_mask():
    """An unmasked off-path split leaks despite correct reconstruction.

    This is an algebraic witness, not an end-to-end C++ execution or a
    security test for the toy PRG. The C++ public-occupancy check omits
    the mask in exactly this situation. Its opened raw sigma is zero;
    the paper's packed correction pair consequently has the form (1, 0).
    """
    for salt in range(64):
        leaves, words, _, calls = setup(
            (0, 1, 4), 4, salt=salt, mask_inactive=False)
        check_invariant(leaves, 4)
        assert calls == 1  # The lower split exists publicly.
        assert words[1] == (1, 0)  # But only an off-path parent splits.

    masked_words = {
        setup((0, 1, 4), 4, salt=salt)[1][1] for salt in range(64)
    }
    assert len(masked_words) > 1


def test_exhaustive_small_supports_and_all_active_indices():
    for size in range(2, 9):
        for support in combinations(range(8), size):
            for alpha in support:
                for salt in (0, 1, 17):
                    leaves, words, _, calls = setup(support, alpha, salt=salt)
                    check_invariant(leaves, alpha)
                    assert calls == size - 2
                    expected_depths = set()

                    def splits(interval):
                        d, children = partition(support, interval)
                        if d:
                            expected_depths.add(d)
                            for child in children:
                                splits(child)

                    splits((1, size))
                    assert set(words) == expected_depths


@pytest.mark.parametrize("modulus", (2, 3, 4, 5, 8, 17))
def test_repeated_payload_expansion(modulus):
    for support in ((0, 1), (0, 1, 4), (1, 2, 7, 10, 14)):
        for alpha in support:
            leaves, _, _, _ = setup(support, alpha, salt=41)
            saved = dict(leaves)
            for q, beta in enumerate(range(modulus)):
                result = expand(leaves, beta, q, modulus)
                assert result == {x: beta if x == alpha else 0 for x in support}
                assert leaves == saved


def test_sign_conversion_for_all_small_tag_sharings():
    for size in range(1, 7):
        for alpha in range(size):
            for left in product((0, 1), repeat=size):
                right = tuple(t ^ (i == alpha) for i, t in enumerate(left))
                c0, c1 = sum(left), sum(right)
                d = (c0 // 2 + c1 // 2 + c1) % 2
                assert d == right[alpha]
                assert c0 - c1 == 1 - 2 * right[alpha]
