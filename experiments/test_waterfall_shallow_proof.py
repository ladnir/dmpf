#!/usr/bin/env python3
"""Exhaustive small-parameter checks for the one-repair proof formula."""

from __future__ import annotations

import itertools
import unittest
from fractions import Fraction

from waterfall_evictions import Configuration, bounded_bfs_repair, waterfall
from waterfall_shallow_proof import exact_one_repair_bound


def exhaustive_rejection_probability(
    t: int, w: int, d: int, stash: int
) -> Fraction:
    configuration = Configuration(t, w, d, stash)
    rejected = 0
    total = d ** (t * w)
    for values in itertools.product(range(d), repeat=t * w):
        candidates = list(values)
        owner, placement, overflow = waterfall(candidates, configuration)
        live = len(overflow)
        if live <= stash:
            failure = False
        elif live >= stash + 2:
            failure = True
        else:
            failure = bounded_bfs_repair(
                candidates,
                configuration,
                owner,
                placement,
                overflow,
                2,
            ) > stash
        rejected += failure
    return Fraction(rejected, total)


class ShallowProofTests(unittest.TestCase):
    def test_exact_formula_matches_exhaustive_enumeration(self) -> None:
        configurations = (
            (3, 2, 2, 0),
            (3, 2, 2, 1),
            (4, 3, 2, 1),
            (4, 2, 3, 1),
            (5, 3, 2, 1),
        )
        for t, w, d, stash in configurations:
            with self.subTest(t=t, w=w, d=d, stash=stash):
                expected, _, _ = exact_one_repair_bound(t, w, d, stash)
                actual = exhaustive_rejection_probability(t, w, d, stash)
                self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
