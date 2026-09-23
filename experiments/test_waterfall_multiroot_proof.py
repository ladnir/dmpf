#!/usr/bin/env python3
"""Exhaustive checks for the all-overflow-roots obstruction formula."""

from __future__ import annotations

import itertools
import unittest
from fractions import Fraction

from waterfall_evictions import Configuration, bounded_bfs_repair, waterfall
from waterfall_multiroot_proof import exact_all_roots_depth2_obstruction


def exhaustive_probability(t: int, w: int, d: int, overflow: int) -> Fraction:
    configuration = Configuration(t, w, d, overflow - 1)
    failures = 0
    total = d ** (t * w)
    for candidates_tuple in itertools.product(range(d), repeat=t * w):
        candidates = list(candidates_tuple)
        owner, placement, unplaced = waterfall(candidates, configuration)
        if len(unplaced) != overflow:
            continue
        failures += bounded_bfs_repair(
            candidates,
            configuration,
            owner,
            placement,
            unplaced,
            2,
        ) > configuration.c
    return Fraction(failures, total)


class MultiRootProofTests(unittest.TestCase):
    def test_exact_formula_matches_exhaustive_enumeration(self) -> None:
        configurations = (
            (3, 2, 2, 1),
            (4, 2, 2, 2),
            (4, 3, 2, 1),
            (5, 2, 2, 2),
        )
        for t, w, d, overflow in configurations:
            with self.subTest(t=t, w=w, d=d, overflow=overflow):
                expected = exact_all_roots_depth2_obstruction(
                    t, w, d, overflow
                )
                actual = exhaustive_probability(t, w, d, overflow)
                self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
