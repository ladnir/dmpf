#!/usr/bin/env python3
"""Exhaustive tests for complete Waterfall reachability repair."""

from __future__ import annotations

import itertools
import unittest

from waterfall_evictions import Configuration, maximum_matching_unplaced, waterfall
from waterfall_reachability import (
    circuit_reachability_repair,
    complete_reachability_repair,
    implicit_owner_dfs_repair,
)


class ReachabilityRepairTests(unittest.TestCase):
    def test_complete_repair_reaches_maximum_matching(self) -> None:
        configurations = (
            Configuration(3, 2, 2, 0),
            Configuration(4, 2, 2, 0),
            Configuration(4, 3, 2, 0),
            Configuration(4, 2, 3, 0),
        )
        for configuration in configurations:
            total_cells = configuration.t * configuration.w
            for values in itertools.product(
                range(configuration.d), repeat=total_cells
            ):
                candidates = list(values)
                owner, placement, _ = waterfall(candidates, configuration)
                actual = complete_reachability_repair(
                    candidates,
                    configuration,
                    owner,
                    placement,
                    configuration.t,
                )
                expected = maximum_matching_unplaced(
                    candidates, configuration
                )
                self.assertEqual(actual, expected)
                circuit_actual = circuit_reachability_repair(
                    candidates,
                    configuration,
                    placement,
                    configuration.t,
                )
                self.assertEqual(circuit_actual, expected)
                dfs_actual = implicit_owner_dfs_repair(
                    candidates,
                    configuration,
                    owner,
                    placement,
                    configuration.t,
                )
                self.assertEqual(dfs_actual, expected)

    def test_literal_circuit_matches_reference_at_every_bound(self) -> None:
        configurations = (
            Configuration(3, 2, 2, 0),
            Configuration(4, 2, 2, 0),
        )
        for configuration in configurations:
            total_cells = configuration.t * configuration.w
            for values in itertools.product(
                range(configuration.d), repeat=total_cells
            ):
                candidates = list(values)
                owner, placement, _ = waterfall(candidates, configuration)
                for bound in range(configuration.t + 1):
                    expected = complete_reachability_repair(
                        candidates,
                        configuration,
                        owner,
                        placement,
                        bound,
                    )
                    actual = circuit_reachability_repair(
                        candidates,
                        configuration,
                        placement,
                        bound,
                    )
                    self.assertEqual(actual, expected)
                    dfs_actual = implicit_owner_dfs_repair(
                        candidates,
                        configuration,
                        owner,
                        placement,
                        bound,
                    )
                    self.assertEqual(dfs_actual, expected)

    def test_augmentation_bound_increases_matching_by_at_most_one_each(self) -> None:
        configuration = Configuration(4, 2, 2, 0)
        for values in itertools.product(range(2), repeat=8):
            candidates = list(values)
            owner, placement, overflow = waterfall(candidates, configuration)
            for bound in range(3):
                remaining = complete_reachability_repair(
                    candidates,
                    configuration,
                    owner,
                    placement,
                    bound,
                )
                self.assertGreaterEqual(remaining, len(overflow) - bound)


if __name__ == "__main__":
    unittest.main()
