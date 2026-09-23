import itertools
import unittest
from fractions import Fraction

from rev_cuckoo_hall import (
    hall_witness_bound,
    hall_witness_bound_sizes,
    hall_occupancy_bound_sizes,
    identical_candidate_obstruction_probability,
    logadd2,
    occupancy_threshold_probability,
)
from waterfall_mpc_cost import Model
from waterfall_space_sweep import (
    Point,
    log2_fraction,
    overflow_distribution,
    overflow_distribution_sizes,
    pareto,
    sweep_asymmetric_reachability,
    tail,
)


class WaterfallSpaceSweepTests(unittest.TestCase):
    def test_overflow_distribution(self) -> None:
        distribution = overflow_distribution(4, 2, 4)
        self.assertEqual(sum(distribution), Fraction(1))
        self.assertEqual(tail(distribution, 0), Fraction(1))
        self.assertEqual(tail(distribution, 5), Fraction(0))

    def test_pareto_empty(self) -> None:
        self.assertEqual(pareto([]), [])

    def test_pareto_uses_ot_count_payload_and_expansion(self) -> None:
        best = Point("basic", 4, 16, 1, 0, -50.0, 42.0, 100, 1_000)
        dominated = Point("basic", 4, 16, 2, 0, -50.0, 42.0, 101, 1_001)
        payload_tradeoff = Point(
            "one repair", 4, 16, 0, 1, -50.0, 42.0, 110, 900
        )
        self.assertEqual(pareto([dominated, payload_tradeoff, best]), [payload_tradeoff, best])

    def test_main_bounded_augmentation_points(self) -> None:
        for sizes, repairs, expected_bits in (
            ((16, 16, 16, 16), 3, 43.36724520107565),
            ((128, 128, 64), 2, 41.13073055719418),
        ):
            hall_log, _ = hall_occupancy_bound_sizes(16, sizes, 0)
            overflow_log = log2_fraction(
                tail(overflow_distribution_sizes(16, sizes), repairs + 1)
            )
            batch_bits = -logadd2(hall_log, overflow_log) - 8
            self.assertAlmostEqual(batch_bits, expected_bits, places=12)

    def test_uniform_helpers_are_asymmetric_special_cases(self) -> None:
        self.assertEqual(
            overflow_distribution(16, 3, 128),
            overflow_distribution_sizes(16, (128, 128, 128)),
        )
        uniform, _ = hall_witness_bound(16, 3, 128, 0)
        heterogeneous, _ = hall_witness_bound_sizes(
            16, (128, 128, 128), 0
        )
        self.assertEqual(uniform, heterogeneous)

    def test_asymmetric_three_n_schedule(self) -> None:
        sizes = (128, 128, 64)
        hall_log, _ = hall_occupancy_bound_sizes(16, sizes, 0)
        overflow_log = log2_fraction(
            tail(overflow_distribution_sizes(16, sizes), 3)
        )
        combined_bits = -logadd2(hall_log, overflow_log) - 8
        self.assertAlmostEqual(combined_bits, 41.13073055719418, places=12)

        points = sweep_asymmetric_reachability(
            Model(), 3, 256, max_columns=384
        )
        best = min(points, key=lambda point: (point.ot_count, point.payload_bits))
        self.assertEqual(best.partition_sizes, sizes)
        self.assertEqual(best.roots, 2)
        self.assertEqual(best.ot_count, 55_358)
        self.assertEqual(best.payload_bits, 3_623_410)

    def test_compact_four_n_schedule(self) -> None:
        points = sweep_asymmetric_reachability(
            Model(), 4, 256, max_columns=256
        )
        best = min(
            points,
            key=lambda point: (
                sum(point.partition_sizes) + point.c,
                point.ot_count,
                point.payload_bits,
            ),
        )
        self.assertEqual(best.partition_sizes, (16, 16, 16, 16))
        self.assertEqual(best.roots, 3)
        self.assertAlmostEqual(best.batch_bits, 43.36724520107565, places=12)
        self.assertEqual(best.ot_count, 39_052)
        self.assertEqual(best.payload_bits, 942_888)

    def test_ot_minimized_four_n_schedule(self) -> None:
        points = sweep_asymmetric_reachability(
            Model(), 4, 256, max_columns=256
        )
        best = min(points, key=lambda point: (point.ot_count, point.payload_bits))
        self.assertEqual(best.partition_sizes, (32, 32, 16, 8))
        self.assertEqual(best.roots, 2)
        self.assertAlmostEqual(best.batch_bits, 41.2014987948858, places=12)
        self.assertEqual(best.ot_count, 38_158)
        self.assertEqual(best.payload_bits, 1_162_710)

    def test_exact_occupancy_threshold_probability(self) -> None:
        brute_force_hits = 0
        assignments = 0
        for values in itertools.product(range(3), repeat=4):
            assignments += 1
            if max(values.count(bin_index) for bin_index in range(3)) >= 3:
                brute_force_hits += 1
        self.assertEqual(
            occupancy_threshold_probability(4, 3, 3),
            Fraction(brute_force_hits, assignments),
        )

    def test_identical_tuple_obstruction_forces_selected_products(self) -> None:
        for sizes, expected_bits in (
            ((128, 128, 32), 38.17030368146004),
            ((128, 128, 64), 41.17029047318425),
            ((16, 16, 8, 16), 39.90764644869978),
            ((16, 16, 16, 16), 43.90744465297357),
        ):
            probability = identical_candidate_obstruction_probability(16, sizes)
            batch_bits = -log2_fraction(probability) - 8
            self.assertAlmostEqual(batch_bits, expected_bits, places=12)


if __name__ == "__main__":
    unittest.main()
