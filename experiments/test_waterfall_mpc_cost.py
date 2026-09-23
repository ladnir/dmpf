import unittest

from waterfall_mpc_cost import (
    Configuration,
    Model,
    configuration_cost,
    direct_hash_cost,
    dpf_initialization_cost,
    field_multiplication_cost,
    karatsuba_multiplication_cost,
    waksman_switch_count,
)


class WaterfallMpcCostTests(unittest.TestCase):
    def test_reference_one_repair_point(self) -> None:
        model = Model()
        cfg = Configuration("one repair", "one-repair", 4, 16, 2)
        cost = configuration_cost(model, cfg)
        self.assertEqual(cost.generator.ot_count, 1_494)
        self.assertEqual(cost.generator.payload_bits, 10_190)
        self.assertEqual(karatsuba_multiplication_cost(20), 153)
        self.assertEqual(field_multiplication_cost(20), 81)
        self.assertEqual(direct_hash_cost(model, cfg), 9_328)
        self.assertEqual(cost.hash_eval.ot_count, 18_272)
        self.assertEqual(cost.hash_eval.payload_bits, 18_656)
        self.assertEqual(cost.sparse_dpf.ot_count, 4_480)
        self.assertEqual(cost.sparse_dpf.payload_bits, 748_800)
        self.assertEqual(waksman_switch_count(cfg.columns), 358)
        self.assertEqual(cost.scatter.ot_count, 716)
        self.assertEqual(cost.scatter.payload_bits, 27_208)
        self.assertEqual(cost.common.ot_count, 6_018)
        self.assertEqual(cost.common.payload_bits, 24_434)
        self.assertEqual(cost.total_ot_count, 30_980)
        self.assertEqual(cost.total_payload_bits, 829_288)

    def test_complete_bfs_is_more_expensive_than_one_repair(self) -> None:
        model = Model()
        one = configuration_cost(
            model, Configuration("one repair", "one-repair", 4, 16, 2)
        )
        complete = configuration_cost(
            model,
            Configuration(
                "complete BFS", "bfs", 4, 16, 0,
                bfs_depth=16, bfs_roots=3,
            ),
        )
        self.assertEqual(complete.generator.ot_count, 323_680)
        self.assertEqual(complete.total_ot_count, 353_122)
        self.assertGreater(complete.total, one.total)

    def test_smaller_d_can_beat_smaller_expansion_coefficient(self) -> None:
        model = Model()
        compact_setup = configuration_cost(
            model, Configuration("one repair", "one-repair", 4, 16, 2)
        )
        smaller_w_plus_c = configuration_cost(
            model, Configuration("one repair", "one-repair", 5, 16, 0)
        )
        self.assertLess(compact_setup.total, smaller_w_plus_c.total)

    def test_compact_expansion_frontier_points(self) -> None:
        model = Model()
        high_bandwidth = configuration_cost(
            model, Configuration("one repair", "one-repair", 4, 64, 0)
        )
        middle = configuration_cost(
            model,
            Configuration(
                "complete reachability", "reachability", 4, 32, 0,
                bfs_roots=2,
            ),
        )
        low_bandwidth = configuration_cost(
            model,
            Configuration(
                "complete reachability", "reachability", 4, 16, 0,
                bfs_roots=3,
            ),
        )
        self.assertEqual(
            (high_bandwidth.total_ot_count, high_bandwidth.total_payload_bits),
            (46_336, 2_915_640),
        )
        self.assertEqual(
            (middle.total_ot_count, middle.total_payload_bits),
            (41_758, 1_643_414),
        )
        self.assertEqual(
            (low_bandwidth.total_ot_count, low_bandwidth.total_payload_bits),
            (39_052, 942_888),
        )

    def test_bounded_eviction_reference_count(self) -> None:
        model = Model()
        eviction = configuration_cost(
            model,
            Configuration(
                "bounded eviction", "eviction", 4, 16, 0,
                bfs_depth=24, bfs_roots=3, proved=False,
            ),
        )
        self.assertEqual(eviction.generator.ot_count, 39_904)
        self.assertEqual(eviction.total_ot_count, 69_346)

    def test_eviction_step_budget_is_global(self) -> None:
        model = Model()
        eight_steps = configuration_cost(
            model,
            Configuration(
                "bounded eviction", "eviction", 4, 16, 0,
                bfs_depth=8, bfs_roots=3, proved=False,
            ),
        )
        self.assertEqual(eight_steps.generator.ot_count, 19_168)
        self.assertEqual(eight_steps.total_ot_count, 48_610)

    def test_complete_reachability_estimates(self) -> None:
        model = Model()
        three_n = configuration_cost(
            model,
            Configuration(
                "3N bounded augmentation", "reachability", 3, 128, 0,
                bfs_roots=2,
            ),
        )
        four_n = configuration_cost(
            model,
            Configuration(
                "4N bounded augmentation", "reachability", 4, 64, 0,
                bfs_roots=1,
            ),
        )
        compact = configuration_cost(
            model,
            Configuration(
                "complete reachability", "reachability", 4, 8, 2,
                bfs_roots=3,
            ),
        )
        smaller_expansion = configuration_cost(
            model,
            Configuration(
                "complete reachability", "reachability", 5, 8, 0,
                bfs_roots=3,
            ),
        )
        self.assertEqual(three_n.generator.ot_count, 6_844)
        self.assertEqual(three_n.total_ot_count, 60_606)
        self.assertEqual(three_n.total_payload_bits, 4_335_634)
        self.assertEqual(four_n.generator.ot_count, 4_846)
        self.assertEqual(four_n.total_ot_count, 49_136)
        self.assertEqual(four_n.total_payload_bits, 2_984_068)
        self.assertEqual(compact.generator.ot_count, 8_970)
        self.assertEqual(compact.total_ot_count, 36_028)
        self.assertEqual(smaller_expansion.generator.ot_count, 10_378)
        self.assertEqual(smaller_expansion.total_ot_count, 37_932)

        target = configuration_cost(
            model,
            Configuration(
                "complete reachability", "reachability", 4, 16, 0,
                bfs_roots=3,
            ),
        )
        self.assertEqual(target.generator.ot_count, 9_610)
        self.assertEqual(target.generator.payload_bits, 117_334)
        self.assertEqual(target.total_ot_count, 39_052)
        self.assertEqual(target.total_payload_bits, 942_888)

    def test_asymmetric_three_n_point(self) -> None:
        model = Model()
        cfg = Configuration(
            "3N asymmetric bounded augmentation",
            "reachability",
            3,
            0,
            0,
            bfs_roots=2,
            partition_sizes=(128, 128, 64),
        )
        cost = configuration_cost(model, cfg)
        self.assertEqual(cfg.columns, 320)
        self.assertEqual(cfg.partition_index_bits, (7, 7, 6))
        self.assertEqual(cost.generator.ot_count, 6_716)
        self.assertEqual(cost.generator.payload_bits, 202_528)
        self.assertEqual(waksman_switch_count(cfg.columns), 2_432)
        self.assertEqual(cost.total_ot_count, 55_358)
        self.assertEqual(cost.total_payload_bits, 3_623_410)

    def test_dpf_reservation_components(self) -> None:
        # Independently enumerate columns and levels, using the implementation's
        # dense-prefix choice. One count here means one OT in each direction.
        for columns in (34, 64, 66, 88, 320):
            for depth in (3, 20):
                model = Model(domain_bits=depth)
                row_depth = (columns - 1).bit_length()
                dense_depth = min(depth, row_depth + 2)
                operations = []
                for _ in range(columns):
                    operations.extend([128] * depth)
                    operations.extend([256] * (depth - dense_depth))
                    operations.append(64)
                for _ in range(model.t):
                    operations.extend([128] * row_depth)
                    operations.append(64)
                cost = dpf_initialization_cost(model, columns)
                self.assertEqual(cost.ot_count, 2 * len(operations))
                self.assertEqual(cost.payload_bits, 2 * sum(operations))

    def test_characteristic_two_omits_only_sign_correlations(self) -> None:
        for columns in (64, 320):
            general = dpf_initialization_cost(Model(), columns)
            binary = dpf_initialization_cost(Model(characteristic_two=True), columns)
            self.assertEqual(general.ot_count - binary.ot_count, 2 * (columns + 16))
            self.assertEqual(general.payload_bits - binary.payload_bits, 128 * (columns + 16))

    def test_asymmetric_four_n_point(self) -> None:
        cfg = Configuration("asymmetric 4N", "reachability", 4, 0, 0,
                            bfs_roots=2, partition_sizes=(32, 32, 16, 8))
        cost = configuration_cost(Model(), cfg)
        self.assertEqual((cost.total_ot_count, cost.total_payload_bits),
                         (38_158, 1_162_710))


if __name__ == "__main__":
    unittest.main()
