from rev_cuckoo_cost import (
    Model,
    bounded_sum_products,
    estimate,
    inverse_cdf_sampler_ots,
    next_power_of_two,
    occupancy_score_ots,
    parse_ints,
    placement_candidate_ots,
    popcount_products,
    public_table_lookup_products,
    rows,
    selector_ots,
    selector_payload_bits,
    waksman_switches,
)


def test_waksman_matches_small_networks() -> None:
    assert [waksman_switches(n) for n in range(1, 9)] == [0, 1, 3, 6, 8, 12, 15, 20]


def test_integer_ranges_are_inclusive() -> None:
    assert parse_ints("12:16:2,19") == [12, 14, 16, 19]


def test_current_geometry_rounds_only_d() -> None:
    model = Model(geometry="current")
    assert model.partition_size(16) == 16
    assert model.partition_size(17) == 32
    assert next_power_of_two(17) == 32
    assert next_power_of_two(1) == 1


def test_k_only_changes_debiasing_preprocessing() -> None:
    model = Model()
    raw = estimate(model, 16, 1)
    debiased = estimate(model, 16, 4)
    assert debiased.selected_dmpf == raw.selected_dmpf
    assert debiased.position_add == raw.position_add
    assert debiased.tensor_per_expansion == raw.tensor_per_expansion
    assert raw.debias == 0
    assert debiased.debias > 0


def test_free_t_is_not_rounded() -> None:
    model = Model(geometry="free")
    assert model.partition_size(17) == 17


def test_explicit_partition_size_overrides_geometry() -> None:
    model = Model(geometry="current", partition_size_override=24)
    assert model.partition_size(17) == 24


def test_unequal_partition_sizes_are_supported() -> None:
    model = Model(partition_sizes_override=(16, 32))
    assert model.partition_sizes(16) == (16, 32)
    assert estimate(model, 16, 1).selected_dmpf > estimate(Model(), 16, 1).selected_dmpf


def test_reference_ratio_compares_noise_choices() -> None:
    records = rows(Model(), [16, 18], [1])
    assert records[0]["ratio_to_reference"] == 1.0
    assert records[1]["ratio_to_reference"] > 1.0


def test_exact_occupancy_score_for_asymmetric_candidate() -> None:
    # 120 comparisons at 7 and 8 OTs, 210 OR-tree ANDs, two 22-AND
    # partition popcounts, and a 4-AND sum of their bounded outputs.
    assert popcount_products(15) == 22
    assert bounded_sum_products((15, 15)) == 4
    assert occupancy_score_ots(16, (4, 5)) == 2_316


def test_weighted_k4_selector_uses_wide_muxes() -> None:
    assert public_table_lookup_products(5) == 15
    # Four 5-bit-address lookups, then three 80-bit comparisons and three
    # wide conditional moves.
    assert selector_ots(4, 16, 2, 80) == 606
    assert selector_payload_bits(4, 16, 2, 80, 160) == 11_520


def test_asymmetric_placement_sampler_is_explicit() -> None:
    assert inverse_cdf_sampler_ots(17, 64) == 722
    assert placement_candidate_ots(16, (16, 16), 64) == 288
    assert placement_candidate_ots(16, (16, 32), 64) == 1_266
