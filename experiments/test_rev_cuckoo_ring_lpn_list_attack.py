import itertools
import math
import random

from rev_cuckoo_regular_noise import (
    LeakageChannel,
    support_lists_one_unknown,
)
from rev_cuckoo_ring_lpn_list_attack import (
    AffinePairGenerator,
    CandidateEntry,
    ScoredPairGenerator,
    add_syndromes,
    pair_beam_candidates,
    pair_surrogate_score,
    partial_collision_factors,
    one_trial,
    paired_bootstrap_gain_interval,
    subtract_syndromes,
)
from rev_cuckoo_ring_lpn_inner_beam import (
    RetentionRow,
    coordinate_retentions,
    fractional_top_retention,
    geometric_retention,
)


def test_syndrome_addition_and_subtraction_round_trip() -> None:
    left = bytes((0, 1, 16, 8))
    right = bytes((3, 16, 2, 10))
    total = add_syndromes(left, right, 17)
    assert subtract_syndromes(total, right, 17) == left


def test_pair_generators_enumerate_cartesian_product_once() -> None:
    entries = tuple(
        CandidateEntry(
            index=index,
            offsets=(index,),
            contribution=bytes((index,)),
            pair_score=float(index),
            exact_score=float(index),
        )
        for index in range(4)
    )
    uniform = AffinePairGenerator(entries, entries, 17, random.Random(1001))
    uniform_pairs = {
        (entry.left_index, entry.right_index) for entry in uniform
    }
    assert len(uniform_pairs) == 16

    scored = ScoredPairGenerator(
        entries, entries, "exact", 17, random.Random(1002)
    )
    scored_pairs = [(entry.left_index, entry.right_index) for entry in scored]
    assert len(set(scored_pairs)) == 16
    assert scored_pairs[0] == (3, 3)


def test_reduced_decoder_recovers_unique_planted_solution() -> None:
    polys = 4
    weight = 2
    block_length = 3
    candidates = tuple(
        itertools.product(range(block_length), repeat=weight)
    )
    channel = LeakageChannel.build("raw", 1, weight, weight)
    results = one_trial(
        trial=0,
        candidates=candidates,
        polys=polys,
        weight=weight,
        block_length=block_length,
        modulus=17,
        channel=channel,
        pair_slope=-0.5,
        seed=1003,
    )
    assert {result.selector for result in results} == {
        "uniform",
        "pair",
        "exact",
    }
    assert all(result.recovered_planted for result in results)
    assert all(result.optimistic_work > 0 for result in results)
    assert all(math.isfinite(result.factor_charged_work) for result in results)


def test_paired_bootstrap_is_zero_for_identical_work() -> None:
    interval = paired_bootstrap_gain_interval(
        [10, 20, 30], [10, 20, 30], 100, random.Random(1004)
    )
    assert interval == (0.0, 0.0)


def test_complete_prefix_score_matches_pair_surrogate() -> None:
    rng = random.Random(1005)
    polys = 2
    weight = 3
    block_length = 4
    known = tuple(rng.randrange(block_length) for _ in range(polys * weight))
    planted = tuple(rng.randrange(block_length) for _ in range(weight))
    channel = LeakageChannel.build("occupancy", 2, weight, weight)
    planted_lists = support_lists_one_unknown(
        known, planted, polys, weight
    )
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng)
        for active in planted_lists
    )
    factors = partial_collision_factors(
        known, descriptors, polys, weight, block_length
    )
    slope = -0.375
    beam = pair_beam_candidates(
        factors,
        weight,
        block_length,
        slope,
        block_length**weight,
        block_length**weight,
        random.Random(1006),
    )
    scores = []
    for entry in beam.candidates:
        active_lists = support_lists_one_unknown(
            known, entry.offsets, polys, weight
        )
        direct = pair_surrogate_score(active_lists, descriptors, slope)
        assert entry.pair_score == direct
        scores.append(entry.pair_score)
    assert scores == sorted(scores, reverse=True)
    assert len({entry.offsets for entry in beam.candidates}) == block_length**weight


def test_fractional_retention_handles_cutoff_ties() -> None:
    scored = [(3.0, (0,)), (2.0, (1,)), (2.0, (2,)), (1.0, (3,))]
    assert fractional_top_retention(scored, (0,), 2) == 1.0
    assert fractional_top_retention(scored, (1,), 2) == 0.5
    assert fractional_top_retention(scored, (2,), 2) == 0.5
    assert fractional_top_retention(scored, (3,), 2) == 0.0


def test_coordinate_retention_uses_geometric_rate() -> None:
    rows = [
        RetentionRow(trial, hidden, "x", retained, 1, 2, 1, 0)
        for trial, values in enumerate(((1.0, 0.0), (1.0, 1.0)))
        for hidden, retained in enumerate(values)
    ]
    rates = coordinate_retentions(rows, 2)
    assert rates == (1.0, 0.5)
    assert geometric_retention(rates) == math.sqrt(0.5)


if __name__ == "__main__":
    test_syndrome_addition_and_subtraction_round_trip()
    test_pair_generators_enumerate_cartesian_product_once()
    test_reduced_decoder_recovers_unique_planted_solution()
    test_paired_bootstrap_is_zero_for_identical_work()
    test_complete_prefix_score_matches_pair_surrogate()
    test_fractional_retention_handles_cutoff_ties()
    test_coordinate_retention_uses_geometric_rate()
