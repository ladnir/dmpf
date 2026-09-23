import math
import random

from rev_cuckoo_full_lpn import placement_count_w2
from rev_cuckoo_goldreich import (
    GoldreichFeatureRoot,
    SharedGoldreichChannel,
    evaluate_solution,
    gf2_rank,
    sample_affine_solution,
)
from rev_cuckoo_regular_noise import LeakageChannel


def test_affine_solution_satisfies_system() -> None:
    rows = (0b10101, 0b01110, 0b11011)
    right = (1, 3, 2)
    solution = sample_affine_solution(rows, right, 5, 2, random.Random(1))
    assert solution is not None
    assert tuple(evaluate_solution(row, solution) for row in rows) == right


def test_affine_solution_rejects_inconsistent_dependency() -> None:
    solution = sample_affine_solution(
        (0b101, 0b101), (0, 1), 3, 1, random.Random(2)
    )
    assert solution is None


def test_goldreich_translation_maps_dummy_to_zero() -> None:
    root = GoldreichFeatureRoot.sample(32, 11, random.Random(3))
    assert root.input_bits == 8
    assert root.intermediate_bits == 16
    assert len(root.partitions[0]) == 32
    # The dummy is not stored, but recomputing translation means its row is 0.
    assert all(feature < 1 << 11 for side in root.partitions for feature in side)


def test_shared_sampler_plants_a_valid_descriptor() -> None:
    base = LeakageChannel.build("raw", 1, 4, 4)
    channel = SharedGoldreichChannel(base, linear_security=8)
    active = (1, 5, 9, 12)
    (descriptor,) = channel.sample_batch((active,), 16, random.Random(4))
    assert placement_count_w2(active, descriptor, 4) > 0
    assert channel.root is not None


def test_shared_sampler_rejects_non_power_of_two_partition() -> None:
    base = LeakageChannel.build("raw", 1, 3, 3)
    try:
        SharedGoldreichChannel(base)
    except ValueError as error:
        assert "power-of-two" in str(error)
    else:
        raise AssertionError("accepted an implementation-invalid partition")


def test_exact_rank_likelihood_matches_ideal_at_full_rank() -> None:
    base = LeakageChannel.build("raw", 1, 4, 4)
    channel = SharedGoldreichChannel(base, linear_security=8, exact_rank_limit=4)
    active = (0, 1, 2, 3)
    (descriptor,) = channel.sample_batch((active,), 8, random.Random(5))
    assert channel.root is not None
    assert gf2_rank([channel.root.partitions[0][item] for item in active]) == 4
    assert gf2_rank([channel.root.partitions[1][item] for item in active]) == 4
    assert math.isclose(
        channel.log_likelihood(active, descriptor),
        base.log_likelihood(active, descriptor),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def self_test() -> None:
    test_affine_solution_satisfies_system()
    test_affine_solution_rejects_inconsistent_dependency()
    test_goldreich_translation_maps_dummy_to_zero()
    test_shared_sampler_plants_a_valid_descriptor()
    test_shared_sampler_rejects_non_power_of_two_partition()
    test_exact_rank_likelihood_matches_ideal_at_full_rank()


if __name__ == "__main__":
    self_test()
