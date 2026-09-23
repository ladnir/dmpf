import math
import random

from rev_cuckoo_posterior_info import exact_posterior_information
from rev_cuckoo_regular_noise import LeakageChannel, support_lists_one_unknown


def test_exact_information_losses_are_ordered() -> None:
    polys = 2
    weight = 3
    block_length = 8
    known = (0, 3, 6, 1, 4, 7)
    candidate = (2, 4, 7)
    channel = LeakageChannel.build("occupancy", 1, weight, 4)
    rng = random.Random(37)
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng)
        for active in support_lists_one_unknown(
            known, candidate, polys, weight
        )
    )
    log_evidence, posterior_mean, information, renyi2 = (
        exact_posterior_information(
            known,
            descriptors,
            weight,
            block_length,
            channel,
        )
    )
    assert math.isfinite(log_evidence)
    assert math.isfinite(posterior_mean)
    assert information >= -1e-12
    assert renyi2 + 1e-12 >= information
    assert renyi2 <= weight * math.log2(block_length) + 1e-12
