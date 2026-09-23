"""Exact toy diagnostic for Reverse-Cuckoo leakage on regular noise.

This experiment uses the implementation-aligned product-list geometry from
``rev_cuckoo_full_lpn.py`` but deliberately omits the syndrome.  For each
planted descriptor batch it enumerates the posterior on the unknown party's
regular offsets and separates:

* full posterior information;
* information in the one-offset marginals;
* pairwise posterior correlations; and
* the gain of a marginally biased, regular Prange-style information-set choice.

The last metric chooses the least likely coordinates from every regular block
and measures the exact posterior probability that all selected coordinates are
error free.  It is only a decoder diagnostic, not a production attack.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
import statistics
from collections.abc import Sequence
from dataclasses import dataclass

from rev_cuckoo_full_lpn import (
    NEG_INF,
    log2_expected_placements,
    logsumexp2,
    plant_descriptor_w2,
    placement_count_w2,
    self_test as full_lpn_self_test,
)
from rev_cuckoo_k_tradeoff import (
    exact_occupancy_race_factors,
    exact_occupancy_weights,
)


Descriptor = tuple[tuple[int, ...], tuple[int, ...]]


def descriptor_deficit(
    active: Sequence[int], descriptor: Descriptor
) -> int:
    left, right = descriptor
    return (
        2 * len(active)
        - len({left[item] for item in active})
        - len({right[item] for item in active})
    )


def uniform_descriptor_w2(
    domain: int, partition_size: int, rng: random.Random
) -> Descriptor:
    return (
        tuple(rng.randrange(partition_size) for _ in range(domain)),
        tuple(rng.randrange(partition_size) for _ in range(domain)),
    )


def exact_feasible_probability(items: int, partition_size: int) -> float:
    """Enumerate the small active-endpoint domain used by this toy study."""

    configuration_count = partition_size ** (2 * items)
    if configuration_count > 2_000_000:
        raise ValueError(
            "feasible-channel enumeration is too large; use a smaller toy "
            "instance or the raw/occupancy channel"
        )
    feasible = 0
    active = tuple(range(items))
    for endpoints in itertools.product(
        range(partition_size), repeat=2 * items
    ):
        descriptor = (endpoints[:items], endpoints[items:])
        feasible += placement_count_w2(
            active, descriptor, partition_size
        ) > 0
    return feasible / configuration_count


@dataclass(frozen=True)
class LeakageChannel:
    mode: str
    selector_k: int
    partition_size: int
    occupancy_rates: dict[int, dict[int, float]]
    race_factors: dict[int, dict[int, float]]
    feasible_probabilities: dict[int, float]

    @classmethod
    def build(
        cls,
        mode: str,
        selector_k: int,
        max_items: int,
        partition_size: int,
    ) -> "LeakageChannel":
        if mode == "occupancy" and selector_k < 1:
            raise ValueError("occupancy selector K must be positive")
        occupancy_rates = {}
        race_factors = {}
        feasible_probabilities = {}
        if mode == "occupancy":
            for items in range(max_items + 1):
                bins = (partition_size, partition_size)
                occupancy_rates[items] = exact_occupancy_weights(items, bins)
                race_factors[items] = exact_occupancy_race_factors(
                    items, bins, selector_k
                )
        elif mode == "feasible":
            for items in range(max_items + 1):
                feasible_probabilities[items] = exact_feasible_probability(
                    items, partition_size
                )
        elif mode != "raw":
            raise ValueError(f"unknown leakage channel: {mode}")
        return cls(
            mode,
            selector_k,
            partition_size,
            occupancy_rates,
            race_factors,
            feasible_probabilities,
        )

    def sample(
        self,
        active: Sequence[int],
        domain: int,
        rng: random.Random,
    ) -> Descriptor:
        if self.mode == "raw":
            return plant_descriptor_w2(
                active, domain, self.partition_size, rng
            )
        if self.mode == "occupancy":
            proposals = [
                plant_descriptor_w2(
                    active, domain, self.partition_size, rng
                )
                for _ in range(self.selector_k)
            ]
            rates = self.occupancy_rates[len(active)]
            weights = [
                rates[descriptor_deficit(active, descriptor)]
                for descriptor in proposals
            ]
            threshold = rng.random() * sum(weights)
            prefix = 0.0
            for descriptor, weight in zip(proposals, weights):
                prefix += weight
                if threshold <= prefix:
                    return descriptor
            return proposals[-1]

        while True:
            descriptor = uniform_descriptor_w2(
                domain, self.partition_size, rng
            )
            if placement_count_w2(
                active, descriptor, self.partition_size
            ):
                return descriptor

    def log_likelihood(
        self, active: Sequence[int], descriptor: Descriptor
    ) -> float:
        placements = placement_count_w2(
            active, descriptor, self.partition_size
        )
        if placements == 0:
            return NEG_INF
        items = len(active)
        if self.mode == "feasible":
            return -math.log2(self.feasible_probabilities[items])

        result = math.log2(placements) - log2_expected_placements(
            items, self.partition_size
        )
        if self.mode == "occupancy":
            deficit = descriptor_deficit(active, descriptor)
            result += math.log2(self.race_factors[items][deficit])
        return result

    def candidate_log_weight(
        self,
        candidate_lists: Sequence[Sequence[int]],
        descriptors: Sequence[Descriptor],
    ) -> float:
        result = 0.0
        for active, descriptor in zip(candidate_lists, descriptors):
            contribution = self.log_likelihood(active, descriptor)
            if contribution == NEG_INF:
                return NEG_INF
            result += contribution
        return result


def support_lists_one_unknown(
    known: Sequence[int], unknown: Sequence[int], polys: int, weight: int
) -> tuple[tuple[int, ...], ...]:
    """Return the P*t factors involving one unknown polynomial.

    The full P-by-P product geometry factors over the unknown-polynomial index.
    For one fixed unknown polynomial, list (a,k) is

        {x[a,u] + y[(k-u) mod t] : u in [t]}.
    """
    lists = []
    for known_poly in range(polys):
        for output_block in range(weight):
            active = {
                known[known_poly * weight + block]
                + unknown[(output_block - block) % weight]
                for block in range(weight)
            }
            lists.append(tuple(sorted(active)))
    return tuple(lists)


def normalize(log_weights: Sequence[float]) -> list[float]:
    log_total = logsumexp2(log_weights)
    return [
        0.0 if value == NEG_INF else 2.0 ** (value - log_total)
        for value in log_weights
    ]


def entropy(probabilities: Sequence[float]) -> float:
    return -sum(p * math.log2(p) for p in probabilities if p)


def marginal_tables(
    candidates: Sequence[tuple[int, ...]],
    probabilities: Sequence[float],
    variable_count: int,
    block_length: int,
) -> list[list[float]]:
    result = [[0.0] * block_length for _ in range(variable_count)]
    for candidate, probability in zip(candidates, probabilities):
        for variable, value in enumerate(candidate):
            result[variable][value] += probability
    return result


def marginal_information(marginals: Sequence[Sequence[float]]) -> float:
    block_length = len(marginals[0])
    return sum(
        math.log2(block_length) - entropy(marginal)
        for marginal in marginals
    )


def pairwise_mutual_information(
    candidates: Sequence[tuple[int, ...]],
    probabilities: Sequence[float],
    marginals: Sequence[Sequence[float]],
    block_length: int,
) -> tuple[float, float]:
    values = []
    for left in range(len(marginals)):
        for right in range(left + 1, len(marginals)):
            joint = [[0.0] * block_length for _ in range(block_length)]
            for candidate, probability in zip(candidates, probabilities):
                joint[candidate[left]][candidate[right]] += probability
            information = 0.0
            for left_value in range(block_length):
                for right_value in range(block_length):
                    probability = joint[left_value][right_value]
                    if probability:
                        information += probability * math.log2(
                            probability
                            / (
                                marginals[left][left_value]
                                * marginals[right][right_value]
                            )
                        )
            values.append(information)
    return statistics.mean(values) if values else 0.0, max(values, default=0.0)


def rank(
    scores: Sequence[float], true_index: int, tolerance: float = 1e-12
) -> float:
    true_score = scores[true_index]
    greater = sum(score > true_score + tolerance for score in scores)
    tied = sum(abs(score - true_score) <= tolerance for score in scores)
    return greater + (tied + 1.0) / 2.0


def marginal_scores(
    candidates: Sequence[tuple[int, ...]],
    marginals: Sequence[Sequence[float]],
) -> list[float]:
    result = []
    for candidate in candidates:
        score = 0.0
        for variable, value in enumerate(candidate):
            probability = marginals[variable][value]
            if probability == 0.0:
                score = NEG_INF
                break
            score += math.log2(probability)
        result.append(score)
    return result


def prange_gain(
    candidates: Sequence[tuple[int, ...]],
    probabilities: Sequence[float],
    marginals: Sequence[Sequence[float]],
    block_length: int,
    selected_per_block: int,
) -> float:
    excluded = [
        set(sorted(range(block_length), key=lambda value: marginal[value])[:selected_per_block])
        for marginal in marginals
    ]
    posterior_success = sum(
        probability
        for candidate, probability in zip(candidates, probabilities)
        if all(value not in excluded[variable] for variable, value in enumerate(candidate))
    )
    uniform_success = (
        (block_length - selected_per_block) / block_length
    ) ** len(marginals)
    return math.log2(posterior_success / uniform_success)


def optimal_pair_prange_gain(
    candidates: Sequence[tuple[int, ...]],
    probabilities: Sequence[float],
    block_length: int,
    selected_per_block: int,
) -> float:
    """Optimize the two-block error-free event over coordinated subsets."""
    if len(candidates[0]) != 2:
        raise ValueError("pair optimizer requires exactly two variables")
    joint = [[0.0] * block_length for _ in range(block_length)]
    for candidate, probability in zip(candidates, probabilities):
        joint[candidate[0]][candidate[1]] += probability

    kept = block_length - selected_per_block
    best = 0.0
    for kept_left in itertools.combinations(range(block_length), kept):
        right_mass = [
            sum(joint[left][right] for left in kept_left)
            for right in range(block_length)
        ]
        mass = sum(sorted(right_mass, reverse=True)[:kept])
        best = max(best, mass)
    uniform_success = (kept / block_length) ** 2
    return math.log2(best / uniform_success)


def one_trial(
    candidates: Sequence[tuple[int, ...]],
    polys: int,
    weight: int,
    block_length: int,
    partition_size: int,
    selected_counts: Sequence[int],
    channel: LeakageChannel,
    rng: random.Random,
) -> dict[str, float]:
    known = tuple(rng.randrange(block_length) for _ in range(polys * weight))
    true_index = rng.randrange(len(candidates))
    true_unknown = candidates[true_index]
    true_lists = support_lists_one_unknown(known, true_unknown, polys, weight)
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng)
        for active in true_lists
    )

    log_weights = []
    for candidate in candidates:
        lists = support_lists_one_unknown(known, candidate, polys, weight)
        log_weights.append(channel.candidate_log_weight(lists, descriptors))
    probabilities = normalize(log_weights)
    marginals = marginal_tables(
        candidates, probabilities, weight, block_length
    )

    full_information = math.log2(len(candidates)) - entropy(probabilities)
    one_variable_information = marginal_information(marginals)
    mean_pair_mi, max_pair_mi = pairwise_mutual_information(
        candidates, probabilities, marginals, block_length
    )
    full_rank = rank(log_weights, true_index)
    marginal_rank = rank(marginal_scores(candidates, marginals), true_index)
    base_rank = (len(candidates) + 1.0) / 2.0

    result = {
        "full_information": full_information,
        "one_variable_information": one_variable_information,
        "correlation_information": full_information - one_variable_information,
        "mean_pair_mi": mean_pair_mi,
        "max_pair_mi": max_pair_mi,
        "full_rank_gain": math.log2(base_rank / full_rank),
        "marginal_rank_gain": math.log2(base_rank / marginal_rank),
    }
    for selected in selected_counts:
        result[f"prange_gain_{selected}"] = prange_gain(
            candidates,
            probabilities,
            marginals,
            block_length,
            selected,
        )
        if weight == 2 and block_length <= 12:
            result[f"joint_prange_gain_{selected}"] = optimal_pair_prange_gain(
                candidates,
                probabilities,
                block_length,
                selected,
            )
    return result


def summarize(name: str, values: Sequence[float]) -> None:
    stderr = statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else 0.0
    print(f"{name}: mean={statistics.mean(values):.6f} se={stderr:.6f}")


def self_test() -> None:
    full_lpn_self_test()
    candidates = tuple(itertools.product(range(2), repeat=2))
    probabilities = [0.25] * 4
    marginals = marginal_tables(candidates, probabilities, 2, 2)
    if abs(marginal_information(marginals)) > 1e-12:
        raise AssertionError("uniform marginals have nonzero information")
    if abs(prange_gain(candidates, probabilities, marginals, 2, 1)) > 1e-12:
        raise AssertionError("uniform posterior changes Prange success")

    active = (0, 1)
    descriptors = [
        (endpoints[:2], endpoints[2:])
        for endpoints in itertools.product(range(2), repeat=4)
    ]
    for mode in ("raw", "occupancy", "feasible"):
        channel = LeakageChannel.build(mode, 4, 2, 2)
        density_mean = sum(
            0.0
            if (log_density := channel.log_likelihood(active, descriptor))
            == NEG_INF
            else 2.0 ** log_density
            for descriptor in descriptors
        ) / len(descriptors)
        if not math.isclose(density_mean, 1.0, rel_tol=1e-9):
            raise AssertionError(
                f"{mode} channel likelihood is not normalized: {density_mean}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=2)
    parser.add_argument("--weight", type=int, default=2)
    parser.add_argument("--block-length", type=int, default=4)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--channel",
        choices=("raw", "occupancy", "feasible"),
        default="raw",
    )
    parser.add_argument("--selector-k", type=int, default=4)
    args = parser.parse_args()

    if min(args.polys, args.weight, args.block_length, args.trials) < 1:
        parser.error("all dimensions and trials must be positive")
    partition_size = args.partition_size or args.weight
    if args.selector_k < 1:
        parser.error("selector K must be positive")
    selected_counts = sorted(
        {
            max(1, args.block_length // 4),
            max(1, args.block_length // 2),
            min(args.block_length - 1, max(1, 3 * args.block_length // 4)),
        }
    )

    self_test()
    candidates = tuple(
        itertools.product(
            range(args.block_length), repeat=args.weight
        )
    )
    channel = LeakageChannel.build(
        args.channel,
        args.selector_k,
        args.weight,
        partition_size,
    )
    rng = random.Random(args.seed)
    results = [
        one_trial(
            candidates,
            args.polys,
            args.weight,
            args.block_length,
            partition_size,
            selected_counts,
            channel,
            rng,
        )
        for _ in range(args.trials)
    ]

    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} d={partition_size} "
        f"unknown_polys=1 lists={args.polys * args.weight} "
        f"candidates={len(candidates)} channel={args.channel} "
        f"selector_k={args.selector_k if args.channel == 'occupancy' else '-'} "
        f"trials={args.trials} seed={args.seed}"
    )
    for metric in results[0]:
        summarize(metric, [result[metric] for result in results])


if __name__ == "__main__":
    main()
