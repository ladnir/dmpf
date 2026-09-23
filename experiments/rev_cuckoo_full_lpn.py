"""Exact toy Ring-LPN game with correlated Reverse-Cuckoo leakage.

This is an information-theoretic enumerator, not a production-parameter
security estimator.  It preserves the implementation's regular-polynomial
support geometry:

    A[a,b,k] = {x[a,u] + y[b,(k-u) mod t] : u in [t]}.

One party's offsets x are included in the adversary view.  The unknown party's
offsets y are enumerated, a binary cyclic-ring syndrome is formed, and an
independent ideal w=2 Reverse-Cuckoo descriptor is planted for every support
list.  Exact placement counts give the likelihood of every syndrome-compatible
candidate.  The reported information densities satisfy, sample by sample,

    i(X; Lambda | known) = i(Y; Lambda | R,known)
                           + i(X; Lambda | R,Y,known).

The binary/unit-coefficient syndrome is only a structural proxy for the field
Ring-LPN game.  It is small enough to validate the likelihood and measurement
method before adding coefficient enumeration or a leakage-aware list decoder.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
import statistics
from dataclasses import dataclass
from typing import Iterable, Sequence


NEG_INF = -math.inf


def logaddexp2(a: float, b: float) -> float:
    if a == NEG_INF:
        return b
    if b == NEG_INF:
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + 2.0 ** (b - a))


def logsumexp2(values: Iterable[float]) -> float:
    total = NEG_INF
    for value in values:
        total = logaddexp2(total, value)
    return total


def mean(values: Sequence[float]) -> float:
    return sum(values) / len(values)


def standard_error(values: Sequence[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values) / math.sqrt(len(values))


def quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def rotate_left_bits(value: int, shift: int, width: int) -> int:
    mask = (1 << width) - 1
    shift %= width
    if shift == 0:
        return value & mask
    return ((value << shift) | (value >> (width - shift))) & mask


def cyclic_product(left: int, right: int, width: int) -> int:
    """Multiply binary polynomials modulo X^width + 1."""
    result = 0
    value = left
    while value:
        bit = value & -value
        result ^= rotate_left_bits(right, bit.bit_length() - 1, width)
        value ^= bit
    return result


def regular_poly(offsets: Sequence[int], block_length: int) -> int:
    result = 0
    for block, offset in enumerate(offsets):
        result |= 1 << (block * block_length + offset)
    return result


def falling_factorial(value: int, count: int) -> int:
    result = 1
    for item in range(value - count + 1, value + 1):
        result *= item
    return result


def log2_expected_placements(active_count: int, partition_size: int) -> float:
    # E_U[Z_A(H)] = (2d)_n / d^n for w=2.
    return math.log2(falling_factorial(2 * partition_size, active_count)) - (
        active_count * math.log2(partition_size)
    )


def placement_count_w2(
    active: Sequence[int], descriptor: tuple[Sequence[int], Sequence[int]], partition_size: int
) -> int:
    """Count injective placements in a bipartite cuckoo multigraph."""
    left, right = descriptor
    vertex_count = 2 * partition_size
    adjacency: list[list[int]] = [[] for _ in range(vertex_count)]
    edges: list[tuple[int, int]] = []
    for item in active:
        edge = (left[item], partition_size + right[item])
        edge_index = len(edges)
        edges.append(edge)
        adjacency[edge[0]].append(edge_index)
        adjacency[edge[1]].append(edge_index)

    seen_vertex = [False] * vertex_count
    seen_edge = [False] * len(edges)
    placements = 1
    for start in range(vertex_count):
        if seen_vertex[start] or not adjacency[start]:
            continue
        stack = [start]
        seen_vertex[start] = True
        component_vertices = 0
        component_edges = 0
        while stack:
            vertex = stack.pop()
            component_vertices += 1
            for edge_index in adjacency[vertex]:
                if not seen_edge[edge_index]:
                    seen_edge[edge_index] = True
                    component_edges += 1
                edge = edges[edge_index]
                other = edge[1] if edge[0] == vertex else edge[0]
                if not seen_vertex[other]:
                    seen_vertex[other] = True
                    stack.append(other)

        if component_edges > component_vertices:
            return 0
        if component_edges == component_vertices:
            placements *= 2
        elif component_edges + 1 == component_vertices:
            placements *= component_vertices
        else:
            raise AssertionError("invalid connected cuckoo component")
    return placements


def plant_descriptor_w2(
    active: Sequence[int], domain: int, partition_size: int, rng: random.Random
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Sample the ideal Reverse-Cuckoo descriptor by planting an injection."""
    if len(active) > 2 * partition_size:
        raise ValueError("more active items than cuckoo bins")
    left = [rng.randrange(partition_size) for _ in range(domain)]
    right = [rng.randrange(partition_size) for _ in range(domain)]
    assigned_bins = rng.sample(range(2 * partition_size), len(active))
    for item, assigned_bin in zip(active, assigned_bins):
        if assigned_bin < partition_size:
            left[item] = assigned_bin
        else:
            right[item] = assigned_bin - partition_size
    return tuple(left), tuple(right)


def support_lists(
    known: Sequence[int], unknown: Sequence[int], polys: int, weight: int
) -> tuple[tuple[int, ...], ...]:
    lists = []
    for known_poly in range(polys):
        for unknown_poly in range(polys):
            for output_block in range(weight):
                active = {
                    known[known_poly * weight + block]
                    + unknown[
                        unknown_poly * weight
                        + ((output_block - block) % weight)
                    ]
                    for block in range(weight)
                }
                lists.append(tuple(sorted(active)))
    return tuple(lists)


def syndrome(
    unknown: Sequence[int], multipliers: Sequence[int], polys: int, weight: int,
    block_length: int
) -> int:
    width = weight * block_length
    result = 0
    for poly in range(polys):
        offsets = unknown[poly * weight : (poly + 1) * weight]
        error = regular_poly(offsets, block_length)
        if poly == 0:
            result ^= error
        else:
            result ^= cyclic_product(multipliers[poly - 1], error, width)
    return result


def candidate_log_weight(
    candidate_lists: Sequence[Sequence[int]],
    descriptors: Sequence[tuple[Sequence[int], Sequence[int]]],
    partition_size: int,
) -> float:
    result = 0.0
    for active, descriptor in zip(candidate_lists, descriptors):
        placements = placement_count_w2(active, descriptor, partition_size)
        if placements == 0:
            return NEG_INF
        result += math.log2(placements) - log2_expected_placements(
            len(active), partition_size
        )
    return result


def posterior_entropy(log_weights: Sequence[float]) -> float:
    log_total = logsumexp2(log_weights)
    entropy = 0.0
    for log_weight in log_weights:
        if log_weight == NEG_INF:
            continue
        probability = 2.0 ** (log_weight - log_total)
        entropy -= probability * (log_weight - log_total)
    return entropy


@dataclass
class TrialResult:
    total_information: float
    distinguishing_information: float
    decoding_information: float
    base_entropy: float
    posterior_entropy: float
    posterior_min_entropy: float
    base_rank: float
    posterior_rank: float
    search_base_rank: float
    search_rank: float


def one_trial(
    candidates: Sequence[tuple[int, ...]], polys: int, weight: int,
    block_length: int, partition_size: int, rng: random.Random
) -> TrialResult:
    width = weight * block_length
    known = tuple(rng.randrange(block_length) for _ in range(polys * weight))
    true_index = rng.randrange(len(candidates))
    true_unknown = candidates[true_index]
    multipliers = tuple(rng.randrange(1 << width) for _ in range(polys - 1))
    true_syndrome = syndrome(
        true_unknown, multipliers, polys, weight, block_length
    )

    true_lists = support_lists(known, true_unknown, polys, weight)
    descriptors = tuple(
        plant_descriptor_w2(active, 2 * block_length, partition_size, rng)
        for active in true_lists
    )

    log_weights = []
    consistent_indices = []
    for index, candidate in enumerate(candidates):
        lists = support_lists(known, candidate, polys, weight)
        log_weights.append(
            candidate_log_weight(lists, descriptors, partition_size)
        )
        if syndrome(candidate, multipliers, polys, weight, block_length) == true_syndrome:
            consistent_indices.append(index)

    true_log_weight = log_weights[true_index]
    if true_log_weight == NEG_INF:
        raise AssertionError("planted support has zero likelihood")

    log_mean_all = logsumexp2(log_weights) - math.log2(len(candidates))
    consistent_weights = [log_weights[index] for index in consistent_indices]
    log_mean_consistent = logsumexp2(consistent_weights) - math.log2(
        len(consistent_indices)
    )
    total_information = true_log_weight - log_mean_all
    decoding_information = true_log_weight - log_mean_consistent
    distinguishing_information = log_mean_consistent - log_mean_all
    if abs(
        total_information - decoding_information - distinguishing_information
    ) > 1e-10:
        raise AssertionError("information-density decomposition failed")

    greater = sum(weight_value > true_log_weight + 1e-12 for weight_value in consistent_weights)
    tied = sum(abs(weight_value - true_log_weight) <= 1e-12 for weight_value in consistent_weights)
    posterior_rank = greater + (tied + 1.0) / 2.0
    search_greater = sum(
        weight_value > true_log_weight + 1e-12 for weight_value in log_weights
    )
    search_tied = sum(
        abs(weight_value - true_log_weight) <= 1e-12
        for weight_value in log_weights
    )
    return TrialResult(
        total_information=total_information,
        distinguishing_information=distinguishing_information,
        decoding_information=decoding_information,
        base_entropy=math.log2(len(consistent_indices)),
        posterior_entropy=posterior_entropy(consistent_weights),
        posterior_min_entropy=(
            logsumexp2(consistent_weights) - max(consistent_weights)
        ),
        base_rank=(len(consistent_indices) + 1.0) / 2.0,
        posterior_rank=posterior_rank,
        search_base_rank=(len(candidates) + 1.0) / 2.0,
        search_rank=search_greater + (search_tied + 1.0) / 2.0,
    )


def self_test() -> None:
    # Three parallel edges have no injective orientation into two vertices.
    descriptor = ((0, 0, 0), (0, 0, 0))
    if placement_count_w2((0, 1, 2), descriptor, 1) != 0:
        raise AssertionError("parallel-edge obstruction was accepted")

    # A path with two edges and three vertices has three placements.
    descriptor = ((0, 1), (0, 0))
    if placement_count_w2((0, 1), descriptor, 2) != 3:
        raise AssertionError("tree placement count is wrong")

    # Two parallel edges form a bipartite two-cycle with two placements.
    descriptor = ((0, 0), (0, 0))
    if placement_count_w2((0, 1), descriptor, 1) != 2:
        raise AssertionError("unicyclic placement count is wrong")

    # Averaging Z over every two-item descriptor recovers E[Z]=(2d)_n/d^n.
    placement_sum = 0
    descriptor_count = 0
    for left in itertools.product(range(2), repeat=2):
        for right in itertools.product(range(2), repeat=2):
            placement_sum += placement_count_w2((0, 1), (left, right), 2)
            descriptor_count += 1
    expected = 2.0 ** log2_expected_placements(2, 2)
    if abs(placement_sum / descriptor_count - expected) > 1e-12:
        raise AssertionError("placement likelihood is not normalized")


def print_metric(name: str, values: Sequence[float]) -> None:
    print(
        f"{name}: mean={mean(values):.6f} se={standard_error(values):.6f} "
        f"p50={quantile(values, 0.5):.6f} "
        f"p90={quantile(values, 0.9):.6f} "
        f"p99={quantile(values, 0.99):.6f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=2)
    parser.add_argument("--weight", type=int, default=2)
    parser.add_argument("--block-length", type=int, default=4)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if args.polys < 1 or args.weight < 1 or args.block_length < 1:
        parser.error("polys, weight, and block length must be positive")
    if args.trials < 1:
        parser.error("trials must be positive")
    partition_size = args.partition_size or args.weight
    if partition_size < 1:
        parser.error("partition size must be positive")

    self_test()
    candidates = tuple(
        itertools.product(
            range(args.block_length), repeat=args.polys * args.weight
        )
    )
    rng = random.Random(args.seed)
    results = [
        one_trial(
            candidates,
            args.polys,
            args.weight,
            args.block_length,
            partition_size,
            rng,
        )
        for _ in range(args.trials)
    ]

    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} d={partition_size} "
        f"lists={args.polys * args.polys * args.weight} "
        f"candidates={len(candidates)} trials={args.trials} seed={args.seed}"
    )
    print_metric(
        "I(X;Lambda|known) density",
        [result.total_information for result in results],
    )
    print_metric(
        "k_dist density",
        [result.distinguishing_information for result in results],
    )
    print_metric(
        "k_dec density", [result.decoding_information for result in results]
    )
    entropy_loss = [
        result.base_entropy - result.posterior_entropy for result in results
    ]
    print_metric("posterior entropy loss", entropy_loss)
    min_entropy_loss = [
        result.base_entropy - result.posterior_min_entropy for result in results
    ]
    print_metric("posterior min-entropy loss", min_entropy_loss)
    attack_gain = [
        math.log2(result.base_rank / result.posterior_rank) for result in results
    ]
    print_metric("posterior-rank gain", attack_gain)
    search_gain = [
        math.log2(result.search_base_rank / result.search_rank)
        for result in results
    ]
    print_metric("support-search mean-log rank gain", search_gain)
    expected_work_gain = -math.log2(
        mean([result.search_rank / result.search_base_rank for result in results])
    )
    print(f"support-search expected-work gain={expected_work_gain:.6f}")
    print(
        "decomposition_check="
        f"{mean([result.total_information for result in results]) - mean([result.distinguishing_information for result in results]) - mean([result.decoding_information for result in results]):.3e}"
    )


if __name__ == "__main__":
    main()
