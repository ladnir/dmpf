#!/usr/bin/env python3
"""Rigorous bounds for rejection-debiased two-choice Reverse Cuckoo.

For a two-choice descriptor, let K be the number of repeated endpoints:

    K = 2t - (# distinct left bins) - (# distinct right bins).

Every feasible cuckoo graph with K <= r and r >= 2 has at least
2^(t-r+1) valid placements.  Rejecting K > r and otherwise accepting an
ordinary placement-first sample with probability 2^(t-r+1) / Z makes the
accepted graph exactly uniform over {Z > 0, K <= r}.  Its distance from a
uniform graph is therefore at most ordinary cuckoo failure plus the occupancy
tail Pr[K > r].

The probability c/Z can be realized without division: conditional on its
graph, the planted placement is uniform among the Z valid placements.  Fix a
canonical ordering and accept exactly the first c placements.

The occupancy probability is evaluated exactly as an integer ratio.  The
ordinary cuckoo-failure term is the rigorous Hall-witness bound implemented in
rev_cuckoo_hall.py.
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass
from fractions import Fraction

from rev_cuckoo_hall import hall_witness_bound


NEG_INF = -math.inf


@dataclass(frozen=True)
class Row:
    items: int
    partition_size: int
    collision_cutoff: int


DEFAULT_ROWS = (
    # These are the smallest implementation-supported powers of two at which
    # the rigorous cuckoo-failure plus occupancy-tail union bound reaches
    # 2^-48 per set, hence 2^-40 over a batch of 256 sets.  Power-of-two d is
    # also lets bins be labeled by F_2^k, matching the characteristic-two
    # solver.  Full binary row rank then gives surjectivity componentwise.
    Row(16, 32_768, 6),
    Row(64, 65_536, 8),
    Row(128, 131_072, 10),
)


def log2_integer(value: int) -> float:
    """Return log2(value) without overflowing a binary64 conversion."""
    if value <= 0:
        return NEG_INF
    shift = max(0, value.bit_length() - 53)
    return math.log2(value >> shift) + shift


def logadd2(left: float, right: float) -> float:
    if left == NEG_INF:
        return right
    if right == NEG_INF:
        return left
    if left < right:
        left, right = right, left
    return left + math.log2(1.0 + math.exp2(right - left))


def endpoint_collision_counts(items: int, bins: int) -> list[int]:
    """Count one-side maps by K_s = items - number of occupied bins."""
    stirling = [[0] * (items + 1) for _ in range(items + 1)]
    stirling[0][0] = 1
    for n in range(1, items + 1):
        for occupied in range(1, n + 1):
            stirling[n][occupied] = (
                stirling[n - 1][occupied - 1]
                + occupied * stirling[n - 1][occupied]
            )

    counts = [0] * items
    falling = 1
    falling_by_occupied = [1] * (items + 1)
    for occupied in range(1, items + 1):
        falling *= bins - occupied + 1
        falling_by_occupied[occupied] = falling
    for repeated in range(items):
        occupied = items - repeated
        counts[repeated] = (
            falling_by_occupied[occupied] * stirling[items][occupied]
        )
    return counts


def two_side_tail_log2(items: int, bins: int, cutoff: int) -> float:
    """Return exact log2 Pr[K_left + K_right > cutoff]."""
    counts = endpoint_collision_counts(items, bins)
    numerator = sum(
        left_count * right_count
        for left_repeated, left_count in enumerate(counts)
        for right_repeated, right_count in enumerate(counts)
        if left_repeated + right_repeated > cutoff
    )
    return log2_integer(numerator) - 2 * items * math.log2(bins)


def analyze(row: Row, batch: int) -> dict[str, float | int]:
    if row.items < 3:
        raise ValueError("the stated placement lower bound requires t >= 3")
    if not 2 <= row.collision_cutoff <= row.items:
        raise ValueError("collision cutoff must lie in [2,t]")
    failure_log2, _ = hall_witness_bound(
        row.items, 2, row.partition_size
    )
    occupancy_log2 = two_side_tail_log2(
        row.items, row.partition_size, row.collision_cutoff
    )
    bad_log2 = logadd2(failure_log2, occupancy_log2)
    bad_probability = math.exp2(bad_log2)
    # The cutoff sampler is uniform on {Z > 0, K <= r}; its TV distance is
    # the probability of the complement, bounded by this union bound.
    tv_log2 = bad_log2
    expected_trials = math.exp2(row.collision_cutoff - 1) / (
        1.0 - bad_probability
    )
    return {
        "items": row.items,
        "partition_size": row.partition_size,
        "cutoff": row.collision_cutoff,
        "failure_log2": failure_log2,
        "occupancy_log2": occupancy_log2,
        "tv_log2": tv_log2,
        "batch_tv_log2": min(0.0, tv_log2 + math.log2(batch)),
        "expected_trials_bound": expected_trials,
    }


def brute_placement_count(edges: tuple[tuple[int, int], ...]) -> int:
    count = 0
    for choices in itertools.product((0, 1), repeat=len(edges)):
        selected = [
            (side, edges[item][side])
            for item, side in enumerate(choices)
        ]
        count += len(set(selected)) == len(selected)
    return count


def valid_placement_choices(
    edges: tuple[tuple[int, int], ...]
) -> list[tuple[int, ...]]:
    """Return valid side-choice vectors in a deterministic canonical order."""
    valid = []
    for choices in itertools.product((0, 1), repeat=len(edges)):
        selected = [
            (side, edges[item][side])
            for item, side in enumerate(choices)
        ]
        if len(set(selected)) == len(selected):
            valid.append(choices)
    return valid


def component_placement_count(
    edges: tuple[tuple[int, int], ...], bins: int
) -> int:
    """Count placements using the tree/unicyclic component formula."""
    vertex_count = 2 * bins
    parent = list(range(vertex_count))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    incident: set[int] = set()
    for left, right in edges:
        left_vertex = left
        right_vertex = bins + right
        incident.add(left_vertex)
        incident.add(right_vertex)
        union(left_vertex, right_vertex)

    vertices_by_root: dict[int, int] = {}
    edges_by_root: dict[int, int] = {}
    for vertex in incident:
        root = find(vertex)
        vertices_by_root[root] = vertices_by_root.get(root, 0) + 1
    for left, _ in edges:
        root = find(left)
        edges_by_root[root] = edges_by_root.get(root, 0) + 1

    placements = 1
    for root, edge_count in edges_by_root.items():
        component_vertices = vertices_by_root[root]
        if edge_count == component_vertices - 1:
            placements *= component_vertices
        elif edge_count == component_vertices:
            placements *= 2
        else:
            return 0
    return placements


def self_test() -> None:
    """Exhaustively verify the graph identities on a small multigraph space."""
    items = 4
    bins = 3
    checked = 0
    graph_counts: list[int] = []
    graph_repetitions: list[int] = []
    for flat in itertools.product(range(bins), repeat=2 * items):
        edges = tuple(
            (flat[2 * item], flat[2 * item + 1])
            for item in range(items)
        )
        brute = brute_placement_count(edges)
        by_components = component_placement_count(edges, bins)
        assert brute == by_components

        distinct_left = len({left for left, _ in edges})
        distinct_right = len({right for _, right in edges})
        repeated = 2 * items - distinct_left - distinct_right
        graph_counts.append(brute)
        graph_repetitions.append(repeated)
        if brute:
            for cutoff in range(2, items + 1):
                if repeated <= cutoff:
                    assert brute >= 1 << (items - cutoff + 1)
        checked += 1

    # A placement-first proposal is uniform on compatible (graph, placement)
    # pairs.  Accepting the first c placements in a canonical order must leave
    # exactly c accepted pairs for every good graph, hence a uniform graph.
    total_placements = sum(graph_counts)
    graph_total = bins ** (2 * items)
    failed_graphs = sum(placements == 0 for placements in graph_counts)
    hall_log2, _ = hall_witness_bound(items, 2, bins)
    actual_failure_log2 = math.log2(Fraction(failed_graphs, graph_total))
    assert hall_log2 + 1e-12 >= actual_failure_log2
    for cutoff in range(2, items + 1):
        threshold = 1 << (items - cutoff + 1)
        enumerated_tail = sum(
            repeated > cutoff for repeated in graph_repetitions
        )
        exact_tail_log2 = two_side_tail_log2(items, bins, cutoff)
        if enumerated_tail:
            assert math.isclose(
                exact_tail_log2,
                math.log2(Fraction(enumerated_tail, graph_total)),
                abs_tol=1e-12,
            )
        else:
            assert exact_tail_log2 == NEG_INF
        good = [
            index
            for index, (placements, repeated) in enumerate(
                zip(graph_counts, graph_repetitions)
            )
            if placements > 0 and repeated <= cutoff
        ]
        accepted_masses = set()
        for index in good:
            # Reconstruct the indexed graph directly from its base-bins
            # representation, matching itertools.product's enumeration.
            digits = [0] * (2 * items)
            value = index
            for position in range(2 * items - 1, -1, -1):
                digits[position] = value % bins
                value //= bins
            edges = tuple(
                (digits[2 * item], digits[2 * item + 1])
                for item in range(items)
            )
            placements = valid_placement_choices(edges)
            assert len(placements) == graph_counts[index]
            accepted_pairs = sum(
                rank < threshold for rank, _ in enumerate(placements)
            )
            assert accepted_pairs == threshold
            accepted_masses.add(Fraction(accepted_pairs, total_placements))
        assert accepted_masses == {Fraction(threshold, total_placements)}

        # TV(U conditioned on good, U) is exactly U[not good].
        conditional_tv = Fraction(graph_total - len(good), graph_total)
        explicit_tv = Fraction(0)
        good_set = set(good)
        for index in range(graph_total):
            conditional = Fraction(1, len(good)) if index in good_set else 0
            explicit_tv += abs(conditional - Fraction(1, graph_total))
        assert explicit_tv / 2 == conditional_tv
    print(f"self-test passed: {checked} labeled two-choice multigraphs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.batch < 1:
        parser.error("--batch must be positive")
    if args.self_test:
        self_test()

    print(
        "t,d,r,log2_cuckoo_failure,log2_occupancy_tail,"
        "log2_tv_bound,log2_batch_tv_bound,expected_trials_bound"
    )
    for row in DEFAULT_ROWS:
        result = analyze(row, args.batch)
        print(
            f"{result['items']},{result['partition_size']},"
            f"{result['cutoff']},{result['failure_log2']:.6f},"
            f"{result['occupancy_log2']:.6f},{result['tv_log2']:.6f},"
            f"{result['batch_tv_log2']:.6f},"
            f"{result['expected_trials_bound']:.6f}"
        )


if __name__ == "__main__":
    main()
