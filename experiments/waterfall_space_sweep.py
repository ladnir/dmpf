"""Sweep proved Waterfall configurations under the OT cost model.

The sweep keeps three proof regimes separate:

* basic Waterfall: exact overflow probability;
* one repair: the exact bound from ``waterfall_shallow_proof``;
* complete reachability: Hall-witness union bound plus an exact cap-overflow
  tail.

Only power-of-two partition sizes are considered, matching the concrete hash
output representation used by ``waterfall_mpc_cost``.
"""

from __future__ import annotations

import argparse
import functools
import itertools
import math
from dataclasses import dataclass
from fractions import Fraction

from rev_cuckoo_hall import (
    NEG_INF,
    hall_occupancy_bound,
    hall_occupancy_bound_sizes,
    logadd2,
)
from waterfall_mpc_cost import Configuration, Model, configuration_cost
from waterfall_shallow_proof import exact_one_repair_bound, exact_transitions


@dataclass(frozen=True)
class Point:
    proof: str
    w: int
    d: int
    c: int
    roots: int
    per_descriptor_log2: float
    batch_bits: float
    ot_count: int
    payload_bits: int
    partition_sizes: tuple[int, ...] = ()

    @property
    def expansion(self) -> int:
        return self.w + self.c

    @property
    def columns(self) -> int:
        return sum(self.partition_sizes or (self.d,) * self.w) + self.c


def log2_fraction(value: Fraction) -> float:
    if value == 0:
        return NEG_INF
    return math.log2(value.numerator) - math.log2(value.denominator)


def overflow_distribution(t: int, w: int, d: int) -> tuple[Fraction, ...]:
    return overflow_distribution_sizes(t, (d,) * w)


@functools.cache
def occupancy_transitions(
    t: int, d: int
) -> tuple[tuple[tuple[int, Fraction], ...], ...]:
    """Cache the exact kernel reused throughout an asymmetric sweep."""

    return tuple(tuple(row) for row in exact_transitions(t, d))


def overflow_distribution_sizes(
    t: int, partition_sizes: tuple[int, ...]
) -> tuple[Fraction, ...]:
    if not partition_sizes or any(d <= 0 for d in partition_sizes):
        raise ValueError("partition sizes must be positive")
    distribution = [Fraction(0)] * (t + 1)
    distribution[t] = Fraction(1)
    for d in partition_sizes:
        transitions = occupancy_transitions(t, d)
        next_distribution = [Fraction(0)] * (t + 1)
        for live, mass in enumerate(distribution):
            if mass == 0:
                continue
            for next_live, probability in transitions[live]:
                next_distribution[next_live] += mass * probability
        distribution = next_distribution
    if sum(distribution) != 1:
        raise AssertionError("overflow distribution does not sum to one")
    return tuple(distribution)


def sweep_asymmetric_reachability(
    model: Model,
    w: int,
    max_d: int,
    max_columns: int | None = None,
) -> list[Point]:
    """Enumerate ordered power-of-two partition schedules for complete repair."""

    if w < 2:
        raise ValueError("w must be at least two")
    target_per_descriptor = -(model.correctness_bits + math.log2(model.batch))
    powers = [1 << exponent for exponent in range(1, max_d.bit_length())]
    powers = [d for d in powers if d <= max_d]
    points: list[Point] = []

    # Hall failure and MPC cost depend only on the multiset of sizes.  Compute
    # them once, then inspect the distinct orders only for the overflow tail.
    for multiset in itertools.combinations_with_replacement(powers, w):
        if max_columns is not None and sum(multiset) > max_columns:
            continue
        hall_log, _ = hall_occupancy_bound_sizes(model.t, multiset)
        if hall_log > target_per_descriptor:
            continue
        incumbent: Point | None = None
        for sizes in set(itertools.permutations(multiset)):
            distribution = overflow_distribution_sizes(model.t, sizes)
            for roots in range(1, model.t + 1):
                cap_log = log2_fraction(tail(distribution, roots + 1))
                failure_log = logadd2(hall_log, cap_log)
                if failure_log > target_per_descriptor:
                    continue
                cfg = Configuration(
                    "asymmetric complete reachability",
                    "reachability",
                    w,
                    0,
                    0,
                    bfs_roots=roots,
                    partition_sizes=sizes,
                )
                cost = configuration_cost(model, cfg)
                point = Point(
                    "asymmetric complete reachability",
                    w,
                    0,
                    0,
                    roots,
                    failure_log,
                    batch_bits(failure_log, model.batch),
                    cost.total_ot_count,
                    cost.total_payload_bits,
                    sizes,
                )
                score = (point.ot_count, point.payload_bits, -point.batch_bits)
                if incumbent is None or score < (
                    incumbent.ot_count,
                    incumbent.payload_bits,
                    -incumbent.batch_bits,
                ):
                    incumbent = point
                break
        if incumbent is not None:
            points.append(incumbent)
    return points


def tail(distribution: tuple[Fraction, ...], threshold: int) -> Fraction:
    return sum(distribution[threshold:], Fraction(0))


def batch_bits(per_descriptor_log2: float, batch: int) -> float:
    return -per_descriptor_log2 - math.log2(batch)


def sweep(model: Model, max_w: int, max_c: int, max_d: int) -> list[Point]:
    points: list[Point] = []
    target_per_descriptor = -(model.correctness_bits + math.log2(model.batch))
    powers = [1 << exponent for exponent in range(1, max_d.bit_length())]
    powers = [d for d in powers if d <= max_d]

    for w in range(2, max_w + 1):
        for d in powers:
            distribution = overflow_distribution(model.t, w, d)

            for c in range(max_c + 1):
                basic_probability = tail(distribution, c + 1)
                basic_log = log2_fraction(basic_probability)
                if basic_log <= target_per_descriptor:
                    cfg = Configuration("basic", "basic", w, d, c)
                    cost = configuration_cost(model, cfg)
                    points.append(
                        Point(
                            "basic", w, d, c, 0, basic_log,
                            batch_bits(basic_log, model.batch),
                            cost.total_ot_count,
                            cost.total_payload_bits,
                        )
                    )

                repair_probability, _, _ = exact_one_repair_bound(
                    model.t, w, d, c
                )
                repair_log = log2_fraction(repair_probability)
                if repair_log <= target_per_descriptor:
                    cfg = Configuration("one repair", "one-repair", w, d, c)
                    cost = configuration_cost(model, cfg)
                    points.append(
                        Point(
                            "one repair", w, d, c, 1, repair_log,
                            batch_bits(repair_log, model.batch),
                            cost.total_ot_count,
                            cost.total_payload_bits,
                        )
                    )

            for c in range(max_c + 1):
                hall_log, _ = hall_occupancy_bound(model.t, w, d, c)
                if hall_log > target_per_descriptor:
                    continue
                for roots in range(1, model.t - c + 1):
                    cap_log = log2_fraction(tail(distribution, c + roots + 1))
                    failure_log = logadd2(hall_log, cap_log)
                    if failure_log > target_per_descriptor:
                        continue
                    cfg = Configuration(
                        "complete reachability", "reachability", w, d, c,
                        bfs_roots=roots,
                    )
                    cost = configuration_cost(model, cfg)
                    points.append(
                        Point(
                            "complete reachability", w, d, c, roots, failure_log,
                            batch_bits(failure_log, model.batch),
                            cost.total_ot_count,
                            cost.total_payload_bits,
                        )
                    )
                    # More root slots increase cost and only strengthen an
                    # already sufficient bound.
                    break

    return points


def pareto(points: list[Point]) -> list[Point]:
    result = []
    for point in points:
        dominated = any(
            other.expansion <= point.expansion
            and other.ot_count <= point.ot_count
            and other.payload_bits <= point.payload_bits
            and (
                other.expansion < point.expansion
                or other.ot_count < point.ot_count
                or other.payload_bits < point.payload_bits
            )
            for other in points
        )
        if not dominated:
            result.append(point)
    return sorted(
        result,
        key=lambda point: (point.expansion, point.ot_count, point.payload_bits),
    )


def show(points: list[Point]) -> None:
    print("proof,w,sizes,c,roots,w+c,m,batch_bits,ot_count,payload_bits,payload_KiB")
    for point in points:
        sizes = point.partition_sizes or (point.d,) * point.w
        print(
            f"{point.proof},{point.w},{'x'.join(map(str, sizes))},{point.c},{point.roots},"
            f"{point.expansion},{point.columns},{point.batch_bits:.4f},"
            f"{point.ot_count},{point.payload_bits},"
            f"{point.payload_bits / 8192:.3f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-w", type=int, default=8)
    parser.add_argument("--max-c", type=int, default=8)
    parser.add_argument("--max-d", type=int, default=512)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--asymmetric-w",
        type=int,
        help="also enumerate ordered heterogeneous schedules of this width",
    )
    args = parser.parse_args()
    model = Model()
    points = sweep(model, args.max_w, args.max_c, args.max_d)
    if args.asymmetric_w is not None:
        points.extend(
            sweep_asymmetric_reachability(
                model,
                args.asymmetric_w,
                args.max_d,
            )
        )

    print("lowest OT count")
    show(
        sorted(
            points,
            key=lambda point: (
                point.ot_count,
                point.payload_bits,
                point.expansion,
            ),
        )[: args.top]
    )
    print("\nPareto frontier: OT count, payload bits, and w+c")
    show(pareto(points))


if __name__ == "__main__":
    main()
