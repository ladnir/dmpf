#!/usr/bin/env python3
"""Exact proof bound for a one-repair Waterfall router.

The proved router has stash capacity ``c``.  It succeeds immediately when
``L_w <= c``.  When ``L_w = c+1``, it applies depth-two BFS to one fixed
overflow row and succeeds iff that row has an augmenting path.  All
``L_w >= c+2`` outcomes are conservatively counted as failures.

For a live-count path ``L_0,...,L_w``, put ``J_s=L_s-L_{s+1}``.  The fixed
overflow row points to one occupied bin in every partition.  Its victim in
partition ``p`` has conditioned-occupied candidates in partitions at most
``p`` and independent, unused candidates in later partitions.  Therefore its
depth-two BFS fails with exact path-conditional probability

    product_{q=1}^{w-1} (J_q/d)^q.

All arithmetic is rational; floating point is used only to display bit counts.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction


def stirling_second_kind(limit: int) -> list[list[int]]:
    values = [[0] * (limit + 1) for _ in range(limit + 1)]
    values[0][0] = 1
    for rows in range(1, limit + 1):
        for blocks in range(1, rows + 1):
            values[rows][blocks] = (
                values[rows - 1][blocks - 1]
                + blocks * values[rows - 1][blocks]
            )
    return values


def exact_transitions(t: int, d: int) -> list[list[tuple[int, Fraction]]]:
    stirling = stirling_second_kind(t)
    transitions: list[list[tuple[int, Fraction]]] = [[] for _ in range(t + 1)]
    transitions[0].append((0, Fraction(1)))
    for live in range(1, t + 1):
        falling = 1
        denominator = d**live
        for occupied in range(1, min(live, d) + 1):
            falling *= d - occupied + 1
            probability = Fraction(
                falling * stirling[live][occupied], denominator
            )
            transitions[live].append((live - occupied, probability))
        if sum(probability for _, probability in transitions[live]) != 1:
            raise AssertionError("occupancy transitions do not sum to one")
    return transitions


def exact_one_repair_bound(
    t: int,
    w: int,
    d: int,
    stash: int,
) -> tuple[Fraction, Fraction, Fraction]:
    """Return (total bound, large-overflow tail, depth-two obstruction)."""
    transitions = exact_transitions(t, d)
    path = [t]
    large_overflow = Fraction(0)
    obstruction = Fraction(0)

    def recurse(partition: int, path_probability: Fraction) -> None:
        nonlocal large_overflow, obstruction
        if partition == w:
            terminal = path[-1]
            if terminal >= stash + 2:
                large_overflow += path_probability
                return
            if terminal != stash + 1:
                return
            occupied = [path[index] - path[index + 1] for index in range(w)]
            obstruction_probability = Fraction(1)
            for later_partition in range(1, w):
                obstruction_probability *= Fraction(
                    occupied[later_partition], d
                ) ** later_partition
            obstruction += path_probability * obstruction_probability
            return

        for next_live, transition_probability in transitions[path[-1]]:
            path.append(next_live)
            recurse(partition + 1, path_probability * transition_probability)
            path.pop()

    recurse(0, Fraction(1))
    return large_overflow + obstruction, large_overflow, obstruction


def bits(probability: Fraction) -> float:
    return math.inf if probability == 0 else -math.log2(float(probability))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument("--w", type=int, default=4)
    parser.add_argument("--d", type=int, default=16)
    parser.add_argument("--stash", type=int, default=2)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument(
        "--sweep-d",
        type=int,
        metavar="MAX_D",
        help="show every d through MAX_D meeting the requested batch bits",
    )
    parser.add_argument("--target-bits", type=float, default=40.0)
    args = parser.parse_args()
    if args.t < 1 or args.w < 2 or args.d < 1 or args.stash < 0:
        raise SystemExit("require t>=1, w>=2, d>=1, and stash>=0")

    total, tail, obstruction = exact_one_repair_bound(
        args.t, args.w, args.d, args.stash
    )
    print(f"parameters: t={args.t} w={args.w} d={args.d} c={args.stash}")
    print(f"Pr[L_w>=c+2]       = {float(tail):.16e} = 2^-{bits(tail):.9f}")
    print(
        "Pr[L_w=c+1 and depth2 fails] "
        f"= {float(obstruction):.16e} = 2^-{bits(obstruction):.9f}"
    )
    print(f"per-descriptor bound = {float(total):.16e} = 2^-{bits(total):.9f}")
    batch_probability = args.batch * total
    print(
        f"{args.batch}-descriptor union bound = {float(batch_probability):.16e} "
        f"= 2^-{bits(batch_probability):.9f}"
    )
    print(f"exact numerator   = {total.numerator}")
    print(f"exact denominator = {total.denominator}")

    if args.sweep_d is not None:
        print("\nqualifying sweep:")
        for d in range(args.d, args.sweep_d + 1):
            candidate, _, _ = exact_one_repair_bound(
                args.t, args.w, d, args.stash
            )
            batch_bits = bits(args.batch * candidate)
            if batch_bits >= args.target_bits:
                print(
                    f"d={d:3d} m={args.w*d+args.stash:4d} "
                    f"per_descriptor={bits(candidate):.6f} "
                    f"batch={batch_bits:.6f}"
                )


if __name__ == "__main__":
    main()
