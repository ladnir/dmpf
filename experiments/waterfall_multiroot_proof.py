#!/usr/bin/env python3
"""Exact bound for trying every overflow root when one repair suffices.

Rows are processed in their public order.  A state records, for each
partition, the number ``J_p`` of bins created so far and the number ``K_p`` of
distinct owners touched by eventual overflow rows.  A non-overflow row placed
in partition ``s`` collides in every ``p<s`` and creates a new bin in ``s``.
An overflow row collides in every partition and either revisits one of the
``K_p`` selected owners or selects a new owner.  This directly captures the
order bias that is lost if overflow rows are incorrectly treated as uniform
samples from the final occupied bins.

After all rows are processed, every distinct victim selected in partition
``p`` is depth-two blocked exactly when each of its unused candidates in
partitions ``q>p`` hits one of the final ``J_q`` occupied bins.  The dynamic
program and this final factor are exact over the rationals.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction

def exact_all_roots_depth2_obstruction(
    t: int,
    w: int,
    d: int,
    overflow: int,
) -> Fraction:
    """Return Pr[L_w=overflow and every overflow root is depth-two blocked]."""
    if not 1 <= overflow <= t:
        raise ValueError("overflow must lie in [1,t]")
    zero = (0,) * w
    # (created bins, selected victims, overflow rows) -> probability mass.
    states: dict[tuple[tuple[int, ...], tuple[int, ...], int], Fraction] = {
        (zero, zero, 0): Fraction(1)
    }

    for _ in range(t):
        next_states: dict[
            tuple[tuple[int, ...], tuple[int, ...], int], Fraction
        ] = {}
        for (created, selected, roots), mass in states.items():
            # The row is placed in its first non-colliding partition.
            collision_prefix = Fraction(1)
            for placement_partition in range(w):
                if created[placement_partition] < d:
                    probability = collision_prefix * Fraction(
                        d - created[placement_partition], d
                    )
                    updated = list(created)
                    updated[placement_partition] += 1
                    key = (tuple(updated), selected, roots)
                    next_states[key] = next_states.get(key, Fraction(0)) + (
                        mass * probability
                    )
                collision_prefix *= Fraction(created[placement_partition], d)

            # Or it collides everywhere and becomes an overflow root.  Roots
            # beyond the requested stratum can be dropped because their count
            # can never decrease.
            if roots >= overflow or collision_prefix == 0:
                continue
            branches: dict[tuple[int, ...], Fraction] = {selected: Fraction(1)}
            for partition in range(w):
                updated_branches: dict[tuple[int, ...], Fraction] = {}
                for selected_state, branch_mass in branches.items():
                    old = selected_state[partition]
                    if old:
                        key = selected_state
                        updated_branches[key] = updated_branches.get(
                            key, Fraction(0)
                        ) + branch_mass * Fraction(old, d)
                    new_choices = created[partition] - old
                    if new_choices:
                        updated = list(selected_state)
                        updated[partition] += 1
                        key = tuple(updated)
                        updated_branches[key] = updated_branches.get(
                            key, Fraction(0)
                        ) + branch_mass * Fraction(new_choices, d)
                branches = updated_branches
            for selected_state, branch_mass in branches.items():
                key = (created, selected_state, roots + 1)
                next_states[key] = next_states.get(key, Fraction(0)) + (
                    mass * branch_mass
                )
        states = next_states

    obstruction = Fraction(0)
    for (created, selected, roots), mass in states.items():
        if roots != overflow:
            continue
        blocked = Fraction(1)
        for victim_partition in range(w - 1):
            for later_partition in range(victim_partition + 1, w):
                blocked *= Fraction(created[later_partition], d) ** selected[
                    victim_partition
                ]
        obstruction += mass * blocked
    return obstruction


def bits(probability: Fraction) -> float:
    return math.inf if probability == 0 else -math.log2(float(probability))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument("--w", type=int, default=4)
    parser.add_argument("--d", type=int, default=8)
    parser.add_argument("--overflow", type=int, nargs="+", default=[3, 4, 5])
    args = parser.parse_args()
    for overflow in args.overflow:
        probability = exact_all_roots_depth2_obstruction(
            args.t, args.w, args.d, overflow
        )
        print(
            f"L={overflow}: Pr[L_w=L and all roots blocked] "
            f"= {float(probability):.16e} = 2^-{bits(probability):.9f}"
        )


if __name__ == "__main__":
    main()
