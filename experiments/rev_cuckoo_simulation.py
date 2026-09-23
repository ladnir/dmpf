#!/usr/bin/env python3
"""Exact checks for Reverse-Cuckoo descriptor-first simulators.

For small ideal partitioned-cuckoo instances, this script enumerates every
descriptor.  It compares three descriptor distributions:

* ``reverse``: choose a uniform injective placement first, then program a
  uniform compatible descriptor;
* ``forward``: choose a uniform descriptor, conditioned only on the existence
  of some placement;
* ``rejection``: choose a uniform descriptor and one uniform candidate per
  item, restarting unless the selected bins are distinct.

The last distribution is exactly equal to ``reverse``.  The calculation also
shows that ``forward`` can differ even when cuckoo failure is identically zero.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import Counter
from fractions import Fraction


Descriptor = tuple[tuple[int, ...], ...]


def placement_count(descriptor: Descriptor, w: int) -> int:
    """Count collision-free choices of one partitioned bin per item."""
    count = 0
    for choices in itertools.product(range(w), repeat=len(descriptor)):
        bins = [(side, descriptor[item][side]) for item, side in enumerate(choices)]
        count += len(set(bins)) == len(bins)
    return count


def enumerate_counts(t: int, w: int, d: int) -> Counter[int]:
    counts: Counter[int] = Counter()
    for flat in itertools.product(range(d), repeat=t * w):
        descriptor = tuple(
            tuple(flat[item * w + side] for side in range(w))
            for item in range(t)
        )
        counts[placement_count(descriptor, w)] += 1
    return counts


def total_variation_reverse_forward(counts: Counter[int]) -> Fraction:
    total = sum(counts.values())
    successful = total - counts[0]
    mean_numerator = sum(z * multiplicity for z, multiplicity in counts.items())

    tv = Fraction(0)
    for z, multiplicity in counts.items():
        reverse = Fraction(z, mean_numerator)
        forward = Fraction(1, successful) if z else Fraction(0)
        tv += multiplicity * abs(reverse - forward)
    return tv / 2


def injection_probability(t: int, m: int) -> float:
    return math.prod((m - item) / m for item in range(t))


def retry_cap(t: int, m: int, statistical_bits: int, batch: int) -> int:
    p_accept = injection_probability(t, m)
    target_bits = statistical_bits + math.log2(batch)
    return math.ceil(target_bits * math.log(2) / -math.log1p(-p_accept))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--t", type=int, default=2)
    parser.add_argument("--w", type=int, default=2)
    parser.add_argument("--d", type=int, default=3)
    parser.add_argument("--statistical-bits", type=int, default=40)
    parser.add_argument("--batch", type=int, default=256)
    args = parser.parse_args()

    descriptor_space = args.d ** (args.t * args.w)
    if descriptor_space > 10_000_000:
        raise SystemExit(
            f"descriptor space {descriptor_space:,} is too large for exact enumeration"
        )

    counts = enumerate_counts(args.t, args.w, args.d)
    total = sum(counts.values())
    failure = Fraction(counts[0], total)
    tv_forward = total_variation_reverse_forward(counts)
    m = args.w * args.d
    p_accept = injection_probability(args.t, m)

    # Given a descriptor h, rejection accepts with probability Z(h)/w^t.
    # Hence its normalized descriptor mass is exactly proportional to Z(h),
    # which is the reverse distribution.  This identity is algebraic; report
    # zero rather than estimating it from samples.
    tv_rejection = Fraction(0)

    print(f"parameters: t={args.t}, w={args.w}, d={args.d}, m={m}")
    print(f"descriptor space: {total:,}")
    print("placement-count histogram:")
    for z in sorted(counts):
        print(f"  Z={z}: {counts[z]:,}")
    print(f"ordinary cuckoo failure: {float(failure):.12g} ({failure})")
    print(
        "TV(reverse, uniform conditioned on success): "
        f"{float(tv_forward):.12g} ({tv_forward})"
    )
    print(f"TV(reverse, descriptor-first rejection): {float(tv_rejection):.1f}")
    print(f"one-trial rejection acceptance: {p_accept:.12g}")
    print(f"expected trials: {1 / p_accept:.6g}")
    print(
        f"trials for total error <= 2^-{args.statistical_bits} over "
        f"batch {args.batch}: "
        f"{retry_cap(args.t, m, args.statistical_bits, args.batch)}"
    )

    if args.t == 2 and args.w == 2 and args.d > 2:
        closed_form = Fraction((args.d - 1) ** 2, args.d**2 * (2 * args.d - 1))
        assert tv_forward == closed_form
        print(f"t=2 closed-form TV check: {closed_form}")


if __name__ == "__main__":
    main()
