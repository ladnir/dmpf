#!/usr/bin/env python3
"""Hall-witness upper bounds for partitioned cuckoo hashing.

Each of n items independently chooses one bin in every partition.  Partition
sizes may be uniform or heterogeneous.  A placement fails exactly when some
k items collectively see at most k-1 bins.  The union bound below enumerates
every such Hall witness.  The dynamic program computes the sum over partition
occupancies without enumerating size tuples, and uses base-2 log space.
"""

import argparse
import math
from fractions import Fraction


NEG_INF = -math.inf


def logadd2(a: float, b: float) -> float:
    if a == NEG_INF:
        return b
    if b == NEG_INF:
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + math.exp2(b - a))


def hall_witness_bound(
    n: int, w: int, d: int, stash: int = 0
) -> tuple[float, list[float]]:
    """Return a union bound for matching deficiency greater than ``stash``.

    Failure with a stash of size c means that some k rows collectively see at
    most k-c-1 candidate positions.  The ordinary perfect-matching bound is the
    special case c=0.
    """
    return hall_witness_bound_sizes(n, (d,) * w, stash)


def hall_witness_bound_sizes(
    n: int, partition_sizes: tuple[int, ...], stash: int = 0
) -> tuple[float, list[float]]:
    """Hall-witness union bound for heterogeneous partition sizes.

    Row ``i`` chooses one independent uniform bin in every partition ``s``,
    whose public size is ``partition_sizes[s]``.  The uniform-size helper
    above is retained as a compatibility wrapper.
    """
    if stash < 0:
        raise ValueError("stash must be nonnegative")
    if not partition_sizes or any(d <= 0 for d in partition_sizes):
        raise ValueError("partition sizes must be positive")
    w = len(partition_sizes)
    total = NEG_INF
    terms: list[float] = []

    for k in range(w + stash + 1, n + 1):
        limit = k - stash
        dp = [NEG_INF] * limit
        dp[0] = 0.0
        for d in partition_sizes:
            # a[r] is the probability contribution from choosing a
            # size-r subset in this partition, summed over all such subsets.
            a = [NEG_INF] * limit
            log_d = math.log2(d)
            for r in range(1, min(d, limit - 1) + 1):
                a[r] = (
                    math.log2(math.comb(d, r))
                    + k * (math.log2(r) - log_d)
                )

            # Convolve the next partition while retaining total size below k.
            next_dp = [NEG_INF] * limit
            for used, prefix in enumerate(dp):
                if prefix == NEG_INF:
                    continue
                for r in range(1, limit - used):
                    next_dp[used + r] = logadd2(
                        next_dp[used + r], prefix + a[r]
                    )
            dp = next_dp

        inner = NEG_INF
        for used in range(w, limit):
            inner = logadd2(inner, dp[used])
        term = math.log2(math.comb(n, k)) + inner
        terms.append(term)
        total = logadd2(total, term)

    return total, terms


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


def hall_occupancy_bound_sizes(
    n: int, partition_sizes: tuple[int, ...], stash: int = 0
) -> tuple[float, list[float]]:
    """Union-bound row sets, but count each set's occupancy event exactly."""

    if stash < 0:
        raise ValueError("stash must be nonnegative")
    if not partition_sizes or any(d <= 0 for d in partition_sizes):
        raise ValueError("partition sizes must be positive")
    w = len(partition_sizes)
    stirling = stirling_second_kind(n)
    total = NEG_INF
    terms: list[float] = []

    for k in range(w + stash + 1, n + 1):
        # A fixed k-row set violates the stash-c Hall condition when the sum
        # of its exact per-partition occupancy counts is below k-stash.
        limit = k - stash
        dp = [NEG_INF] * limit
        dp[0] = 0.0
        for d in partition_sizes:
            occupancy = [NEG_INF] * limit
            log_falling = 0.0
            for occupied in range(1, min(d, k, limit - 1) + 1):
                log_falling += math.log2(d - occupied + 1)
                occupancy[occupied] = (
                    log_falling
                    + math.log2(stirling[k][occupied])
                    - k * math.log2(d)
                )

            next_dp = [NEG_INF] * limit
            for used, prefix in enumerate(dp):
                if prefix == NEG_INF:
                    continue
                for occupied in range(1, limit - used):
                    next_dp[used + occupied] = logadd2(
                        next_dp[used + occupied],
                        prefix + occupancy[occupied],
                    )
            dp = next_dp

        inner = NEG_INF
        for used in range(w, limit):
            inner = logadd2(inner, dp[used])
        term = math.log2(math.comb(n, k)) + inner
        terms.append(term)
        total = logadd2(total, term)

    return total, terms


def hall_occupancy_bound(
    n: int, w: int, d: int, stash: int = 0
) -> tuple[float, list[float]]:
    """Uniform-partition wrapper for the exact fixed-set occupancy bound."""

    return hall_occupancy_bound_sizes(n, (d,) * w, stash)


def occupancy_threshold_probability(
    balls: int, bins: int, threshold: int
) -> Fraction:
    """Exact probability that some bin receives at least ``threshold`` balls."""

    if balls < 0 or bins <= 0 or threshold <= 0:
        raise ValueError("require balls>=0, bins>0, and threshold>0")
    if threshold > balls:
        return Fraction(0)

    counts = [0] * threshold
    safe_assignments = 0
    balls_factorial = math.factorial(balls)

    def enumerate_counts(size: int, remaining: int) -> None:
        nonlocal safe_assignments
        if size == threshold:
            if remaining != 0:
                return
            occupied_bins = sum(counts[1:])
            falling_bins = 1
            for offset in range(occupied_bins):
                falling_bins *= bins - offset
            denominator = 1
            for occupancy in range(1, threshold):
                denominator *= (
                    math.factorial(counts[occupancy])
                    * math.factorial(occupancy) ** counts[occupancy]
                )
            safe_assignments += (
                falling_bins * balls_factorial // denominator
            )
            return

        for count in range(remaining // size + 1):
            counts[size] = count
            enumerate_counts(size + 1, remaining - size * count)

    enumerate_counts(1, balls)
    return Fraction(bins**balls - safe_assignments, bins**balls)


def identical_candidate_obstruction_probability(
    rows: int, partition_sizes: tuple[int, ...]
) -> Fraction:
    """Probability that w+1 rows share one full candidate tuple.

    Such rows collectively touch only w main bins and therefore form a Hall
    witness even if the placement algorithm performs unlimited repairs.
    """

    if not partition_sizes or any(d <= 0 for d in partition_sizes):
        raise ValueError("partition sizes must be positive")
    return occupancy_threshold_probability(
        rows,
        math.prod(partition_sizes),
        len(partition_sizes) + 1,
    )


def next_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def recommend(
    n: int, w: int, target_bits: float, batch: int, stash: int = 0
) -> tuple[int, float]:
    # A union bound over the batch requires the per-set failure probability
    # to be at most 2^-target_bits / batch.
    required = target_bits + math.log2(batch)
    d = next_power_of_two(max(1, math.ceil(n / w)))
    while True:
        bound, _ = hall_witness_bound(n, w, d, stash)
        if bound <= -required:
            return d, bound
        d *= 2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, nargs="+", default=[16, 64, 128])
    parser.add_argument("--w", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--target", type=float, default=40.0)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--stash", type=int, default=0)
    parser.add_argument("--d", type=int)
    args = parser.parse_args()

    if args.batch < 1:
        parser.error("--batch must be positive")
    if args.stash < 0:
        parser.error("--stash must be nonnegative")
    for w in args.w:
        if w < 2:
            parser.error("--w must be at least two")
        for n in args.n:
            if n <= w:
                parser.error("each --n must exceed w")
            if args.d is None:
                d, bound = recommend(n, w, args.target, args.batch, args.stash)
            else:
                d = args.d
                bound, _ = hall_witness_bound(n, w, d, args.stash)
            batch_bound = min(0.0, bound + math.log2(args.batch))
            print(
                f"n={n} w={w} d={d} stash={args.stash} total_bins={w*d} "
                f"log2_per_set_bound={bound:.6f} "
                f"log2_batch_bound={batch_bound:.6f}"
            )


if __name__ == "__main__":
    main()
