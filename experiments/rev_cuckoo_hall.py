#!/usr/bin/env python3
"""Hall-witness upper bounds for partitioned cuckoo hashing.

Each of n items independently chooses one bin in each of w partitions of
size d.  A placement fails exactly when some k items collectively see at
most k-1 bins.  The union bound below enumerates every such Hall witness.
The dynamic program computes the sum over partition occupancies without
enumerating w-tuples, and uses base-2 log space throughout.
"""

import argparse
import math


NEG_INF = -math.inf


def logadd2(a: float, b: float) -> float:
    if a == NEG_INF:
        return b
    if b == NEG_INF:
        return a
    if a < b:
        a, b = b, a
    return a + math.log2(1.0 + math.exp2(b - a))


def hall_witness_bound(n: int, w: int, d: int) -> tuple[float, list[float]]:
    """Return log2 of the full Hall-witness union bound and its k terms."""
    total = NEG_INF
    terms: list[float] = []
    log_d = math.log2(d)

    for k in range(w + 1, n + 1):
        # a[r] is the probability contribution from choosing a particular
        # size-r bin subset in one partition, summed over all such subsets.
        a = [NEG_INF] * k
        for r in range(1, min(d, k - 1) + 1):
            a[r] = (
                math.log2(math.comb(d, r))
                + k * (math.log2(r) - log_d)
            )

        # Convolve w partitions, retaining total subset size below k.
        dp = [NEG_INF] * k
        dp[0] = 0.0
        for _ in range(w):
            next_dp = [NEG_INF] * k
            for used, prefix in enumerate(dp):
                if prefix == NEG_INF:
                    continue
                for r in range(1, k - used):
                    next_dp[used + r] = logadd2(
                        next_dp[used + r], prefix + a[r]
                    )
            dp = next_dp

        inner = NEG_INF
        for used in range(w, k):
            inner = logadd2(inner, dp[used])
        term = math.log2(math.comb(n, k)) + inner
        terms.append(term)
        total = logadd2(total, term)

    return total, terms


def next_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def recommend(n: int, w: int, target_bits: float, batch: int) -> tuple[int, float]:
    # A union bound over the batch requires the per-set failure probability
    # to be at most 2^-target_bits / batch.
    required = target_bits + math.log2(batch)
    d = next_power_of_two(max(1, math.ceil(n / w)))
    while True:
        bound, _ = hall_witness_bound(n, w, d)
        if bound <= -required:
            return d, bound
        d *= 2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, nargs="+", default=[16, 64, 128])
    parser.add_argument("--w", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--target", type=float, default=40.0)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--d", type=int)
    args = parser.parse_args()

    if args.batch < 1:
        parser.error("--batch must be positive")
    for w in args.w:
        if w < 2:
            parser.error("--w must be at least two")
        for n in args.n:
            if n <= w:
                parser.error("each --n must exceed w")
            if args.d is None:
                d, bound = recommend(n, w, args.target, args.batch)
            else:
                d = args.d
                bound, _ = hall_witness_bound(n, w, d)
            batch_bound = min(0.0, bound + math.log2(args.batch))
            print(
                f"n={n} w={w} d={d} total_bins={w*d} "
                f"log2_per_set_bound={bound:.6f} "
                f"log2_batch_bound={batch_bound:.6f}"
            )


if __name__ == "__main__":
    main()
