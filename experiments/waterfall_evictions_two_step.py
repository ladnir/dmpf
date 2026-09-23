#!/usr/bin/env python3
"""Exact failure probabilities for biased two-step Waterfall repair."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

from waterfall_evictions import Configuration, exact_two_step_biased_failure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument("--w", type=int, default=3)
    parser.add_argument("--d", type=int, nargs="+", default=[32, 40, 48, 52, 56, 64, 74, 128])
    parser.add_argument("--c", type=int, default=2)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    if args.t < 1 or args.w < 2 or args.c < 0 or args.batch < 1:
        raise SystemExit("require t>=1, w>=2, c>=0, and batch>=1")
    if any(d < 1 for d in args.d):
        raise SystemExit("every --d must be positive")

    rows = []
    for d in args.d:
        configuration = Configuration(args.t, args.w, d, args.c)
        probability = exact_two_step_biased_failure(configuration)
        log2_probability = -math.inf if probability == 0.0 else math.log2(probability)
        rows.append(
            {
                "t": args.t,
                "w": args.w,
                "d": d,
                "c": args.c,
                "output_columns": args.w * d + args.c,
                "failure_probability": probability,
                "failure_log2": log2_probability,
                "batch": args.batch,
                "batch_failure_log2": min(
                    0.0, log2_probability + math.log2(args.batch)
                ),
                "batch_security_bits": max(
                    0.0, -log2_probability - math.log2(args.batch)
                ),
                "power_of_two_d": d & (d - 1) == 0,
            }
        )

    stream = sys.stdout if args.csv is None else args.csv.open("w", newline="")
    try:
        writer = csv.DictWriter(
            stream,
            fieldnames=tuple(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if args.csv is not None:
            stream.close()


if __name__ == "__main__":
    main()
