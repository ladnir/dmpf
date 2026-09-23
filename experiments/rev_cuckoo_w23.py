#!/usr/bin/env python3
"""Empirical/extrapolated reverse-cuckoo parameters for w=2 and w=3.

The observable tail is fit as lambda = a log2(d) + b log2(n) + c.  That fit
is diagnostic only: its slope is too steep to extrapolate all the way to a
40-bit failure target.  The deep tail is instead anchored by the smallest
possible obstruction, w+1 items with all w choices equal:

    p_fail ~ binom(n, w + 1) / d ** (w * w).

The experiments show the measured ratio converging toward this term as d
grows.  Parameter estimates below invert this asymptotic expression while
treating d as continuous; integer rounding is displayed separately.

For w=3 the final selector uses a two-component model.  The threshold term is
fit as 2^-(a_n * e + b_n), while the local term explicitly includes the
four-, five-, and six-item Hall obstructions.  The two probabilities are added
before the target is inverted.
"""

import argparse
import csv
import math
from pathlib import Path


def solve3(matrix, vector):
    aug = [list(row) + [rhs] for row, rhs in zip(matrix, vector)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda row: abs(aug[row][col]))
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [value / scale for value in aug[col]]
        for row in range(3):
            if row == col:
                continue
            scale = aug[row][col]
            aug[row] = [x - scale * y for x, y in zip(aug[row], aug[col])]
    return [aug[row][3] for row in range(3)]


def fit(rows):
    xs = [(math.log2(row["d"]), math.log2(row["n"]), 1.0) for row in rows]
    ys = [row["lambda"] for row in rows]
    xtx = [[sum(x[i] * x[j] for x in xs) for j in range(3)] for i in range(3)]
    xty = [sum(x[i] * y for x, y in zip(xs, ys)) for i in range(3)]
    beta = solve3(xtx, xty)
    residuals = [y - sum(b * x for b, x in zip(beta, x)) for x, y in zip(xs, ys)]
    rmse = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    return beta, rmse, max(abs(value) for value in residuals)


def fit_line(points):
    mean_x = sum(x for x, _ in points) / len(points)
    mean_y = sum(y for _, y in points) / len(points)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points)
    slope /= sum((x - mean_x) ** 2 for x, _ in points)
    intercept = mean_y - slope * mean_x
    residuals = [y - (slope * x + intercept) for x, y in points]
    rmse = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    return slope, intercept, rmse, max(abs(value) for value in residuals)


def local_terms(n, d):
    # Four items have the same three partitioned choices.
    four = math.comb(n, 4) / d**9

    # Five items collectively use two bins in one partition and one bin in
    # each other partition. The factor 30 counts onto maps to the two bins.
    five = 0.0
    if n >= 5:
        five = 45 * math.comb(n, 5) * (d - 1) / d**12

    # Six items collectively use either (3,1,1) or (2,2,1) bins across the
    # partitions. These events can overlap smaller obstructions, but retaining
    # them accurately models the measured crossover and is conservative there.
    six = 0.0
    if n >= 6:
        configurations = (
            3 * math.comb(d, 3) * d * d * 540
            + 3 * math.comb(d, 2) ** 2 * d * 3844
        )
        six = math.comb(n, 6) * configurations / d**18
    return four, five, six


def local_probability(n, d):
    return sum(local_terms(n, d))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).with_name("rev_cuckoo_w23.csv"),
    )
    parser.add_argument("--target", type=float, default=40.0)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--min-failures", type=int, default=20)
    parser.add_argument("--min-lambda", type=float, default=8.0)
    args = parser.parse_args()

    groups = {2: [], 3: []}
    with args.csv.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            row = {key: int(value) for key, value in raw.items()}
            if row["failures"]:
                row["lambda"] = -math.log2(row["failures"] / row["trials"])
            groups[row["w"]].append(row)

    aggregate3 = {}
    for row in groups[3]:
        key = (row["n"], row["d"])
        trials, failures = aggregate3.get(key, (0, 0))
        aggregate3[key] = (trials + row["trials"], failures + row["failures"])

    print("observable-tail fit: lambda = a log2(d) + b log2(n) + c")
    print("w,points,a,b,c,rmse,max_abs_residual,asymptotic_d_slope")
    for w, all_rows in groups.items():
        rows = [
            row
            for row in all_rows
            if row["failures"] >= args.min_failures
            and row["lambda"] >= args.min_lambda
        ]
        beta, rmse, max_residual = fit(rows)
        print(
            f"{w},{len(rows)},{beta[0]:.6f},{beta[1]:.6f},{beta[2]:.6f},"
            f"{rmse:.6f},{max_residual:.6f},{w * w}"
        )

    print("\nmeasured / leading-obstruction probability at n=16")
    print("w,d,failures,trials,lambda,ratio")
    for w, minimum_d in ((2, 40), (3, 8)):
        for row in groups[w]:
            if row["n"] != 16 or row["d"] < minimum_d or not row["failures"]:
                continue
            leading = math.comb(row["n"], w + 1) / row["d"] ** (w * w)
            observed = row["failures"] / row["trials"]
            print(
                f"{w},{row['d']},{row['failures']},{row['trials']},"
                f"{row['lambda']:.6f},{observed / leading:.6f}"
            )

    batch_target = args.target + math.log2(args.batch)
    print("\ncontinuous asymptotic d estimates")
    print("w,n,d_single,ceil_single,d_batch,ceil_batch")
    for w in (2, 3):
        exponent = w * w
        for n in (16, 64, 128):
            coefficient = math.comb(n, w + 1)
            single = (coefficient * 2.0**args.target) ** (1.0 / exponent)
            batch = (coefficient * 2.0**batch_target) ** (1.0 / exponent)
            print(
                f"{w},{n},{single:.6f},{math.ceil(single)},"
                f"{batch:.6f},{math.ceil(batch)}"
            )

    print("\nthree-choice local-obstruction comparison")
    print(
        "n,normal_bins_single,normal_e_single,partition_d_single,partition_e_single,"
        "normal_bins_batch,normal_e_batch,partition_d_batch,partition_e_batch"
    )
    for n in (16, 64, 128, 256):
        coefficient = math.comb(n, 4)

        def normal_bins(security):
            bins = n
            while coefficient / math.comb(bins, 3) ** 3 > 2.0**-security:
                bins += 1
            return bins

        normal_single = normal_bins(args.target)
        normal_batch = normal_bins(batch_target)
        partition_single = (coefficient * 2.0**args.target) ** (1.0 / 9.0)
        partition_batch = (coefficient * 2.0**batch_target) ** (1.0 / 9.0)
        print(
            f"{n},{normal_single},{normal_single / n:.6f},"
            f"{partition_single:.6f},{3 * partition_single / n:.6f},"
            f"{normal_batch},{normal_batch / n:.6f},"
            f"{partition_batch:.6f},{3 * partition_batch / n:.6f}"
        )

    print("\ntwo-component partitioned w=3 model")
    print("n,points,a_threshold,b_threshold,rmse,max_abs_residual")
    threshold_fits = {}
    for n in (16, 32, 64, 128):
        points = []
        for (row_n, d), (trials, failures) in aggregate3.items():
            if row_n != n or failures < args.min_failures:
                continue
            observed = failures / trials
            local = local_probability(n, d)
            # Keep the fit in the threshold-dominated region. The remaining
            # points are validation data for the crossover.
            if observed < 2.5 * local:
                continue
            points.append((3 * d / n, -math.log2(observed - local)))
        a, b, rmse, max_residual = fit_line(points)
        threshold_fits[n] = (a, b)
        print(f"{n},{len(points)},{a:.6f},{b:.6f},{rmse:.6f},{max_residual:.6f}")

    print("\nn=16 held-out crossover validation")
    print("d,trials,failures,observed_lambda,model_lambda,error_bits")
    a, b = threshold_fits[16]
    for d in (10, 11, 12, 13, 14):
        trials, failures = aggregate3[(16, d)]
        observed = failures / trials
        threshold = 2.0 ** -(a * (3 * d / 16) + b)
        model = threshold + local_probability(16, d)
        print(
            f"{d},{trials},{failures},{-math.log2(observed):.6f},"
            f"{-math.log2(model):.6f},{math.log2(model / observed):.6f}"
        )

    print("\nw=3 selected integer parameters")
    print("target,n,d,expansion,lambda_total,lambda_threshold,lambda_local")
    targets = (
        ("single", args.target),
        ("batch", batch_target),
        ("batch-plus-2", batch_target + 2),
    )
    for name, security in targets:
        for n in (16, 64, 128):
            a, b = threshold_fits[n]
            d = math.ceil(n / 3)
            while True:
                threshold = 2.0 ** -(a * (3 * d / n) + b)
                local = local_probability(n, d)
                if threshold + local <= 2.0**-security:
                    break
                d += 1
            print(
                f"{name},{n},{d},{3 * d / n:.6f},"
                f"{-math.log2(threshold + local):.6f},"
                f"{-math.log2(threshold):.6f},{-math.log2(local):.6f}"
            )


if __name__ == "__main__":
    main()
