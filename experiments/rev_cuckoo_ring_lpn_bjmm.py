"""Optimistic q-ary two-level representation-ISD cost estimator.

This generalizes the regular MMT/BJMM enumeration of Esser--Santini to a
regular error whose nonzero coefficients range over F_q^*.  It models one
complete syndrome and deliberately omits polynomial factors, elimination,
memory-access costs, and integral rounding of representation filters.

For a regular target of weight p, a representation x_1 + x_2 may split the
target coefficient in s active blocks and introduce cancelling pairs in u
inactive blocks.  Both addends have weight p_x exactly when

    2 p_x = p + s + 2 u.

The same count is applied recursively when x_i is represented as a sum of two
vectors of weight p_y.  The resulting estimate is attack-favorable; a result
above regular Prange is therefore a useful screen against implementing the
full decoder.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import math
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class RepresentationEstimate:
    algorithm: str
    p: int
    p_x: int
    p_y: int
    ell: int
    ell_x: float
    ell_y: float
    leaf_log2: float
    first_merge_log2: float
    second_merge_log2: float
    final_merge_log2: float
    iteration_log2: float
    free_final_work_log2: float
    inverse_success_log2: float
    expected_work_log2: float
    memory_log2: float


def log2_binomial(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
    ) / math.log(2.0)


def log2_add(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    maximum = max(left, right)
    return maximum + math.log2(
        2.0 ** (left - maximum) + 2.0 ** (right - maximum)
    )


def representation_log2(
    target_weight: int,
    addend_weight: int,
    blocks: int,
    block_width: float,
    modulus: int,
) -> float:
    """Log2 number of ordered equal-weight regular representations."""

    if not 0 <= target_weight <= blocks:
        return -math.inf
    if not 0 <= addend_weight <= blocks:
        return -math.inf
    if block_width <= 0.0 or modulus < 2:
        return -math.inf

    total = -math.inf
    for split_active in range(target_weight + 1):
        if split_active and modulus == 2:
            continue
        cancellation_numerator = (
            2 * addend_weight - target_weight - split_active
        )
        if cancellation_numerator < 0 or cancellation_numerator % 2:
            continue
        cancel_inactive = cancellation_numerator // 2
        if cancel_inactive > blocks - target_weight:
            continue
        unsplit_active = target_weight - split_active
        if unsplit_active % 2:
            continue

        term = (
            log2_binomial(target_weight, split_active)
            + log2_binomial(
                unsplit_active, unsplit_active // 2
            )
            + log2_binomial(
                blocks - target_weight, cancel_inactive
            )
            + cancel_inactive
            * (math.log2(block_width) + math.log2(modulus - 1))
        )
        if split_active:
            term += split_active * math.log2(modulus - 2)
        total = log2_add(total, term)
    return total


def support_success_log2(
    blocks: int,
    block_length: int,
    information_width: float,
    p: int,
) -> float:
    fraction = information_width / block_length
    if not 0.0 < fraction < 1.0:
        return -math.inf
    return (
        log2_binomial(blocks, p)
        + p * math.log2(fraction)
        + (blocks - p) * math.log2(1.0 - fraction)
    )


def estimate_candidate(
    *,
    polys: int,
    weight: int,
    block_length: int,
    modulus: int,
    p: int,
    p_x: int,
    p_y: int,
    ell: int,
) -> RepresentationEstimate | None:
    """Return the optimistic two-level representation cost, if admissible."""

    blocks = polys * weight
    if blocks % 2 or p_y % 2:
        return None
    information_width = (
        (polys - 1) * weight * block_length + ell
    ) / blocks
    if information_width >= block_length:
        return None
    log2_rx = representation_log2(
        p, p_x, blocks, information_width, modulus
    )
    log2_ry = representation_log2(
        p_x, p_y, blocks, information_width, modulus
    )
    return estimate_from_representation_counts(
        blocks=blocks,
        block_length=block_length,
        modulus=modulus,
        information_width=information_width,
        p=p,
        p_x=p_x,
        p_y=p_y,
        ell=ell,
        log2_rx=log2_rx,
        log2_ry=log2_ry,
    )


def estimate_from_representation_counts(
    *,
    blocks: int,
    block_length: int,
    modulus: int,
    information_width: float,
    p: int,
    p_x: int,
    p_y: int,
    ell: int,
    log2_rx: float,
    log2_ry: float,
) -> RepresentationEstimate | None:
    if log2_rx == -math.inf or log2_ry == -math.inf:
        return None
    if blocks % 2 or p_y % 2:
        return None
    log2_q = math.log2(modulus)
    ell_x = log2_rx / log2_q
    ell_y = log2_ry / log2_q
    tolerance = 1e-12
    if ell_y > ell_x + tolerance or ell_x > ell + tolerance:
        return None

    coordinate_log2 = (
        math.log2(information_width) + math.log2(modulus - 1)
    )
    leaf_log2 = (
        log2_binomial(blocks // 2, p_y // 2)
        + (p_y // 2) * coordinate_log2
    )
    first_merge_log2 = (
        log2_binomial(blocks, p_y)
        + p_y * coordinate_log2
        - ell_y * log2_q
    )
    second_merge_log2 = (
        2 * log2_binomial(blocks, p_y)
        + 2 * p_y * coordinate_log2
        - (ell_x + ell_y) * log2_q
    )
    final_merge_log2 = (
        2 * log2_binomial(blocks, p_x)
        + 2 * p_x * coordinate_log2
        - (ell + ell_x) * log2_q
    )
    iteration_log2 = max(
        0.0,
        leaf_log2,
        first_merge_log2,
        second_merge_log2,
        final_merge_log2,
    )
    success_log2 = support_success_log2(
        blocks, block_length, information_width, p
    )
    if success_log2 == -math.inf:
        return None
    final_x_list_log2 = (
        log2_binomial(blocks, p_x)
        + p_x * coordinate_log2
        - ell_x * log2_q
    )
    memory_log2 = max(
        0.0, leaf_log2, first_merge_log2, final_x_list_log2
    )
    inverse_success_log2 = -success_log2
    return RepresentationEstimate(
        algorithm="BJMM",
        p=p,
        p_x=p_x,
        p_y=p_y,
        ell=ell,
        ell_x=ell_x,
        ell_y=ell_y,
        leaf_log2=leaf_log2,
        first_merge_log2=first_merge_log2,
        second_merge_log2=second_merge_log2,
        final_merge_log2=final_merge_log2,
        iteration_log2=iteration_log2,
        free_final_work_log2=(
            max(
                0.0,
                leaf_log2,
                first_merge_log2,
                second_merge_log2,
            )
            + inverse_success_log2
        ),
        inverse_success_log2=inverse_success_log2,
        expected_work_log2=iteration_log2 + inverse_success_log2,
        memory_log2=memory_log2,
    )


def estimate_mmt_from_representation_count(
    *,
    blocks: int,
    block_length: int,
    modulus: int,
    information_width: float,
    p: int,
    p_x: int,
    ell: int,
    log2_rx: float,
) -> RepresentationEstimate | None:
    """Return the optimistic one-level regular MMT cost."""

    if log2_rx == -math.inf or blocks % 2 or p_x % 2:
        return None
    log2_q = math.log2(modulus)
    ell_x = log2_rx / log2_q
    if ell_x > ell + 1e-12:
        return None
    coordinate_log2 = (
        math.log2(information_width) + math.log2(modulus - 1)
    )
    leaf_log2 = (
        log2_binomial(blocks // 2, p_x // 2)
        + (p_x // 2) * coordinate_log2
    )
    x_list_log2 = (
        log2_binomial(blocks, p_x)
        + p_x * coordinate_log2
        - ell_x * log2_q
    )
    final_merge_log2 = (
        2 * log2_binomial(blocks, p_x)
        + 2 * p_x * coordinate_log2
        - (ell + ell_x) * log2_q
    )
    iteration_log2 = max(
        0.0, leaf_log2, x_list_log2, final_merge_log2
    )
    success_log2 = support_success_log2(
        blocks, block_length, information_width, p
    )
    if success_log2 == -math.inf:
        return None
    inverse_success_log2 = -success_log2
    return RepresentationEstimate(
        algorithm="MMT",
        p=p,
        p_x=p_x,
        p_y=-1,
        ell=ell,
        ell_x=ell_x,
        ell_y=0.0,
        leaf_log2=leaf_log2,
        first_merge_log2=x_list_log2,
        second_merge_log2=-math.inf,
        final_merge_log2=final_merge_log2,
        iteration_log2=iteration_log2,
        free_final_work_log2=(
            max(0.0, leaf_log2, x_list_log2)
            + inverse_success_log2
        ),
        inverse_success_log2=inverse_success_log2,
        expected_work_log2=iteration_log2 + inverse_success_log2,
        memory_log2=max(0.0, leaf_log2, x_list_log2),
    )


def search(
    *,
    polys: int,
    weight: int,
    block_length: int,
    modulus: int,
    ell_max: int,
    retain: int = 1000,
) -> list[RepresentationEstimate]:
    """Return the global leading rows and the optimum for every p."""

    if retain < 1:
        raise ValueError("retain must be positive")
    blocks = polys * weight
    leading = []
    best_by_algorithm_p = {}
    best_free_final = {}
    serial = 0

    def retain_row(row: RepresentationEstimate) -> None:
        nonlocal serial
        key = (row.algorithm, row.p)
        best = best_by_algorithm_p.get(key)
        if best is None or row.expected_work_log2 < best.expected_work_log2:
            best_by_algorithm_p[key] = row
        free_best = best_free_final.get(row.algorithm)
        if (
            free_best is None
            or row.free_final_work_log2 < free_best.free_final_work_log2
        ):
            best_free_final[row.algorithm] = row
        entry = (-row.expected_work_log2, serial, row)
        serial += 1
        if len(leading) < retain:
            heapq.heappush(leading, entry)
        elif row.expected_work_log2 < -leading[0][0]:
            heapq.heapreplace(leading, entry)

    for ell in range(ell_max + 1):
        information_width = (
            (polys - 1) * weight * block_length + ell
        ) / blocks
        if information_width >= block_length:
            break
        representation_logs = tuple(
            tuple(
                representation_log2(
                    target,
                    addend,
                    blocks,
                    information_width,
                    modulus,
                )
                for addend in range(blocks + 1)
            )
            for target in range(blocks + 1)
        )
        for p in range(blocks + 1):
            for p_x in range(blocks + 1):
                mmt = estimate_mmt_from_representation_count(
                    blocks=blocks,
                    block_length=block_length,
                    modulus=modulus,
                    information_width=information_width,
                    p=p,
                    p_x=p_x,
                    ell=ell,
                    log2_rx=representation_logs[p][p_x],
                )
                if mmt is not None:
                    retain_row(mmt)
                for p_y in range(0, blocks + 1, 2):
                    row = estimate_from_representation_counts(
                        blocks=blocks,
                        block_length=block_length,
                        modulus=modulus,
                        information_width=information_width,
                        p=p,
                        p_x=p_x,
                        p_y=p_y,
                        ell=ell,
                        log2_rx=representation_logs[p][p_x],
                        log2_ry=representation_logs[p_x][p_y],
                    )
                    if row is not None:
                        retain_row(row)
    retained = {
        (row.algorithm, row.p, row.p_x, row.p_y, row.ell): row
        for row in (
            *best_by_algorithm_p.values(),
            *best_free_final.values(),
            *(entry[2] for entry in leading),
        )
    }
    return sorted(
        retained.values(), key=lambda row: row.expected_work_log2
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=4096)
    parser.add_argument("--modulus", type=int, default=65537)
    parser.add_argument("--ell-max", type=int, default=128)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--output-limit", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys < 2 or args.block_length % args.polys:
        parser.error("P must divide L and be at least two")
    if args.weight < 1 or args.modulus < 2 or args.ell_max < 0:
        parser.error("weight and modulus must be positive")
    if args.top < 1 or args.output_limit < 1:
        parser.error("top and output-limit must be positive")

    rows = search(
        polys=args.polys,
        weight=args.weight,
        block_length=args.block_length,
        modulus=args.modulus,
        ell_max=args.ell_max,
        retain=max(args.top, args.output_limit),
    )
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} modulus={args.modulus} "
        f"ell_max={args.ell_max} retained={len(rows)}"
    )
    selected = [row for row in rows if row.p > 0]
    if selected:
        best = selected[0]
        print(
            f"best_nontrivial_algorithm={best.algorithm} "
            f"p={best.p} p_x={best.p_x} "
            f"p_y={best.p_y} ell={best.ell} "
            f"work={best.expected_work_log2:.6f}"
        )
        free_best = min(
            selected, key=lambda row: row.free_final_work_log2
        )
        print(
            f"best_with_free_final_algorithm={free_best.algorithm} "
            f"p={free_best.p} p_x={free_best.p_x} "
            f"p_y={free_best.p_y} ell={free_best.ell} "
            f"work={free_best.free_final_work_log2:.6f}"
        )
    else:
        selected = rows
    for rank, row in enumerate(selected[: args.top], 1):
        print(
            f"rank={rank} algorithm={row.algorithm} p={row.p} "
            f"p_x={row.p_x} p_y={row.p_y} "
            f"ell={row.ell} ell_x={row.ell_x:.6f} "
            f"ell_y={row.ell_y:.6f} leaf={row.leaf_log2:.6f} "
            f"merge1={row.first_merge_log2:.6f} "
            f"merge2={row.second_merge_log2:.6f} "
            f"merge3={row.final_merge_log2:.6f} "
            f"iteration={row.iteration_log2:.6f} "
            f"free_final_work={row.free_final_work_log2:.6f} "
            f"inverse_success={row.inverse_success_log2:.6f} "
            f"work={row.expected_work_log2:.6f} "
            f"memory={row.memory_log2:.6f}"
        )
    if args.output:
        best_by_algorithm_p = {}
        for row in rows:
            best_by_algorithm_p.setdefault((row.algorithm, row.p), row)
        output_rows = {
            (row.algorithm, row.p, row.p_x, row.p_y, row.ell): row
            for row in (
                *best_by_algorithm_p.values(),
                *selected[: args.output_limit],
            )
        }
        ordered_output = sorted(
            output_rows.values(),
            key=lambda row: row.expected_work_log2,
        )
        fields = tuple(RepresentationEstimate.__dataclass_fields__)
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(asdict(row) for row in ordered_output)
        print(f"wrote {args.output} rows={len(ordered_output)}")


if __name__ == "__main__":
    main()
