"""Rank audit for leakage-guided regular Ring-LPN information sets.

The leakage selector is independent of the public Ring-LPN multipliers.  This
experiment samples its actual fixed-cardinality information sets and compares
the rank of the selected structured matrix against both ordinary uniform
regular Prange and the selector's beta-zero stratified distribution.
"""

from __future__ import annotations

import argparse
import csv
import random
import time
from pathlib import Path

import numpy as np

from rev_cuckoo_basic_ring_lpn_isd import (
    sample_balanced_stratified_subset,
)
from rev_cuckoo_basic_ring_lpn_isd_production import (
    active_product_lists,
    factor_descriptors,
    sample_poly,
)
from rev_cuckoo_full_lpn_field import FieldModel
from rev_cuckoo_regular_isd import (
    collision_factors,
    fit_pair_slope,
    pair_mean_field_distributions,
)
from rev_cuckoo_regular_noise import LeakageChannel
from rev_cuckoo_ring_lpn_sequential_isd import (
    sample_sequential_polynomial_selection,
)


Selection = tuple[frozenset[int], ...]


def sample_beliefs(
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    pair_rounds: int,
    rng: random.Random,
) -> tuple[
    tuple[tuple[tuple[float, ...], ...], ...],
    tuple[tuple, ...],
]:
    """Sample the complete product channel and return one belief table per polynomial."""

    known = tuple(
        sample_poly(model.weight, model.block_length, model.modulus, rng)
        for _ in range(model.polys)
    )
    hidden = tuple(
        sample_poly(model.weight, model.block_length, model.modulus, rng)
        for _ in range(model.polys)
    )
    descriptors = tuple(
        channel.sample(active, 2 * model.block_length, rng)
        for active in active_product_lists(model, known, hidden)
    )
    known_positions = tuple(
        position for state in known for position in state.positions
    )
    kept = model.block_length // model.polys
    result = []
    factor_sets = []
    for hidden_index in range(model.polys):
        hidden_descriptors = factor_descriptors(
            descriptors, hidden_index, model.polys, model.weight
        )
        factors = collision_factors(
            known_positions,
            hidden_descriptors,
            polys=model.polys,
            weight=model.weight,
            block_length=model.block_length,
        )
        factor_sets.append(factors)
        result.append(
            pair_mean_field_distributions(
                known_positions,
                hidden_descriptors,
                polys=model.polys,
                weight=model.weight,
                block_length=model.block_length,
                kept=kept,
                slope=pair_slope,
                restarts=0,
                rounds=pair_rounds,
                damping=0.5,
                rng=rng,
                factors=factors,
            )
        )
    return tuple(result), tuple(factor_sets)


def uniform_selection(
    polys: int,
    weight: int,
    block_length: int,
    kept: int,
    rng: random.Random,
) -> tuple[Selection, ...]:
    return tuple(
        tuple(
            frozenset(rng.sample(range(block_length), kept))
            for _ in range(weight)
        )
        for _ in range(polys)
    )


def stratified_selection(
    beliefs: tuple[tuple[tuple[float, ...], ...], ...],
    kept: int,
    beta: float,
    rng: random.Random,
) -> tuple[Selection, ...]:
    return tuple(
        tuple(
            sample_balanced_stratified_subset(
                distribution, kept, beta, rng
            )
            for distribution in polynomial
        )
        for polynomial in beliefs
    )


def structured_selected_matrix(
    selections: tuple[Selection, ...],
    multipliers: tuple[np.ndarray, ...],
    block_length: int,
    modulus: int,
) -> np.ndarray:
    """Construct ``[I | M_r1 | ...]`` restricted to selected columns."""

    polys = len(selections)
    weight = len(selections[0])
    width = weight * block_length
    kept_columns = sum(len(block) for polynomial in selections for block in polynomial)
    if kept_columns != width:
        raise ValueError("selected matrix is not square")
    matrix = np.zeros((width, width), dtype=np.uint64)
    rows = np.arange(width, dtype=np.int64)
    output_column = 0
    for polynomial, selection in enumerate(selections):
        columns = np.fromiter(
            (
                block * block_length + offset
                for block, offsets in enumerate(selection)
                for offset in sorted(offsets)
            ),
            dtype=np.int64,
        )
        count = len(columns)
        if polynomial == 0:
            matrix[columns, output_column + np.arange(count)] = 1
        else:
            multiplier = multipliers[polynomial - 1]
            indices = (rows[:, None] - columns[None, :]) % width
            values = multiplier[indices]
            wrapped = rows[:, None] < columns[None, :]
            values = np.where(wrapped, (modulus - values) % modulus, values)
            matrix[:, output_column : output_column + count] = values
        output_column += count
    return matrix


def matrix_rank_mod(matrix: np.ndarray, modulus: int) -> int:
    """Vectorized Gaussian rank over a small prime field."""

    matrix = matrix.copy()
    rows, columns = matrix.shape
    rank = 0
    for column in range(columns):
        candidates = np.flatnonzero(matrix[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            matrix[[rank, pivot]] = matrix[[pivot, rank]]
        inverse = pow(int(matrix[rank, column]), -1, modulus)
        matrix[rank, column:] = matrix[rank, column:] * inverse % modulus
        if rank + 1 < rows:
            factors = matrix[rank + 1 :, column].copy()
            active = np.flatnonzero(factors)
            if len(active):
                target = rank + 1 + active
                matrix[target, column:] = (
                    matrix[target, column:]
                    + modulus * modulus
                    - factors[active, None] * matrix[rank, column:]
                ) % modulus
        rank += 1
        if rank == rows:
            break
    return rank


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=16)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--product-modulus", type=int, default=65537)
    parser.add_argument("--rank-modulus", type=int, default=257)
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--trials", type=int, default=16)
    parser.add_argument("--samples-per-trial", type=int, default=2)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys < 2 or args.block_length % args.polys:
        parser.error("P must be at least two and divide the block length")
    if args.rank_modulus > 65537:
        parser.error("the vectorized audit requires a rank modulus at most 65537")
    partition_size = args.partition_size or args.weight
    model = FieldModel(
        polys=args.polys,
        weight=args.weight,
        block_length=args.block_length,
        modulus=args.product_modulus,
        poly_states=(),
        full_states=(),
        position_candidates=(),
    )
    channel = LeakageChannel.build(
        "occupancy", args.selector_k, args.weight, partition_size
    )
    pair_slope, pair_r_squared = fit_pair_slope(
        channel,
        args.weight,
        args.pair_fit_samples,
        random.Random(args.seed ^ 0xA5A5A5A5),
    )
    kept = args.block_length // args.polys
    width = args.weight * args.block_length
    rows = []
    start = time.perf_counter()
    for trial in range(args.trials):
        beliefs, factor_sets = sample_beliefs(
            model,
            channel,
            pair_slope,
            args.pair_rounds,
            random.Random(args.seed + trial),
        )
        for sample in range(args.samples_per_trial):
            sample_rng = random.Random(
                args.seed ^ (trial << 32) ^ sample ^ 0xC0FFEE
            )
            selections = {
                "uniform": uniform_selection(
                    args.polys,
                    args.weight,
                    args.block_length,
                    kept,
                    sample_rng,
                ),
                "stratified_beta_0": stratified_selection(
                    beliefs, kept, 0.0, sample_rng
                ),
                "guided": stratified_selection(
                    beliefs, kept, args.beta, sample_rng
                ),
                "sequential": tuple(
                    sample_sequential_polynomial_selection(
                        beliefs[polynomial],
                        factor_sets[polynomial],
                        kept,
                        args.beta,
                        pair_slope,
                        sample_rng,
                    )
                    for polynomial in range(args.polys)
                ),
            }
            multipliers = tuple(
                np.fromiter(
                    (
                        sample_rng.randrange(args.rank_modulus)
                        for _ in range(width)
                    ),
                    dtype=np.uint64,
                    count=width,
                )
                for _ in range(args.polys - 1)
            )
            for selector, selection in selections.items():
                rank = matrix_rank_mod(
                    structured_selected_matrix(
                        selection,
                        multipliers,
                        args.block_length,
                        args.rank_modulus,
                    ),
                    args.rank_modulus,
                )
                rows.append(
                    {
                        "trial": trial,
                        "sample": sample,
                        "selector": selector,
                        "rank": rank,
                        "deficiency": width - rank,
                        "full_rank": int(rank == width),
                    }
                )
    seconds = time.perf_counter() - start
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} width={width} kept={kept} "
        f"K={args.selector_k} beta={args.beta} "
        f"rank_modulus={args.rank_modulus} trials={args.trials} "
        f"samples_per_trial={args.samples_per_trial} "
        f"pair_slope={pair_slope} pair_r_squared={pair_r_squared} "
        f"seconds={seconds}"
    )
    for selector in (
        "uniform",
        "stratified_beta_0",
        "guided",
        "sequential",
    ):
        selected = [row for row in rows if row["selector"] == selector]
        full = sum(row["full_rank"] for row in selected)
        deficiency = sum(row["deficiency"] for row in selected)
        print(
            f"selector={selector} full_rank={full}/{len(selected)} "
            f"failure_rate={(len(selected) - full) / len(selected)} "
            f"total_deficiency={deficiency}"
        )
    if args.output:
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
