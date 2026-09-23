"""Syndrome-projected bucket ordering for the reduced four-list decoder.

This is the first global projection experiment.  It does not guess a planted
projection value.  Instead, it groups the two left-polynomial and two
right-polynomial sums by a public projection of the Ring-LPN syndrome.  The
attack computes the posterior mass of every compatible bucket from the four
per-polynomial leakage scores and searches buckets in decreasing mass.

All candidate and bucket construction is charged.  Consequently a projection
can help only by avoiding enough low-mass pair buckets to repay the posterior
mass convolution.  The experiment still enumerates each per-polynomial list;
its purpose is to determine whether a global projected merge is structurally
useful before combining it with a bounded inner generator.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import random
import statistics
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from rev_cuckoo_full_lpn import NEG_INF, logaddexp2
from rev_cuckoo_regular_isd import fit_pair_slope
from rev_cuckoo_regular_noise import LeakageChannel, support_lists_one_unknown
from rev_cuckoo_ring_lpn_list_attack import (
    Candidate,
    CandidateEntry,
    PairEntry,
    Syndrome,
    add_syndromes,
    build_candidate_lists,
    decode_lists,
    descriptors_for_hidden,
    polynomial_contribution,
    subtract_syndromes,
)


Projection = tuple[int, ...]


@dataclass(frozen=True)
class ProjectedDecodeResult:
    selector: str
    buckets_visited: int
    pair_entries_materialized: int
    mass_convolution_terms: int
    recovered_planted: bool


@dataclass(frozen=True)
class ProjectedTrialRow:
    trial: int
    selector: str
    projection_width: int
    buckets_visited: int
    candidate_scores: int
    leakage_factor_evaluations: int
    mass_convolution_terms: int
    pair_entries_materialized: int
    charged_work: int
    recovered_planted: int


def project(syndrome: Syndrome, width: int) -> Projection:
    return tuple(syndrome[:width])


def add_projection(
    left: Projection, right: Projection, modulus: int
) -> Projection:
    return tuple((a + b) % modulus for a, b in zip(left, right))


def subtract_projection(
    left: Projection, right: Projection, modulus: int
) -> Projection:
    return tuple((a - b) % modulus for a, b in zip(left, right))


def projection_buckets(
    entries: tuple[CandidateEntry, ...], width: int
) -> dict[Projection, tuple[CandidateEntry, ...]]:
    mutable: dict[Projection, list[CandidateEntry]] = {}
    for entry in entries:
        mutable.setdefault(project(entry.contribution, width), []).append(entry)
    return {key: tuple(values) for key, values in mutable.items()}


def projection_log_masses(
    buckets: dict[Projection, tuple[CandidateEntry, ...]],
) -> dict[Projection, float]:
    result = {}
    for key, entries in buckets.items():
        total = NEG_INF
        for entry in entries:
            total = logaddexp2(total, entry.exact_score)
        result[key] = total
    return result


def pair_projection_masses(
    left: dict[Projection, float],
    right: dict[Projection, float],
    modulus: int,
) -> tuple[dict[Projection, float], int]:
    result: dict[Projection, float] = {}
    terms = 0
    for left_key, left_mass in left.items():
        for right_key, right_mass in right.items():
            key = add_projection(left_key, right_key, modulus)
            result[key] = logaddexp2(
                result.get(key, NEG_INF), left_mass + right_mass
            )
            terms += 1
    return result, terms


def dense_projection_masses(
    buckets: dict[Projection, tuple[CandidateEntry, ...]],
    width: int,
    modulus: int,
) -> np.ndarray:
    """Return normalized posterior mass on F_q^width."""

    shape = (modulus,) * width
    result = np.zeros(shape, dtype=np.float64)
    maximum = max(
        entry.exact_score
        for entries in buckets.values()
        for entry in entries
        if math.isfinite(entry.exact_score)
    )
    for key, entries in buckets.items():
        result[key] = sum(
            2.0 ** (entry.exact_score - maximum)
            for entry in entries
            if math.isfinite(entry.exact_score)
        )
    return result


def cyclic_mass_convolution(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Circular convolution on the additive group F_q^r."""

    transformed = np.fft.fftn(left) * np.fft.fftn(right)
    result = np.fft.ifftn(transformed).real
    result[result < 0.0] = 0.0
    return result


def materialize_pair_bucket(
    left: dict[Projection, tuple[CandidateEntry, ...]],
    right: dict[Projection, tuple[CandidateEntry, ...]],
    target: Projection,
    selector: str,
    modulus: int,
    rng: random.Random,
) -> list[PairEntry]:
    scored = []
    for left_key, left_entries in left.items():
        right_entries = right.get(
            subtract_projection(target, left_key, modulus), ()
        )
        for left_entry in left_entries:
            for right_entry in right_entries:
                scored.append(
                    (
                        left_entry.exact_score + right_entry.exact_score,
                        rng.random(),
                        PairEntry(
                            add_syndromes(
                                left_entry.contribution,
                                right_entry.contribution,
                                modulus,
                            ),
                            left_entry.index,
                            right_entry.index,
                        ),
                    )
                )
    if selector == "exact_projected":
        scored.sort(key=lambda item: (item[0], item[1]), reverse=True)
    else:
        rng.shuffle(scored)
    return [entry for _, _, entry in scored]


def search_pair_buckets(
    candidate_lists: tuple[tuple[CandidateEntry, ...], ...],
    target: Syndrome,
    true_indices: tuple[int, ...],
    projection_width: int,
    selector: str,
    modulus: int,
    rng: random.Random,
) -> ProjectedDecodeResult:
    if selector not in ("uniform_projected", "exact_projected"):
        raise ValueError(f"unknown projected selector: {selector}")
    buckets = tuple(
        projection_buckets(entries, projection_width)
        for entries in candidate_lists
    )
    target_projection = project(target, projection_width)
    mass_convolution_terms = 0
    if selector == "exact_projected":
        masses = tuple(
            dense_projection_masses(
                values, projection_width, modulus
            )
            for values in buckets
        )
        left_mass = cyclic_mass_convolution(masses[0], masses[1])
        right_mass = cyclic_mass_convolution(masses[2], masses[3])
        group_size = modulus**projection_width
        mass_convolution_terms = 2 * math.ceil(
            group_size * math.log2(group_size)
        )
        order = [
            (
                left_mass[key]
                * right_mass[
                    subtract_projection(target_projection, key, modulus)
                ],
                rng.random(),
                key,
            )
            for key in itertools.product(
                range(modulus), repeat=projection_width
            )
            if left_mass[key] > 0.0
            and right_mass[
                subtract_projection(target_projection, key, modulus)
            ]
            > 0.0
        ]
        order.sort(key=lambda item: (item[0], item[1]), reverse=True)
        projection_order = [key for _, _, key in order]
    else:
        projection_order = list(
            itertools.product(range(modulus), repeat=projection_width)
        )
        rng.shuffle(projection_order)

    pair_entries_materialized = 0
    for bucket_index, left_projection in enumerate(projection_order, start=1):
        right_projection = subtract_projection(
            target_projection, left_projection, modulus
        )
        left_pairs = materialize_pair_bucket(
            buckets[0],
            buckets[1],
            left_projection,
            selector,
            modulus,
            rng,
        )
        right_pairs = materialize_pair_bucket(
            buckets[2],
            buckets[3],
            right_projection,
            selector,
            modulus,
            rng,
        )
        pair_entries_materialized += len(left_pairs) + len(right_pairs)
        left_table: dict[Syndrome, PairEntry] = {}
        right_table: dict[Syndrome, PairEntry] = {}
        for position in range(max(len(left_pairs), len(right_pairs))):
            if position < len(left_pairs):
                left = left_pairs[position]
                needed = subtract_syndromes(target, left.contribution, modulus)
                right = right_table.get(needed)
                if right is not None:
                    recovered = (
                        left.left_index,
                        left.right_index,
                        right.left_index,
                        right.right_index,
                    ) == true_indices
                    return ProjectedDecodeResult(
                        selector,
                        bucket_index,
                        pair_entries_materialized,
                        mass_convolution_terms,
                        recovered,
                    )
                left_table.setdefault(left.contribution, left)
            if position < len(right_pairs):
                right = right_pairs[position]
                needed = subtract_syndromes(target, right.contribution, modulus)
                left = left_table.get(needed)
                if left is not None:
                    recovered = (
                        left.left_index,
                        left.right_index,
                        right.left_index,
                        right.right_index,
                    ) == true_indices
                    return ProjectedDecodeResult(
                        selector,
                        bucket_index,
                        pair_entries_materialized,
                        mass_convolution_terms,
                        recovered,
                    )
                right_table.setdefault(right.contribution, right)
    raise RuntimeError(
        f"{selector} decoder exhausted all compatible projection buckets"
    )


def one_trial(
    trial: int,
    candidates: tuple[Candidate, ...],
    polys: int,
    weight: int,
    block_length: int,
    modulus: int,
    projection_width: int,
    channel: LeakageChannel,
    pair_slope: float,
    seed: int,
) -> tuple[ProjectedTrialRow, ...]:
    rng = random.Random(seed)
    known = tuple(rng.randrange(block_length) for _ in range(polys * weight))
    true_indices = tuple(rng.randrange(len(candidates)) for _ in range(polys))
    true_candidates = tuple(candidates[index] for index in true_indices)
    ring_width = weight * block_length
    multipliers = tuple(
        tuple(rng.randrange(modulus) for _ in range(ring_width))
        for _ in range(polys - 1)
    )
    lists_by_hidden = tuple(
        support_lists_one_unknown(known, hidden, polys, weight)
        for hidden in true_candidates
    )
    true_lists = tuple(
        lists_by_hidden[hidden][known_poly * weight + output]
        for known_poly in range(polys)
        for hidden in range(polys)
        for output in range(weight)
    )
    descriptors = tuple(
        channel.sample(active, 2 * block_length, rng) for active in true_lists
    )
    candidate_lists = build_candidate_lists(
        candidates,
        known,
        descriptors,
        multipliers,
        polys,
        weight,
        block_length,
        modulus,
        channel,
        pair_slope,
    )
    target = bytes([0] * ring_width)
    for polynomial, candidate in enumerate(true_candidates):
        target = add_syndromes(
            target,
            polynomial_contribution(
                candidate,
                polynomial,
                multipliers,
                block_length,
                modulus,
            ),
            modulus,
        )

    rows = []
    candidate_scores = polys * len(candidates)
    leakage_factors = candidate_scores * polys * weight
    for selector_index, selector in enumerate(
        ("uniform_projected", "exact_projected")
    ):
        decoded = search_pair_buckets(
            candidate_lists,
            target,
            true_indices,
            projection_width,
            selector,
            modulus,
            random.Random(seed ^ (selector_index << 40) ^ 0x50524F4A),
        )
        factors = leakage_factors if selector == "exact_projected" else 0
        rows.append(
            ProjectedTrialRow(
                trial,
                selector,
                projection_width,
                decoded.buckets_visited,
                candidate_scores,
                factors,
                decoded.mass_convolution_terms,
                decoded.pair_entries_materialized,
                candidate_scores
                + factors
                + decoded.mass_convolution_terms
                + decoded.pair_entries_materialized,
                int(decoded.recovered_planted),
            )
        )

    for selector_index, selector in enumerate(("uniform", "exact"), start=2):
        decoded = decode_lists(
            candidate_lists,
            target,
            true_indices,
            selector,
            modulus,
            random.Random(seed ^ (selector_index << 40) ^ 0x50524F4A),
        )
        factors = leakage_factors if selector == "exact" else 0
        rows.append(
            ProjectedTrialRow(
                trial,
                f"{selector}_flat",
                projection_width,
                0,
                candidate_scores,
                factors,
                0,
                decoded.merge_entries,
                candidate_scores + factors + decoded.merge_entries,
                int(decoded.recovered_planted),
            )
        )
    return tuple(rows)


def geometric_mean(values: list[int]) -> float:
    return 2.0 ** statistics.fmean(math.log2(value) for value in values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=4)
    parser.add_argument("--block-length", type=int, default=4)
    parser.add_argument("--modulus", type=int, default=17)
    parser.add_argument("--projection-width", type=int, default=1)
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--trials", type=int, default=64)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys != 4:
        parser.error("the projected decoder requires P=4")
    if not 2 <= args.modulus <= 256:
        parser.error("the byte-packed syndrome modulus must lie in [2,256]")
    if not 1 <= args.projection_width <= args.weight * args.block_length:
        parser.error("projection width must lie in [1,tL]")
    if min(args.weight, args.block_length, args.trials) < 1:
        parser.error("dimensions and trials must be positive")
    partition_size = args.partition_size or args.weight
    channel = LeakageChannel.build(
        "occupancy", args.selector_k, args.weight, partition_size
    )
    pair_slope, pair_r_squared = fit_pair_slope(
        channel,
        args.weight,
        args.pair_fit_samples,
        random.Random(args.seed ^ 0xA5A5A5A5),
    )
    candidates = tuple(
        itertools.product(range(args.block_length), repeat=args.weight)
    )
    start = time.perf_counter()
    rows = [
        row
        for trial in range(args.trials)
        for row in one_trial(
            trial,
            candidates,
            args.polys,
            args.weight,
            args.block_length,
            args.modulus,
            args.projection_width,
            channel,
            pair_slope,
            args.seed + trial,
        )
    ]
    seconds = time.perf_counter() - start
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} modulus={args.modulus} "
        f"projection_width={args.projection_width} K={args.selector_k} "
        f"candidates_per_poly={len(candidates)} pair_slope={pair_slope} "
        f"pair_r_squared={pair_r_squared} trials={args.trials} "
        f"seconds={seconds} seed={args.seed}"
    )
    by_selector = {
        selector: [row for row in rows if row.selector == selector]
        for selector in (
            "uniform_flat",
            "exact_flat",
            "uniform_projected",
            "exact_projected",
        )
    }
    baseline = statistics.fmean(
        row.charged_work for row in by_selector["uniform_flat"]
    )
    for selector, selected in by_selector.items():
        mean_work = statistics.fmean(row.charged_work for row in selected)
        print(
            f"selector={selector} charged_work_mean={mean_work} "
            f"charged_work_geomean="
            f"{geometric_mean([row.charged_work for row in selected])} "
            f"gain_vs_uniform_flat_bits={math.log2(baseline / mean_work)} "
            f"buckets_mean={statistics.fmean(row.buckets_visited for row in selected)} "
            f"pair_entries_mean="
            f"{statistics.fmean(row.pair_entries_materialized for row in selected)} "
            f"mass_terms_mean="
            f"{statistics.fmean(row.mass_convolution_terms for row in selected)} "
            f"planted_recovery={sum(row.recovered_planted for row in selected)}"
            f"/{len(selected)}"
        )
    if args.output:
        fields = tuple(ProjectedTrialRow.__dataclass_fields__)
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {field: getattr(row, field) for field in fields}
                )
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
