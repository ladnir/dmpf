"""Blockwise inner-list search for the leakage-weighted Ring-LPN decoder.

The complete four-list decoder enumerates all L^t supports of every hidden
polynomial.  This experiment replaces that step with a bounded beam.  Prefixes
are ranked by the monotone pair-collision upper bound, and the completed beam
is optionally reranked with the exact Reverse-Cuckoo matching likelihood.

The primary statistic is planted-support retention at a fixed final list size
B.  A uniform list retains one polynomial with probability B/L^t.  For the
work comparison we increase the uniform list width until it has the same
retention, then compare its quadratic outer merge cost with the charged beam
cost.  The fixed-B four-polynomial probability ratio is reported separately;
it is not treated as a restart exponent because the leakage instance is fixed.
This is an inner-list diagnostic, not a production Ring-LPN decoder.
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

from rev_cuckoo_regular_isd import fit_pair_slope
from rev_cuckoo_regular_noise import LeakageChannel, support_lists_one_unknown
from rev_cuckoo_ring_lpn_list_attack import (
    descriptors_for_hidden,
    pair_beam_candidates,
    partial_collision_factors,
)


Candidate = tuple[int, ...]


@dataclass(frozen=True)
class RetentionRow:
    trial: int
    hidden: int
    selector: str
    retained: float
    list_size: int
    candidate_space: int
    generated_children: int
    exact_factor_evaluations: int


def candidate_from_index(index: int, weight: int, block_length: int) -> Candidate:
    values = [0] * weight
    for variable in range(weight - 1, -1, -1):
        index, values[variable] = divmod(index, block_length)
    return tuple(values)


def exact_score(
    candidate: Candidate,
    known: tuple[int, ...],
    hidden_descriptors,
    polys: int,
    weight: int,
    channel: LeakageChannel,
) -> float:
    return channel.candidate_log_weight(
        support_lists_one_unknown(known, candidate, polys, weight),
        hidden_descriptors,
    )


def fractional_top_retention(
    scored_candidates: list[tuple[float, Candidate]],
    true_candidate: Candidate,
    list_size: int,
) -> float:
    """Expected inclusion under uniform tie breaking at the score cutoff."""

    true_scores = [
        score
        for score, candidate in scored_candidates
        if candidate == true_candidate
    ]
    if not true_scores:
        return 0.0
    true_score = true_scores[0]
    above = sum(score > true_score for score, _ in scored_candidates)
    tied = sum(score == true_score for score, _ in scored_candidates)
    if above >= list_size:
        return 0.0
    return min(1.0, (list_size - above) / tied)


def one_trial(
    trial: int,
    polys: int,
    weight: int,
    block_length: int,
    channel: LeakageChannel,
    pair_slope: float,
    list_size: int,
    oversample: int,
    beam_width: int,
    oracle_limit: int,
    seed: int,
) -> tuple[RetentionRow, ...]:
    rng = random.Random(seed)
    candidate_space = block_length**weight
    known = tuple(rng.randrange(block_length) for _ in range(polys * weight))
    true_candidates = tuple(
        tuple(rng.randrange(block_length) for _ in range(weight))
        for _ in range(polys)
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

    rows = []
    generated_count = min(candidate_space, max(list_size, list_size * oversample))
    for hidden, true_candidate in enumerate(true_candidates):
        hidden_descriptors = descriptors_for_hidden(
            descriptors, hidden, polys, weight
        )
        factors = partial_collision_factors(
            known,
            hidden_descriptors,
            polys,
            weight,
            block_length,
        )
        beam = pair_beam_candidates(
            factors,
            weight,
            block_length,
            pair_slope,
            max(beam_width, generated_count),
            generated_count,
            random.Random(seed ^ (hidden << 32) ^ 0x4245414D),
        )
        pair_scored = [
            (entry.pair_score, entry.offsets) for entry in beam.candidates
        ]
        rows.append(
            RetentionRow(
                trial,
                hidden,
                "pair_beam",
                fractional_top_retention(
                    pair_scored, true_candidate, list_size
                ),
                min(list_size, len(pair_scored)),
                candidate_space,
                beam.generated_children,
                0,
            )
        )

        scored = [
            (
                exact_score(
                    entry.offsets,
                    known,
                    hidden_descriptors,
                    polys,
                    weight,
                    channel,
                ),
                entry.offsets,
            )
            for entry in beam.candidates
        ]
        scored.sort(key=lambda item: item[0], reverse=True)
        rows.append(
            RetentionRow(
                trial,
                hidden,
                "exact_rerank",
                fractional_top_retention(
                    scored, true_candidate, list_size
                ),
                min(list_size, len(scored)),
                candidate_space,
                beam.generated_children,
                len(beam.candidates) * polys * weight,
            )
        )

        if candidate_space <= oracle_limit:
            oracle = []
            for index in range(candidate_space):
                candidate = candidate_from_index(index, weight, block_length)
                oracle.append(
                    (
                        exact_score(
                            candidate,
                            known,
                            hidden_descriptors,
                            polys,
                            weight,
                            channel,
                        ),
                        candidate,
                    )
                )
            oracle.sort(key=lambda item: item[0], reverse=True)
            rows.append(
                RetentionRow(
                    trial,
                    hidden,
                    "oracle_exact",
                    fractional_top_retention(
                        oracle, true_candidate, list_size
                    ),
                    min(list_size, len(oracle)),
                    candidate_space,
                    candidate_space,
                    candidate_space * polys * weight,
                )
            )
    return tuple(rows)


def equal_success_gain_bits(
    retention: float,
    candidate_space: int,
    polys: int,
    selected_attempt_work: float,
) -> float:
    if retention <= 0.0:
        return -math.inf
    equivalent_uniform_width = retention * candidate_space
    equivalent_uniform_work = (
        polys * equivalent_uniform_width
        + 2 * equivalent_uniform_width**2
    )
    return math.log2(equivalent_uniform_work / selected_attempt_work)


def coordinate_retentions(
    rows: list[RetentionRow],
    polys: int,
    selected_trials: list[int] | None = None,
) -> tuple[float, ...]:
    if selected_trials is None:
        selected_trials = sorted({row.trial for row in rows})
    by_coordinate = {
        (row.trial, row.hidden): row.retained for row in rows
    }
    return tuple(
        statistics.fmean(
            by_coordinate[trial, hidden] for trial in selected_trials
        )
        for hidden in range(polys)
    )


def geometric_retention(retentions: tuple[float, ...]) -> float:
    if any(retention <= 0.0 for retention in retentions):
        return 0.0
    return math.exp(statistics.fmean(math.log(value) for value in retentions))


def bootstrap_interval(
    rows: list[RetentionRow],
    trials: int,
    candidate_space: int,
    polys: int,
    selected_attempt_work: float,
    samples: int,
    rng: random.Random,
) -> tuple[float, float]:
    gains = []
    for _ in range(samples):
        selected_trials = [rng.randrange(trials) for _ in range(trials)]
        retained = geometric_retention(
            coordinate_retentions(rows, polys, selected_trials)
        )
        gains.append(
            equal_success_gain_bits(
                retained,
                candidate_space,
                polys,
                selected_attempt_work,
            )
        )
    gains.sort()
    return (
        gains[math.floor(0.025 * (samples - 1))],
        gains[math.ceil(0.975 * (samples - 1))],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=6)
    parser.add_argument("--block-length", type=int, default=6)
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--list-size", type=int, default=1024)
    parser.add_argument("--oversample", type=int, default=4)
    parser.add_argument("--beam-width", type=int, default=4096)
    parser.add_argument("--trials", type=int, default=64)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--bootstrap-samples", type=int, default=5000)
    parser.add_argument("--oracle-limit", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys != 4:
        parser.error("the current outer decoder model requires P=4")
    if min(
        args.weight,
        args.block_length,
        args.selector_k,
        args.list_size,
        args.oversample,
        args.beam_width,
        args.trials,
        args.bootstrap_samples,
    ) < 1:
        parser.error("all dimensions, budgets, and sample counts must be positive")
    candidate_space = args.block_length**args.weight
    if args.list_size > candidate_space:
        parser.error("list size cannot exceed L^t")
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
    start = time.perf_counter()
    rows = [
        row
        for trial in range(args.trials)
        for row in one_trial(
            trial,
            args.polys,
            args.weight,
            args.block_length,
            channel,
            pair_slope,
            args.list_size,
            args.oversample,
            args.beam_width,
            args.oracle_limit,
            args.seed + trial,
        )
    ]
    seconds = time.perf_counter() - start
    uniform_retention = args.list_size / candidate_space
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} K={args.selector_k} "
        f"candidate_space={candidate_space} list_size={args.list_size} "
        f"oversample={args.oversample} beam_width={args.beam_width} "
        f"uniform_retention={uniform_retention} pair_slope={pair_slope} "
        f"pair_r_squared={pair_r_squared} trials={args.trials} "
        f"seconds={seconds} seed={args.seed}"
    )
    for selector in ("pair_beam", "exact_rerank", "oracle_exact"):
        selected = [row for row in rows if row.selector == selector]
        if not selected:
            continue
        coordinate_rates = coordinate_retentions(selected, args.polys)
        retention = geometric_retention(coordinate_rates)
        generated_children = statistics.fmean(
            row.generated_children for row in selected
        )
        exact_factors = statistics.fmean(
            row.exact_factor_evaluations for row in selected
        )
        selected_attempt_work = (
            args.polys * (generated_children + exact_factors)
            + 2 * args.list_size**2
        )
        raw_gain = (
            -math.inf
            if retention == 0.0
            else args.polys * math.log2(retention / uniform_retention)
        )
        equal_gain = equal_success_gain_bits(
            retention,
            candidate_space,
            args.polys,
            selected_attempt_work,
        )
        interval = bootstrap_interval(
            selected,
            args.trials,
            candidate_space,
            args.polys,
            selected_attempt_work,
            args.bootstrap_samples,
            random.Random(args.seed ^ len(selector) ^ 0x424F4F54),
        )
        print(
            f"selector={selector} retained_mass={sum(row.retained for row in selected)}"
            f"/{len(selected)} coordinate_retentions={coordinate_rates} "
            f"equivalent_uniform_retention={retention} "
            f"raw_four_poly_gain_bits={raw_gain} "
            f"generated_children_per_poly={generated_children} "
            f"exact_factor_evaluations_per_poly={exact_factors} "
            f"attempt_work={selected_attempt_work} "
            f"equivalent_uniform_list_size={retention * candidate_space} "
            f"equal_success_charged_gain_bits={equal_gain} "
            f"equal_success_charged_gain_ci95=[{interval[0]},{interval[1]}]"
        )
    if args.output:
        fields = tuple(RetentionRow.__dataclass_fields__)
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
