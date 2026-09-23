"""Matching-aware dual-check selectors after sparse-factor reduction.

This experiment asks whether complete Reverse-Cuckoo matching factors retain
useful structure after reducing the Ring-LPN errors modulo ``X^n + c``.  It
uses the existing pair-proposal SMC sampler, which targets the aggressively
leakage-weighted heuristic distribution ``q_pair * W_match``.  Two concrete selectors
are evaluated.

The marginal selector projects the final particles to residue-occupancy
scores and feeds those scores to the global product-weighted ``n-1``-zero
selector.  The planted-support selector preserves the complete support of a
particle: it plants one predicted reduced support per error polynomial and
fills the remaining zero coordinates uniformly.  Its success probability is
averaged exactly over the empirical particle populations and the random fill;
no rare support-hit Monte Carlo is used.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import random
import statistics
import time
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from rev_cuckoo_basic_ring_lpn_isd import BALANCED_BETAS
from rev_cuckoo_basic_ring_lpn_isd_production import factor_descriptors
from rev_cuckoo_full_lpn_field import FieldModel
from rev_cuckoo_regular_isd import fit_pair_slope
from rev_cuckoo_regular_noise import LeakageChannel
from rev_cuckoo_ring_lpn_factor_attack import (
    product_weighted_support_probability,
    projected_occupancy_scores,
)
from rev_cuckoo_ring_lpn_matching_isd import pair_proposal_target_particles
from rev_cuckoo_ring_lpn_sequential_isd import sample_trial_instance


Candidate = tuple[int, ...]


@dataclass(frozen=True)
class MatchingFactorTrial:
    trial: int
    reduced_weight: int
    uniform_completion_probability: float
    pair_completion_probability: float
    matching_completion_probability: float
    marginal_success_probabilities: tuple[tuple[float, ...], ...]
    minimum_smc_ess: float
    minimum_smc_unique: float
    smc_acceptance: float


def projected_candidate_support(
    candidate: Sequence[int], block_length: int, factor_degree: int
) -> frozenset[int]:
    """Project a complete regular-error candidate to reduced coordinates."""

    return frozenset(
        (block * block_length + offset) % factor_degree
        for block, offset in enumerate(candidate)
    )


def particle_occupancy_scores(
    candidates: Sequence[Candidate],
    block_length: int,
    factor_degree: int,
) -> tuple[float, ...]:
    """Empirical reduced-coordinate occupancy marginals of particles."""

    if not candidates:
        raise ValueError("particle population is empty")
    counts = [0] * factor_degree
    for candidate in candidates:
        for residue in projected_candidate_support(
            candidate, block_length, factor_degree
        ):
            counts[residue] += 1
    scale = 1.0 / len(candidates)
    return tuple(count * scale for count in counts)


def completion_probability(
    width: int,
    kept: int,
    planted_size: int,
    true_size: int,
    planted_hits: int,
) -> float:
    """Probability that uniform completion covers the remaining true support."""

    missing = true_size - planted_hits
    remaining_width = width - planted_size
    remaining_kept = kept - planted_size
    if missing < 0 or remaining_kept < missing or remaining_kept > remaining_width:
        return 0.0
    return math.comb(
        remaining_width - missing, remaining_kept - missing
    ) / math.comb(remaining_width, remaining_kept)


def planted_completion_probability(
    populations: Sequence[Sequence[Candidate]],
    true_positions: Sequence[Sequence[int]],
    block_length: int,
    factor_degree: int,
    kept: int,
) -> float:
    """Average planted-support selector success over empirical populations."""

    if len(populations) != len(true_positions):
        raise ValueError("one particle population is required per polynomial")
    state_probabilities: dict[tuple[int, int], float] = {(0, 0): 1.0}
    true_size = 0
    for population, positions in zip(populations, true_positions):
        if not population:
            raise ValueError("particle population is empty")
        truth = projected_candidate_support(
            positions, block_length, factor_degree
        )
        true_size += len(truth)
        histogram = Counter(
            (
                len(support),
                len(support & truth),
            )
            for support in (
                projected_candidate_support(
                    candidate, block_length, factor_degree
                )
                for candidate in population
            )
        )
        following: dict[tuple[int, int], float] = {}
        inverse_population = 1.0 / len(population)
        for (old_size, old_hits), old_probability in state_probabilities.items():
            for (size, hits), count in histogram.items():
                key = (old_size + size, old_hits + hits)
                following[key] = following.get(key, 0.0) + (
                    old_probability * count * inverse_population
                )
        state_probabilities = following

    width = len(populations) * factor_degree
    return sum(
        probability
        * completion_probability(
            width, kept, planted_size, true_size, planted_hits
        )
        for (planted_size, planted_hits), probability in state_probabilities.items()
    )


def sample_pair_population(
    distributions: Sequence[Sequence[float]],
    population_size: int,
    rng: random.Random,
) -> tuple[Candidate, ...]:
    """Sample complete candidates independently from pair marginals."""

    cumulative = []
    totals = []
    for distribution in distributions:
        running = 0.0
        row = []
        for probability in distribution:
            running += probability
            row.append(running)
        if running <= 0.0:
            raise ValueError("pair distribution has zero mass")
        cumulative.append(tuple(row))
        totals.append(running)
    return tuple(
        tuple(
            bisect.bisect_left(row, rng.random() * total)
            for row, total in zip(cumulative, totals)
        )
        for _ in range(population_size)
    )


def uniform_zero_set_probability(
    width: int, kept: int, support_size: int
) -> float:
    if support_size > kept:
        return 0.0
    return math.comb(width - support_size, kept - support_size) / math.comb(
        width, kept
    )


def log2_mean_inverse(probabilities: Sequence[float]) -> float:
    """Expected-work exponent, allowing a zero-success trial."""

    if any(probability <= 0.0 for probability in probabilities):
        return math.inf
    inverse_logs = [-math.log2(probability) for probability in probabilities]
    maximum = max(inverse_logs)
    return maximum + math.log2(
        sum(2.0 ** (value - maximum) for value in inverse_logs)
        / len(inverse_logs)
    )


def log2_inverse_mean(probabilities: Sequence[float]) -> float:
    """Average-signal exponent, allowing complete numerical underflow."""

    mean = sum(probabilities) / len(probabilities)
    return -math.log2(mean) if mean > 0.0 else math.inf


def one_trial(
    trial: int,
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    factor_degree: int,
    betas: Sequence[float],
    mixes: Sequence[float],
    pair_rounds: int,
    posterior_particles: int,
    posterior_moves: int,
    pair_control_particles: int,
    seed: int,
) -> MatchingFactorTrial:
    instance = sample_trial_instance(
        model,
        channel,
        pair_slope,
        pair_rounds,
        random.Random(seed),
    )
    pair_scores = []
    particle_scores = []
    matching_populations = []
    pair_populations = []
    diagnostics = []
    global_support: set[int] = set()
    for hidden_index, (positions, distributions) in enumerate(
        zip(instance.true_positions, instance.distributions)
    ):
        descriptors = factor_descriptors(
            instance.descriptors,
            hidden_index,
            model.polys,
            model.weight,
        )
        posterior_rng = random.Random(
            seed ^ (hidden_index << 32) ^ 0x4D41544348
        )
        diagnostic, particles = pair_proposal_target_particles(
            instance.known_positions,
            descriptors,
            distributions,
            model.weight,
            channel,
            posterior_particles,
            posterior_moves,
            posterior_rng,
        )
        if not particles:
            raise RuntimeError("matching SMC particle population collapsed")
        population = tuple(particle.candidate for particle in particles)
        matching_populations.append(population)
        diagnostics.append(diagnostic)
        pair_populations.append(
            sample_pair_population(
                distributions,
                pair_control_particles,
                random.Random(seed ^ (hidden_index << 32) ^ 0x50414952),
            )
        )
        pair_scores.append(
            projected_occupancy_scores(
                distributions, model.block_length, factor_degree
            )
        )
        particle_scores.append(
            particle_occupancy_scores(
                population, model.block_length, factor_degree
            )
        )
        base = hidden_index * factor_degree
        global_support.update(
            base + residue
            for residue in projected_candidate_support(
                positions, model.block_length, factor_degree
            )
        )

    marginal_successes = []
    for mix in mixes:
        scores = tuple(
            (1.0 - mix) * pair_score + mix * particle_score
            for pair_row, particle_row in zip(pair_scores, particle_scores)
            for pair_score, particle_score in zip(pair_row, particle_row)
        )
        marginal_successes.append(
            tuple(
                product_weighted_support_probability(
                    scores,
                    frozenset(global_support),
                    factor_degree - 1,
                    beta,
                )
                for beta in betas
            )
        )

    width = model.polys * factor_degree
    kept = factor_degree - 1
    return MatchingFactorTrial(
        trial=trial,
        reduced_weight=len(global_support),
        uniform_completion_probability=uniform_zero_set_probability(
            width, kept, len(global_support)
        ),
        pair_completion_probability=planted_completion_probability(
            pair_populations,
            instance.true_positions,
            model.block_length,
            factor_degree,
            kept,
        ),
        matching_completion_probability=planted_completion_probability(
            matching_populations,
            instance.true_positions,
            model.block_length,
            factor_degree,
            kept,
        ),
        marginal_success_probabilities=tuple(marginal_successes),
        minimum_smc_ess=min(
            diagnostic.minimum_ess_fraction for diagnostic in diagnostics
        ),
        minimum_smc_unique=min(
            diagnostic.minimum_unique_fraction for diagnostic in diagnostics
        ),
        smc_acceptance=statistics.fmean(
            diagnostic.mutation_acceptance for diagnostic in diagnostics
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=4096)
    parser.add_argument("--factor-degree", type=int, default=128)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--modulus", type=int, default=65537)
    parser.add_argument("--selector-k", type=int, default=1)
    parser.add_argument("--trials", type=int, default=2)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--posterior-particles", type=int, default=128)
    parser.add_argument("--posterior-moves", type=int, default=2)
    parser.add_argument("--pair-control-particles", type=int, default=4096)
    parser.add_argument(
        "--mixes", type=float, nargs="+", default=(0.0, 0.25, 0.5, 0.75)
    )
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys < 2:
        parser.error("the number of polynomials must be at least two")
    if args.factor_degree - 1 < args.polys * args.weight:
        parser.error("the dual zero set cannot contain the unreduced support")
    if min(
        args.trials,
        args.posterior_particles,
        args.pair_control_particles,
    ) < 1:
        parser.error("trial and particle counts must be positive")
    if any(not 0.0 <= mix < 1.0 for mix in args.mixes):
        parser.error("mixes must lie in [0,1)")

    partition_size = args.partition_size or args.weight
    model = FieldModel(
        polys=args.polys,
        weight=args.weight,
        block_length=args.block_length,
        modulus=args.modulus,
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
    start = time.perf_counter()
    trials = tuple(
        one_trial(
            trial,
            model,
            channel,
            pair_slope,
            args.factor_degree,
            BALANCED_BETAS,
            args.mixes,
            args.pair_rounds,
            args.posterior_particles,
            args.posterior_moves,
            args.pair_control_particles,
            args.seed + trial,
        )
        for trial in range(args.trials)
    )
    seconds = time.perf_counter() - start
    check_bits = math.log2(args.factor_degree)

    def work_bits(field: str) -> float:
        return log2_mean_inverse(
            [getattr(trial, field) for trial in trials]
        ) + check_bits

    completion_bits = {
        name: work_bits(f"{name}_completion_probability")
        for name in ("uniform", "pair", "matching")
    }
    uniform_bits = completion_bits["uniform"]
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} factor_degree={args.factor_degree} "
        f"K={args.selector_k} trials={args.trials} "
        f"posterior_particles={args.posterior_particles} "
        f"posterior_moves={args.posterior_moves} "
        f"pair_slope={pair_slope} pair_r_squared={pair_r_squared} "
        f"uniform_completion_bits={uniform_bits} "
        f"pair_completion_bits={completion_bits['pair']} "
        f"pair_completion_gain_bits={uniform_bits - completion_bits['pair']} "
        f"matching_completion_bits={completion_bits['matching']} "
        f"matching_completion_gain_bits={uniform_bits - completion_bits['matching']} "
        f"seconds={seconds} seed={args.seed}"
    )
    best = None
    for mix_index, mix in enumerate(args.mixes):
        for beta_index, beta in enumerate(BALANCED_BETAS):
            probabilities = [
                trial.marginal_success_probabilities[mix_index][beta_index]
                for trial in trials
            ]
            bits = log2_mean_inverse(probabilities) + check_bits
            signal_bits = log2_inverse_mean(probabilities) + check_bits
            gain = uniform_bits - bits
            signal_gain = (
                log2_inverse_mean(
                    [trial.uniform_completion_probability for trial in trials]
                )
                + check_bits
                - signal_bits
            )
            candidate = (gain, mix, beta, bits, signal_gain, signal_bits)
            if best is None or candidate[0] > best[0]:
                best = candidate
            print(
                f"mix={mix} beta={beta} marginal_bits={bits} "
                f"marginal_gain_bits={gain} marginal_signal_bits={signal_bits} "
                f"marginal_signal_gain_bits={signal_gain}"
            )
    assert best is not None
    print(
        f"best_marginal_mix={best[1]} best_marginal_beta={best[2]} "
        f"best_marginal_bits={best[3]} best_marginal_gain_bits={best[0]} "
        f"best_marginal_signal_bits={best[5]} "
        f"best_marginal_signal_gain_bits={best[4]}"
    )
    for trial in trials:
        print(
            f"trial={trial.trial} reduced_weight={trial.reduced_weight} "
            f"pair_completion_log_ratio="
            f"{math.log2(trial.pair_completion_probability / trial.uniform_completion_probability)} "
            f"matching_completion_log_ratio="
            f"{math.log2(trial.matching_completion_probability / trial.uniform_completion_probability)} "
            f"minimum_smc_ess={trial.minimum_smc_ess} "
            f"minimum_smc_unique={trial.minimum_smc_unique} "
            f"smc_acceptance={trial.smc_acceptance}"
        )

    if args.output:
        fields = (
            "trial",
            "reduced_weight",
            "uniform_completion_probability",
            "pair_completion_probability",
            "matching_completion_probability",
            "minimum_smc_ess",
            "minimum_smc_unique",
            "smc_acceptance",
        ) + tuple(
            f"marginal_mix_{mix}_beta_{beta}"
            for mix in args.mixes
            for beta in BALANCED_BETAS
        )
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for trial in trials:
                row = {
                    field: getattr(trial, field)
                    for field in fields[:8]
                }
                for mix_index, mix in enumerate(args.mixes):
                    for beta_index, beta in enumerate(BALANCED_BETAS):
                        row[f"marginal_mix_{mix}_beta_{beta}"] = (
                            trial.marginal_success_probabilities[mix_index][
                                beta_index
                            ]
                        )
                writer.writerow(row)
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
