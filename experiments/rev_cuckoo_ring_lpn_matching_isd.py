"""Global matching-particle selector for repeatable regular Prange.

The pair-only sequential attack reacts to collision labels selected in earlier
blocks, but never evaluates a complete cuckoo matching factor.  This variant
first samples either the exact matching-likelihood posterior with the existing
factorwise SMC engine or the explicit attack distribution q_pair W_match with
a pair-marginal proposal.  During each Prange iteration, the particles supply
a conditional histogram for the next block.  The histogram is blended with
the full-support pair distribution and passed to the same fixed-cardinality
stratified sampler.  After a subset is chosen, the particle population is
conditioned and resampled, so later blocks follow the mode selected by the
earlier subsets.

The attack never sees the planted support.  As in the pair-only experiment, a
second particle filter conditions the attack on the planted coordinates only
to evaluate an approximately 2^-128 success probability.
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

from rev_cuckoo_basic_ring_lpn_isd import (
    best_weighted_expected_work_gain,
    best_weighted_half_success_gain,
    sample_balanced_stratified_subset_conditioned,
)
from rev_cuckoo_basic_ring_lpn_isd_production import factor_descriptors
from rev_cuckoo_full_lpn_field import FieldModel
from rev_cuckoo_regular_isd import CollisionFactor, fit_pair_slope
from rev_cuckoo_regular_isd_smc import (
    Particle,
    SmcResult,
    effective_sample_fraction,
    evaluate_processed_factors,
    factor_log_likelihood,
    normalized_weights,
    smc_target_particles,
    systematic_indices,
)
from rev_cuckoo_regular_isd_mc import logmeanexp2
from rev_cuckoo_regular_noise import LeakageChannel
from rev_cuckoo_ring_lpn_sequential_isd import (
    add_distribution_masses,
    add_subset_masses,
    corrected_distribution,
    entropy,
    independent_polynomial_success,
    label_arrays,
    sequential_polynomial_success,
    systematic_resample,
    sample_trial_instance,
)


Candidate = tuple[int, ...]


@dataclass(frozen=True)
class TrialResult:
    independent_success: float
    pair_success: float
    matching_success: float
    minimum_smc_ess: float
    minimum_smc_unique: float
    smc_acceptance: float


def matching_distribution(
    pair_distribution: np.ndarray,
    contexts: tuple[Candidate, ...],
    variable: int,
    block_length: int,
    posterior_mix: float,
) -> np.ndarray:
    """Blend the pair message with a global-posterior particle histogram."""

    if not 0.0 <= posterior_mix < 1.0:
        raise ValueError("posterior mix must lie in [0,1)")
    if not contexts or posterior_mix == 0.0:
        return pair_distribution
    counts = np.bincount(
        np.fromiter(
            (context[variable] for context in contexts),
            dtype=np.int64,
            count=len(contexts),
        ),
        minlength=block_length,
    ).astype(np.float64)
    posterior = counts / len(contexts)
    result = (1.0 - posterior_mix) * pair_distribution
    result += posterior_mix * posterior
    return result / np.sum(result)


def condition_contexts(
    contexts: tuple[Candidate, ...],
    variable: int,
    subset: frozenset[int],
    rng: random.Random,
) -> tuple[Candidate, ...]:
    """Condition and bootstrap-resample an unweighted posterior population."""

    survivors = tuple(
        context for context in contexts if context[variable] in subset
    )
    if not survivors:
        return ()
    return tuple(rng.choices(survivors, k=len(contexts)))


def matching_sequential_polynomial_success(
    true_positions: tuple[int, ...],
    distributions: tuple[tuple[float, ...], ...],
    factors: tuple[CollisionFactor, ...],
    posterior_contexts: tuple[Candidate, ...],
    kept: int,
    beta: float,
    slope: float,
    posterior_mix: float,
    particles: int,
    rng: random.Random,
) -> float:
    """Evaluate one matching-guided hidden-polynomial support sampler."""

    base = np.asarray(distributions, dtype=np.float64)
    block_length = base.shape[1]
    labels, label_count = label_arrays(factors)
    order = tuple(
        sorted(
            range(len(distributions)),
            key=lambda index: entropy(base[index]),
        )
    )
    baseline_masses = np.zeros(
        (labels.shape[0], label_count), dtype=np.float64
    )
    particle_masses = [np.zeros_like(baseline_masses) for _ in range(particles)]
    context_populations = [posterior_contexts for _ in range(particles)]
    log_success = 0.0
    for variable in order:
        children_masses = []
        children_contexts = []
        inclusion_probabilities = []
        for masses, contexts in zip(particle_masses, context_populations):
            pair_distribution = corrected_distribution(
                base[variable],
                labels,
                variable,
                masses,
                baseline_masses,
                slope=slope,
            )
            distribution = matching_distribution(
                pair_distribution,
                contexts,
                variable,
                block_length,
                posterior_mix,
            )
            subset, inclusion = sample_balanced_stratified_subset_conditioned(
                distribution,
                kept,
                beta,
                true_positions[variable],
                rng,
            )
            child_masses = masses.copy()
            add_subset_masses(
                child_masses,
                labels,
                variable,
                subset,
                kept,
                label_count,
            )
            children_masses.append(child_masses)
            children_contexts.append(
                condition_contexts(contexts, variable, subset, rng)
            )
            inclusion_probabilities.append(inclusion)
        mean_inclusion = sum(inclusion_probabilities) / particles
        log_success += math.log2(mean_inclusion)
        indices = systematic_resample(
            inclusion_probabilities, particles, rng
        )
        particle_masses = [children_masses[index].copy() for index in indices]
        context_populations = [children_contexts[index] for index in indices]
        add_distribution_masses(
            baseline_masses,
            labels,
            variable,
            base[variable],
            label_count,
        )
    return 2.0 ** log_success


def exact_posterior_contexts(
    known: tuple[int, ...],
    descriptors,
    weight: int,
    block_length: int,
    channel: LeakageChannel,
    count: int,
    rng: random.Random,
) -> tuple[Candidate, ...]:
    """Sample contexts from an exactly enumerated reduced posterior."""

    candidate_count = block_length ** weight
    if candidate_count > 1_000_000:
        raise ValueError(
            "exact posterior is restricted to at most one million candidates"
        )
    candidates = tuple(
        itertools.product(range(block_length), repeat=weight)
    )
    log_weights = []
    for candidate in candidates:
        values = [
            factor_log_likelihood(
                known,
                candidate,
                factor,
                descriptors,
                weight,
                channel,
            )
            for factor in range(len(descriptors))
        ]
        log_weights.append(
            -math.inf
            if any(not math.isfinite(value) for value in values)
            else sum(values)
        )
    maximum = max(log_weights)
    if not math.isfinite(maximum):
        raise RuntimeError("reduced exact posterior has no feasible candidate")
    weights = [
        0.0 if not math.isfinite(value) else 2.0 ** (value - maximum)
        for value in log_weights
    ]
    return tuple(rng.choices(candidates, weights=weights, k=count))


def pair_proposal_rejuvenate(
    particles: list[Particle],
    distributions: tuple[tuple[float, ...], ...],
    known: tuple[int, ...],
    processed: tuple[int, ...],
    descriptors,
    weight: int,
    channel: LeakageChannel,
    moves: int,
    rng: random.Random,
) -> tuple[list[Particle], int, int]:
    """Rejuvenate target q_pair times W using q_pair proposals.

    Replacing one coordinate from its pair marginal is an independence
    proposal.  Its forward/reverse ratio cancels the corresponding q_pair
    factor in the target, so acceptance depends only on the exact processed
    matching likelihoods.
    """

    result = []
    accepted = 0
    attempted = 0
    domain = range(len(distributions[0]))
    for source in particles:
        particle = Particle(
            source.candidate,
            list(source.factor_logs),
            source.total_log,
        )
        for _ in range(moves):
            variable = rng.randrange(weight)
            value = rng.choices(domain, weights=distributions[variable], k=1)[0]
            if value == particle.candidate[variable]:
                continue
            proposal_values = list(particle.candidate)
            proposal_values[variable] = value
            proposal = tuple(proposal_values)
            attempted += 1
            factor_logs, total_log = evaluate_processed_factors(
                known,
                proposal,
                processed,
                descriptors,
                weight,
                channel,
            )
            if not math.isfinite(total_log):
                continue
            delta = total_log - particle.total_log
            if delta >= 0.0 or rng.random() < 2.0 ** delta:
                particle = Particle(proposal, factor_logs, total_log)
                accepted += 1
        result.append(particle)
    return result, accepted, attempted


def pair_proposal_target_particles(
    known: tuple[int, ...],
    descriptors,
    distributions: tuple[tuple[float, ...], ...],
    weight: int,
    channel: LeakageChannel,
    particle_count: int,
    moves_per_stage: int,
    rng: random.Random,
) -> tuple[SmcResult, list[Particle]]:
    """Sample the explicit attack target q_pair times W_match."""

    domain = range(len(distributions[0]))
    particles = [
        Particle(
            tuple(
                rng.choices(domain, weights=distribution, k=1)[0]
                for distribution in distributions
            ),
            [],
            0.0,
        )
        for _ in range(particle_count)
    ]
    factors = list(range(len(descriptors)))
    rng.shuffle(factors)
    processed = []
    log_mean_weight = 0.0
    ess_fractions = []
    unique_fractions = []
    total_accepted = 0
    total_attempted = 0
    for factor in factors:
        incremental_logs = []
        for particle in particles:
            value = factor_log_likelihood(
                known,
                particle.candidate,
                factor,
                descriptors,
                weight,
                channel,
            )
            particle.factor_logs.append(value)
            if math.isfinite(particle.total_log) and math.isfinite(value):
                particle.total_log += value
            else:
                particle.total_log = -math.inf
            incremental_logs.append(value)
        weights = normalized_weights(incremental_logs)
        if not weights:
            return SmcResult(-math.inf, 0.0, 0.0, 0.0, 0.0), []
        log_mean_weight += logmeanexp2(incremental_logs)
        ess_fractions.append(effective_sample_fraction(weights))
        indices = systematic_indices(weights, rng)
        particles = [
            Particle(
                particles[index].candidate,
                list(particles[index].factor_logs),
                particles[index].total_log,
            )
            for index in indices
        ]
        unique_fractions.append(
            len({particle.candidate for particle in particles}) / particle_count
        )
        processed.append(factor)
        particles, accepted, attempted = pair_proposal_rejuvenate(
            particles,
            distributions,
            known,
            tuple(processed),
            descriptors,
            weight,
            channel,
            moves_per_stage,
            rng,
        )
        total_accepted += accepted
        total_attempted += attempted
    return (
        SmcResult(
            log_mean_weight,
            min(ess_fractions),
            statistics.fmean(ess_fractions),
            min(unique_fractions),
            total_accepted / total_attempted if total_attempted else 0.0,
        ),
        particles,
    )


def one_trial(
    model: FieldModel,
    channel: LeakageChannel,
    pair_slope: float,
    pair_rounds: int,
    beta: float,
    posterior_mix: float,
    posterior_particles: int,
    posterior_moves: int,
    exact_posterior: bool,
    posterior_proposal: str,
    selection_particles: int,
    selection_repeats: int,
    seed: int,
) -> TrialResult:
    instance = sample_trial_instance(
        model,
        channel,
        pair_slope,
        pair_rounds,
        random.Random(seed),
    )
    kept = model.block_length // model.polys
    full_blocks = tuple(
        tuple(range(model.block_length)) for _ in range(model.weight)
    )
    posterior_populations = []
    diagnostics: list[SmcResult] = []
    for hidden_index in range(model.polys):
        descriptors = factor_descriptors(
            instance.descriptors,
            hidden_index,
            model.polys,
            model.weight,
        )
        posterior_rng = random.Random(
            seed ^ (hidden_index << 32) ^ 0x4D41544348
        )
        if exact_posterior:
            contexts = exact_posterior_contexts(
                instance.known_positions,
                descriptors,
                model.weight,
                model.block_length,
                channel,
                posterior_particles,
                posterior_rng,
            )
            diagnostic = SmcResult(0.0, 1.0, 1.0, 1.0, 1.0)
            particles = ()
        elif posterior_proposal == "pair":
            diagnostic, particles = pair_proposal_target_particles(
                instance.known_positions,
                descriptors,
                instance.distributions[hidden_index],
                model.weight,
                channel,
                posterior_particles,
                posterior_moves,
                posterior_rng,
            )
            contexts = tuple(particle.candidate for particle in particles)
        else:
            diagnostic, particles = smc_target_particles(
                instance.known_positions,
                descriptors,
                full_blocks,
                model.weight,
                channel,
                posterior_particles,
                posterior_moves,
                posterior_rng,
            )
            contexts = tuple(particle.candidate for particle in particles)
        diagnostics.append(diagnostic)
        posterior_populations.append(contexts)

    independent = 1.0
    for positions, distributions in zip(
        instance.true_positions, instance.distributions
    ):
        independent *= independent_polynomial_success(
            positions, distributions, kept, beta
        )

    pair_repeats = []
    matching_repeats = []
    for repeat in range(selection_repeats):
        pair = 1.0
        matching = 1.0
        for hidden_index in range(model.polys):
            repeat_seed = seed ^ (repeat << 40) ^ (hidden_index << 24)
            pair *= sequential_polynomial_success(
                instance.true_positions[hidden_index],
                instance.distributions[hidden_index],
                instance.factor_sets[hidden_index],
                kept,
                beta,
                pair_slope,
                selection_particles,
                random.Random(repeat_seed ^ 0x50414952),
            )
            matching *= matching_sequential_polynomial_success(
                instance.true_positions[hidden_index],
                instance.distributions[hidden_index],
                instance.factor_sets[hidden_index],
                posterior_populations[hidden_index],
                kept,
                beta,
                pair_slope,
                posterior_mix,
                selection_particles,
                random.Random(repeat_seed ^ 0x4D41544348),
            )
        pair_repeats.append(pair)
        matching_repeats.append(matching)

    finite_diagnostics = [
        diagnostic
        for diagnostic in diagnostics
        if math.isfinite(diagnostic.log_mean_weight)
    ]
    return TrialResult(
        independent_success=independent,
        pair_success=sum(pair_repeats) / len(pair_repeats),
        matching_success=sum(matching_repeats) / len(matching_repeats),
        minimum_smc_ess=min(
            (diagnostic.minimum_ess_fraction for diagnostic in finite_diagnostics),
            default=0.0,
        ),
        minimum_smc_unique=min(
            (diagnostic.minimum_unique_fraction for diagnostic in finite_diagnostics),
            default=0.0,
        ),
        smc_acceptance=(
            sum(diagnostic.mutation_acceptance for diagnostic in finite_diagnostics)
            / len(finite_diagnostics)
            if finite_diagnostics
            else 0.0
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=16)
    parser.add_argument("--block-length", type=int, default=4096)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--modulus", type=int, default=65537)
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--posterior-mix", type=float, default=0.5)
    parser.add_argument("--posterior-particles", type=int, default=128)
    parser.add_argument("--posterior-moves", type=int, default=2)
    parser.add_argument("--exact-posterior", action="store_true")
    parser.add_argument(
        "--posterior-proposal",
        choices=("uniform", "pair"),
        default="uniform",
    )
    parser.add_argument("--selection-particles", type=int, default=8)
    parser.add_argument("--selection-repeats", type=int, default=2)
    parser.add_argument("--trials", type=int, default=2)
    parser.add_argument("--pair-rounds", type=int, default=8)
    parser.add_argument("--pair-fit-samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.polys < 2 or args.block_length % args.polys:
        parser.error("P must be at least two and divide the block length")
    if not 0.0 <= args.posterior_mix < 1.0:
        parser.error("posterior mix must lie in [0,1)")
    if min(
        args.posterior_particles,
        args.selection_particles,
        args.selection_repeats,
        args.trials,
    ) < 1:
        parser.error("particle, repeat, and trial counts must be positive")
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
    results = [
        one_trial(
            model,
            channel,
            pair_slope,
            args.pair_rounds,
            args.beta,
            args.posterior_mix,
            args.posterior_particles,
            args.posterior_moves,
            args.exact_posterior,
            args.posterior_proposal,
            args.selection_particles,
            args.selection_repeats,
            args.seed + trial,
        )
        for trial in range(args.trials)
    ]
    seconds = time.perf_counter() - start
    baseline_bits = args.polys * args.weight * math.log2(args.polys)
    rows = {
        "independent": [(result.independent_success,) for result in results],
        "pair": [(result.pair_success,) for result in results],
        "matching": [(result.matching_success,) for result in results],
    }
    gains = {
        name: best_weighted_expected_work_gain(values, baseline_bits, (args.beta,))[0]
        for name, values in rows.items()
    }
    half_gains = {
        name: best_weighted_half_success_gain(values, baseline_bits, (args.beta,))[0]
        for name, values in rows.items()
    }
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} K={args.selector_k} "
        f"beta={args.beta} posterior_mix={args.posterior_mix} "
        f"posterior_particles={args.posterior_particles} "
        f"posterior_moves={args.posterior_moves} "
        f"exact_posterior={args.exact_posterior} "
        f"posterior_proposal={args.posterior_proposal} "
        f"selection_particles={args.selection_particles} "
        f"selection_repeats={args.selection_repeats} trials={args.trials} "
        f"pair_slope={pair_slope} pair_r_squared={pair_r_squared} "
        f"independent_expected_work_gain_bits={gains['independent']} "
        f"pair_expected_work_gain_bits={gains['pair']} "
        f"matching_expected_work_gain_bits={gains['matching']} "
        f"independent_half_success_gain_bits={half_gains['independent']} "
        f"pair_half_success_gain_bits={half_gains['pair']} "
        f"matching_half_success_gain_bits={half_gains['matching']} "
        f"seconds={seconds} seed={args.seed}"
    )
    for trial, result in enumerate(results):
        print(
            f"trial={trial} "
            f"independent_log_gain={math.log2(result.independent_success) + baseline_bits} "
            f"pair_log_gain={math.log2(result.pair_success) + baseline_bits} "
            f"matching_log_gain={math.log2(result.matching_success) + baseline_bits} "
            f"minimum_smc_ess={result.minimum_smc_ess} "
            f"minimum_smc_unique={result.minimum_smc_unique} "
            f"smc_acceptance={result.smc_acceptance}"
        )
    if args.output:
        fields = tuple(TrialResult.__dataclass_fields__) + ("trial",)
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for trial, result in enumerate(results):
                row = {
                    field: getattr(result, field)
                    for field in TrialResult.__dataclass_fields__
                }
                row["trial"] = trial
                writer.writerow(row)
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
