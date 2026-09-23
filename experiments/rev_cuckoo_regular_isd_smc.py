"""Sequential-Monte-Carlo scorer for leakage-aware regular Prange.

The exact posterior rectangle mass is

    P^(-t) E[W(Y) | Y uniform in R] / E[W(Y) | Y uniform],

where W is the product of the finite-K Reverse-Cuckoo list likelihoods.  Naive
importance sampling has very poor ESS at t=16 because most wrong candidates
eventually violate a hard matching constraint.  This program estimates each
expectation by introducing the P*t list factors sequentially, resampling after
every factor, and rejuvenating particles with exact single-coordinate
Metropolis moves.

This is a research diagnostic.  Its SMC variance and particle diversity must
be reported with every gain estimate; a numerical gain without those
diagnostics is not evidence of a decoding attack.

The optional joint-conditional selector constructs a rectangle from complete
batch likelihoods rather than factorwise messages.  Posterior particles supply
contexts for the other offsets; every value of the current offset is scored by
the product of all exact list likelihoods.  Rectangle scoring uses independent
SMC populations to prevent selector overfitting.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
import statistics
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

from rev_cuckoo_regular_isd import (
    fit_pair_slope,
    pair_mean_field_distributions,
    pair_refined_rectangle_selection,
)
from rev_cuckoo_regular_isd_mc import (
    logmeanexp2,
    sample_uniform_candidate,
)
from rev_cuckoo_goldreich import SharedGoldreichChannel
from rev_cuckoo_regular_noise import LeakageChannel, support_lists_one_unknown


Candidate = tuple[int, ...]
Descriptor = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass
class Particle:
    candidate: Candidate
    factor_logs: list[float]
    total_log: float


@dataclass(frozen=True)
class SmcResult:
    log_mean_weight: float
    minimum_ess_fraction: float
    mean_ess_fraction: float
    minimum_unique_fraction: float
    mutation_acceptance: float


@dataclass(frozen=True)
class SplitResult:
    gain_bits: float
    posterior: SmcResult
    minimum_survival_fraction: float
    minimum_unique_fraction: float
    mutation_acceptance: float


def factor_active_set(
    known: Sequence[int],
    candidate: Sequence[int],
    known_poly: int,
    output_block: int,
    weight: int,
) -> tuple[int, ...]:
    active = {
        known[known_poly * weight + item]
        + candidate[(output_block - item) % weight]
        for item in range(weight)
    }
    return tuple(active)


def factor_log_likelihood(
    known: Sequence[int],
    candidate: Sequence[int],
    factor: int,
    descriptors: Sequence[Descriptor],
    weight: int,
    channel: LeakageChannel,
) -> float:
    known_poly, output_block = divmod(factor, weight)
    active = factor_active_set(
        known, candidate, known_poly, output_block, weight
    )
    return channel.log_likelihood(active, descriptors[factor])


def matching_refined_rectangle_selection(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    initial: Sequence[Sequence[int]],
    weight: int,
    block_length: int,
    channel: LeakageChannel,
    contexts: int,
    rounds: int,
    rng: random.Random,
) -> tuple[frozenset[int], ...]:
    """Refine a rectangle with sampled exact matching-factor messages."""

    selection = [set(block) for block in initial]
    kept = len(selection[0])
    infeasible_penalty = -64.0
    for _ in range(rounds):
        changed = False
        for variable in range(weight):
            scores = [0.0] * block_length
            for factor, descriptor in enumerate(descriptors):
                known_poly, output_block = divmod(factor, weight)
                item = (output_block - variable) % weight
                variable_base = known[known_poly * weight + item]
                context_scores = []
                for _ in range(contexts):
                    candidate = [
                        rng.choice(tuple(selection[other]))
                        for other in range(weight)
                    ]
                    active_base = {
                        known[
                            known_poly * weight
                            + (output_block - other) % weight
                        ]
                        + candidate[other]
                        for other in range(weight)
                        if other != variable
                    }
                    memo: dict[tuple[int, ...], float] = {}
                    values = []
                    for value in range(block_length):
                        address = variable_base + value
                        if address in active_base:
                            key = (1,)
                            active = active_base
                        else:
                            key = (
                                0,
                                descriptor[0][address],
                                descriptor[1][address],
                            )
                            active = active_base | {address}
                        likelihood = memo.get(key)
                        if likelihood is None:
                            likelihood = channel.log_likelihood(
                                tuple(active), descriptor
                            )
                            memo[key] = likelihood
                        values.append(likelihood)
                    context_scores.append(values)
                for value in range(block_length):
                    likelihood = logmeanexp2(
                        [row[value] for row in context_scores]
                    )
                    scores[value] += (
                        likelihood
                        if math.isfinite(likelihood)
                        else infeasible_penalty
                    )
            following = set(
                sorted(
                    range(block_length),
                    key=lambda value: scores[value],
                    reverse=True,
                )[:kept]
            )
            if following != selection[variable]:
                selection[variable] = following
                changed = True
        if not changed:
            break
    return tuple(frozenset(block) for block in selection)


def conditional_value_logs(
    known: Sequence[int],
    context: Candidate,
    variable: int,
    descriptors: Sequence[Descriptor],
    weight: int,
    block_length: int,
    channel: LeakageChannel,
) -> list[float]:
    """Exact full-batch log likelihood for all values of one offset."""

    totals = [0.0] * block_length
    possible = [True] * block_length
    for factor, descriptor in enumerate(descriptors):
        known_poly, output_block = divmod(factor, weight)
        item = (output_block - variable) % weight
        variable_base = known[known_poly * weight + item]
        active_base = {
            known[
                known_poly * weight
                + (output_block - other) % weight
            ]
            + context[other]
            for other in range(weight)
            if other != variable
        }
        memo: dict[tuple[int, ...], float] = {}
        for value in range(block_length):
            if not possible[value]:
                continue
            address = variable_base + value
            if address in active_base:
                key = (1,)
                active = active_base
            else:
                key = (
                    0,
                    descriptor[0][address],
                    descriptor[1][address],
                )
                active = active_base | {address}
            likelihood = memo.get(key)
            if likelihood is None:
                likelihood = channel.log_likelihood(
                    tuple(active), descriptor
                )
                memo[key] = likelihood
            if math.isfinite(likelihood):
                totals[value] += likelihood
            else:
                possible[value] = False
                totals[value] = -math.inf
    return totals


def normalized_log_distribution(log_values: Sequence[float]) -> list[float]:
    finite = [value for value in log_values if math.isfinite(value)]
    if not finite:
        return [0.0] * len(log_values)
    maximum = max(finite)
    weights = [
        0.0 if not math.isfinite(value) else 2.0 ** (value - maximum)
        for value in log_values
    ]
    total = sum(weights)
    return [weight / total for weight in weights]


def joint_conditional_rectangle_selection(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    weight: int,
    block_length: int,
    kept: int,
    channel: LeakageChannel,
    particles_count: int,
    moves_per_stage: int,
    restriction_moves: int,
    contexts: int,
    rng: random.Random,
) -> tuple[tuple[frozenset[int], ...], SmcResult]:
    """Greedily build a rectangle from exact full-posterior conditionals.

    Posterior particles supply contexts for the other offsets.  For every
    context, all values of the current offset are scored by the product of all
    descriptor likelihoods and normalized before contexts are averaged.  The
    construction conditions and rejuvenates its training population after
    each selected block.  Scoring must use an independent SMC population.
    """

    if contexts < 1:
        raise ValueError("joint selector contexts must be positive")
    full_blocks = tuple(
        tuple(range(block_length)) for _ in range(weight)
    )
    diagnostic, particles = smc_target_particles(
        known,
        descriptors,
        full_blocks,
        weight,
        channel,
        particles_count,
        moves_per_stage,
        rng,
    )
    if not particles:
        return tuple(frozenset(block) for block in full_blocks), diagnostic

    allowed_blocks = [tuple(block) for block in full_blocks]
    selection = [frozenset(block) for block in full_blocks]
    all_factors = tuple(range(len(descriptors)))
    variables = list(range(weight))
    rng.shuffle(variables)
    for variable in variables:
        unique = list({particle.candidate for particle in particles})
        sources = unique if len(unique) <= contexts else rng.sample(unique, contexts)
        scores = [0.0] * block_length
        used_contexts = 0
        for context in sources:
            distribution = normalized_log_distribution(
                conditional_value_logs(
                    known,
                    context,
                    variable,
                    descriptors,
                    weight,
                    block_length,
                    channel,
                )
            )
            if any(distribution):
                used_contexts += 1
                for value, probability in enumerate(distribution):
                    scores[value] += probability
        if used_contexts == 0:
            continue
        selected = frozenset(
            sorted(
                range(block_length),
                key=lambda value: scores[value],
                reverse=True,
            )[:kept]
        )
        selection[variable] = selected
        allowed_blocks[variable] = tuple(selected)

        survivors = [
            particle
            for particle in particles
            if particle.candidate[variable] in selected
        ]
        if not survivors:
            refreshed, particles = smc_target_particles(
                known,
                descriptors,
                allowed_blocks,
                weight,
                channel,
                particles_count,
                moves_per_stage,
                rng,
            )
            diagnostic = refreshed
            if not particles:
                break
            continue
        indices = systematic_indices(
            [1.0 / len(survivors)] * len(survivors),
            rng,
            particles_count,
        )
        particles = [
            Particle(
                survivors[index].candidate,
                list(survivors[index].factor_logs),
                survivors[index].total_log,
            )
            for index in indices
        ]
        particles, _, _ = rejuvenate(
            particles,
            allowed_blocks,
            known,
            all_factors,
            descriptors,
            weight,
            channel,
            restriction_moves,
            rng,
        )
    return tuple(selection), diagnostic


def diverse_posterior_anchors(
    particles: Sequence[Particle], count: int
) -> tuple[Candidate, ...]:
    """Select frequent, Hamming-diverse posterior modes."""

    if count < 1:
        raise ValueError("portfolio size must be positive")
    frequencies = Counter(particle.candidate for particle in particles)
    if not frequencies:
        return ()
    ordered = sorted(
        frequencies,
        key=lambda candidate: (frequencies[candidate], candidate),
        reverse=True,
    )
    anchors = [ordered[0]]
    remaining = set(ordered[1:])
    while remaining and len(anchors) < count:
        following = max(
            remaining,
            key=lambda candidate: (
                min(
                    sum(left != right for left, right in zip(candidate, anchor))
                    for anchor in anchors
                ),
                frequencies[candidate],
                candidate,
            ),
        )
        anchors.append(following)
        remaining.remove(following)
    return tuple(anchors)


def anchored_rectangle_selection(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    anchor: Candidate,
    weight: int,
    block_length: int,
    kept: int,
    channel: LeakageChannel,
) -> tuple[frozenset[int], ...]:
    """Build one mode-specific rectangle around a posterior anchor."""

    selection = []
    for variable in range(weight):
        distribution = normalized_log_distribution(
            conditional_value_logs(
                known,
                anchor,
                variable,
                descriptors,
                weight,
                block_length,
                channel,
            )
        )
        selection.append(
            frozenset(
                sorted(
                    range(block_length),
                    key=lambda value: distribution[value],
                    reverse=True,
                )[:kept]
            )
        )
    return tuple(selection)


def contextual_rectangle_portfolio(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    weight: int,
    block_length: int,
    kept: int,
    channel: LeakageChannel,
    portfolio_size: int,
    particles_count: int,
    moves_per_stage: int,
    rng: random.Random,
) -> tuple[tuple[tuple[frozenset[int], ...], ...], SmcResult]:
    """Construct a charged portfolio of diverse mode-specific rectangles."""

    full_blocks = tuple(
        tuple(range(block_length)) for _ in range(weight)
    )
    diagnostic, particles = smc_target_particles(
        known,
        descriptors,
        full_blocks,
        weight,
        channel,
        particles_count,
        moves_per_stage,
        rng,
    )
    anchors = diverse_posterior_anchors(particles, portfolio_size)
    rectangles = tuple(
        anchored_rectangle_selection(
            known,
            descriptors,
            anchor,
            weight,
            block_length,
            kept,
            channel,
        )
        for anchor in anchors
    )
    return rectangles, diagnostic


def posterior_single_offset_set(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    weight: int,
    block_length: int,
    channel: LeakageChannel,
    particles_count: int,
    moves_per_stage: int,
    contexts: int,
    guess_width: int,
    rng: random.Random,
) -> tuple[
    tuple[frozenset[int], ...],
    tuple[int, tuple[int, ...]],
    float,
    SmcResult,
    SmcResult,
]:
    """Select and independently score a posterior offset-value set."""

    if contexts < 1:
        raise ValueError("guess contexts must be positive")
    if not 1 <= guess_width <= block_length:
        raise ValueError("guess width must lie in the offset domain")
    full_blocks = tuple(
        tuple(range(block_length)) for _ in range(weight)
    )
    diagnostic, particles = smc_target_particles(
        known,
        descriptors,
        full_blocks,
        weight,
        channel,
        particles_count,
        moves_per_stage,
        rng,
    )
    if not particles:
        raise RuntimeError("guess selector posterior is empty")
    candidates = [particle.candidate for particle in particles]
    best_variable = 0
    best_values: tuple[int, ...] = ()
    best_count = -1
    for variable in range(weight):
        top = Counter(
            candidate[variable] for candidate in candidates
        ).most_common(guess_width)
        observed_values = tuple(value for value, _ in top)
        count = sum(frequency for _, frequency in top)
        if count > best_count:
            best_variable = variable
            best_values = observed_values
            best_count = count
    if len(best_values) < guess_width:
        selected = set(best_values)
        fillers = rng.sample(
            [value for value in range(block_length) if value not in selected],
            guess_width - len(best_values),
        )
        best_values += tuple(fillers)
    selection = [frozenset(block) for block in full_blocks]
    selection[best_variable] = frozenset(best_values)
    evaluation_diagnostic, evaluation_particles = smc_target_particles(
        known,
        descriptors,
        full_blocks,
        weight,
        channel,
        particles_count,
        moves_per_stage,
        rng,
    )
    if not evaluation_particles:
        raise RuntimeError("guess evaluator posterior is empty")
    evaluation_candidates = [
        particle.candidate for particle in evaluation_particles
    ]
    evaluation_sources = (
        rng.sample(evaluation_candidates, contexts)
        if contexts <= len(evaluation_candidates)
        else rng.choices(evaluation_candidates, k=contexts)
    )
    evaluated_probability = 0.0
    for context in evaluation_sources:
        distribution = normalized_log_distribution(
            conditional_value_logs(
                known,
                context,
                best_variable,
                descriptors,
                weight,
                block_length,
                channel,
            )
        )
        evaluated_probability += sum(
            distribution[value] for value in best_values
        )
    evaluated_probability /= len(evaluation_sources)
    return (
        tuple(selection),
        (best_variable, best_values),
        evaluated_probability,
        diagnostic,
        evaluation_diagnostic,
    )


def normalized_weights(log_weights: Sequence[float]) -> list[float]:
    maximum = max(log_weights)
    if not math.isfinite(maximum):
        return []
    weights = [
        0.0 if not math.isfinite(value) else 2.0 ** (value - maximum)
        for value in log_weights
    ]
    total = sum(weights)
    return [weight / total for weight in weights]


def effective_sample_fraction(weights: Sequence[float]) -> float:
    return 1.0 / (len(weights) * sum(weight * weight for weight in weights))


def systematic_indices(
    weights: Sequence[float], rng: random.Random, output_count: int | None = None
) -> list[int]:
    count = output_count or len(weights)
    start = rng.random() / count
    indices = []
    cumulative = weights[0]
    source = 0
    for target in range(count):
        point = start + target / count
        while point > cumulative and source + 1 < count:
            source += 1
            cumulative += weights[source]
        indices.append(source)
    return indices


def sample_from_blocks(
    blocks: Sequence[Sequence[int]], rng: random.Random
) -> Candidate:
    return tuple(rng.choice(block) for block in blocks)


def propose_coordinate(
    candidate: Candidate,
    blocks: Sequence[Sequence[int]],
    rng: random.Random,
) -> Candidate:
    variable = rng.randrange(len(candidate))
    following = list(candidate)
    following[variable] = rng.choice(blocks[variable])
    return tuple(following)


def evaluate_processed_factors(
    known: Sequence[int],
    candidate: Candidate,
    processed: Sequence[int],
    descriptors: Sequence[Descriptor],
    weight: int,
    channel: LeakageChannel,
) -> tuple[list[float], float]:
    factor_logs = [
        factor_log_likelihood(
            known, candidate, factor, descriptors, weight, channel
        )
        for factor in processed
    ]
    if any(not math.isfinite(value) for value in factor_logs):
        return factor_logs, -math.inf
    return factor_logs, sum(factor_logs)


def rejuvenate(
    particles: Sequence[Particle],
    blocks: Sequence[Sequence[int]],
    known: Sequence[int],
    processed: Sequence[int],
    descriptors: Sequence[Descriptor],
    weight: int,
    channel: LeakageChannel,
    moves: int,
    rng: random.Random,
) -> tuple[list[Particle], int, int]:
    result = []
    accepted = 0
    attempted = 0
    for source in particles:
        particle = Particle(
            source.candidate,
            list(source.factor_logs),
            source.total_log,
        )
        for _ in range(moves):
            proposal = propose_coordinate(particle.candidate, blocks, rng)
            if proposal == particle.candidate:
                continue
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


def evaluate_factor_powers(
    known: Sequence[int],
    candidate: Candidate,
    powers: Sequence[float],
    descriptors: Sequence[Descriptor],
    weight: int,
    channel: LeakageChannel,
) -> tuple[list[float], float]:
    """Evaluate a product with one nonnegative exponent per factor."""

    factor_logs = [
        factor_log_likelihood(
            known, candidate, factor, descriptors, weight, channel
        )
        for factor in range(len(descriptors))
    ]
    if any(
        power > 0.0 and not math.isfinite(value)
        for power, value in zip(powers, factor_logs)
    ):
        return factor_logs, -math.inf
    return factor_logs, sum(
        power * value for power, value in zip(powers, factor_logs)
    )


def rejuvenate_powered(
    particles: Sequence[Particle],
    blocks: Sequence[Sequence[int]],
    known: Sequence[int],
    powers: Sequence[float],
    descriptors: Sequence[Descriptor],
    weight: int,
    channel: LeakageChannel,
    moves: int,
    rng: random.Random,
) -> tuple[list[Particle], int, int]:
    """Metropolis rejuvenation for a factorwise-powered target."""

    result = []
    accepted = 0
    attempted = 0
    for source in particles:
        particle = Particle(
            source.candidate,
            list(source.factor_logs),
            source.total_log,
        )
        for _ in range(moves):
            proposal = propose_coordinate(particle.candidate, blocks, rng)
            if proposal == particle.candidate:
                continue
            attempted += 1
            factor_logs, total_log = evaluate_factor_powers(
                known,
                proposal,
                powers,
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


def smc_target_particles(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    blocks: Sequence[Sequence[int]],
    weight: int,
    channel: LeakageChannel,
    particles_count: int,
    moves_per_stage: int,
    rng: random.Random,
) -> tuple[SmcResult, list[Particle]]:
    factors = list(range(len(descriptors)))
    rng.shuffle(factors)
    particles = [
        Particle(sample_from_blocks(blocks, rng), [], 0.0)
        for _ in range(particles_count)
    ]
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
            len({particle.candidate for particle in particles})
            / particles_count
        )
        processed.append(factor)
        particles, accepted, attempted = rejuvenate(
            particles,
            blocks,
            known,
            processed,
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
            statistics.mean(ess_fractions),
            min(unique_fractions),
            total_accepted / total_attempted if total_attempted else 0.0,
        ),
        particles,
    )


def smc_continue_to_squared_target(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    blocks: Sequence[Sequence[int]],
    weight: int,
    channel: LeakageChannel,
    source_particles: Sequence[Particle],
    moves_per_stage: int,
    steps_per_factor: int,
    rng: random.Random,
) -> tuple[SmcResult, list[Particle]]:
    """Continue an unweighted W-posterior population to target W squared.

    A second copy of every descriptor likelihood is introduced one factor at
    a time.  The returned normalizing constant is log2(Z_2 / Z_1), not Z_2.
    Keeping this as a continuation avoids a much sharper prior-to-W^2 jump.
    """

    if not source_particles:
        return SmcResult(-math.inf, 0.0, 0.0, 0.0, 0.0), []
    if steps_per_factor < 1:
        raise ValueError("power steps must be positive")
    particles = [
        Particle(
            particle.candidate,
            list(particle.factor_logs),
            particle.total_log,
        )
        for particle in source_particles
    ]
    powers = [1.0] * len(descriptors)
    increment = 1.0 / steps_per_factor
    log_mean_increment = 0.0
    ess_fractions = []
    unique_fractions = []
    total_accepted = 0
    total_attempted = 0

    for _ in range(steps_per_factor):
        factors = list(range(len(descriptors)))
        rng.shuffle(factors)
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
                weighted_value = increment * value
                if (
                    math.isfinite(particle.total_log)
                    and math.isfinite(weighted_value)
                ):
                    particle.total_log += weighted_value
                else:
                    particle.total_log = -math.inf
                incremental_logs.append(weighted_value)

            weights = normalized_weights(incremental_logs)
            if not weights:
                return SmcResult(-math.inf, 0.0, 0.0, 0.0, 0.0), []
            log_mean_increment += logmeanexp2(incremental_logs)
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
                len({particle.candidate for particle in particles})
                / len(particles)
            )
            powers[factor] += increment
            particles, accepted, attempted = rejuvenate_powered(
                particles,
                blocks,
                known,
                powers,
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
            log_mean_increment,
            min(ess_fractions),
            statistics.mean(ess_fractions),
            min(unique_fractions),
            total_accepted / total_attempted if total_attempted else 0.0,
        ),
        particles,
    )


def smc_mean_weight(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    blocks: Sequence[Sequence[int]],
    weight: int,
    channel: LeakageChannel,
    particles_count: int,
    moves_per_stage: int,
    rng: random.Random,
) -> SmcResult:
    result, _ = smc_target_particles(
        known,
        descriptors,
        blocks,
        weight,
        channel,
        particles_count,
        moves_per_stage,
        rng,
    )
    return result


def smc_rectangle_gain(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    full_blocks: Sequence[Sequence[int]],
    selection: Sequence[Sequence[int]],
    weight: int,
    channel: LeakageChannel,
    particles_count: int,
    moves_per_stage: int,
    restriction_moves: int,
    refresh_threshold: float,
    rng: random.Random,
) -> SplitResult:
    posterior, particles = smc_target_particles(
        known,
        descriptors,
        full_blocks,
        weight,
        channel,
        particles_count,
        moves_per_stage,
        rng,
    )
    if not particles:
        return SplitResult(-math.inf, posterior, 0.0, 0.0, 0.0)

    variables = [
        variable
        for variable in range(weight)
        if len(selection[variable]) < len(full_blocks[variable])
    ]
    rng.shuffle(variables)
    allowed_blocks = [tuple(block) for block in full_blocks]
    all_factors = tuple(range(len(descriptors)))
    gain_bits = 0.0
    survival_fractions = []
    unique_fractions = []
    total_accepted = 0
    total_attempted = 0
    posterior_runs = [posterior]

    for variable in variables:
        current_unique = (
            len({particle.candidate for particle in particles})
            / particles_count
        )
        if current_unique < refresh_threshold:
            refreshed, particles = smc_target_particles(
                known,
                descriptors,
                allowed_blocks,
                weight,
                channel,
                particles_count,
                moves_per_stage,
                rng,
            )
            posterior_runs.append(refreshed)
            if not particles:
                return SplitResult(
                    -math.inf,
                    aggregate_posterior_diagnostics(posterior_runs),
                    0.0,
                    min(unique_fractions, default=0.0),
                    total_accepted / total_attempted
                    if total_attempted
                    else 0.0,
                )
        selected = set(selection[variable])
        survivors = [
            index
            for index, particle in enumerate(particles)
            if particle.candidate[variable] in selected
        ]
        if not survivors:
            refreshed, particles = smc_target_particles(
                known,
                descriptors,
                allowed_blocks,
                weight,
                channel,
                particles_count,
                moves_per_stage,
                rng,
            )
            posterior_runs.append(refreshed)
            survivors = [
                index
                for index, particle in enumerate(particles)
                if particle.candidate[variable] in selected
            ]
        survival = len(survivors) / particles_count
        survival_fractions.append(survival)
        if not survivors:
            return SplitResult(
                -math.inf,
                aggregate_posterior_diagnostics(posterior_runs),
                0.0,
                min(unique_fractions, default=0.0),
                total_accepted / total_attempted if total_attempted else 0.0,
            )
        baseline = len(selected) / len(full_blocks[variable])
        gain_bits += math.log2(survival / baseline)

        survivor_weights = [1.0 / len(survivors)] * len(survivors)
        source_indices = systematic_indices(
            survivor_weights, rng, particles_count
        )
        particles = [
            Particle(
                particles[survivors[index]].candidate,
                list(particles[survivors[index]].factor_logs),
                particles[survivors[index]].total_log,
            )
            for index in source_indices
        ]
        unique_fractions.append(
            len({particle.candidate for particle in particles})
            / particles_count
        )
        allowed_blocks[variable] = tuple(selection[variable])
        particles, accepted, attempted = rejuvenate(
            particles,
            allowed_blocks,
            known,
            all_factors,
            descriptors,
            weight,
            channel,
            restriction_moves,
            rng,
        )
        total_accepted += accepted
        total_attempted += attempted

    return SplitResult(
        gain_bits,
        aggregate_posterior_diagnostics(posterior_runs),
        min(survival_fractions),
        min(unique_fractions),
        total_accepted / total_attempted if total_attempted else 0.0,
    )


def aggregate_posterior_diagnostics(
    results: Sequence[SmcResult],
) -> SmcResult:
    return SmcResult(
        results[-1].log_mean_weight,
        min(result.minimum_ess_fraction for result in results),
        statistics.mean(result.mean_ess_fraction for result in results),
        min(result.minimum_unique_fraction for result in results),
        statistics.mean(result.mutation_acceptance for result in results),
    )


def combine_replicates(results: Sequence[SmcResult]) -> tuple[float, float]:
    values = [result.log_mean_weight for result in results]
    if not any(math.isfinite(value) for value in values):
        return -math.inf, math.inf
    estimate = logmeanexp2(values)
    if len(values) < 2:
        return estimate, 0.0
    relative = [
        0.0 if not math.isfinite(value) else 2.0 ** (value - estimate)
        for value in values
    ]
    relative_se = statistics.stdev(relative) / math.sqrt(len(relative))
    log_se = math.log2(1.0 + relative_se)
    return estimate, log_se


def combine_split_replicates(
    results: Sequence[SplitResult],
) -> tuple[float, float]:
    values = [result.gain_bits for result in results]
    if not any(math.isfinite(value) for value in values):
        return -math.inf, math.inf
    estimate = logmeanexp2(values)
    if len(values) < 2:
        return estimate, 0.0
    relative = [
        0.0 if not math.isfinite(value) else 2.0 ** (value - estimate)
        for value in values
    ]
    relative_se = statistics.stdev(relative) / math.sqrt(len(relative))
    return estimate, math.log2(1.0 + relative_se)


def score_space(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    blocks: Sequence[Sequence[int]],
    weight: int,
    channel: LeakageChannel,
    particles: int,
    moves: int,
    replicates: int,
    rng: random.Random,
) -> tuple[float, float, list[SmcResult]]:
    results = [
        smc_mean_weight(
            known,
            descriptors,
            blocks,
            weight,
            channel,
            particles,
            moves,
            random.Random(rng.getrandbits(128)),
        )
        for _ in range(replicates)
    ]
    estimate, log_se = combine_replicates(results)
    return estimate, log_se, results


def summarize_diagnostics(results: Sequence[SmcResult]) -> str:
    return (
        f"min_ess={min(result.minimum_ess_fraction for result in results):.4f} "
        f"mean_ess={statistics.mean(result.mean_ess_fraction for result in results):.4f} "
        f"min_unique={min(result.minimum_unique_fraction for result in results):.4f} "
        f"accept={statistics.mean(result.mutation_acceptance for result in results):.4f}"
    )


def score_rectangle_splitting(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    full_blocks: Sequence[Sequence[int]],
    selection: Sequence[Sequence[int]],
    weight: int,
    channel: LeakageChannel,
    particles: int,
    moves: int,
    restriction_moves: int,
    refresh_threshold: float,
    replicates: int,
    rng: random.Random,
) -> tuple[float, float, list[SplitResult]]:
    results = [
        smc_rectangle_gain(
            known,
            descriptors,
            full_blocks,
            selection,
            weight,
            channel,
            particles,
            moves,
            restriction_moves,
            refresh_threshold,
            random.Random(rng.getrandbits(128)),
        )
        for _ in range(replicates)
    ]
    estimate, log_se = combine_split_replicates(results)
    return estimate, log_se, results


def summarize_split_diagnostics(results: Sequence[SplitResult]) -> str:
    return (
        f"posterior_min_ess="
        f"{min(result.posterior.minimum_ess_fraction for result in results):.4f} "
        f"posterior_mean_ess="
        f"{statistics.mean(result.posterior.mean_ess_fraction for result in results):.4f} "
        f"min_survival="
        f"{min(result.minimum_survival_fraction for result in results):.4f} "
        f"min_unique="
        f"{min(result.minimum_unique_fraction for result in results):.4f} "
        f"accept="
        f"{statistics.mean(result.mutation_acceptance for result in results):.4f}"
    )


def exact_mean_weight(
    known: Sequence[int],
    descriptors: Sequence[Descriptor],
    blocks: Sequence[Sequence[int]],
    weight: int,
    channel: LeakageChannel,
) -> float:
    log_weights = []
    for candidate in itertools.product(*blocks):
        factor_logs = [
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
            if any(not math.isfinite(value) for value in factor_logs)
            else sum(factor_logs)
        )
    return logmeanexp2(log_weights)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--weight", type=int, default=3)
    parser.add_argument("--block-length", type=int, default=8)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--trials", type=int, default=4)
    parser.add_argument("--start-trial", type=int, default=0)
    parser.add_argument("--hidden-polys", type=int)
    parser.add_argument("--particles", type=int, default=256)
    parser.add_argument("--moves", type=int, default=2)
    parser.add_argument("--restriction-moves", type=int, default=8)
    parser.add_argument("--refresh-threshold", type=float, default=0.1)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--replay-replicates", type=int)
    parser.add_argument(
        "--score-mode",
        choices=("splitting", "partition"),
        default="splitting",
    )
    parser.add_argument("--pair-restarts", type=int, default=4)
    parser.add_argument("--pair-rounds", type=int, default=24)
    parser.add_argument("--pair-damping", type=float, default=0.5)
    parser.add_argument("--rectangle-restarts", type=int, default=4)
    parser.add_argument("--rectangle-rounds", type=int, default=16)
    parser.add_argument("--matching-contexts", type=int, default=4)
    parser.add_argument("--matching-rounds", type=int, default=0)
    parser.add_argument("--joint-contexts", type=int, default=0)
    parser.add_argument("--joint-particles", type=int)
    parser.add_argument("--joint-moves", type=int)
    parser.add_argument("--joint-restriction-moves", type=int)
    parser.add_argument("--pair-slope", type=float)
    parser.add_argument("--pair-fit-samples", type=int, default=20000)
    parser.add_argument(
        "--channel", choices=("raw", "occupancy"), default="occupancy"
    )
    parser.add_argument("--selector-k", type=int, default=4)
    parser.add_argument(
        "--hash-family",
        choices=("ideal", "goldreich"),
        default="ideal",
    )
    parser.add_argument("--linear-security", type=int, default=10)
    parser.add_argument("--exact-rank-limit", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--exact-check", action="store_true")
    args = parser.parse_args()

    if min(
        args.polys,
        args.weight,
        args.block_length,
        args.trials,
        args.particles,
        args.replicates,
        args.pair_rounds,
        args.rectangle_rounds,
    ) < 1:
        parser.error("all size and sampling parameters must be positive")
    if args.moves < 0:
        parser.error("moves must be nonnegative")
    if args.rectangle_restarts < 0:
        parser.error("rectangle-restarts must be nonnegative")
    if args.matching_contexts < 1:
        parser.error("matching-contexts must be positive")
    if args.matching_rounds < 0:
        parser.error("matching-rounds must be nonnegative")
    if args.joint_contexts < 0:
        parser.error("joint-contexts must be nonnegative")
    joint_particles = args.joint_particles or args.particles
    joint_moves = args.joint_moves if args.joint_moves is not None else args.moves
    joint_restriction_moves = (
        args.joint_restriction_moves
        if args.joint_restriction_moves is not None
        else args.restriction_moves
    )
    if joint_particles < 1:
        parser.error("joint-particles must be positive")
    if min(joint_moves, joint_restriction_moves) < 0:
        parser.error("joint selector moves must be nonnegative")
    if args.restriction_moves < 0:
        parser.error("restriction-moves must be nonnegative")
    if not 0.0 <= args.refresh_threshold <= 1.0:
        parser.error("refresh-threshold must lie between zero and one")
    replay_replicates = args.replay_replicates or args.replicates
    if replay_replicates < 1:
        parser.error("replay-replicates must be positive")
    if not 0 <= args.start_trial < args.trials:
        parser.error("start-trial must lie between zero and trials-1")
    if args.block_length % args.polys:
        parser.error("block length must be divisible by the polynomial count")
    if args.linear_security < 0:
        parser.error("linear-security must be nonnegative")
    if args.exact_rank_limit < 0:
        parser.error("exact-rank-limit must be nonnegative")
    hidden_polys = args.hidden_polys or args.polys
    if not 1 <= hidden_polys <= args.polys:
        parser.error("hidden-polys must lie between one and polys")

    partition_size = args.partition_size or args.weight
    base_channel = LeakageChannel.build(
        args.channel,
        args.selector_k,
        args.weight,
        partition_size,
    )
    channel = (
        SharedGoldreichChannel(
            base_channel,
            linear_security=args.linear_security,
            exact_rank_limit=args.exact_rank_limit,
        )
        if args.hash_family == "goldreich"
        else base_channel
    )
    fit_rng = random.Random(args.seed ^ 0x4B1D)
    if args.pair_slope is None:
        pair_slope, pair_r2 = fit_pair_slope(
            channel, args.weight, args.pair_fit_samples, fit_rng
        )
    else:
        pair_slope = args.pair_slope
        pair_r2 = math.nan

    rng = random.Random(args.seed)
    kept = args.block_length // args.polys
    trial_gains = []
    for trial in range(args.trials):
        known = tuple(
            rng.randrange(args.block_length)
            for _ in range(args.polys * args.weight)
        )
        trial_gain = 0.0
        for hidden_poly in range(hidden_polys):
            true_candidate = sample_uniform_candidate(
                args.weight, args.block_length, rng
            )
            true_lists = support_lists_one_unknown(
                known, true_candidate, args.polys, args.weight
            )
            if isinstance(channel, SharedGoldreichChannel):
                descriptors = channel.sample_batch(
                    true_lists, 2 * args.block_length, rng
                )
            else:
                descriptors = tuple(
                    channel.sample(active, 2 * args.block_length, rng)
                    for active in true_lists
                )
            distributions = pair_mean_field_distributions(
                known,
                descriptors,
                args.polys,
                args.weight,
                args.block_length,
                kept,
                pair_slope,
                args.pair_restarts,
                args.pair_rounds,
                args.pair_damping,
                rng,
            )
            selection = pair_refined_rectangle_selection(
                known,
                descriptors,
                args.polys,
                args.weight,
                args.block_length,
                kept,
                pair_slope,
                distributions,
                args.rectangle_restarts,
                args.rectangle_rounds,
                rng,
            )
            if args.matching_rounds:
                selection = matching_refined_rectangle_selection(
                    known,
                    descriptors,
                    selection,
                    args.weight,
                    args.block_length,
                    channel,
                    args.matching_contexts,
                    args.matching_rounds,
                    rng,
                )
            pre_joint_selection = selection
            selector_diagnostic = ""
            if args.joint_contexts:
                selection, joint_diagnostic = (
                    joint_conditional_rectangle_selection(
                        known,
                        descriptors,
                        args.weight,
                        args.block_length,
                        kept,
                        channel,
                        joint_particles,
                        joint_moves,
                        joint_restriction_moves,
                        args.joint_contexts,
                        random.Random(rng.getrandbits(128)),
                    )
                )
                selector_diagnostic = (
                    f" joint_selector_min_ess="
                    f"{joint_diagnostic.minimum_ess_fraction:.4f}"
                    f" joint_selector_mean_ess="
                    f"{joint_diagnostic.mean_ess_fraction:.4f}"
                )
            full_blocks = tuple(
                tuple(range(args.block_length)) for _ in range(args.weight)
            )
            rectangle_blocks = tuple(tuple(block) for block in selection)
            if trial < args.start_trial:
                seed_groups = 1 if args.score_mode == "splitting" else 2
                for _ in range(seed_groups * replay_replicates):
                    rng.getrandbits(128)
                continue
            if args.score_mode == "splitting":
                gain, gain_se, split_results = score_rectangle_splitting(
                    known,
                    descriptors,
                    full_blocks,
                    rectangle_blocks,
                    args.weight,
                    channel,
                    args.particles,
                    args.moves,
                    args.restriction_moves,
                    args.refresh_threshold,
                    args.replicates,
                    rng,
                )
                diagnostic_text = (
                    f"gain_se={gain_se:.4f} "
                    f"split[{summarize_split_diagnostics(split_results)}]"
                )
            else:
                full, full_se, full_results = score_space(
                    known,
                    descriptors,
                    full_blocks,
                    args.weight,
                    channel,
                    args.particles,
                    args.moves,
                    args.replicates,
                    rng,
                )
                rectangle, rectangle_se, rectangle_results = score_space(
                    known,
                    descriptors,
                    rectangle_blocks,
                    args.weight,
                    channel,
                    args.particles,
                    args.moves,
                    args.replicates,
                    rng,
                )
                gain = rectangle - full
                diagnostic_text = (
                    f"full_se={full_se:.4f} rectangle_se={rectangle_se:.4f} "
                    f"full[{summarize_diagnostics(full_results)}] "
                    f"rectangle[{summarize_diagnostics(rectangle_results)}]"
                )
            trial_gain += gain
            exact_text = ""
            if args.exact_check:
                exact_full = exact_mean_weight(
                    known,
                    descriptors,
                    full_blocks,
                    args.weight,
                    channel,
                )
                exact_rectangle = exact_mean_weight(
                    known,
                    descriptors,
                    rectangle_blocks,
                    args.weight,
                    channel,
                )
                exact_gain = exact_rectangle - exact_full
                pre_joint_text = ""
                if args.joint_contexts:
                    pre_joint_rectangle = exact_mean_weight(
                        known,
                        descriptors,
                        tuple(tuple(block) for block in pre_joint_selection),
                        args.weight,
                        channel,
                    )
                    pre_joint_gain = pre_joint_rectangle - exact_full
                    pre_joint_text = (
                        f" pre_joint_exact_gain={pre_joint_gain:.6f}"
                    )
                exact_text = (
                    f" exact_gain={exact_gain:.6f} "
                    f"error={gain - exact_gain:.6f}"
                    f"{pre_joint_text}"
                )
            print(
                f"trial={trial} poly={hidden_poly} gain={gain:.6f} "
                f"{diagnostic_text}"
                f"{selector_diagnostic}"
                f"{exact_text}"
            )
        if trial >= args.start_trial:
            trial_gains.append(trial_gain)
            print(f"trial={trial} total_gain={trial_gain:.6f}")

    print(
        f"polys={args.polys} weight={args.weight} block_length={args.block_length} "
        f"d={partition_size} channel={args.channel} selector_k={args.selector_k} "
        f"hash_family={args.hash_family} "
        f"linear_security={args.linear_security} "
        f"exact_rank_limit={args.exact_rank_limit} "
        f"pair_slope={pair_slope:.9f} pair_r2={pair_r2:.9f} "
        f"hidden_polys={hidden_polys} "
        f"score_mode={args.score_mode} "
        f"particles={args.particles} moves={args.moves} "
        f"restriction_moves={args.restriction_moves} "
        f"joint_contexts={args.joint_contexts} "
        f"joint_particles={joint_particles} "
        f"joint_moves={joint_moves} "
        f"joint_restriction_moves={joint_restriction_moves} "
        f"refresh_threshold={args.refresh_threshold:.4f} "
        f"replicates={args.replicates} trials={args.trials} seed={args.seed}"
    )
    if isinstance(channel, SharedGoldreichChannel):
        print(
            f"goldreich_feature_bits={channel.feature_bits} "
            f"sampled_systems={channel.sampled_systems} "
            f"rank_deficient_systems={channel.rank_deficient_systems} "
            f"rank_failures={channel.rank_failures} "
            f"approximate_rank_queries={channel.approximate_rank_queries}"
        )
    print(
        f"gain_bits_logmean={logmeanexp2(trial_gains):.6f} "
        f"gain_bits_meanlog={statistics.mean(trial_gains):.6f}"
    )


if __name__ == "__main__":
    main()
