"""Posterior information diagnostic for Reverse-Cuckoo leakage.

For one sampled transcript, estimate

    KL(p(x | leakage) || uniform)
      = E_post[log2 W(x)] - log2(E_uniform[W(x)]).

The quantity is posterior entropy loss, not an attack exponent.  Independent
SMC replicates expose evidence instability and mode collapse.  Small instances
can be checked by exact enumeration.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
import statistics

from rev_cuckoo_goldreich import SharedGoldreichChannel
from rev_cuckoo_regular_isd_mc import logmeanexp2, sample_uniform_candidate
from rev_cuckoo_regular_isd_smc import (
    factor_log_likelihood,
    smc_continue_to_squared_target,
    smc_target_particles,
)
from rev_cuckoo_regular_noise import LeakageChannel, support_lists_one_unknown


def exact_posterior_information(
    known: tuple[int, ...],
    descriptors: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
    weight: int,
    block_length: int,
    channel: LeakageChannel,
) -> tuple[float, float, float, float]:
    """Return exact evidence, posterior mean, KL, and Renyi-2 loss."""

    logs = []
    for candidate in itertools.product(range(block_length), repeat=weight):
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
        logs.append(
            -math.inf
            if any(not math.isfinite(value) for value in values)
            else sum(values)
        )
    log_evidence = logmeanexp2(logs)
    finite = [value for value in logs if math.isfinite(value)]
    maximum = max(finite)
    weights = [
        0.0 if not math.isfinite(value) else 2.0 ** (value - maximum)
        for value in logs
    ]
    total = sum(weights)
    posterior_mean_log_weight = sum(
        weight * value
        for weight, value in zip(weights, logs)
        if weight
    ) / total
    return (
        log_evidence,
        posterior_mean_log_weight,
        posterior_mean_log_weight - log_evidence,
        logmeanexp2(
            [
                2.0 * value if math.isfinite(value) else -math.inf
                for value in logs
            ]
        )
        - 2.0 * log_evidence,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=2)
    parser.add_argument("--weight", type=int, default=3)
    parser.add_argument("--block-length", type=int, default=8)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--particles", type=int, default=1024)
    parser.add_argument("--moves", type=int, default=4)
    parser.add_argument("--power-moves", type=int)
    parser.add_argument("--power-steps", type=int, default=1)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--replicate-offset", type=int, default=0)
    parser.add_argument(
        "--channel", choices=("raw", "occupancy"), default="occupancy"
    )
    parser.add_argument("--selector-k", type=int)
    parser.add_argument(
        "--hash-family", choices=("ideal", "goldreich"), default="goldreich"
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
        args.particles,
        args.replicates,
    ) < 1:
        parser.error("dimensions and population sizes must be positive")
    if args.replicate_offset < 0:
        parser.error("replicate offset must be nonnegative")
    if args.moves < 0:
        parser.error("moves must be nonnegative")
    power_moves = args.moves if args.power_moves is None else args.power_moves
    if power_moves < 0:
        parser.error("power moves must be nonnegative")
    if args.power_steps < 1:
        parser.error("power steps must be positive")
    if args.exact_check and args.block_length ** args.weight > 2_000_000:
        parser.error("exact candidate space is too large")
    selector_k = (
        args.selector_k
        if args.selector_k is not None
        else (4 if args.channel == "occupancy" else 1)
    )
    if args.channel == "raw" and selector_k != 1:
        parser.error("selector K is undefined for the raw channel")

    partition_size = args.partition_size or args.weight
    base = LeakageChannel.build(
        args.channel,
        selector_k,
        args.weight,
        partition_size,
    )
    channel = (
        SharedGoldreichChannel(
            base,
            linear_security=args.linear_security,
            exact_rank_limit=args.exact_rank_limit,
        )
        if args.hash_family == "goldreich"
        else base
    )
    rng = random.Random(args.seed)
    known = tuple(
        rng.randrange(args.block_length)
        for _ in range(args.polys * args.weight)
    )
    true_candidate = sample_uniform_candidate(
        args.weight, args.block_length, rng
    )
    active_sets = support_lists_one_unknown(
        known, true_candidate, args.polys, args.weight
    )
    domain = 2 * args.block_length
    if isinstance(channel, SharedGoldreichChannel):
        descriptors = channel.sample_batch(active_sets, domain, rng)
    else:
        descriptors = tuple(
            channel.sample(active, domain, rng) for active in active_sets
        )
    blocks = tuple(
        tuple(range(args.block_length)) for _ in range(args.weight)
    )

    for _ in range(args.replicate_offset):
        rng.getrandbits(128)
        rng.getrandbits(128)

    estimates = []
    renyi2_particle_estimates = []
    renyi2_tempered_estimates = []
    prior_entropy = args.weight * math.log2(args.block_length)
    for local_replicate in range(args.replicates):
        replicate = args.replicate_offset + local_replicate
        diagnostic, particles = smc_target_particles(
            known,
            descriptors,
            blocks,
            args.weight,
            channel,
            args.particles,
            args.moves,
            random.Random(rng.getrandbits(128)),
        )
        if not particles:
            print(f"replicate={replicate} failed=1")
            continue
        posterior_mean_log_weight = statistics.fmean(
            particle.total_log for particle in particles
        )
        information = posterior_mean_log_weight - diagnostic.log_mean_weight
        renyi2_particle_information = (
            logmeanexp2([particle.total_log for particle in particles])
            - diagnostic.log_mean_weight
        )
        power_diagnostic, power_particles = smc_continue_to_squared_target(
            known,
            descriptors,
            blocks,
            args.weight,
            channel,
            particles,
            power_moves,
            args.power_steps,
            random.Random(rng.getrandbits(128)),
        )
        if not power_particles:
            print(f"replicate={replicate} power_failed=1")
            continue
        renyi2_tempered_information = (
            power_diagnostic.log_mean_weight - diagnostic.log_mean_weight
        )
        estimates.append(information)
        renyi2_particle_estimates.append(renyi2_particle_information)
        renyi2_tempered_estimates.append(renyi2_tempered_information)
        final_unique = len({particle.candidate for particle in particles}) / len(
            particles
        )
        print(
            f"replicate={replicate} information={information:.6f} "
            f"renyi2_particle={renyi2_particle_information:.6f} "
            f"renyi2_tempered={renyi2_tempered_information:.6f} "
            f"shannon_remaining={prior_entropy - information:.6f} "
            f"collision_remaining="
            f"{prior_entropy - renyi2_tempered_information:.6f} "
            f"log_evidence={diagnostic.log_mean_weight:.6f} "
            f"posterior_mean_log_weight={posterior_mean_log_weight:.6f} "
            f"min_ess={diagnostic.minimum_ess_fraction:.4f} "
            f"mean_ess={diagnostic.mean_ess_fraction:.4f} "
            f"min_unique={diagnostic.minimum_unique_fraction:.4f} "
            f"final_unique={final_unique:.4f} "
            f"accept={diagnostic.mutation_acceptance:.4f} "
            f"power_min_ess={power_diagnostic.minimum_ess_fraction:.4f} "
            f"power_mean_ess={power_diagnostic.mean_ess_fraction:.4f} "
            f"power_min_unique={power_diagnostic.minimum_unique_fraction:.4f} "
            f"power_accept={power_diagnostic.mutation_acceptance:.4f}"
        )
    if not estimates:
        raise RuntimeError("all posterior SMC replicates failed")
    spread = max(estimates) - min(estimates)
    print(
        f"information_mean={statistics.fmean(estimates):.6f} "
        f"information_spread={spread:.6f} "
        f"renyi2_particle_mean="
        f"{statistics.fmean(renyi2_particle_estimates):.6f} "
        f"renyi2_tempered_mean="
        f"{statistics.fmean(renyi2_tempered_estimates):.6f} "
        f"renyi2_tempered_spread="
        f"{max(renyi2_tempered_estimates) - min(renyi2_tempered_estimates):.6f} "
        f"shannon_remaining_mean="
        f"{prior_entropy - statistics.fmean(estimates):.6f} "
        f"collision_remaining_mean="
        f"{prior_entropy - statistics.fmean(renyi2_tempered_estimates):.6f} "
        f"successful={len(estimates)}"
    )
    if args.exact_check:
        (
            log_evidence,
            posterior_mean,
            information,
            renyi2_information,
        ) = exact_posterior_information(
            known, descriptors, args.weight, args.block_length, channel
        )
        print(
            f"exact_information={information:.6f} "
            f"exact_renyi2_information={renyi2_information:.6f} "
            f"exact_shannon_remaining={prior_entropy - information:.6f} "
            f"exact_collision_remaining="
            f"{prior_entropy - renyi2_information:.6f} "
            f"exact_log_evidence={log_evidence:.6f} "
            f"exact_posterior_mean_log_weight={posterior_mean:.6f}"
        )
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} d={partition_size} "
        f"channel={args.channel} "
        f"selector_k={selector_k if args.channel == 'occupancy' else '-'} "
        f"hash_family={args.hash_family} particles={args.particles} "
        f"moves={args.moves} power_moves={power_moves} "
        f"power_steps={args.power_steps} "
        f"replicates={args.replicates} seed={args.seed}"
    )


if __name__ == "__main__":
    main()
