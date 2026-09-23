"""Charged posterior attacks for multimodal Reverse-Cuckoo leakage.

The rectangle mode builds and charges a portfolio of mode-specific rectangles.
The offset mode selects an exact-width value set from a posterior histogram
and scores it with an independent posterior population and exact full-batch
conditionals.  Independent scoring prevents a maximum over noisy training
estimates from being reported as attack gain.
"""

from __future__ import annotations

import argparse
import math
import random

from rev_cuckoo_goldreich import SharedGoldreichChannel
from rev_cuckoo_regular_isd_mc import logmeanexp2, sample_uniform_candidate
from rev_cuckoo_regular_isd_smc import (
    contextual_rectangle_portfolio,
    exact_mean_weight,
    posterior_single_offset_set,
    score_rectangle_splitting,
    summarize_split_diagnostics,
)
from rev_cuckoo_regular_noise import LeakageChannel, support_lists_one_unknown


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=2)
    parser.add_argument("--weight", type=int, default=3)
    parser.add_argument("--block-length", type=int, default=8)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument("--portfolio-size", type=int, default=4)
    parser.add_argument(
        "--offset-set",
        "--guess-one",
        dest="offset_set",
        action="store_true",
        help="run the one-offset value-set attack",
    )
    parser.add_argument("--guess-contexts", type=int, default=16)
    parser.add_argument("--guess-width", type=int, default=1)
    parser.add_argument("--selector-particles", type=int, default=256)
    parser.add_argument("--selector-moves", type=int, default=2)
    parser.add_argument("--particles", type=int, default=256)
    parser.add_argument("--moves", type=int, default=2)
    parser.add_argument("--restriction-moves", type=int, default=8)
    parser.add_argument("--refresh-threshold", type=float, default=0.2)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument(
        "--channel", choices=("raw", "occupancy"), default="occupancy"
    )
    parser.add_argument("--selector-k", type=int)
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
        args.portfolio_size,
        args.guess_contexts,
        args.selector_particles,
        args.particles,
        args.replicates,
    ) < 1:
        parser.error("all dimensions and population sizes must be positive")
    if min(
        args.selector_moves, args.moves, args.restriction_moves
    ) < 0:
        parser.error("move counts must be nonnegative")
    if args.block_length % args.polys:
        parser.error("block length must be divisible by polynomial count")
    if not 1 <= args.guess_width <= args.block_length:
        parser.error("guess width must lie in the offset domain")
    if not 0.0 <= args.refresh_threshold <= 1.0:
        parser.error("refresh-threshold must lie between zero and one")
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

    kept = args.block_length // args.polys
    guess_text = ""
    guess_gain = math.nan
    if args.offset_set:
        (
            selection,
            guess,
            guess_probability,
            selector_diagnostic,
            evaluator_diagnostic,
        ) = posterior_single_offset_set(
            known,
            descriptors,
            args.weight,
            args.block_length,
            channel,
            args.selector_particles,
            args.selector_moves,
            args.guess_contexts,
            args.guess_width,
            random.Random(rng.getrandbits(128)),
        )
        rectangles = (selection,)
        guess_gain = math.log2(
            args.block_length * guess_probability / args.guess_width
        )
        guess_text = (
            f" guessed_variable={guess[0]}"
            f" guessed_width={len(guess[1])}"
            f" guessed_probability={guess_probability:.8g}"
            f" evaluator_min_ess="
            f"{evaluator_diagnostic.minimum_ess_fraction:.4f}"
            f" evaluator_mean_ess="
            f"{evaluator_diagnostic.mean_ess_fraction:.4f}"
            f" selector_unique="
            f"{selector_diagnostic.minimum_unique_fraction:.4f}"
            f" evaluator_unique="
            f"{evaluator_diagnostic.minimum_unique_fraction:.4f}"
            f" selector_accept="
            f"{selector_diagnostic.mutation_acceptance:.4f}"
            f" evaluator_accept="
            f"{evaluator_diagnostic.mutation_acceptance:.4f}"
        )
    else:
        rectangles, selector_diagnostic = contextual_rectangle_portfolio(
            known,
            descriptors,
            args.weight,
            args.block_length,
            kept,
            channel,
            args.portfolio_size,
            args.selector_particles,
            args.selector_moves,
            random.Random(rng.getrandbits(128)),
        )
    if not rectangles:
        raise RuntimeError("portfolio selector produced no rectangle")
    full_blocks = tuple(
        tuple(range(args.block_length)) for _ in range(args.weight)
    )
    exact_full = (
        exact_mean_weight(
            known, descriptors, full_blocks, args.weight, channel
        )
        if args.exact_check
        else math.nan
    )

    estimated_gains = []
    exact_gains = []
    for index, selection in enumerate(rectangles):
        rectangle = tuple(tuple(block) for block in selection)
        if args.offset_set:
            estimated_gains.append(guess_gain)
            exact_text = ""
            if args.exact_check:
                exact_gain = (
                    exact_mean_weight(
                        known,
                        descriptors,
                        rectangle,
                        args.weight,
                        channel,
                    )
                    - exact_full
                )
                exact_gains.append(exact_gain)
                exact_text = f" exact_gain={exact_gain:.6f}"
            print(
                f"guess={index} gain={guess_gain:.6f} "
                f"posterior_probability={guess_probability:.8g}"
                f"{exact_text}"
            )
            continue
        gain, gain_se, results = score_rectangle_splitting(
            known,
            descriptors,
            full_blocks,
            rectangle,
            args.weight,
            channel,
            args.particles,
            args.moves,
            args.restriction_moves,
            args.refresh_threshold,
            args.replicates,
            rng,
        )
        estimated_gains.append(gain)
        exact_text = ""
        if args.exact_check:
            exact_gain = (
                exact_mean_weight(
                    known, descriptors, rectangle, args.weight, channel
                )
                - exact_full
            )
            exact_gains.append(exact_gain)
            exact_text = f" exact_gain={exact_gain:.6f}"
        print(
            f"rectangle={index} gain={gain:.6f} gain_se={gain_se:.4f} "
            f"split[{summarize_split_diagnostics(results)}]"
            f"{exact_text}"
        )

    charged_gain = logmeanexp2(estimated_gains)
    print(
        f"portfolio_size={len(rectangles)} charged_gain={charged_gain:.6f} "
        f"selector_min_ess={selector_diagnostic.minimum_ess_fraction:.4f} "
        f"selector_mean_ess={selector_diagnostic.mean_ess_fraction:.4f}"
        f"{guess_text}"
    )
    if args.exact_check:
        print(f"exact_charged_gain={logmeanexp2(exact_gains):.6f}")
    print(
        f"polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} d={partition_size} "
        f"channel={args.channel} "
        f"selector_k={selector_k if args.channel == 'occupancy' else '-'} "
        f"hash_family={args.hash_family} particles={args.particles} "
        f"selector_particles={args.selector_particles} "
        f"replicates={args.replicates} seed={args.seed}"
    )
    if isinstance(channel, SharedGoldreichChannel):
        print(
            f"goldreich_feature_bits={channel.feature_bits} "
            f"sampled_systems={channel.sampled_systems} "
            f"rank_deficient_systems={channel.rank_deficient_systems} "
            f"rank_failures={channel.rank_failures} "
            f"approximate_rank_queries={channel.approximate_rank_queries}"
        )


if __name__ == "__main__":
    main()
