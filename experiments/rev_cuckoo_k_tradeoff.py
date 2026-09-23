"""Join finite-K Reverse-Cuckoo leakage estimates with the OT cost model.

The uniform descriptor distribution U is sampled as a labeled two-choice
multigraph.  Its placement count Z and endpoint occupancy deficit D determine
the distributions considered here:

* raw: ordinary Reverse Cuckoo, with density proportional to Z;
* exact-K: draw K raw proposals and select with weight 1/Z;
* occupancy-K: select with a public lookup-table weight 1 / E_U[Z | D];
* moments-K: use 1 / E_U[Z | D,C_2], where C_2 is the total number of
  same-partition colliding endpoint pairs.
* feasible: the limiting ideal U conditioned on Z > 0.

For finite K, selection is simulated only over the small empirical (Z,D)
histogram.  Rao-Blackwellization adds each proposal's conditional selection
probability instead of sampling the final categorical choice.

The reported search-gain diagnostic is the same independent-list likelihood
ranking model used in rev_cuckoo_leakage.md.  It is not a reduction in the
bit-security of Ring-LPN or a leakage-aware decoding attack.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence

from rev_cuckoo_cost import Model as CostModel
from rev_cuckoo_cost import estimate as estimate_cost


Feature = tuple[int, int, int]


@dataclass(frozen=True)
class Leakage:
    kl_bits: float
    tv: float
    pair_chernoff_bits: float
    search_gain_bits: float


def falling_factorial(value: int, count: int) -> int:
    return math.prod(range(value - count + 1, value + 1))


def expected_placements(items: int, bins: tuple[int, int]) -> float:
    left_bins, right_bins = bins
    return sum(
        math.comb(items, left_count)
        * falling_factorial(left_bins, left_count)
        * falling_factorial(right_bins, items - left_count)
        / left_bins**left_count
        / right_bins ** (items - left_count)
        for left_count in range(items + 1)
        if left_count <= left_bins and items - left_count <= right_bins
    )


def graph_feature(
    items: int, bins: tuple[int, int], rng: random.Random
) -> Feature:
    """Return (Z(H), endpoint occupancy deficit, colliding endpoint pairs)."""

    left_bins, right_bins = bins
    parent = list(range(left_bins + right_bins))
    incident = 0
    edges: list[tuple[int, int]] = []
    left_occupied = 0
    right_occupied = 0
    left_loads = [0] * left_bins
    right_loads = [0] * right_bins

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for _ in range(items):
        left = rng.randrange(left_bins)
        right_bin = rng.randrange(right_bins)
        right = left_bins + right_bin
        edges.append((left, right))
        left_occupied |= 1 << left
        right_occupied |= 1 << right_bin
        left_loads[left] += 1
        right_loads[right_bin] += 1
        incident |= (1 << left) | (1 << right)
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    vertex_counts: dict[int, int] = defaultdict(int)
    edge_counts: dict[int, int] = defaultdict(int)
    bits = incident
    while bits:
        bit = bits & -bits
        vertex = bit.bit_length() - 1
        vertex_counts[find(vertex)] += 1
        bits ^= bit
    for left, _ in edges:
        edge_counts[find(left)] += 1

    placements = 1
    for root, edge_count in edge_counts.items():
        vertices = vertex_counts[root]
        if edge_count == vertices - 1:
            placements *= vertices
        elif edge_count == vertices:
            placements *= 2
        else:
            placements = 0
            break

    deficit = 2 * items - left_occupied.bit_count() - right_occupied.bit_count()
    collision_pairs = sum(load * (load - 1) // 2 for load in left_loads)
    collision_pairs += sum(load * (load - 1) // 2 for load in right_loads)
    return placements, deficit, collision_pairs


def sample_uniform_histogram(
    items: int, bins: tuple[int, int], trials: int, rng: random.Random
) -> Counter[Feature]:
    return Counter(graph_feature(items, bins, rng) for _ in range(trials))


def normalized(values: Mapping[Feature, float]) -> dict[Feature, float]:
    total = sum(values.values())
    if total <= 0:
        raise ValueError("distribution has zero mass")
    return {feature: mass / total for feature, mass in values.items() if mass > 0}


def proposal_distribution(
    uniform: Mapping[Feature, float]
) -> dict[Feature, float]:
    return normalized(
        {feature: mass * feature[0] for feature, mass in uniform.items() if feature[0]}
    )


def conditional_inverse_mean_weights(
    uniform: Mapping[Feature, float],
    key: Callable[[Feature], object],
) -> dict[object, float]:
    mass: dict[object, float] = defaultdict(float)
    placement_mass: dict[object, float] = defaultdict(float)
    for feature, probability in uniform.items():
        placements = feature[0]
        if placements:
            score = key(feature)
            mass[score] += probability
            placement_mass[score] += probability * placements
    return {
        score: mass[score] / placement_mass[score]
        for score in mass
    }


def occupancy_weights(uniform: Mapping[Feature, float]) -> dict[object, float]:
    return conditional_inverse_mean_weights(uniform, lambda feature: feature[1])


def occupied_bin_distribution(
    bins: int,
    random_draws: int,
    initially_occupied: int = 0,
) -> list[float]:
    """Distribution of occupied bins after uniform draws.

    ``initially_occupied`` bins are fixed and distinct.  Each subsequent draw
    is independent and uniform in the complete bin set.
    """

    if not 0 <= initially_occupied <= bins:
        raise ValueError("bad initial occupancy")
    if random_draws < 0:
        raise ValueError("random draws must be nonnegative")
    probability = [0.0] * (bins + 1)
    probability[initially_occupied] = 1.0
    for _ in range(random_draws):
        following = [0.0] * (bins + 1)
        for occupied, mass in enumerate(probability):
            if not mass:
                continue
            following[occupied] += mass * occupied / bins
            if occupied < bins:
                following[occupied + 1] += mass * (bins - occupied) / bins
        probability = following
    return probability


def exact_deficit_masses(
    items: int,
    bins: tuple[int, int],
) -> tuple[dict[int, float], dict[int, float]]:
    """Return ``(Pr_U[D], E_U[Z 1_D])`` exactly up to floating arithmetic."""

    left_bins, right_bins = bins
    left_uniform = occupied_bin_distribution(left_bins, items)
    right_uniform = occupied_bin_distribution(right_bins, items)
    deficit_mass: dict[int, float] = defaultdict(float)
    for left_occupied, left_mass in enumerate(left_uniform):
        for right_occupied, right_mass in enumerate(right_uniform):
            deficit = 2 * items - left_occupied - right_occupied
            deficit_mass[deficit] += left_mass * right_mass

    placement_mass: dict[int, float] = defaultdict(float)
    for left_count in range(items + 1):
        right_count = items - left_count
        if left_count > left_bins or right_count > right_bins:
            continue
        compatible_placements = (
            math.comb(items, left_count)
            * falling_factorial(left_bins, left_count)
            * falling_factorial(right_bins, right_count)
        )
        compatibility_probability = (
            left_bins ** (-left_count)
            * right_bins ** (-right_count)
        )
        coefficient = compatible_placements * compatibility_probability

        # Conditioned on this fixed placement, the planted endpoints are
        # distinct.  Only the endpoints on the opposite side remain random.
        left_conditioned = occupied_bin_distribution(
            left_bins,
            right_count,
            left_count,
        )
        right_conditioned = occupied_bin_distribution(
            right_bins,
            left_count,
            right_count,
        )
        for left_occupied, left_mass in enumerate(left_conditioned):
            for right_occupied, right_mass in enumerate(right_conditioned):
                deficit = 2 * items - left_occupied - right_occupied
                placement_mass[deficit] += (
                    coefficient * left_mass * right_mass
                )
    return dict(deficit_mass), dict(placement_mass)


def exact_occupancy_weights(
    items: int,
    bins: tuple[int, int],
) -> dict[int, float]:
    """Compute ``1/E_U[Z|D]`` without a sampled graph histogram."""

    deficit_mass, placement_mass = exact_deficit_masses(items, bins)
    return {
        deficit: mass / placement_mass[deficit]
        for deficit, mass in deficit_mass.items()
        if mass and placement_mass.get(deficit, 0.0)
    }


def _adaptive_simpson(
    function: Callable[[float], float],
    left: float,
    right: float,
    tolerance: float,
    whole: float,
    depth: int,
) -> float:
    midpoint = (left + right) / 2.0
    left_midpoint = (left + midpoint) / 2.0
    right_midpoint = (midpoint + right) / 2.0
    left_area = (midpoint - left) * (
        function(left) + 4.0 * function(left_midpoint) + function(midpoint)
    ) / 6.0
    right_area = (right - midpoint) * (
        function(midpoint) + 4.0 * function(right_midpoint) + function(right)
    ) / 6.0
    refined = left_area + right_area
    if depth == 0 or abs(refined - whole) <= 15.0 * tolerance:
        return refined + (refined - whole) / 15.0
    return _adaptive_simpson(
        function,
        left,
        midpoint,
        tolerance / 2.0,
        left_area,
        depth - 1,
    ) + _adaptive_simpson(
        function,
        midpoint,
        right,
        tolerance / 2.0,
        right_area,
        depth - 1,
    )


def finite_race_factors(
    proposal: Mapping[object, float],
    candidates: int,
    rate: Mapping[object, float],
    tolerance: float = 1e-11,
) -> dict[object, float]:
    """Return the finite-race likelihood multiplier for every score.

    Draw ``candidates`` independent scores from ``proposal`` and select one
    with probability proportional to ``rate[score]``.  If the selected score
    has distribution ``Q``, the returned value is ``Q(score)/proposal(score)``.

    For a fixed rate ``a`` the multiplier is

        K * integral_0^1 (sum_s p_s u^(a_s/a))^(K-1) du.

    This follows by representing ``a/(a+S)`` as a Laplace integral.  The
    bounded integral is evaluated deterministically with adaptive Simpson
    quadrature, avoiding Monte Carlo noise in downstream posterior studies.
    """

    if candidates < 1:
        raise ValueError("candidate count must be positive")
    positive = {
        score: probability
        for score, probability in proposal.items()
        if probability > 0.0
    }
    total = sum(positive.values())
    if total <= 0.0:
        raise ValueError("proposal distribution has zero mass")
    probabilities = {
        score: probability / total for score, probability in positive.items()
    }
    if candidates == 1:
        return {score: 1.0 for score in probabilities}
    if any(rate.get(score, 0.0) <= 0.0 for score in probabilities):
        raise ValueError("selection rates must be positive on proposal support")

    result: dict[object, float] = {}
    for score in probabilities:
        reference_rate = rate[score]

        def integrand(unit: float) -> float:
            if unit == 0.0:
                return 0.0
            transform = sum(
                probability * unit ** (rate[other] / reference_rate)
                for other, probability in probabilities.items()
            )
            return transform ** (candidates - 1)

        midpoint = 0.5
        whole = (
            integrand(0.0)
            + 4.0 * integrand(midpoint)
            + integrand(1.0)
        ) / 6.0
        result[score] = candidates * _adaptive_simpson(
            integrand,
            0.0,
            1.0,
            tolerance,
            whole,
            24,
        )
    return result


def exact_occupancy_race_factors(
    items: int,
    bins: tuple[int, int],
    candidates: int,
) -> dict[int, float]:
    """Return ``Q_K(D)/P_raw(D)`` for the occupancy-weighted race."""

    deficit_mass, placement_mass = exact_deficit_masses(items, bins)
    total_placements = sum(placement_mass.values())
    proposal = {
        deficit: mass / total_placements
        for deficit, mass in placement_mass.items()
        if mass > 0.0
    }
    rates = {
        deficit: deficit_mass[deficit] / placement_mass[deficit]
        for deficit in proposal
    }
    return finite_race_factors(proposal, candidates, rates)


def cumulative_distribution(
    distribution: Mapping[Feature, float]
) -> tuple[list[Feature], list[float]]:
    features = sorted(distribution)
    cumulative = []
    total = 0.0
    for feature in features:
        total += distribution[feature]
        cumulative.append(total)
    cumulative[-1] = 1.0
    return features, cumulative


def draw_feature(
    features: Sequence[Feature], cumulative: Sequence[float], rng: random.Random
) -> Feature:
    return features[bisect.bisect_left(cumulative, rng.random())]


def selected_distribution(
    proposal: Mapping[Feature, float],
    candidates: int,
    batches: int,
    weight: Callable[[Feature], float],
    rng: random.Random,
) -> dict[Feature, float]:
    if candidates == 1:
        return dict(proposal)
    features, cumulative = cumulative_distribution(proposal)
    selected: dict[Feature, float] = defaultdict(float)
    for _ in range(batches):
        batch = [draw_feature(features, cumulative, rng) for _ in range(candidates)]
        weights = [weight(feature) for feature in batch]
        total = sum(weights)
        for feature, proposal_weight in zip(batch, weights):
            selected[feature] += proposal_weight / total
    return normalized(selected)


def leakage_metrics(
    uniform: Mapping[Feature, float],
    selected: Mapping[Feature, float],
    descriptor_batch: int,
) -> Leakage:
    kl = 0.0
    tv = 0.0
    bhattacharyya = 0.0
    midpoint_terms: list[tuple[float, float]] = []
    for feature, uniform_mass in uniform.items():
        selected_mass = selected.get(feature, 0.0)
        tv += abs(selected_mass - uniform_mass)
        if selected_mass:
            ratio = selected_mass / uniform_mass
            kl += selected_mass * math.log2(ratio)
            midpoint_mass = math.sqrt(selected_mass * uniform_mass)
            bhattacharyya += midpoint_mass
            midpoint_terms.append((midpoint_mass, math.log(ratio)))

    pair_chernoff = -2.0 * math.log2(bhattacharyya)
    midpoint_mean = sum(mass * log_ratio for mass, log_ratio in midpoint_terms)
    midpoint_mean /= bhattacharyya
    midpoint_variance = sum(
        mass * (log_ratio - midpoint_mean) ** 2
        for mass, log_ratio in midpoint_terms
    ) / bhattacharyya
    if midpoint_variance < 1e-15:
        search_gain = descriptor_batch * pair_chernoff
    else:
        search_gain = (
            descriptor_batch * pair_chernoff
            + 0.5 * math.log2(math.pi * descriptor_batch * midpoint_variance)
            - 1.0
        )
    return Leakage(kl, tv / 2.0, pair_chernoff, search_gain)


def evaluate(
    items: int,
    bins: tuple[int, int],
    graph_trials: int,
    selector_batches: int,
    candidates: Sequence[int],
    descriptor_batch: int,
    seed: int,
    expansions: int,
) -> tuple[list[dict[str, object]], dict[str, float]]:
    graph_rng = random.Random(seed)
    histogram = sample_uniform_histogram(items, bins, graph_trials, graph_rng)
    uniform = {feature: count / graph_trials for feature, count in histogram.items()}
    proposal = proposal_distribution(uniform)
    occupancy = exact_occupancy_weights(items, bins)
    moments = conditional_inverse_mean_weights(
        uniform, lambda feature: (feature[1], feature[2])
    )
    feasible_probability = sum(
        probability for feature, probability in uniform.items() if feature[0]
    )
    feasible = normalized(
        {feature: probability for feature, probability in uniform.items() if feature[0]}
    )

    empirical_mean = sum(
        feature[0] * probability for feature, probability in uniform.items()
    )
    exact_mean = expected_placements(items, bins)

    rows: list[dict[str, object]] = []
    reference_cost = estimate_cost(
        CostModel(n=1 << 20, expansions=expansions), 16, 1
    ).total(expansions)

    def append_row(
        name: str,
        selector: str,
        k: int,
        distribution: Mapping[Feature, float],
        score_model: str | None,
    ) -> None:
        leakage = leakage_metrics(uniform, distribution, descriptor_batch)
        if score_model is None:
            cost_ratio = 1.0
            debias_ots = 0
            total_ots = estimate_cost(
                CostModel(
                    n=1 << 20,
                    expansions=expansions,
                    partition_sizes_override=bins,
                ),
                items,
                1,
            ).total(expansions)
        else:
            model = CostModel(
                n=1 << 20,
                expansions=expansions,
                geometry="free",
                score=score_model,
                partition_sizes_override=bins,
            )
            raw_cost = estimate_cost(model, items, 1).total(expansions)
            cost = estimate_cost(model, items, k)
            total_ots = cost.total(expansions)
            cost_ratio = total_ots / raw_cost
            debias_ots = cost.debias
        rows.append(
            {
                "name": name,
                "selector": selector,
                "t": items,
                "d": "/".join(str(size) for size in bins),
                "K": k,
                "kl_bits_per_list": leakage.kl_bits,
                "tv_per_list": leakage.tv,
                "pair_chernoff_per_list": leakage.pair_chernoff_bits,
                "search_gain_bits": leakage.search_gain_bits,
                "cost_ratio": cost_ratio,
                "cost_ratio_to_t16_d16": total_ots / reference_cost,
                "debias_ots": debias_ots,
                "total_ots": total_ots,
            }
        )

    append_row("raw", "raw", 1, proposal, None)
    for k in candidates:
        if k <= 1:
            continue
        exact_rng = random.Random(seed ^ (0x9E3779B9 * k))
        exact = selected_distribution(
            proposal,
            k,
            selector_batches,
            lambda feature: 1.0 / feature[0],
            exact_rng,
        )
        append_row(f"exact-K{k}", "inverse-Z", k, exact, "functional")

        occupancy_rng = random.Random(seed ^ (0x85EBCA6B * k))
        by_occupancy = selected_distribution(
            proposal,
            k,
            selector_batches,
            lambda feature: occupancy[feature[1]],
            occupancy_rng,
        )
        append_row(
            f"occupancy-K{k}",
            "inverse-E[Z|D]",
            k,
            by_occupancy,
            "occupancy",
        )

        moments_rng = random.Random(seed ^ (0xC2B2AE35 * k))
        by_moments = selected_distribution(
            proposal,
            k,
            selector_batches,
            lambda feature: moments[(feature[1], feature[2])],
            moments_rng,
        )
        append_row(
            f"moments-K{k}",
            "inverse-E[Z|D,C2]",
            k,
            by_moments,
            "occupancy",
        )
    append_row("feasible-limit", "U|Z>0", 0, feasible, None)
    diagnostics = {
        "feasible_probability": feasible_probability,
        "failure_probability": 1.0 - feasible_probability,
        "empirical_mean_Z": empirical_mean,
        "exact_mean_Z": exact_mean,
        "mean_Z_relative_error": empirical_mean / exact_mean - 1.0,
        "feature_categories": float(len(histogram)),
    }
    return rows, diagnostics


def print_rows(rows: Sequence[Mapping[str, object]]) -> None:
    print(
        "name             K   KL/list   TV/list   Cpair/list   search-gain"
        "  local-x   base-x   debias-OTs"
    )
    for row in rows:
        print(
            f"{row['name']:<16} {row['K']:>2}"
            f" {row['kl_bits_per_list']:9.5f}"
            f" {row['tv_per_list']:9.5f}"
            f" {row['pair_chernoff_per_list']:12.6f}"
            f" {row['search_gain_bits']:13.3f}"
            f" {row['cost_ratio']:8.2f}x"
            f" {row['cost_ratio_to_t16_d16']:7.2f}x"
            f" {row['debias_ots']:12,d}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=16)
    parser.add_argument(
        "--d",
        default="16",
        help="one size for equal partitions, or two comma-separated sizes",
    )
    parser.add_argument("--graph-trials", type=int, default=1_000_000)
    parser.add_argument("--selector-batches", type=int, default=1_000_000)
    parser.add_argument("--K", default="2,4,8")
    parser.add_argument(
        "--descriptor-batch",
        type=int,
        default=0,
        help="number of leaked lists; 0 uses polys^2*t = 16t",
    )
    parser.add_argument("--expansions", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0xC0FFEE)
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    parsed_bins = tuple(int(value) for value in args.d.split(","))
    if len(parsed_bins) == 1:
        bins = (parsed_bins[0], parsed_bins[0])
    elif len(parsed_bins) == 2:
        bins = parsed_bins
    else:
        parser.error("--d must contain one or two sizes")
    if min(args.t, *bins, args.graph_trials, args.selector_batches) < 1:
        parser.error("t, d, and trial counts must be positive")
    candidates = sorted({int(value) for value in args.K.split(",")})
    descriptor_batch = args.descriptor_batch or 16 * args.t

    rows, diagnostics = evaluate(
        args.t,
        bins,
        args.graph_trials,
        args.selector_batches,
        candidates,
        descriptor_batch,
        args.seed,
        args.expansions,
    )
    print(
        "diagnostics: "
        + " ".join(f"{key}={value:.9g}" for key, value in diagnostics.items())
    )
    print_rows(rows)
    if args.csv:
        with args.csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
