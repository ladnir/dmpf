"""Exact coefficient-aware toy Ring-LPN game with Reverse-Cuckoo leakage.

Unlike rev_cuckoo_full_lpn.py's binary unit-coefficient structural model, this
enumerator assigns every regular sparse position a nonzero coefficient in a
small prime field and forms the syndrome in F_q[X]/(X^N+1).  It remains a toy
model: descriptors are ideal and conditionally independent, and parameters
must be small enough to enumerate all unknown positions and coefficients.

For ``polys=2``, the public sample has the basic Ring-LPN shape
``(r, f + r*e)`` with regular ``e`` and ``f``.  The Reverse-Cuckoo leakage is
scored only after conditioning on that full syndrome.  The raw channel is the
original planted descriptor; the occupancy channel samples ``K`` independent
planted descriptors and selects one using the finite-K occupancy rule.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
from dataclasses import dataclass
from typing import Sequence

from rev_cuckoo_full_lpn import (
    NEG_INF,
    TrialResult,
    logsumexp2,
    mean,
    posterior_entropy,
    print_metric,
    self_test as cuckoo_self_test,
    support_lists,
)
from rev_cuckoo_regular_noise import LeakageChannel


def add_polys(left: Sequence[int], right: Sequence[int], modulus: int) -> tuple[int, ...]:
    return tuple((a + b) % modulus for a, b in zip(left, right))


def negacyclic_product(
    left: Sequence[int], right: Sequence[int], modulus: int
) -> tuple[int, ...]:
    width = len(left)
    result = [0] * width
    for left_index, left_value in enumerate(left):
        if left_value == 0:
            continue
        for right_index, right_value in enumerate(right):
            if right_value == 0:
                continue
            index = left_index + right_index
            product = left_value * right_value
            if index >= width:
                result[index - width] -= product
            else:
                result[index] += product
    return tuple(value % modulus for value in result)


def regular_poly(
    positions: Sequence[int], coefficients: Sequence[int], block_length: int,
    modulus: int
) -> tuple[int, ...]:
    result = [0] * (len(positions) * block_length)
    for block, (position, coefficient) in enumerate(zip(positions, coefficients)):
        result[block * block_length + position] = coefficient % modulus
    return tuple(result)


@dataclass(frozen=True)
class PolyState:
    positions: tuple[int, ...]
    coefficients: tuple[int, ...]
    polynomial: tuple[int, ...]


@dataclass(frozen=True)
class FieldModel:
    polys: int
    weight: int
    block_length: int
    modulus: int
    poly_states: tuple[PolyState, ...]
    full_states: tuple[tuple[int, ...], ...]
    position_candidates: tuple[tuple[int, ...], ...]


def make_model(
    polys: int,
    weight: int,
    block_length: int,
    modulus: int,
    enumerate_full_states: bool = True,
) -> FieldModel:
    poly_states = []
    for positions in itertools.product(range(block_length), repeat=weight):
        for coefficients in itertools.product(range(1, modulus), repeat=weight):
            poly_states.append(
                PolyState(
                    positions=positions,
                    coefficients=coefficients,
                    polynomial=regular_poly(
                        positions, coefficients, block_length, modulus
                    ),
                )
            )
    states = tuple(poly_states)
    return FieldModel(
        polys=polys,
        weight=weight,
        block_length=block_length,
        modulus=modulus,
        poly_states=states,
        full_states=(
            tuple(itertools.product(range(len(states)), repeat=polys))
            if enumerate_full_states
            else ()
        ),
        position_candidates=(
            tuple(
                itertools.product(
                    range(block_length), repeat=polys * weight
                )
            )
            if enumerate_full_states
            else ()
        ),
    )


def flatten_positions(model: FieldModel, state: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        position
        for state_index in state
        for position in model.poly_states[state_index].positions
    )


def product_active_list(
    model: FieldModel,
    known_poly: PolyState,
    hidden_poly: PolyState,
    output_block: int,
) -> tuple[int, ...]:
    """Return nonzero local addresses after product-row deduplication."""

    coefficients_by_address: dict[int, int] = {}
    for known_block in range(model.weight):
        hidden_block = (output_block - known_block) % model.weight
        address = (
            known_poly.positions[known_block]
            + hidden_poly.positions[hidden_block]
        )
        coefficient = (
            known_poly.coefficients[known_block]
            * hidden_poly.coefficients[hidden_block]
        )
        if known_block + hidden_block >= model.weight:
            coefficient = -coefficient
        coefficients_by_address[address] = (
            coefficients_by_address.get(address, 0) + coefficient
        ) % model.modulus
    return tuple(
        sorted(
            address
            for address, coefficient in coefficients_by_address.items()
            if coefficient
        )
    )


def product_raw_addresses(
    model: FieldModel,
    known_poly: PolyState,
    hidden_poly: PolyState,
    output_block: int,
) -> tuple[int, ...]:
    """Return local addresses after address deduplication but before values."""

    return tuple(
        sorted(
            {
                known_poly.positions[known_block]
                + hidden_poly.positions[(output_block - known_block) % model.weight]
                for known_block in range(model.weight)
            }
        )
    )


def product_active_lists(
    model: FieldModel,
    known_states: Sequence[PolyState],
    hidden_state: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    """Return the deduplicated nonzero rows of every sparse product DMPF.

    Product terms are grouped exactly as in ``RingLpnTriple::genDpf``: list
    ``(a,b,k)`` contains the ``t`` pairs whose input block indices sum to
    ``k`` modulo ``t``.  Its local address is the sum of the two within-block
    offsets in ``[0,2L)``.  Terms wrapping in the block index are negated for
    reduction modulo ``X^N+1``.  Equal local addresses are deduplicated and
    their coefficients are added; a zero sum becomes an inactive row.
    """

    if len(known_states) != model.polys or len(hidden_state) != model.polys:
        raise ValueError("expected one state for every polynomial")
    lists = []
    for known_poly in known_states:
        for hidden_state_index in hidden_state:
            hidden_poly = model.poly_states[hidden_state_index]
            for output_block in range(model.weight):
                lists.append(
                    product_active_list(
                        model, known_poly, hidden_poly, output_block
                    )
                )
    return tuple(lists)


def product_active_lists_one_hidden(
    model: FieldModel,
    known_states: Sequence[PolyState],
    hidden_poly: PolyState,
) -> tuple[tuple[int, ...], ...]:
    """Return the ``P*t`` product lists involving one hidden polynomial."""

    if len(known_states) != model.polys:
        raise ValueError("expected one known state for every polynomial")
    lists = []
    for known_poly in known_states:
        for output_block in range(model.weight):
            lists.append(
                product_active_list(
                    model, known_poly, hidden_poly, output_block
                )
            )
    return tuple(lists)


def random_ring_element(width: int, modulus: int, rng: random.Random) -> tuple[int, ...]:
    return tuple(rng.randrange(modulus) for _ in range(width))


def contribution_tables(
    model: FieldModel, multipliers: Sequence[Sequence[int]]
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    tables = [tuple(state.polynomial for state in model.poly_states)]
    for poly_index in range(1, model.polys):
        multiplier = multipliers[poly_index - 1]
        tables.append(
            tuple(
                negacyclic_product(multiplier, state.polynomial, model.modulus)
                for state in model.poly_states
            )
        )
    return tuple(tables)


def state_syndrome(
    model: FieldModel, state: Sequence[int],
    contributions: Sequence[Sequence[Sequence[int]]]
) -> tuple[int, ...]:
    result = (0,) * (model.weight * model.block_length)
    for poly_index, state_index in enumerate(state):
        result = add_polys(result, contributions[poly_index][state_index], model.modulus)
    return result


def field_trial(
    model: FieldModel,
    channel: LeakageChannel,
    rng: random.Random,
    descriptor_rng: random.Random | None = None,
    product_values: bool = False,
) -> TrialResult:
    leakage_rng = descriptor_rng if descriptor_rng is not None else rng
    if product_values:
        known_states = tuple(
            model.poly_states[rng.randrange(len(model.poly_states))]
            for _ in range(model.polys)
        )
        known = tuple(
            position
            for state in known_states
            for position in state.positions
        )
    else:
        known_states = ()
        known = tuple(
            rng.randrange(model.block_length)
            for _ in range(model.polys * model.weight)
        )
    true_state_index = rng.randrange(len(model.full_states))
    true_state = model.full_states[true_state_index]
    true_positions = flatten_positions(model, true_state)
    width = model.weight * model.block_length
    multipliers = tuple(
        random_ring_element(width, model.modulus, rng)
        for _ in range(model.polys - 1)
    )
    contributions = contribution_tables(model, multipliers)
    true_syndrome = state_syndrome(model, true_state, contributions)

    if product_values:
        true_lists = product_active_lists(model, known_states, true_state)
    else:
        true_lists = support_lists(
            known, true_positions, model.polys, model.weight
        )
    descriptors = tuple(
        channel.sample(active, 2 * model.block_length, leakage_rng)
        for active in true_lists
    )

    if product_values:
        factor_weights = []
        for hidden_poly_index in range(model.polys):
            factor_descriptors = tuple(
                descriptors[
                    (known_poly_index * model.polys + hidden_poly_index)
                    * model.weight
                    + output_block
                ]
                for known_poly_index in range(model.polys)
                for output_block in range(model.weight)
            )
            factor_weights.append(
                tuple(
                    channel.candidate_log_weight(
                        product_active_lists_one_hidden(
                            model, known_states, hidden_poly
                        ),
                        factor_descriptors,
                    )
                    for hidden_poly in model.poly_states
                )
            )
        state_log_weights = {
            state: sum(
                factor_weights[poly_index][state_index]
                for poly_index, state_index in enumerate(state)
            )
            for state in model.full_states
        }
        weights_by_position: dict[tuple[int, ...], list[float]] = {}
        for state, log_weight in state_log_weights.items():
            positions = flatten_positions(model, state)
            weights_by_position.setdefault(positions, []).append(log_weight)
        position_log_weights = {
            positions: logsumexp2(log_weights) - math.log2(len(log_weights))
            for positions, log_weights in weights_by_position.items()
        }
        true_log_weight = state_log_weights[true_state]
    else:
        position_log_weights = {}
        for positions in model.position_candidates:
            lists = support_lists(known, positions, model.polys, model.weight)
            position_log_weights[positions] = channel.candidate_log_weight(
                lists, descriptors
            )
        state_log_weights = {
            state: position_log_weights[flatten_positions(model, state)]
            for state in model.full_states
        }
        true_log_weight = position_log_weights[true_positions]
    if true_log_weight == NEG_INF:
        raise AssertionError("planted support has zero likelihood")

    log_mean_all = logsumexp2(state_log_weights.values()) - math.log2(
        len(state_log_weights)
    )
    consistent_weights = []
    for state in model.full_states:
        if state_syndrome(model, state, contributions) == true_syndrome:
            consistent_weights.append(
                state_log_weights[state]
            )
    if not consistent_weights:
        raise AssertionError("true state is not syndrome-compatible")

    log_mean_consistent = logsumexp2(consistent_weights) - math.log2(
        len(consistent_weights)
    )
    total_information = true_log_weight - log_mean_all
    decoding_information = true_log_weight - log_mean_consistent
    distinguishing_information = log_mean_consistent - log_mean_all
    if abs(
        total_information - decoding_information - distinguishing_information
    ) > 1e-10:
        raise AssertionError("information-density decomposition failed")

    greater = sum(value > true_log_weight + 1e-12 for value in consistent_weights)
    tied = sum(abs(value - true_log_weight) <= 1e-12 for value in consistent_weights)
    position_weights = tuple(position_log_weights.values())
    search_greater = sum(
        value > true_log_weight + 1e-12 for value in position_weights
    )
    search_tied = sum(
        abs(value - true_log_weight) <= 1e-12 for value in position_weights
    )
    log_total_consistent = logsumexp2(consistent_weights)
    return TrialResult(
        total_information=total_information,
        distinguishing_information=distinguishing_information,
        decoding_information=decoding_information,
        base_entropy=math.log2(len(consistent_weights)),
        posterior_entropy=posterior_entropy(consistent_weights),
        posterior_min_entropy=log_total_consistent - max(consistent_weights),
        base_rank=(len(consistent_weights) + 1.0) / 2.0,
        posterior_rank=greater + (tied + 1.0) / 2.0,
        search_base_rank=(len(position_weights) + 1.0) / 2.0,
        search_rank=search_greater + (search_tied + 1.0) / 2.0,
    )


def self_test() -> None:
    cuckoo_self_test()
    # X^(N-1) * X = -1 modulo X^N+1.
    left = (0, 0, 0, 1)
    right = (0, 1, 0, 0)
    if negacyclic_product(left, right, 5) != (4, 0, 0, 0):
        raise AssertionError("negacyclic wrap sign is wrong")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polys", type=int, default=2)
    parser.add_argument("--weight", type=int, default=2)
    parser.add_argument("--block-length", type=int, default=3)
    parser.add_argument("--modulus", type=int, default=3)
    parser.add_argument("--partition-size", type=int)
    parser.add_argument(
        "--channel", choices=("raw", "occupancy", "feasible"), default="raw"
    )
    parser.add_argument("--selector-k", type=int, default=1)
    parser.add_argument(
        "--leakage-input",
        choices=("support", "product-values"),
        default="support",
        help="whether deduplication may deactivate coefficient cancellations",
    )
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    partition_size = args.partition_size or args.weight
    if min(
        args.polys,
        args.weight,
        args.block_length,
        args.modulus - 1,
        partition_size,
        args.selector_k,
        args.trials,
    ) < 1:
        parser.error("all dimensions, trials, and modulus must be positive")

    self_test()
    model = make_model(
        args.polys, args.weight, args.block_length, args.modulus
    )
    channel = LeakageChannel.build(
        args.channel, args.selector_k, args.weight, partition_size
    )
    rng = random.Random(args.seed)
    results = [
        field_trial(
            model,
            channel,
            rng,
            product_values=args.leakage_input == "product-values",
        )
        for _ in range(args.trials)
    ]
    print(
        f"field={args.modulus} polys={args.polys} weight={args.weight} "
        f"block_length={args.block_length} d={partition_size} "
        f"channel={args.channel} K={args.selector_k} "
        f"leakage_input={args.leakage_input} "
        f"lists={args.polys * args.polys * args.weight} "
        f"position_candidates={len(model.position_candidates)} "
        f"full_candidates={len(model.full_states)} trials={args.trials} "
        f"seed={args.seed}"
    )
    print_metric(
        "I(X;Lambda|known) density",
        [result.total_information for result in results],
    )
    print_metric(
        "k_dist density",
        [result.distinguishing_information for result in results],
    )
    print_metric(
        "k_dec density", [result.decoding_information for result in results]
    )
    print_metric(
        "posterior entropy loss",
        [result.base_entropy - result.posterior_entropy for result in results],
    )
    print_metric(
        "posterior min-entropy loss",
        [result.base_entropy - result.posterior_min_entropy for result in results],
    )
    print_metric(
        "posterior-rank gain",
        [
            math.log2(result.base_rank / result.posterior_rank)
            for result in results
        ],
    )
    decoding_expected_work_gain = -math.log2(
        mean([result.posterior_rank / result.base_rank for result in results])
    )
    print(f"decoding expected-work gain={decoding_expected_work_gain:.6f}")
    print_metric(
        "support-search mean-log rank gain",
        [
            math.log2(result.search_base_rank / result.search_rank)
            for result in results
        ],
    )
    expected_work_gain = -math.log2(
        mean([result.search_rank / result.search_base_rank for result in results])
    )
    print(f"support-search expected-work gain={expected_work_gain:.6f}")
    print(
        "decomposition_check="
        f"{mean([r.total_information for r in results]) - mean([r.distinguishing_information for r in results]) - mean([r.decoding_information for r in results]):.3e}"
    )


if __name__ == "__main__":
    main()
