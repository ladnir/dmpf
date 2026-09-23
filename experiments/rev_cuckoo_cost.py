"""Preliminary nonlinear-cost model for LPN-optimized Reverse Cuckoo.

The model counts total 1-out-of-2 OTs across both directions, following the
paper convention that one shared-bit multiplication costs two OTs.  It mirrors
the current libOTe formulas for deduplication, the shared Waksman permutation,
Goldreich hash evaluation, BinarySolver, SparseDpf, and conditional negation.

The proposed K-candidate debiaser has not been implemented.  The occupancy
score nevertheless has an implementation-level count: HybEquality computes
all endpoint equalities, the existing Dedup OR tree produces one duplicate
indicator per noninitial endpoint, and a balanced adder tree sums those bits.
The other graph scores remain deliberately visible estimates:

* ``functional`` uses the planted-witness functional-graph representation and
  assumes a pointer-jumping/component-count circuit with O(t^2 log t) shared
  products.
* ``closure`` uses a conservative t^3 transitive-closure circuit.

The default ``free`` geometry permits every integer t and uses d=t.  Passing
``--geometry current`` overlays the current implementation restriction
d=2^ceil(log2(t)); the current RingLpnTriple also restricts t itself to a power
of two, but the table intentionally shows the intermediate integer values.

This is an OT/nonlinear-work model, not a runtime or complete communication
model.  XORs, openings, PRG work, local full-domain expansion, framing, and OT
extension setup are omitted.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from typing import Iterable, Sequence


def ceil_log2(value: int) -> int:
    if value <= 0:
        raise ValueError("value must be positive")
    return max(1, (value - 1).bit_length())


def next_power_of_two(value: int) -> int:
    if value <= 0:
        raise ValueError("value must be positive")
    return 1 << (value - 1).bit_length()


def waksman_switches(size: int) -> int:
    """Match WaksmanPermute::switchCount in libOTe."""

    if size <= 1:
        return 0
    if size == 2:
        return 1
    lower = size // 2
    upper = (size + 1) // 2
    return lower + waksman_switches(upper) + waksman_switches(lower) + lower


@dataclass(frozen=True)
class Model:
    n: int = 1 << 20
    polys: int = 4
    w: int = 2
    linear_sec: int = 10
    field_bits: int = 64
    expansions: int = 10
    geometry: str = "free"
    score: str = "functional"
    selector_key_bits: int = 80
    selector_record_bits_override: int | None = None
    placement_precision_bits: int = 64
    partition_size_override: int | None = None
    partition_sizes_override: tuple[int, ...] | None = None

    def partition_sizes(self, t: int) -> tuple[int, ...]:
        if self.partition_sizes_override is not None:
            if len(self.partition_sizes_override) != self.w:
                raise ValueError("partition sizes must match w")
            if any(size < 1 for size in self.partition_sizes_override):
                raise ValueError("partition sizes must be positive")
            return self.partition_sizes_override
        return (self.partition_size(t),) * self.w

    def partition_size(self, t: int) -> int:
        if self.partition_size_override is not None:
            if self.partition_size_override < 1:
                raise ValueError("partition-size override must be positive")
            return self.partition_size_override
        if self.geometry == "free":
            return t
        if self.geometry == "current":
            return next_power_of_two(t)
        raise ValueError(f"unknown geometry: {self.geometry}")

    def selector_record_bits(self, t: int) -> int:
        if self.selector_record_bits_override is not None:
            if self.selector_record_bits_override < 0:
                raise ValueError("selector record bits must be nonnegative")
            return self.selector_record_bits_override
        # Both endpoint labels plus one planted-side bit per active row.
        return t * (sum(ceil_log2(d) for d in self.partition_sizes(t)) + 1)


@dataclass(frozen=True)
class Cost:
    dedup: int
    selected_permutation: int
    hash_eval: int
    solver: int
    sparse_dpf: int
    negate: int
    position_add: int
    extra_proposals: int
    graph_scores: int
    selector: int
    selector_payload_bits: int
    tensor_per_expansion: int

    @property
    def selected_dmpf(self) -> int:
        return (
            self.dedup
            + self.selected_permutation
            + self.hash_eval
            + self.solver
            + self.sparse_dpf
            + self.negate
        )

    @property
    def debias(self) -> int:
        return self.extra_proposals + self.graph_scores + self.selector

    @property
    def one_time(self) -> int:
        return self.selected_dmpf + self.position_add + self.debias

    def total(self, expansions: int) -> int:
        return self.one_time + expansions * self.tensor_per_expansion

    def amortized(self, expansions: int) -> float:
        return self.one_time / expansions + self.tensor_per_expansion


def dedup_ots(t: int, key_bits: int, batches: int) -> int:
    pairs = math.comb(t, 2)
    equality_ots = pairs * (key_bits + ceil_log2(key_bits + 1))
    shared_products = 2 * pairs + t
    return batches * (equality_ots + 2 * shared_products)


def binary_solver_products(rows: int, columns: int, label_bits: int) -> int:
    one_hot = rows * (2 * columns - 3)
    matrix_vector = rows * columns
    matrix_update = rows * rows
    random_solution = columns
    right_multiply = min(columns * rows, columns * label_bits)
    correction_multiply = min(columns * rows, rows * label_bits)
    return (
        one_hot
        + matrix_vector
        + matrix_update
        + random_solution
        + right_multiply
        + correction_multiply
    )


def occupancy_products(t: int, label_bits: Sequence[int]) -> int:
    """Legacy product estimate for the exact repeated-endpoint count."""

    equality = math.comb(t, 2) * sum(max(0, bits - 1) for bits in label_bits)
    duplicate_or = len(label_bits) * (t - 1) * (t - 2) // 2
    return equality + duplicate_or


def bounded_sum_products(maxima: Sequence[int]) -> int:
    """ANDs in a balanced sum tree with public bounds on every input."""

    if any(maximum < 0 for maximum in maxima):
        raise ValueError("bounds must be nonnegative")
    bounds = list(maxima)
    products = 0
    while len(bounds) > 1:
        next_bounds: list[int] = []
        for i in range(0, len(bounds) - 1, 2):
            maximum = bounds[i] + bounds[i + 1]
            products += ceil_log2(maximum + 1) - 1
            next_bounds.append(maximum)
        if len(bounds) & 1:
            next_bounds.append(bounds[-1])
        bounds = next_bounds
    return products


def popcount_products(bit_count: int) -> int:
    """ANDs in the range-aware balanced popcount used by DistinctCount."""

    if bit_count < 0:
        raise ValueError("bit count must be nonnegative")
    return bounded_sum_products([1] * bit_count)


def occupancy_score_ots(t: int, label_bits: Sequence[int]) -> int:
    """Exact OT count for the occupancy deficit of one proposal graph."""

    pairs = math.comb(t, 2)
    equality = pairs * sum(
        bits + ceil_log2(bits + 1) for bits in label_bits
    )
    duplicate_or_products = len(label_bits) * (t - 1) * (t - 2) // 2
    partition_popcounts = len(label_bits) * popcount_products(t - 1)
    combine_partitions = bounded_sum_products([t - 1] * len(label_bits))
    return equality + 2 * (
        duplicate_or_products + partition_popcounts + combine_partitions
    )


def graph_score_products(t: int, label_bits: Sequence[int], method: str) -> int:
    occupied = occupancy_products(t, label_bits)
    if method == "occupancy":
        return occupied
    if method == "functional":
        # The planted placement turns every item into a partial pointer to the
        # item occupying its unplanted endpoint.  The coefficient 4 is an
        # explicit preliminary allowance for pointer jumping, root/cycle
        # flags, component sizes, and the product of component factors.
        return occupied + 4 * t * t * ceil_log2(t)
    if method == "closure":
        return occupied + t**3 + 2 * t * t * ceil_log2(t)
    raise ValueError(f"unknown graph score: {method}")


def graph_score_ots(t: int, label_bits: Sequence[int], method: str) -> int:
    if method == "occupancy":
        return occupancy_score_ots(t, label_bits)
    return 2 * graph_score_products(t, label_bits, method)


def public_table_lookup_products(address_bits: int) -> int:
    """Wide products for a public table indexed by secret-shared bits.

    The first mux level selects between public strings and is linear.  The
    remaining binary mux tree has 2^(q-1)-1 shared-bit by shared-string
    products.  Each such product costs two OTs independent of string width.
    """

    if address_bits < 1:
        return 0
    return (1 << (address_bits - 1)) - 1


def inverse_cdf_sampler_ots(outcomes: int, precision_bits: int) -> int:
    """Conservative OT count for a secret fixed-distribution sample."""

    if outcomes < 2 or precision_bits < 1:
        raise ValueError("require at least two outcomes and one precision bit")
    depth = ceil_log2(outcomes)
    comparisons = 2 * depth * precision_bits
    # At level ell, the next public CDF threshold is a 2^ell-entry table
    # indexed by the preceding branch bits.  The first mux layer is linear.
    threshold_muxes = 2 * sum(
        public_table_lookup_products(prefix_bits)
        for prefix_bits in range(1, depth)
    )
    unary_decoder = 2 * ((1 << depth) - 2)
    return comparisons + threshold_muxes + unary_decoder


def placement_candidate_ots(
    t: int,
    partition_sizes: Sequence[int],
    precision_bits: int,
) -> int:
    """OTs to sample one hidden planted placement.

    Equal partitions use the ordinary full-slot Waksman shuffle.  For the
    asymmetric two-partition construction, sample the left row count by an
    inverse CDF, decode its prefix, and use independent hidden permutations
    of the rows and both bin sets.  One wide endpoint mux per row applies the
    secret side choice.
    """

    if len(set(partition_sizes)) == 1:
        return 2 * waksman_switches(sum(partition_sizes))
    if len(partition_sizes) != 2:
        raise ValueError("asymmetric placement model currently supports w=2")
    permutation_switches = waksman_switches(t) + sum(
        waksman_switches(size) for size in partition_sizes
    )
    side_count = inverse_cdf_sampler_ots(t + 1, precision_bits)
    endpoint_muxes = 2 * t
    return 2 * permutation_switches + side_count + endpoint_muxes


def selector_ots(
    k: int,
    t: int,
    partitions: int,
    key_bits: int,
) -> int:
    if k <= 1:
        return 0
    if key_bits < 1:
        raise ValueError("selector key bits must be positive")

    # For every candidate, public randomness supplies a table of exponential
    # race keys indexed by D.  D is in [0, partitions*(t-1)].  A binary mux
    # tree performs the secret lookup.  A K-way tournament then uses one
    # key comparison and one wide conditional move per loser.  cryptoTools'
    # size-optimized L-bit unsigned comparison uses L ANDs.
    max_deficit = partitions * (t - 1)
    address_bits = ceil_log2(max_deficit + 1)
    lookups = 2 * k * public_table_lookup_products(address_bits)
    tournament = (k - 1) * (2 * key_bits + 2)
    return lookups + tournament


def selector_payload_bits(
    k: int,
    t: int,
    partitions: int,
    key_bits: int,
    record_bits: int,
) -> int:
    """Correlated-OT message payload used by the weighted selector.

    This excludes OT-extension setup and the proposal-generation payload.  A
    bit AND has two one-bit OT messages.  A wide conditional move has two OT
    messages as wide as its string operand.
    """

    if k <= 1:
        return 0
    if min(key_bits, record_bits) < 0:
        raise ValueError("bit widths must be nonnegative")
    max_deficit = partitions * (t - 1)
    address_bits = ceil_log2(max_deficit + 1)
    lookup_products = k * public_table_lookup_products(address_bits)
    lookup_payload = 2 * lookup_products * key_bits
    comparison_payload = (k - 1) * 2 * key_bits
    winner_payload = (k - 1) * 2 * (key_bits + record_bits)
    return lookup_payload + comparison_payload + winner_payload


def estimate(model: Model, t: int, k: int) -> Cost:
    if t < 2 or k < 1:
        raise ValueError("require t >= 2 and K >= 1")

    partition_sizes = model.partition_sizes(t)
    f = sum(partition_sizes)
    batches = model.polys * model.polys * t
    key_bits = ceil_log2((2 * model.n) // t + 1)
    label_bits = tuple(ceil_log2(d) for d in partition_sizes)
    dpf_domain = math.ceil(2 * model.n / t)

    dedup = dedup_ots(t, key_bits, batches)
    permutation_per_candidate = placement_candidate_ots(
        t,
        partition_sizes,
        model.placement_precision_bits,
    ) * batches
    hash_products = batches * sum(
        d * (math.ceil((d + model.linear_sec) / 8) * 8)
        for d in partition_sizes
    )
    hash_eval = 2 * hash_products
    solver = (
        2
        * batches
        * sum(
            binary_solver_products(d, d + model.linear_sec, bits)
            for d, bits in zip(partition_sizes, label_bits)
        )
    )
    sparse_dpf = 2 * ceil_log2(dpf_domain) * batches * f
    negate = 2 * batches * f

    total_noise = model.polys * t
    position_add = 2 * total_noise * total_noise * ceil_log2(dpf_domain)
    tensor_per_expansion = total_noise * model.field_bits

    if k == 1:
        graph_scores = 0
    else:
        graph_scores = k * batches * graph_score_ots(t, label_bits, model.score)
    selector = batches * selector_ots(
        k,
        t,
        len(partition_sizes),
        model.selector_key_bits,
    )
    selector_payload = batches * selector_payload_bits(
        k,
        t,
        len(partition_sizes),
        model.selector_key_bits,
        model.selector_record_bits(t),
    )

    return Cost(
        dedup=dedup,
        selected_permutation=permutation_per_candidate,
        hash_eval=hash_eval,
        solver=solver,
        sparse_dpf=sparse_dpf,
        negate=negate,
        position_add=position_add,
        extra_proposals=(k - 1) * permutation_per_candidate,
        graph_scores=graph_scores,
        selector=selector,
        selector_payload_bits=selector_payload,
        tensor_per_expansion=tensor_per_expansion,
    )


def parse_ints(spec: str) -> list[int]:
    values: list[int] = []
    for term in spec.split(","):
        if ":" not in term:
            values.append(int(term))
            continue
        parts = [int(part) for part in term.split(":")]
        if len(parts) == 2:
            start, stop = parts
            step = 1
        elif len(parts) == 3:
            start, stop, step = parts
        else:
            raise ValueError(f"bad integer range: {term}")
        values.extend(range(start, stop + (1 if step > 0 else -1), step))
    return sorted(set(values))


def rows(
    model: Model,
    ts: Iterable[int],
    ks: Iterable[int],
    reference_t: int = 16,
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    reference_total = estimate(model, reference_t, 1).total(model.expansions)
    for t in ts:
        baseline = estimate(model, t, 1)
        baseline_total = baseline.total(model.expansions)
        for k in ks:
            cost = estimate(model, t, k)
            result.append(
                {
                    "t": t,
                    "d": "/".join(str(size) for size in model.partition_sizes(t)),
                    "K": k,
                    "noise_weight": model.polys * t,
                    "baseline_isd_bits": model.polys * t * math.log2(model.polys),
                    "descriptors": model.polys * model.polys * t,
                    "product_points": model.polys * model.polys * t * t,
                    "dedup_ots": cost.dedup,
                    "permutation_ots": cost.selected_permutation,
                    "hash_eval_ots": cost.hash_eval,
                    "solver_ots": cost.solver,
                    "sparse_dpf_ots": cost.sparse_dpf,
                    "negate_ots": cost.negate,
                    "selected_dmpf_ots": cost.selected_dmpf,
                    "position_add_ots": cost.position_add,
                    "extra_proposal_ots": cost.extra_proposals,
                    "graph_score_ots": cost.graph_scores,
                    "selector_ots": cost.selector,
                    "selector_payload_bits": cost.selector_payload_bits,
                    "debias_ots": cost.debias,
                    "one_time_ots": cost.one_time,
                    "tensor_ots_per_expansion": cost.tensor_per_expansion,
                    "amortized_ots": cost.amortized(model.expansions),
                    "total_ots": cost.total(model.expansions),
                    "ratio_to_raw_same_t": cost.total(model.expansions) / baseline_total,
                    "ratio_to_reference": cost.total(model.expansions) / reference_total,
                }
            )
    return result


def print_table(records: list[dict[str, object]]) -> None:
    header = (
        " t   d   K  noise  ISD   desc   points     selected      debias"
        "     one-time   amort/exp  raw-x  ref-x"
    )
    print(header)
    for row in records:
        print(
            f"{row['t']:2d} {row['d']:>5} {row['K']:3d}"
            f" {row['noise_weight']:6d} {row['baseline_isd_bits']:5.0f}"
            f" {row['descriptors']:6d} {row['product_points']:8d}"
            f" {row['selected_dmpf_ots']:12,d}"
            f" {row['debias_ots']:11,d}"
            f" {row['one_time_ots']:12,d}"
            f" {row['amortized_ots']:11,.0f}"
            f" {row['ratio_to_raw_same_t']:5.2f}x"
            f" {row['ratio_to_reference']:5.2f}x"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", default="12:24:2", help="values/ranges, e.g. 12:24:2,32")
    parser.add_argument("--K", default="1,2,4,8", help="candidate counts")
    parser.add_argument("--n", type=int, default=1 << 20)
    parser.add_argument("--polys", type=int, default=4)
    parser.add_argument("--linear-sec", type=int, default=10)
    parser.add_argument("--field-bits", type=int, default=64)
    parser.add_argument("--selector-key-bits", type=int, default=80)
    parser.add_argument("--selector-record-bits", type=int)
    parser.add_argument("--placement-precision-bits", type=int, default=64)
    parser.add_argument("--expansions", type=int, default=10)
    parser.add_argument("--reference-t", type=int, default=16)
    parser.add_argument("--geometry", choices=("free", "current"), default="free")
    parser.add_argument(
        "--d",
        help="explicit size for every partition, or comma-separated sizes",
    )
    parser.add_argument(
        "--score",
        choices=("occupancy", "functional", "closure"),
        default="functional",
    )
    parser.add_argument("--csv", action="store_true")
    args = parser.parse_args()

    explicit_sizes = None
    explicit_size = None
    if args.d:
        parsed_sizes = tuple(int(value) for value in args.d.split(","))
        if len(parsed_sizes) == 1:
            explicit_size = parsed_sizes[0]
        else:
            explicit_sizes = parsed_sizes

    model = Model(
        n=args.n,
        polys=args.polys,
        linear_sec=args.linear_sec,
        field_bits=args.field_bits,
        selector_key_bits=args.selector_key_bits,
        selector_record_bits_override=args.selector_record_bits,
        placement_precision_bits=args.placement_precision_bits,
        expansions=args.expansions,
        geometry=args.geometry,
        score=args.score,
        partition_size_override=explicit_size,
        partition_sizes_override=explicit_sizes,
    )
    records = rows(
        model,
        parse_ints(args.t),
        parse_ints(args.K),
        reference_t=args.reference_t,
    )
    if args.csv:
        writer = csv.DictWriter(sys.stdout, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    else:
        print_table(records)


if __name__ == "__main__":
    main()
