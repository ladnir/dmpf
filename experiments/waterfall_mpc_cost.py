"""Concrete OT model for Waterfall Cuckoo.

A product of a shared bit and a shared ``ell``-bit string costs two
``ell``-bit OTs.  In particular, a Boolean AND is the ``ell = 1`` case; a wide
conditional move is *not* expanded into ``ell`` separate ANDs. Shared DPF
correction selections also cost two OTs; privately controlled permutation
switches cost one OT each.

The two reported coordinates are the number of OTs and the sum of their
payload lengths.  XORs, public linear maps, openings, local PRG calls, circuit
depth, and OT-extension setup are outside this model.  The payload coordinate
is therefore not complete wire communication.
"""

from __future__ import annotations

import argparse
import functools
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Model:
    t: int = 16
    batch: int = 256
    correctness_bits: int = 40
    domain_bits: int = 20
    payload_bits: int = 64
    kappa: int = 128
    characteristic_two: bool = False

    @property
    def address_bits(self) -> int:
        # The wire domain contains [0,N) together with the dummy label N.
        return self.domain_bits + 1

    @property
    def scatter_record_bits(self) -> int:
        # Forward permutation: one t-bit placement column.  Inverse
        # permutation: activity, address, and one-hot row identity.  Payloads
        # use the cached row-to-column point DPFs instead.
        return self.t + 1 + self.address_bits


@dataclass(frozen=True)
class Configuration:
    name: str
    generator: str
    w: int
    d: int
    c: int
    bfs_depth: int = 0
    bfs_roots: int = 0
    proved: bool = True
    partition_sizes: tuple[int, ...] = ()

    @property
    def sizes(self) -> tuple[int, ...]:
        sizes = self.partition_sizes or (self.d,) * self.w
        if len(sizes) != self.w:
            raise ValueError("partition_sizes must contain exactly w entries")
        if any(d <= 0 or d & (d - 1) for d in sizes):
            raise ValueError("partition sizes must be positive powers of two")
        return sizes

    @property
    def partition_index_bits(self) -> tuple[int, ...]:
        return tuple(max(1, d.bit_length() - 1) for d in self.sizes)

    @property
    def index_bits(self) -> int:
        bits = self.partition_index_bits
        if any(k != bits[0] for k in bits[1:]):
            raise ValueError("this component requires equal partition sizes")
        return bits[0]

    @property
    def columns(self) -> int:
        return sum(self.sizes) + self.c


@dataclass(frozen=True)
class MultiplicationProfile:
    """Shared-bit multiplication count and sum of right-operand widths."""

    count: int
    width_sum: int

    @property
    def ot_count(self) -> int:
        return 2 * self.count

    @property
    def payload_bits(self) -> int:
        return 2 * self.width_sum


@dataclass(frozen=True)
class OtCost:
    ot_count: int
    payload_bits: int


@dataclass(frozen=True)
class Cost:
    common: OtCost
    generator: OtCost
    hash_eval: OtCost
    sparse_dpf: OtCost
    scatter: OtCost

    @property
    def total_ot_count(self) -> int:
        return sum(
            component.ot_count
            for component in (
                self.common,
                self.generator,
                self.hash_eval,
                self.sparse_dpf,
                self.scatter,
            )
        )

    @property
    def total_payload_bits(self) -> int:
        return sum(
            component.payload_bits
            for component in (
                self.common,
                self.generator,
                self.hash_eval,
                self.sparse_dpf,
                self.scatter,
            )
        )

    @property
    def total(self) -> int:
        """Primary optimization objective: total number of OTs."""

        return self.total_ot_count


def dedup_and_activity_cost(model: Model) -> int:
    """Sum of multiplication widths for the configuration-independent path."""

    t = model.t
    a = model.address_bits
    g = model.payload_bits
    equalities = math.comb(t, 2) * (a - 1)
    first_occurrence = (t - 1) * (t - 2) // 2
    address_muxes = t * a
    value_masks = (t * (t + 1) // 2) * g
    control_masks = t
    scatter_activity = t * ((a - 1) + a)
    return (
        equalities
        + first_occurrence
        + address_muxes
        + value_masks
        + control_masks
        + scatter_activity
    )


def dedup_and_activity_profile(model: Model) -> MultiplicationProfile:
    """Pack address and payload masks as bit-times-string products."""

    t = model.t
    a = model.address_bits
    scalar_products = (
        math.comb(t, 2) * (a - 1)
        + (t - 1) * (t - 2) // 2
        + t
        + t * (a - 1)
    )
    wide_products = 2 * t + t * (t + 1) // 2
    profile = MultiplicationProfile(
        count=scalar_products + wide_products,
        width_sum=dedup_and_activity_cost(model),
    )
    return profile


def karatsuba_multiplication_cost(bits: int) -> int:
    """Bilinear AND count for recursive polynomial-basis Karatsuba."""

    if bits <= 0:
        raise ValueError("field width must be positive")
    if bits == 1:
        return 1
    low = bits // 2
    high = bits - low
    return karatsuba_multiplication_cost(low) + 2 * karatsuba_multiplication_cost(high)


def field_multiplication_cost(bits: int) -> int:
    """Bilinear AND count of the selected exact field multiplier.

    The 20-bit implementation converts linearly to GF(16)^5, evaluates at
    nine GF(16) points, and uses a rank-nine GF(16) multiplier at each point.
    All basis conversion, interpolation, and reduction operations are linear.
    Other widths retain the recursive polynomial-basis Karatsuba circuit.
    """

    if bits == 20:
        return 9 * 9
    return karatsuba_multiplication_cost(bits)


def direct_hash_cost(model: Model, cfg: Configuration) -> int:
    """Exact polynomial hash evaluation plus dummy-candidate selection.

    Squaring in a fixed binary field basis is linear.  To obtain every power
    below t, only the odd powers x^3,x^5,... require multiplication; even
    powers are squares of earlier powers.  The powers are shared across all
    partitions.  Field products use the selected exact bilinear circuit, and
    reduction modulo the public irreducible polynomial is linear.
    """

    ell = model.domain_bits
    nonlinear_powers = max(0, (model.t - 2) // 2)
    powers = model.t * nonlinear_powers * field_multiplication_cost(ell)
    dummy_selection = model.t * sum(cfg.partition_index_bits)
    return powers + dummy_selection


def direct_hash_profile(model: Model, cfg: Configuration) -> MultiplicationProfile:
    """Field products are bit products; each dummy mux is one k-bit product."""

    ell = model.domain_bits
    nonlinear_powers = max(0, (model.t - 2) // 2)
    field_bit_products = (
        model.t
        * nonlinear_powers
        * field_multiplication_cost(ell)
    )
    dummy_muxes = cfg.w * model.t
    return MultiplicationProfile(
        count=field_bit_products + dummy_muxes,
        width_sum=direct_hash_cost(model, cfg),
    )


def waterfall_collision_cost(t: int, w: int, k: int) -> int:
    # A k-bit conditioned equality and the live-row OR/guard cost k+1 gates
    # per unordered pair and partition.
    return w * (k + 1) * math.comb(t, 2)


def basic_generator_cost(model: Model, cfg: Configuration) -> int:
    """Sum of operand widths for the occupancy-scan implementation."""

    t, c = model.t, cfg.c
    # Per row and partition: read the d-bit occupancy array, multiply the
    # resulting occupied bit by the live bit, and generate the winning one-hot
    # row.  The two mux trees each have aggregate width d-1.
    main_placement = t * sum(2 * d - 1 for d in cfg.sizes)
    stash_and_reject = t * (c + 1)
    return main_placement + stash_and_reject


def basic_generator_profile(model: Model, cfg: Configuration) -> MultiplicationProfile:
    """Scan secret occupancy arrays and emit the winning one-hot rows."""

    t, c = model.t, cfg.c
    scan_products = t * sum(2 * k + 1 for k in cfg.partition_index_bits)
    stash_products = t
    return MultiplicationProfile(
        count=scan_products + stash_products,
        width_sum=basic_generator_cost(model, cfg),
    )


def one_repair_generator_cost(model: Model, cfg: Configuration) -> int:
    """Sum of operand widths for the OT-minimized one-repair evaluator."""

    t, w, d, c, k = model.t, cfg.w, cfg.d, cfg.c, cfg.index_bits
    h = max(1, (t - 1).bit_length())
    log_w = max(1, (w - 1).bit_length())
    pairs = math.comb(w, 2)
    return (
        # Occupancy read, winning one-hot row, cached raw one-hot row, and
        # final gating by the repaired row-to-partition matching.
        w * t * (4 * d - 3)
        # Initial overflow classification and final stash compaction.
        + t * (2 * c + 3)
        # Read all root candidates, then their h-bit owners.
        + w * k * (t - 1)
        + w * h * (d - 1)
        # Read each victim's later candidate and the corresponding occupancy.
        + pairs * k * (t - 1)
        + pairs * (d - 1)
        # Gate and select the first feasible path.
        + 2 * pairs
        # Select the victim identifier and route the two w-bit matching deltas.
        + (1 << log_w) * h
        + 2 * w * (t - 1)
        + 1
    )


def one_repair_generator_profile(
    model: Model, cfg: Configuration
) -> MultiplicationProfile:
    """OT-minimized occupancy/owner scan with one bounded repair."""

    t, w, k = model.t, cfg.w, cfg.index_bits
    h = max(1, (t - 1).bit_length())
    log_w = max(1, (w - 1).bit_length())
    pairs = math.comb(w, 2)
    return MultiplicationProfile(
        count=(
            # Occupancy read, live/raw decoder generation, and final gating.
            w * t * (2 * k + 2)
            # Initial and final saturated prefix scans.
            + 2 * t
            # Root-row read, owner lookups, victim-row reads, occupancy reads.
            + h
            + w * k
            + (w - 1) * h
            + pairs * k
            # Gate the feasibility vector and select its first live entry.
            + 1
            + pairs
            # Select the victim identifier and route root/victim row updates.
            + log_w
            + 1
            + 2 * h
            # Final rejection gate.
            + 1
        ),
        width_sum=one_repair_generator_cost(model, cfg),
    )


def eviction_generator_cost(model: Model, cfg: Configuration) -> int:
    """Fixed-circuit upper count for a global bounded token-eviction walk."""

    if cfg.bfs_depth <= 0 or cfg.bfs_roots <= 0:
        raise ValueError("eviction configurations use bfs_depth as the step bound")
    t, w, d, c, k = model.t, cfg.w, cfg.d, cfg.c, cfg.index_bits
    collisions = waterfall_collision_cost(t, w, k)
    final_dense_placement = w * t * (2 * d - 2)
    stash_and_reject = t * (c + 1)
    root_selection = t * (cfg.bfs_roots + 1)

    # Select the token's w labels, scan each against t owners, select one
    # candidate, and update the token and row-to-partition matching.
    one_step = 2 * w * t * k + 2 * w * t + 2 * w
    return (
        collisions
        + final_dense_placement
        + stash_and_reject
        + root_selection
        + cfg.bfs_depth * one_step
    )


def bfs_generator_cost(model: Model, cfg: Configuration) -> int:
    """Fixed-circuit upper count for deterministic bounded BFS.

    The circuit stores row-to-row parent edges.  Candidate equalities from the
    initial Waterfall are reused.  Each BFS level scans every row, partition,
    and possible owner, selects the first parent, and stores enough information
    to backtrack the selected augmenting path.
    """

    if cfg.bfs_depth <= 0 or cfg.bfs_roots <= 0:
        raise ValueError("BFS configurations require positive depth and root bounds")
    t, w, d, c, k = model.t, cfg.w, cfg.d, cfg.c, cfg.index_bits

    collisions = waterfall_collision_cost(t, w, k)
    # Retain raw candidate decoders and gate them once by the final matching.
    final_dense_placement = w * t * (2 * d - 2)
    stash_and_reject = t * (c + 1)
    root_selection = t * (cfg.bfs_roots + 1)

    owner_graph = w * t * t
    one_level_search = 3 * w * t * t + 2 * w * t + t
    one_root = owner_graph + cfg.bfs_depth * one_level_search
    return (
        collisions
        + final_dense_placement
        + stash_and_reject
        + root_selection
        + cfg.bfs_roots * one_root
    )


def reachability_generator_cost(model: Model, cfg: Configuration) -> int:
    """Width sum for complete repair by implicit-owner multi-source DFS."""

    if cfg.bfs_roots <= 0:
        raise ValueError("reachability configurations require a positive root bound")
    t, w, c = model.t, cfg.w, cfg.c
    sizes = cfg.sizes
    h = max(1, (t - 1).bit_length())

    # The scan caches both raw and winning candidate decoders.  The final
    # matching gates the raw decoder once, after every repair is complete.
    base = t * sum(4 * d - 3 for d in sizes) + t * (c + 1)

    # One DFS expansion selects a frontier row, reads its w owner records,
    # tests free candidates, inserts previously unseen owners, and records
    # their parent row and partition.  Common controls are packed as wide
    # shared-bit products.
    forward_step = (
        (t - 1) * (2 + w * (h + 4))
        + h * t
        + 4 * w
    )
    forward_search = t * forward_step

    # Backtracking reads one parent record per possible hop, emits the new
    # row-to-partition assignments, and clears every old path assignment.
    backtrack = (
        (w - 1) * (t - 1)
        + (t - 1) * ((h + 2 * w + 1) * (t - 1) + 1)
        + t * w
    )

    one_augmentation = (
        # Materialize occupancy and h-bit owner arrays from the current match.
        t * sum(d - 1 for d in sizes)
        # Read the occupancy/owner record at every candidate address.
        + t * (h + 1) * sum(d - 1 for d in sizes)
        + forward_search
        + backtrack
    )
    return base + cfg.bfs_roots * one_augmentation


def reachability_generator_profile(
    model: Model, cfg: Configuration
) -> MultiplicationProfile:
    """OT-minimized implicit-owner DFS with fixed public loop bounds."""

    if cfg.bfs_roots <= 0:
        raise ValueError("reachability configurations require a positive root bound")
    t, w = model.t, cfg.w
    index_bits = cfg.partition_index_bits
    r = cfg.bfs_roots
    h = max(1, (t - 1).bit_length())

    base = t * sum(2 * k + 2 for k in index_bits) + t
    forward_step = t + 2 * h + w * h + 3 * w + 2
    backtrack = (w - 1) + (t - 1) * (2 * h + w) + t
    one_augmentation = (
        w * t
        + t * sum(index_bits)
        + t * forward_step
        + backtrack
    )
    return MultiplicationProfile(
        count=base + r * one_augmentation,
        width_sum=reachability_generator_cost(model, cfg),
    )


def generator_cost(model: Model, cfg: Configuration) -> int:
    if cfg.generator == "basic":
        return basic_generator_cost(model, cfg)
    if cfg.generator == "one-repair":
        return one_repair_generator_cost(model, cfg)
    if cfg.generator == "eviction":
        return eviction_generator_cost(model, cfg)
    if cfg.generator == "bfs":
        return bfs_generator_cost(model, cfg)
    if cfg.generator == "reachability":
        return reachability_generator_cost(model, cfg)
    raise ValueError(f"unknown generator: {cfg.generator}")


def generator_profile(model: Model, cfg: Configuration) -> MultiplicationProfile:
    if cfg.generator == "basic":
        return basic_generator_profile(model, cfg)
    if cfg.generator == "one-repair":
        return one_repair_generator_profile(model, cfg)
    if cfg.generator == "reachability":
        return reachability_generator_profile(model, cfg)

    # The experimental eviction and bounded-BFS circuits have not yet been
    # scheduled for wide products.  Counting every product as width one gives
    # a conservative OT-count upper bound while preserving their bit payload.
    width_sum = generator_cost(model, cfg)
    return MultiplicationProfile(count=width_sum, width_sum=width_sum)


@functools.cache
def waksman_switch_count(size: int) -> int:
    """Number of hidden controls in the implementation's Waksman network."""

    if size <= 1:
        return 0
    if size == 2:
        return 1
    lower = size // 2
    upper = size - lower
    return 2 * lower + waksman_switch_count(lower) + waksman_switch_count(upper)


def dpf_initialization_cost(model: Model, columns: int) -> OtCost:
    """Reservation bound for column and row-to-column reusable DPFs.

    Match WaterfallDmpf's dense-prefix choice. Column correction selections
    use two kappa-bit OTs per level. Residual sparse levels additionally mask
    both children using two 2*kappa-bit OTs. Row-to-column trees are dense.
    Non-characteristic-two leaf updates cache one sign product per DPF.
    Empty levels can reduce consumption, but not this reservation bound.
    """
    if columns < 1:
        raise ValueError("there must be at least one column")
    row_depth = (columns - 1).bit_length()
    dense_depth = min(model.domain_bits, row_depth + 2)
    correction_levels = columns * model.domain_bits + model.t * row_depth
    masked_levels = columns * (model.domain_bits - dense_depth)
    sign_products = 0 if model.characteristic_two else columns + model.t
    return OtCost(
        ot_count=2 * (correction_levels + masked_levels + sign_products),
        payload_bits=(2 * model.kappa * correction_levels
                      + 4 * model.kappa * masked_levels
                      + 2 * model.payload_bits * sign_products),
    )


def configuration_cost(model: Model, cfg: Configuration) -> Cost:
    common_profile = dedup_and_activity_profile(model)
    generator = generator_profile(model, cfg)
    hash_eval = direct_hash_profile(model, cfg)

    sparse_dpf = dpf_initialization_cost(model, cfg.columns)

    # Two serial Waksman passes are used.  Each party privately controls one
    # pass, so a switch costs one OT rather than the two OTs needed for a shared
    # control.  The inverse reuses the correlated OT with domain separation.
    switches = waksman_switch_count(cfg.columns)
    scatter = OtCost(
        ot_count=2 * switches,
        payload_bits=2 * switches * model.scatter_record_bits,
    )
    return Cost(
        common=OtCost(common_profile.ot_count, common_profile.payload_bits),
        generator=OtCost(generator.ot_count, generator.payload_bits),
        hash_eval=OtCost(hash_eval.ot_count, hash_eval.payload_bits),
        sparse_dpf=sparse_dpf,
        scatter=scatter,
    )


def default_configurations(model: Model) -> tuple[Configuration, ...]:
    return (
        Configuration("basic", "basic", 2, 512, 3),
        Configuration("basic", "basic", 3, 128, 2),
        Configuration("basic", "basic", 4, 64, 1),
        Configuration("basic", "basic", 5, 64, 0),
        Configuration("one repair", "one-repair", 3, 64, 2),
        Configuration("one repair", "one-repair", 4, 32, 1),
        Configuration("one repair", "one-repair", 5, 16, 0),
        Configuration("one repair", "one-repair", 3, 32, 3),
        Configuration("one repair", "one-repair", 4, 16, 2),
        Configuration(
            "bounded eviction", "eviction", 4, 16, 0,
            bfs_depth=8, bfs_roots=3, proved=False,
        ),
        Configuration(
            "bounded BFS", "bfs", 4, 16, 0,
            bfs_depth=4, bfs_roots=3, proved=False,
        ),
        Configuration(
            "complete BFS", "bfs", 4, 16, 0,
            bfs_depth=model.t, bfs_roots=3,
        ),
    )


def print_table(model: Model) -> None:
    header = (
        "generator",
        "w",
        "d",
        "c",
        "m",
        "G_ots",
        "hash_ots",
        "DPF_ots",
        "scatter_ots",
        "total_ots",
        "payload_bits",
        "proved",
    )
    print(",".join(header))
    for cfg in default_configurations(model):
        cost = configuration_cost(model, cfg)
        print(
            ",".join(
                map(
                    str,
                    (
                        cfg.name,
                        cfg.w,
                        cfg.d,
                        cfg.c,
                        cfg.columns,
                        cost.generator.ot_count,
                        cost.hash_eval.ot_count,
                        cost.sparse_dpf.ot_count,
                        cost.scatter.ot_count,
                        cost.total_ot_count,
                        cost.total_payload_bits,
                        "yes" if cfg.proved else "no",
                    ),
                )
            )
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain-bits", type=int, default=20)
    parser.add_argument("--payload-bits", type=int, default=64)
    parser.add_argument("--kappa", type=int, default=128)
    args = parser.parse_args()
    print_table(
        Model(
            domain_bits=args.domain_bits,
            payload_bits=args.payload_bits,
            kappa=args.kappa,
        )
    )


if __name__ == "__main__":
    main()
