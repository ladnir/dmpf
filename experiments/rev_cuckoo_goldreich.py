"""Implementation-shaped shared-Goldreich Reverse-Cuckoo leakage channel.

The C++ implementation samples one public Goldreich root per cuckoo partition
and reuses those roots for every set in a batch.  A set reveals fresh affine-
fiber solutions under the shared feature matrices.  This module reproduces
that distribution while exposing the same descriptor-table interface used by
the ideal-channel attack experiments.

For a candidate whose translated feature rows are independent, its likelihood
is exactly the usual Reverse-Cuckoo placement likelihood.  At reduced sizes we
can also enumerate placements and weight each one by the ranks of its two
assigned feature submatrices.  That is the exact coefficient-vector likelihood
and includes feature-rank leakage rather than silently treating it as an abort.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from dataclasses import dataclass

from rev_cuckoo_full_lpn import NEG_INF, falling_factorial
from rev_cuckoo_regular_noise import Descriptor, LeakageChannel, descriptor_deficit


def _parity_inner_product(left: int, right: int) -> int:
    return (left & right).bit_count() & 1


def gf2_rank(rows: Sequence[int]) -> int:
    """Rank of integer-packed binary rows."""

    basis: dict[int, int] = {}
    for source in rows:
        row = source
        while row:
            pivot = row.bit_length() - 1
            previous = basis.get(pivot)
            if previous is None:
                basis[pivot] = row
                break
            row ^= previous
    return len(basis)


def sample_affine_solution(
    rows: Sequence[int],
    right_hand_sides: Sequence[int],
    columns: int,
    output_bits: int,
    rng: random.Random,
) -> tuple[int, ...] | None:
    """Uniformly sample S with ``rows[i] * S = right_hand_sides[i]``.

    The returned tuple stores one ``output_bits``-wide coefficient for every
    feature column.  ``None`` denotes an inconsistent dependent equation.
    """

    if len(rows) != len(right_hand_sides):
        raise ValueError("row and right-hand-side counts differ")
    mask = (1 << output_bits) - 1
    basis: dict[int, tuple[int, int]] = {}
    for source, source_rhs in zip(rows, right_hand_sides):
        row = source
        rhs = source_rhs & mask
        while row:
            pivot = row.bit_length() - 1
            previous = basis.get(pivot)
            if previous is None:
                basis[pivot] = (row, rhs)
                break
            row ^= previous[0]
            rhs ^= previous[1]
        if row == 0 and rhs:
            return None

    solution = [rng.getrandbits(output_bits) for _ in range(columns)]
    for pivot in sorted(basis):
        row, rhs = basis[pivot]
        lower = row & ((1 << pivot) - 1)
        value = rhs
        while lower:
            bit = lower & -lower
            value ^= solution[bit.bit_length() - 1]
            lower ^= bit
        solution[pivot] = value
    return tuple(solution)


def evaluate_solution(feature: int, solution: Sequence[int]) -> int:
    result = 0
    while feature:
        bit = feature & -feature
        result ^= solution[bit.bit_length() - 1]
        feature ^= bit
    return result


@dataclass(frozen=True)
class GoldreichFeatureRoot:
    domain: int
    input_bits: int
    intermediate_bits: int
    feature_bits: int
    partitions: tuple[tuple[int, ...], tuple[int, ...]]

    @classmethod
    def sample(
        cls,
        domain: int,
        feature_bits: int,
        rng: random.Random,
    ) -> "GoldreichFeatureRoot":
        if domain < 1 or feature_bits < 1:
            raise ValueError("domain and feature width must be positive")
        # This matches the byte-rounded dimensions in GoldreichHash::init.
        input_bits = 8 * max(1, (math.ceil(math.log2(domain + 1)) + 7) // 8)
        intermediate_bits = 8 * ((feature_bits + 7) // 8)
        input_mask_bits = (1 << input_bits) - 1
        combined_bits = input_bits + intermediate_bits

        partitions = []
        for _ in range(2):
            linear0 = tuple(
                rng.getrandbits(input_bits) for _ in range(intermediate_bits)
            )
            linear1 = tuple(
                rng.getrandbits(input_bits) for _ in range(intermediate_bits)
            )
            compression = tuple(
                rng.getrandbits(combined_bits) for _ in range(feature_bits)
            )

            def raw_feature(value: int) -> int:
                x = value & input_mask_bits
                left = sum(
                    _parity_inner_product(x, row) << bit
                    for bit, row in enumerate(linear0)
                )
                right = sum(
                    _parity_inner_product(x, row) << bit
                    for bit, row in enumerate(linear1)
                )
                nonlinear = left & right
                combined = x | (nonlinear << input_bits)
                return sum(
                    _parity_inner_product(combined, row) << bit
                    for bit, row in enumerate(compression)
                )

            dummy = raw_feature(domain)
            partitions.append(
                tuple(raw_feature(value) ^ dummy for value in range(domain))
            )

        return cls(
            domain,
            input_bits,
            intermediate_bits,
            feature_bits,
            (partitions[0], partitions[1]),
        )


class FeatureRankError(RuntimeError):
    pass


@dataclass(frozen=True)
class _GoldreichProposal:
    solutions: tuple[tuple[int, ...], tuple[int, ...]]


@dataclass
class SharedGoldreichChannel:
    """Shared-feature wrapper around an ideal finite-K leakage channel."""

    base: LeakageChannel
    linear_security: int = 10
    exact_rank_limit: int = 10
    root: GoldreichFeatureRoot | None = None
    sampled_systems: int = 0
    rank_deficient_systems: int = 0
    rank_failures: int = 0
    approximate_rank_queries: int = 0

    def __post_init__(self) -> None:
        if self.partition_size & (self.partition_size - 1):
            raise ValueError(
                "the implemented Goldreich channel requires a power-of-two "
                "partition size"
            )

    @property
    def mode(self) -> str:
        return self.base.mode

    @property
    def selector_k(self) -> int:
        return self.base.selector_k

    @property
    def partition_size(self) -> int:
        return self.base.partition_size

    @property
    def feature_bits(self) -> int:
        return self.partition_size + self.linear_security

    def begin_batch(self, domain: int, rng: random.Random) -> None:
        self.root = GoldreichFeatureRoot.sample(domain, self.feature_bits, rng)

    def _require_root(self, domain: int | None = None) -> GoldreichFeatureRoot:
        if self.root is None:
            raise RuntimeError("begin_batch must be called before sampling")
        if domain is not None and self.root.domain != domain:
            raise ValueError("descriptor domain differs from the shared root")
        return self.root

    def _sample_proposal(
        self,
        active: Sequence[int],
        domain: int,
        rng: random.Random,
    ) -> _GoldreichProposal:
        root = self._require_root(domain)
        if len(active) > 2 * self.partition_size:
            raise ValueError("more active items than cuckoo bins")
        assigned = rng.sample(range(2 * self.partition_size), len(active))
        solutions = []
        for side in range(2):
            rows = []
            labels = []
            for item, position in zip(active, assigned):
                if position // self.partition_size == side:
                    rows.append(root.partitions[side][item])
                    labels.append(position % self.partition_size)
            self.sampled_systems += 1
            if gf2_rank(rows) != len(rows):
                self.rank_deficient_systems += 1
            solution = sample_affine_solution(
                rows,
                labels,
                self.feature_bits,
                math.ceil(math.log2(self.partition_size)),
                rng,
            )
            if solution is None:
                self.rank_failures += 1
                raise FeatureRankError(
                    "the planted Goldreich feature equations are inconsistent"
                )
            solutions.append(solution)
        return _GoldreichProposal((solutions[0], solutions[1]))

    def _materialize(self, proposal: _GoldreichProposal) -> Descriptor:
        root = self._require_root()
        tables = tuple(
            tuple(
                evaluate_solution(feature, proposal.solutions[side])
                for feature in root.partitions[side]
            )
            for side in range(2)
        )
        return tables[0], tables[1]

    def _proposal_deficit(
        self, active: Sequence[int], proposal: _GoldreichProposal
    ) -> int:
        root = self._require_root()
        labels = (
            {
                evaluate_solution(
                    root.partitions[side][item], proposal.solutions[side]
                )
                for item in active
            }
            for side in range(2)
        )
        occupied = sum(len(side_labels) for side_labels in labels)
        return 2 * len(active) - occupied

    def sample(
        self,
        active: Sequence[int],
        domain: int,
        rng: random.Random,
    ) -> Descriptor:
        if self.base.mode == "raw":
            return self._materialize(
                self._sample_proposal(active, domain, rng)
            )
        if self.base.mode != "occupancy":
            raise ValueError("shared Goldreich sampling supports raw/occupancy")
        proposals = [
            self._sample_proposal(active, domain, rng)
            for _ in range(self.base.selector_k)
        ]
        rates = self.base.occupancy_rates[len(active)]
        weights = [
            rates[self._proposal_deficit(active, proposal)]
            for proposal in proposals
        ]
        threshold = rng.random() * sum(weights)
        prefix = 0.0
        for proposal, weight in zip(proposals, weights):
            prefix += weight
            if threshold <= prefix:
                return self._materialize(proposal)
        return self._materialize(proposals[-1])

    def sample_batch(
        self,
        active_sets: Sequence[Sequence[int]],
        domain: int,
        rng: random.Random,
        retries: int = 8,
    ) -> tuple[Descriptor, ...]:
        """Sample one shared root and all descriptor solutions atomically."""

        for _ in range(retries):
            self.begin_batch(domain, rng)
            try:
                return tuple(
                    self.sample(active, domain, rng) for active in active_sets
                )
            except FeatureRankError:
                continue
        raise FeatureRankError("shared Goldreich batch exceeded rank retries")

    def _exact_raw_log_likelihood(
        self, active: Sequence[int], descriptor: Descriptor
    ) -> float:
        root = self._require_root()
        items = tuple(active)
        if len(items) > self.exact_rank_limit:
            self.approximate_rank_queries += 1
            return self.base.log_likelihood(active, descriptor)

        left, right = descriptor
        rank_weight = 0
        for orientation in range(1 << len(items)):
            used_left = 0
            used_right = 0
            left_rows = []
            right_rows = []
            valid = True
            for index, item in enumerate(items):
                if orientation >> index & 1:
                    label = right[item]
                    bit = 1 << label
                    if used_right & bit:
                        valid = False
                        break
                    used_right |= bit
                    right_rows.append(root.partitions[1][item])
                else:
                    label = left[item]
                    bit = 1 << label
                    if used_left & bit:
                        valid = False
                        break
                    used_left |= bit
                    left_rows.append(root.partitions[0][item])
            if valid:
                exponent = math.ceil(math.log2(self.partition_size)) * (
                    gf2_rank(left_rows) + gf2_rank(right_rows)
                )
                rank_weight += 1 << exponent
        if rank_weight == 0:
            return NEG_INF
        return math.log2(rank_weight) - math.log2(
            falling_factorial(2 * self.partition_size, len(items))
        )

    def log_likelihood(
        self, active: Sequence[int], descriptor: Descriptor
    ) -> float:
        if self.root is None:
            # Pair-slope fitting samples ideal public endpoint tables.  The
            # good-rank Goldreich likelihood has the same slope.
            return self.base.log_likelihood(active, descriptor)
        raw = self._exact_raw_log_likelihood(active, descriptor)
        if raw == NEG_INF or self.base.mode == "raw":
            return raw
        # Under the good-rank event the endpoint/deficit law is exactly ideal.
        # The exact rank enumerator above corrects the raw density; the finite-K
        # race factor remains the ideal one.  This is the only approximation.
        deficit = descriptor_deficit(active, descriptor)
        return raw + math.log2(
            self.base.race_factors[len(active)][deficit]
        )

    def candidate_log_weight(
        self,
        candidate_lists: Sequence[Sequence[int]],
        descriptors: Sequence[Descriptor],
    ) -> float:
        result = 0.0
        for active, descriptor in zip(candidate_lists, descriptors):
            contribution = self.log_likelihood(active, descriptor)
            if contribution == NEG_INF:
                return NEG_INF
            result += contribution
        return result
