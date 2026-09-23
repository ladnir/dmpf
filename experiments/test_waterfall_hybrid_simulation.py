"""Finite checks of the Waterfall hybrid's distributional identities.

These tests concern ideal-call outputs and the explicit opening. They do not
test a concrete MPC/DPF realization or inspect an ideal subprotocol's state.
The output-fiber checks preserve the former uniform-output argument as
historical algebra; the current functionality instead permits corrupt share
selection. Its aggregation checks are in test_reusable_dpf_sharing_contract.
"""

from collections import Counter, defaultdict
from itertools import permutations, product

import pytest


def compose(left, right):
    """Return left after right, as in sigma_1 composed with sigma_0."""
    return tuple(left[right[i]] for i in range(len(right)))


def complete_outputs(domain, sparse_sets, output, free, modulus):
    """Reference fiber sampler with deterministic first-containing pivots."""
    pivots = {x: next(j for j, subset in enumerate(sparse_sets) if x in subset)
              for x in domain}
    entries = tuple((j, x) for j, subset in enumerate(sparse_sets) for x in subset)
    free_entries = tuple(entry for entry in entries if entry[0] != pivots[entry[1]])
    values = dict(zip(free_entries, free, strict=True))
    for x, target in zip(domain, output, strict=True):
        values[pivots[x], x] = (
            target - sum(value for (j, address), value in values.items() if address == x)
        ) % modulus
    return tuple(values[entry] for entry in entries)


def test_opening_uniform_for_either_fixed_party_permutation():
    size, rows = 4, 2
    all_permutations = tuple(permutations(range(size)))
    injections = tuple(permutations(range(size), rows))
    for party in (0, 1):
        for known in all_permutations:
            for placement in injections:
                counts = Counter()
                for hidden in all_permutations:
                    combined = compose(hidden, known) if party == 0 else compose(known, hidden)
                    counts[tuple(combined[column] for column in placement)] += 1
                assert counts == Counter({injection: 2 for injection in injections})


def test_opening_message_completes_the_corrupted_share():
    # Matrices are flattened bit strings; the opening message is V XOR O_p.
    target = 0b01001000
    for local_share in range(1 << 8):
        honest_message = target ^ local_share
        assert local_share ^ honest_message == target


def test_unmatched_logical_labels_survive_column_permutation():
    partial = (0, None, 2)
    for permutation in permutations(range(4)):
        opened = tuple(None if column is None else permutation[column] for column in partial)
        assert tuple(i for i, column in enumerate(opened) if column is None) == (1,)


@pytest.mark.parametrize('modulus', (2, 3, 4))
def test_column_sampler_matches_every_conditional_fiber(modulus):
    domain = (0, 1, 2)
    sparse_sets = ((0, 1), (2,), (0, 1, 2), ())
    entries = tuple((j, x) for j, subset in enumerate(sparse_sets) for x in subset)
    fibers = defaultdict(set)
    for values in product(range(modulus), repeat=len(entries)):
        output = tuple(sum(value for (_, address), value in zip(entries, values)
                           if address == x) % modulus for x in domain)
        fibers[output].add(values)
    free_count = len(entries) - len(domain)
    for output in product(range(modulus), repeat=len(domain)):
        sampled = Counter(complete_outputs(domain, sparse_sets, output, free, modulus)
                          for free in product(range(modulus), repeat=free_count))
        assert set(sampled) == fibers[output]
        assert set(sampled.values()) == {1}
        assert len(sampled) == modulus ** free_count


def test_completion_without_free_entries_and_with_nonconsecutive_addresses():
    assert complete_outputs((2, 7), ((), (2,), (7,)), (1, 2), (), 3) == (1, 2)


def test_fresh_completion_after_an_adaptively_selected_output():
    domain, sparse_sets = (0, 1), ((0, 1), (0, 1))
    for first_free in product(range(3), repeat=2):
        first = complete_outputs(domain, sparse_sets, (0, 1), first_free, 3)
        next_output = ((first[0] + first[3]) % 3, (first[1] + 2 * first[2]) % 3)
        second = [complete_outputs(domain, sparse_sets, next_output, free, 3)
                  for free in product(range(3), repeat=2)]
        assert len(set(second)) == 9
        assert all(((values[0] + values[2]) % 3, (values[1] + values[3]) % 3)
                   == next_output for values in second)
