"""Conditional algebra for an output-completing DPF contract.

The honest active conversion is idealized as a uniform group element.
These checks do not establish that the actual cached-seed construction
admits this hybrid, nor prove privacy or general composition.
"""

from collections import Counter
from itertools import product

import pytest

from test_reusable_dpf_output_boundary import corrected_output, equal_tag_pair


def real_expansion(local, tags, alpha, beta, honest_pad, party, modulus):
    """Complete an idealized leaf state and perform the real scalar update."""
    other = tuple(honest_pad if x == alpha else -value
                  for x, value in enumerate(local))
    other_tags = tuple(tag ^ (x == alpha) for x, tag in enumerate(tags))
    sigma = (1 - 2 * party) * (tags[alpha] - other_tags[alpha])
    gamma = sigma * (beta - local[alpha] - honest_pad) % modulus
    corrupt = corrected_output(local, tags, gamma, party, modulus)
    honest = corrected_output(other, other_tags, gamma, 1 - party, modulus)
    return gamma, corrupt, honest


def simulated_expansion(local, tags, gamma, party, modulus):
    # Deliberately has no alpha, beta, or honest pad as an input.
    return gamma, corrected_output(local, tags, gamma, party, modulus)


@pytest.mark.parametrize("modulus", (2, 3, 4, 5, 8, 17))
def test_uniform_hidden_pad_matches_input_only_local_simulator(modulus):
    local = (2, 4, 7)
    for party in (0, 1):
        for tags in product((0, 1), repeat=3):
            ideal = Counter(simulated_expansion(local, tags, gamma, party, modulus)
                            for gamma in range(modulus))
            i, j = equal_tag_pair(tags)
            for alpha in range(3):
                for beta in range(modulus):
                    real = Counter()
                    for honest_pad in range(modulus):
                        gamma, corrupt, honest = real_expansion(
                            local, tags, alpha, beta, honest_pad, party, modulus)
                        real[gamma, corrupt] += 1
                        assert tuple((a + b) % modulus for a, b in
                                     zip(corrupt, honest, strict=True)) == tuple(
                                         beta if x == alpha else 0 for x in range(3))
                        assert (corrupt[i] - corrupt[j]) % modulus == (
                            local[i] - local[j]) % modulus
                    assert real == ideal


@pytest.mark.parametrize("modulus", (2, 3, 5))
def test_two_round_adaptive_payloads_match_joint_completed_outputs(modulus):
    # The second payload depends on the previous honest output as well as
    # the corrupt view. Neither simulator call needs to know that payload.
    tags = (0, 1, 0)
    for party in (0, 1):
        for alpha in range(3):
            real, ideal = Counter(), Counter()
            for randomness in product(range(modulus), repeat=2):
                for simulated, distribution in ((False, real), (True, ideal)):
                    history = []
                    beta = 1 % modulus
                    for round_number, random_element in enumerate(randomness):
                        local = tuple((3 * x + round_number) % modulus
                                      for x in range(3))
                        if simulated:
                            gamma, corrupt = simulated_expansion(
                                local, tags, random_element, party, modulus)
                            honest = tuple(((beta if x == alpha else 0) - value)
                                           % modulus for x, value in enumerate(corrupt))
                        else:
                            gamma, corrupt, honest = real_expansion(
                                local, tags, alpha, beta, random_element, party, modulus)
                        history.append((gamma, corrupt, honest))
                        beta = (honest[0] + 2 * corrupt[1] + gamma) % modulus
                    distribution[tuple(history)] += 1
            assert real == ideal


@pytest.mark.parametrize("modulus", (2, 3, 4, 17))
def test_overlapping_sparse_map_aggregation_commutes_with_completion(modulus):
    domains = ((0, 1, 2), (0, 2), (1, 2))
    functions = ({0: 1}, {0: 2}, {})  # Repeated target and an empty contribution.
    for mask in range(modulus):
        target = [0] * 3
        corrupt_sum = [0] * 3
        honest_sum = [0] * 3
        for domain, function in zip(domains, functions, strict=True):
            for x in domain:
                # Strongly correlated, nonuniform shares are intentional.
                corrupt = (mask + x) % modulus
                honest = (function.get(x, 0) - corrupt) % modulus
                target[x] += function.get(x, 0)
                corrupt_sum[x] += corrupt
                honest_sum[x] += honest
        assert tuple(value % modulus for value in honest_sum) == tuple(
            (value - share) % modulus for value, share in
            zip(target, corrupt_sum, strict=True))


@pytest.mark.parametrize("prime", (2, 3, 5, 17))
def test_ole_reverse_completion(prime):
    for corrupt_x, corrupt_z in product(range(prime), repeat=2):
        completed = [(u, (corrupt_x * u - corrupt_z) % prime)
                     for u in range(prime)]
        assert len(set(completed)) == prime
        for honest_x, honest_z in completed:
            assert (corrupt_z + honest_z) % prime == corrupt_x * honest_x % prime


def test_private_zero_sharing_does_not_promise_uniform_masks():
    # A fixed zero function with no private inputs has trivial input privacy.
    # Always returning (0,0) satisfies completion, but not uniform sharing.
    modulus = 5
    deterministic_shares = Counter({(0, 0): modulus})
    uniform_shares = Counter((r, -r % modulus) for r in range(modulus))
    assert all((a + b) % modulus == 0 for a, b in deterministic_shares)
    assert deterministic_shares != uniform_shares
