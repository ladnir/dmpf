#!/usr/bin/env python3
"""Reference implementations of complete Waterfall repair.

The exact-distance routines provide a simple shortest-path correctness oracle.
The multi-source routine emulates the optimized fixed circuit: every unmatched
row is an initial root, rows are expanded in public order, and first-discovery
parents define the selected augmenting path.  Both increase the matching size
by exactly one whenever an augmenting path exists.
"""

from __future__ import annotations

from waterfall_evictions import Configuration


def reachability_distances(
    candidates: list[int],
    configuration: Configuration,
    owner: list[int],
    placement: list[int],
) -> list[int]:
    """Return exact row distance to a free candidate, or -1 if unreachable."""
    t, w, d = configuration.t, configuration.w, configuration.d
    distance = [-1] * t
    for row in range(t):
        for partition in range(w):
            slot = partition * d + candidates[row * w + partition]
            if owner[slot] == -1:
                distance[row] = 0
                break

    # A simple augmenting path visits at most t rows, hence has at most t-1
    # occupied row-to-owner edges before its terminal free candidate.
    for layer in range(1, t):
        changed = False
        previous = distance.copy()
        for row in range(t):
            if previous[row] != -1:
                continue
            base = row * w
            for partition in range(w):
                slot = partition * d + candidates[base + partition]
                victim = owner[slot]
                if victim == row:
                    continue
                if victim != -1 and previous[victim] != -1:
                    distance[row] = layer
                    changed = True
                    break
        if not changed:
            break
    return distance


def augment_first_reachable_root(
    candidates: list[int],
    configuration: Configuration,
    owner: list[int],
    placement: list[int],
) -> bool:
    """Apply one canonical shortest augmentation and report whether it exists."""
    t, w, d = configuration.t, configuration.w, configuration.d
    distance = reachability_distances(
        candidates, configuration, owner, placement
    )
    root = next(
        (
            row
            for row in range(t)
            if placement[row] == -1 and distance[row] != -1
        ),
        -1,
    )
    if root == -1:
        return False

    path_rows = [root]
    path_slots: list[int] = []
    row = root
    while distance[row] > 0:
        next_row = -1
        next_slot = -1
        base = row * w
        for partition in range(w):
            slot = partition * d + candidates[base + partition]
            victim = owner[slot]
            if victim != -1 and distance[victim] == distance[row] - 1:
                next_row = victim
                next_slot = slot
                break
        if next_row == -1:
            raise AssertionError("finite distance has no decreasing edge")
        path_slots.append(next_slot)
        path_rows.append(next_row)
        row = next_row

    terminal_slot = next(
        (
            partition * d + candidates[row * w + partition]
            for partition in range(w)
            if owner[partition * d + candidates[row * w + partition]] == -1
        ),
        -1,
    )
    if terminal_slot == -1:
        raise AssertionError("distance-zero row has no free candidate")
    path_slots.append(terminal_slot)

    # Apply from the terminal backwards so every destination is free when used.
    for path_row, slot in reversed(list(zip(path_rows, path_slots))):
        old_slot = placement[path_row]
        if owner[slot] != -1:
            raise AssertionError("augmenting path destination is not free")
        owner[slot] = path_row
        placement[path_row] = slot
        if old_slot != -1:
            owner[old_slot] = -1
    return True


def complete_reachability_repair(
    candidates: list[int],
    configuration: Configuration,
    initial_owner: list[int],
    initial_placement: list[int],
    augmentation_bound: int,
) -> int:
    """Return the unplaced count after bounded complete augmentations."""
    if augmentation_bound < 0:
        raise ValueError("augmentation_bound must be nonnegative")
    owner = initial_owner.copy()
    placement = initial_placement.copy()
    for _ in range(augmentation_bound):
        if sum(slot == -1 for slot in placement) <= configuration.c:
            break
        if not augment_first_reachable_root(
            candidates, configuration, owner, placement
        ):
            break
    return sum(slot == -1 for slot in placement)


def augment_from_all_roots_dfs(
    candidates: list[int],
    configuration: Configuration,
    owner: list[int],
    placement: list[int],
) -> bool:
    """Apply one complete multi-source reachability augmentation.

    Every unmatched row is an initial root.  Rows are expanded in public row
    order, candidates are inspected in public partition order, and each owner
    is visited once.  The parent forest therefore has at most ``t`` rows and
    can be traversed by a fixed ``t``-step circuit.
    """
    t, w, d = configuration.t, configuration.w, configuration.d
    visited = [slot == -1 for slot in placement]
    expanded = [False] * t
    parent_row = [-1] * t
    parent_slot = [-1] * t
    terminal_row = -1
    terminal_slot = -1

    for _ in range(t):
        current = next(
            (
                row
                for row in range(t)
                if visited[row] and not expanded[row]
            ),
            -1,
        )
        if current == -1:
            break
        expanded[current] = True
        for partition in range(w):
            slot = partition * d + candidates[current * w + partition]
            victim = owner[slot]
            if victim == -1:
                terminal_row = current
                terminal_slot = slot
                break
            if not visited[victim]:
                visited[victim] = True
                parent_row[victim] = current
                parent_slot[victim] = slot
        if terminal_row != -1:
            break

    if terminal_row == -1:
        return False

    assignments = {terminal_row: terminal_slot}
    child = terminal_row
    while parent_row[child] != -1:
        parent = parent_row[child]
        assignments[parent] = parent_slot[child]
        child = parent

    for row in assignments:
        old_slot = placement[row]
        if old_slot != -1:
            if owner[old_slot] != row:
                raise AssertionError("placement and owner arrays disagree")
            owner[old_slot] = -1
    for row, slot in assignments.items():
        if owner[slot] != -1:
            raise AssertionError("augmentation destinations are not distinct")
        owner[slot] = row
        placement[row] = slot
    return True


def implicit_owner_dfs_repair(
    candidates: list[int],
    configuration: Configuration,
    initial_owner: list[int],
    initial_placement: list[int],
    augmentation_bound: int,
) -> int:
    """Return the unplaced count after bounded multi-source repairs."""
    if augmentation_bound < 0:
        raise ValueError("augmentation_bound must be nonnegative")
    owner = initial_owner.copy()
    placement = initial_placement.copy()
    for _ in range(augmentation_bound):
        if sum(slot == -1 for slot in placement) <= configuration.c:
            break
        if not augment_from_all_roots_dfs(
            candidates, configuration, owner, placement
        ):
            break
    return sum(slot == -1 for slot in placement)


def circuit_reachability_repair(
    candidates: list[int],
    configuration: Configuration,
    initial_placement: list[int],
    augmentation_bound: int,
) -> int:
    """Emulate the Boolean equations in the paper literally.

    This deliberately keeps the dense ``e``, ``G``, ``rho``, ``delta``,
    ``x``, and ``p`` wires.  It is an executable specification for auditing
    the fixed circuit, not the implementation to use in a hot path.
    """
    if augmentation_bound < 0:
        raise ValueError("augmentation_bound must be nonnegative")
    t, w, d = configuration.t, configuration.w, configuration.d
    if len(candidates) != t * w or len(initial_placement) != t:
        raise ValueError("candidate or placement dimensions do not match")

    mu = [[0] * w for _ in range(t)]
    for row, slot in enumerate(initial_placement):
        if slot == -1:
            continue
        partition, bin_index = divmod(slot, d)
        if not 0 <= partition < w:
            raise ValueError("initial placement lies outside the main table")
        if candidates[row * w + partition] != bin_index:
            raise ValueError("initial placement is not a candidate of its row")
        mu[row][partition] = 1

    def check_valid_matching(state: list[list[int]]) -> None:
        occupied: set[tuple[int, int]] = set()
        for row in range(t):
            if sum(state[row]) > 1:
                raise AssertionError("a row occupies multiple positions")
            for partition in range(w):
                if not state[row][partition]:
                    continue
                position = (partition, candidates[row * w + partition])
                if position in occupied:
                    raise AssertionError("two rows occupy the same position")
                occupied.add(position)

    check_valid_matching(mu)
    for _ in range(augmentation_bound):
        q = [1 ^ (sum(mu[row]) & 1) for row in range(t)]
        e = [
            [
                [
                    int(
                        candidates[row * w + partition]
                        == candidates[owner * w + partition]
                    )
                    & mu[owner][partition]
                    for owner in range(t)
                ]
                for partition in range(w)
            ]
            for row in range(t)
        ]
        free = [
            [1 ^ (sum(e[row][partition]) & 1) for partition in range(w)]
            for row in range(t)
        ]
        graph = [
            [
                int(
                    row != owner
                    and any(e[row][partition][owner] for partition in range(w))
                )
                for owner in range(t)
            ]
            for row in range(t)
        ]

        rho = [[int(any(free[row])) for row in range(t)]]
        delta = [rho[0].copy()]
        for _layer in range(1, t):
            previous = rho[-1]
            transition = [
                int(
                    any(
                        graph[row][owner] and previous[owner]
                        for owner in range(t)
                        if owner != row
                    )
                )
                for row in range(t)
            ]
            exact = [
                transition[row] & (1 ^ previous[row]) for row in range(t)
            ]
            delta.append(exact)
            rho.append(
                [previous[row] ^ exact[row] for row in range(t)]
            )

        eligible = [q[row] & rho[-1][row] for row in range(t)]
        first_eligible = next(
            (row for row, bit in enumerate(eligible) if bit), -1
        )
        g = [int(row == first_eligible) for row in range(t)]
        x = [
            [g[row] & delta[layer][row] for row in range(t)]
            for layer in range(t)
        ]
        p = [[0] * w for _ in range(t)]

        for layer in range(t - 1, 0, -1):
            alpha = sum(x[layer]) & 1
            selected_pair: tuple[int, int] | None = None
            for partition in range(w):
                selected_label = sum(
                    x[layer][row] * candidates[row * w + partition]
                    for row in range(t)
                )
                for owner in range(t):
                    v = (
                        alpha
                        & delta[layer - 1][owner]
                        & mu[owner][partition]
                        & int(
                            selected_label
                            == candidates[owner * w + partition]
                        )
                    )
                    if v and selected_pair is None:
                        selected_pair = (partition, owner)
            if selected_pair is not None:
                partition, owner = selected_pair
                p[layer][partition] = 1
                x[layer - 1][owner] ^= 1

        selected_terminal_partition = next(
            (
                partition
                for partition in range(w)
                if any(
                    x[0][row] and free[row][partition]
                    for row in range(t)
                )
            ),
            -1,
        )
        if selected_terminal_partition != -1:
            p[0][selected_terminal_partition] = 1

        next_mu = [[0] * w for _ in range(t)]
        for row in range(t):
            y = sum(x[layer][row] for layer in range(t)) & 1
            for partition in range(w):
                new_position = 0
                for layer in range(t):
                    new_position ^= x[layer][row] & p[layer][partition]
                next_mu[row][partition] = (
                    mu[row][partition]
                    ^ (mu[row][partition] & y)
                    ^ new_position
                )
        mu = next_mu
        check_valid_matching(mu)

    return sum(1 ^ (sum(mu[row]) & 1) for row in range(t))
