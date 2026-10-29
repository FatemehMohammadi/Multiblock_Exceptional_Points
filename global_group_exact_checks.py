"""Independent finite permutation-group checks for the proved hierarchy.

Standard-library only. These checks enumerate concrete finite groups; they
are supplementary checks, not a substitute for the all-m/all-r proofs.
"""
from collections import deque
from itertools import product
from math import factorial
import json


def compose(p, q):
    """Composition p after q, with zero-based images."""
    return tuple(p[q[j]] for j in range(len(p)))


def transposition(n, i, j):
    p = list(range(n))
    p[i], p[j] = p[j], p[i]
    return tuple(p)


def order(generators):
    n = len(generators[0])
    identity = tuple(range(n))
    seen = {identity}
    todo = deque([identity])
    while todo:
        p = todo.popleft()
        for g in generators:
            q = compose(g, p)
            if q not in seen:
                seen.add(q)
                todo.append(q)
    return len(seen)


def synchronous_shift(r, m):
    return tuple(j * m + (k + 1) % m for j in range(r) for k in range(m))


def main():
    records = []
    for m in range(2, 7):
        n = 2 * m
        a = synchronous_shift(2, m)
        b = transposition(n, 0, m)
        central = order([a])
        middle = order([a, b])
        assert central == m
        assert middle == 2**m * m
        rec = {"m": m, "central_order": central, "middle_order": middle}
        if m <= 3:
            c = compose(a, b)
            d = transposition(n, 0, c[0])
            full = order([a, b, d])
            assert full == factorial(n)
            rec["global_order"] = full
        records.append(rec)

    
    
    phase_records = []
    for r, m, edges in [
        (3, 2, [(0, 1), (1, 2)]),
        (3, 3, [(0, 1), (1, 2)]),
        (4, 2, [(0, 1), (0, 2), (0, 3)]),
        (4, 2, [(0, 1), (1, 2), (2, 3)]),
    ]:
        a = synchronous_shift(r, m)
        expected = factorial(r) ** m * m
        count = 0
        for phases in product(range(m), repeat=2 * (r - 1)):
            generators = [a]
            for j, (u, v) in enumerate(edges):
                generators.append(
                    transposition(r * m, u * m + phases[2*j],
                                  v * m + phases[2*j+1]))
            actual = order(generators)
            assert actual == expected, (r, m, edges, phases, actual)
            count += 1
        phase_records.append({"r": r, "m": m, "edges": edges,
                              "phase_assignments_checked": count,
                              "group_order": expected})
    print(json.dumps({"status": "passed", "rank_two": records,
                      "higher_rank_tree_phases": phase_records}, indent=2))


if __name__ == "__main__":
    main()
