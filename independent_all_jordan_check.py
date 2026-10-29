"""Independent exact full determinant checks; Python standard library only."""
from fractions import Fraction
from pathlib import Path
import json


def partitions(n, cap=None):
    if not n:
        yield ()
        return
    for a in range(min(n, cap or n), 0, -1):
        for rest in partitions(n-a, a):
            yield (a,) + rest


def multiply(p, q):
    out = {}
    for (a, b), x in p.items():
        for (c, d), y in q.items():
            key = (a+c, b+d)
            out[key] = out.get(key, 0) + x*y
    return {k: v for k, v in out.items() if v}


def full_determinant(entries):
    n = len(entries)
    out = {}

    def visit(i, cols, product, inversions):
        if i == n:
            sign = (-1)**inversions
            for k, v in product.items():
                out[k] = out.get(k, 0) + sign*v
            return
        for j in range(n):
            if j not in cols and entries[i][j]:
                visit(i+1, cols + (j,), multiply(product, entries[i][j]),
                      inversions + sum(c > j for c in cols))

    visit(0, (), {(0, 0): 1}, 0)
    return {k: v for k, v in out.items() if v}


def exact_rank(a):
    x = [[Fraction(v) for v in row] for row in a]
    n, rank = len(x), 0
    for col in range(n):
        pivot = next((j for j in range(rank, n) if x[j][col]), None)
        if pivot is None:
            continue
        x[pivot], x[rank] = x[rank], x[pivot]
        z = x[rank][col]
        x[rank] = [v/z for v in x[rank]]
        for j in range(rank+1, n):
            z = x[j][col]
            x[j] = [v-z*w for v, w in zip(x[j], x[rank])]
        rank += 1
    return rank


def main():
    checks = []
    for n in range(2, 10):
        for ps in partitions(n):
            if max(ps) == 1:
                continue
            a = [[0]*n for _ in range(n)]
            b = [[0]*n for _ in range(n)]
            starts, ends, k = [], [], 0
            for p in ps:
                starts.append(k)
                ends.append(k+p-1)
                for j in range(1, p):
                    a[k+j-1][k+j] = 1
                k += p
            g = len(ps)
            for i in range(g):
                b[ends[(i+1) % g]][starts[i]] = 1
            i = next(j for j, p in enumerate(ps) if p >= 2)
            b[ends[(i+1) % g]][starts[i]+1] = 1
            entries = [[{} for _ in range(n)] for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i == j:
                        entries[i][j][(1, 0)] = 1
                    if a[i][j]:
                        entries[i][j][(0, 0)] = -a[i][j]
                    if b[i][j]:
                        entries[i][j][(0, 1)] = -b[i][j]
            got = full_determinant(entries)
            expected = {(n, 0): 1, (0, g): -1, (1, g): -1}
            if got != expected:
                raise ValueError((ps, got, expected))
            rank = exact_rank(b)
            nonzeros = sum(bool(v) for row in b for v in row)
            if rank != g or nonzeros != g+1:
                raise ValueError((ps, rank, nonzeros))
            gram = [[sum(b[i][k]*b[j][k] for k in range(n))
                     for j in range(n)] for i in range(n)]
            if any(gram[i][j] for i in range(n) for j in range(n) if i != j):
                raise ValueError((ps, 'non-diagonal row Gram matrix'))
            if sorted(gram[i][i] for i in range(n)) != [0]*(n-g)+[1]*(g-1)+[2]:
                raise ValueError((ps, 'incorrect row Gram eigenvalues'))
            checks.append({'N': n, 'partition': ps, 'rank': rank,
                           'nonzero_entries': nonzeros})
    report = {'status': 'passed', 'cases': len(checks),
              'arithmetic': 'exact integer polynomials and rational row reduction',
              'scope': 'all nonzero nilpotent Jordan partitions in dimensions 2 through 9',
              'checks': checks}
    path = Path(__file__).with_name('all_jordan_symbolic_audit.json')
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(f'Independent full bivariate determinant, rank and norm checks: {len(checks)} passed')


if __name__ == '__main__':
    main()
