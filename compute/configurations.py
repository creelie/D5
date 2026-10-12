"""Checks the statements of Section 6 of the paper about the five extremal configurations.

For m = 0..4 let Lambda_m and b_m be as in Table 1.1 (unscaled, minimal norm 2) and let a_1..a_5, b be
the vectors of Table 6.1.  The script checks, in exact arithmetic:

  1. the Gram matrix of sqrt2*(a_1, ..., a_5, b) is the equality point E_m of data/equality_points.json;
  2. a_1..a_5 lie in Lambda_m and |det(a_1..a_5)| = 4 = covol(Lambda_m), so they form a basis;
     b - b_m lies in Lambda_m;
  3. 2 b_m lies in Lambda_m exactly for m = 0, 3, 4, and then Lambda_m u (Lambda_m + b_m) = D_5;
  4. |Min Lambda_m| = 16, 20, 24, 24, 20 and every point of P_m = Lambda_m u (Lambda_m + b_m) has
     40 neighbours at squared distance 2 (none closer);
  5. the neighbours of 0 in P_2 are the 40 points of Leech's configuration L_5, reflected in x_5 = 0;
  6. the map x -> -(x_4, x_2, x_3, x_1, x_5) maps Lambda_1 onto the lattice of the Q_5 uniform packing of
     Cohn and Rajagopal (Table 4.2 of their paper) and b_1 onto their translation vector (0,0,0,-1,-1);
  7. D_5 has 31 sublattices of index 2; they are the kernels of x_i, x_i + x_j and (1/2) sum e_i x_i
     modulo 2 (5 + 10 + 16), with 24, 16 and 20 minimal vectors, of the types of Lambda_3, Lambda_0
     and Lambda_4.

Usage: python3 configurations.py ../data/equality_points.json
"""
import sys, json, itertools
from fractions import Fraction as F

def fr(s):
    return F(s)

def dot(u, v):
    return sum(x * y for x, y in zip(u, v))

def det(M):
    M = [row[:] for row in M]
    n, d = len(M), F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d

def inverse(M):
    n = len(M)
    A = [row[:] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]

def is_int(x):
    return x.denominator == 1

half = F(1, 2)

def in_L0(x):
    return all(map(is_int, x)) and (x[0] + x[1] + x[2]) % 2 == 0 and (x[3] + x[4]) % 2 == 0

def in_L1(x):   # A_4 + Z (4/5) 1
    k = sum(x) / 4
    if not is_int(k):
        return False
    return all(is_int(c - k * F(4, 5)) for c in x)

def in_L2(x):   # (D_4 x 0) + Z (h + 2 e_5), h = (1/2, 1/2, 1/2, 1/2, 0)
    k = x[4] / 2
    if not is_int(k):
        return False
    y = [x[i] - k * half for i in range(4)]
    return all(map(is_int, y)) and sum(y) % 2 == 0

def in_L3(x):
    return all(map(is_int, x)) and (x[0] + x[1] + x[2] + x[3]) % 2 == 0 and x[4] % 2 == 0

def in_L4(x):
    return all(map(is_int, x)) and sum(x) % 4 == 0

def in_D5(x):
    return all(map(is_int, x)) and sum(x) % 2 == 0

IN = [in_L0, in_L1, in_L2, in_L3, in_L4]
BM = [[1, 0, 0, 1, 0], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 0, 0, 0, 1], [1, 1, 0, 0, 0]]
BM = [[F(c) for c in v] for v in BM]

# Table 6.1: a_1, ..., a_5, b
TAB = {
    0: ["0 1 -1 0 0", "1 0 1 0 0", "1 -1 0 0 0", "0 0 0 1 -1", "0 0 0 1 1", "0 1 0 -1 0"],
    1: ["1/5 -9/5 -4/5 -4/5 -4/5", "-1 0 0 0 1", "-1 1 0 0 0", "0 -1 0 1 0", "0 -1 1 0 0",
        "1/5 -4/5 -4/5 -4/5 1/5"],
    2: ["1 0 1 0 0", "0 1 0 1 0", "0 0 -1 1 0", "0 0 -1 -1 0", "-1/2 -1/2 -1/2 -1/2 -2", "1 1 2 1 1"],
    3: ["0 1 0 -1 0", "1 0 -1 0 0", "1 -1 0 0 0", "0 0 1 1 0", "0 0 0 0 2", "0 1 -1 -1 -1"],
    4: ["-1 -1 -1 0 -1", "0 1 0 0 -1", "0 1 0 -1 0", "1 -1 0 0 0", "0 0 1 -1 0", "-1 0 -1 1 -1"],
}
TAB = {m: [[fr(c) for c in row.split()] for row in rows] for m, rows in TAB.items()}

def points_near(a, centre, bound):
    """All x = sum n_i a_i with |x - centre|^2 <= bound, by an exact coefficient box."""
    G = [[dot(u, v) for v in a] for u in a]
    Gi = inverse(G)
    # coordinates of the centre in the basis a: c = G^{-1} (a_i . centre)
    rhs = [dot(u, centre) for u in a]
    c = [sum(Gi[i][j] * rhs[j] for j in range(5)) for i in range(5)]
    # |x - centre|^2 = G[n - c] >= (n_i - c_i)^2 / (G^{-1})_{ii}
    rng = []
    for i in range(5):
        r2 = bound * Gi[i][i]
        r = 0
        while F(r + 1) ** 2 <= r2:
            r += 1
        lo = (c[i] - r - 1).__floor__()
        hi = (c[i] + r + 1).__ceil__()
        rng.append(range(lo, hi + 1))
    out = []
    for n in itertools.product(*rng):
        x = [sum(n[i] * a[i][k] for i in range(5)) for k in range(5)]
        d = [x[k] - centre[k] for k in range(5)]
        if dot(d, d) <= bound:
            out.append(x)
    return out

def main(path):
    E = json.load(open(path))
    zero = [F(0)] * 5
    nmins = []
    for m in range(5):
        ws = TAB[m]
        a, b = ws[:5], ws[5]
        G = [[2 * dot(u, v) for v in ws] for u in ws]
        assert G == [[fr(x) for x in row] for row in E[m]['J']], m
        assert all(IN[m](v) for v in a), m
        assert abs(det(a)) == 4, m
        assert IN[m]([b[k] - BM[m][k] for k in range(5)]), m
        lattice = IN[m]([2 * c for c in BM[m]])
        assert lattice == (m in (0, 3, 4)), m
        if lattice:
            # Lambda_m in D_5, b_m in D_5 \ Lambda_m, covol 4 = 2 covol(D_5): the union is D_5
            assert all(in_D5(v) for v in a) and in_D5(BM[m]) and not IN[m](BM[m])
        # minimal vectors of Lambda_m and neighbours of 0 in P_m
        near0 = [x for x in points_near(a, zero, F(2)) if x != zero]
        assert all(dot(x, x) == 2 for x in near0)
        nmin = len(near0)
        nb = [[c + bc for c, bc in zip(x, BM[m])] for x in points_near(a, [-c for c in BM[m]], F(2))]
        assert all(dot(x, x) == 2 for x in nb)
        kiss = nmin + len(nb)
        # neighbours of b_m: points of Lambda_m near b_m, and b_m + (Min Lambda_m)
        nbb = points_near(a, BM[m], F(2))
        assert all(dot([x[k] - BM[m][k] for k in range(5)], [x[k] - BM[m][k] for k in range(5)]) == 2
                   for x in nbb)
        assert len(nbb) + nmin == kiss
        nmins.append(nmin)
        print(f"m={m}: Gram = E_{m}, basis of Lambda_{m}, b = b_{m} mod Lambda_{m}; "
              f"{'lattice D5' if lattice else 'not a lattice'}; |Min Lambda| = {nmin}; kissing number {kiss}")
        assert kiss == 40
        if m == 2:
            K = {tuple(x) for x in near0} | {tuple(x) for x in nb}
            L5 = set()
            for i, j in itertools.combinations(range(4), 2):
                for s, t in itertools.product((1, -1), repeat=2):
                    v = [F(0)] * 5
                    v[i], v[j] = F(s), F(t)
                    L5.add(tuple(v))
            for i in range(4):
                for s in (1, -1):
                    v = [F(0)] * 5
                    v[i], v[4] = F(s), F(1)
                    L5.add(tuple(v))
            for sg in itertools.product((1, -1), repeat=4):
                if sg.count(-1) % 2 == 1:
                    L5.add(tuple([half * s for s in sg] + [F(-1)]))
            assert K == L5
            print("      the neighbours of 0 in P_2 form L_5 (layers x_5 = 1, 0, -1 with 8, 24, 8 points)")
        if m == 1:
            phi = lambda x: [-x[3], -x[1], -x[2], -x[0], -x[4]]
            cr = [[1, -1, 0, 0, 0], [0, 1, -1, 0, 0], [0, 0, 1, -1, 0], [0, 0, 0, 1, -1],
                  [F(4, 5)] * 5]
            cr = [[F(c) for c in v] for v in cr]
            assert all(in_L1(v) for v in cr) and abs(det(cr)) == 4
            assert all(in_L1(phi(v)) for v in a)
            assert phi(BM[1]) == [0, 0, 0, -1, -1]
            print("      x -> -(x4,x2,x3,x1,x5) maps (Lambda_1, b_1) to the Q5 uniform packing of Cohn-Rajagopal")
    assert nmins == [16, 20, 24, 24, 20]
    # index-2 sublattices of D_5
    D = [[1, -1, 0, 0, 0], [0, 1, -1, 0, 0], [0, 0, 1, -1, 0], [0, 0, 0, 1, -1], [0, 0, 0, 1, 1]]
    D = [[F(c) for c in v] for v in D]
    assert abs(det(D)) == 2
    roots = [x for x in points_near(D, zero, F(2)) if x != zero]
    assert len(roots) == 40
    forms = []
    for i in range(5):
        forms.append(('x_i', lambda x, i=i: x[i] % 2))
    for i, j in itertools.combinations(range(5), 2):
        forms.append(('x_i+x_j', lambda x, i=i, j=j: (x[i] + x[j]) % 2))
    for sg in itertools.product((1, -1), repeat=4):
        e = (1,) + sg
        forms.append(('sum', lambda x, e=e: (sum(c * s for c, s in zip(x, e)) / 2) % 2))
    vals = set()
    count = {}
    for name, f in forms:
        v = tuple(int(f(x)) for x in D)
        assert any(v)
        vals.add(v)
        n = sum(1 for r in roots if f(r) == 0)
        count[(name, n)] = count.get((name, n), 0) + 1
    assert len(vals) == 31 == 2 ** 5 - 1
    assert count == {('x_i', 24): 5, ('x_i+x_j', 16): 10, ('sum', 20): 16}
    print("D5: 31 sublattices of index 2 = 5 (x_i, 24 minimal vectors) + 10 (x_i+x_j, 16) + 16 (sum, 20)")
    print("configurations: all checks passed")

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '../data/equality_points.json')
