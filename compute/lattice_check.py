"""For a configuration with 2b in Lambda, the set Lambda u (Lambda + b) is the lattice L = Lambda + Zb.
Check that L has determinant 128, minimum 4 and 40 minimal vectors forming a root system of type D5
(5 of them with Gram matrix 2 x Cartan(D5) generate L)."""
import sys, json, itertools
from fractions import Fraction as F
from forms import blocks, det, solve, short_vectors, qform

CARTAN_D5 = [[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, -1, -1], [0, 0, -1, 2, 0], [0, 0, -1, 0, 2]]

def main(path):
    data = json.load(open(path))
    for c, C in enumerate(data['classes']):
        J = [[F(x) for x in row] for row in C['J']]
        Q, r, s = blocks(J)
        q = solve(Q, r)           # b = sum q_i a_i
        if not all((2 * x).denominator == 1 for x in q):
            print(f"class {c}: not a lattice (2b not in Lambda)")
            continue
        # basis of L = Z^5 + Z q in coordinates of the a_i: q has entries in (1/2)Z
        # L has index 2 over Z^5; its Gram in the a-basis is Q, on rational coordinates
        # minimal vectors of L: x in Z^5 or x in q + Z^5 with Q[x] = 4
        mins = [tuple(x) for x in short_vectors(Q, [F(0)] * 5, F(4)) if any(x) and qform(Q, x) == 4]
        mins += [tuple(F(xi) + qi for xi, qi in zip(x, q)) for x in short_vectors(Q, q, F(4))
                 if qform(Q, [F(xi) + qi for xi, qi in zip(x, q)]) == 4]
        allv = short_vectors(Q, [F(0)] * 5, F(4 - 1, 1))
        below = [x for x in allv if any(x)]
        belowc = [x for x in short_vectors(Q, q, F(3))]
        mins = list(set(mins) | set(tuple(-v for v in m) for m in mins))
        detL = det(Q) / 4
        def bil(x, y):
            return sum(Q[i][j] * x[i] * y[j] for i in range(5) for j in range(5))
        # greedy search for a simple system
        def search(chosen):
            k = len(chosen)
            if k == 5:
                return chosen
            for v in mins:
                if all(bil(v, chosen[j]) == 2 * CARTAN_D5[k][j] for j in range(k)):
                    res = search(chosen + [v])
                    if res:
                        return res
            return None
        S = search([])
        G = [[bil(a, b) for b in S] for a in S]
        detS = det(G)
        print(f"class {c}: lattice, det L = {detL}, vectors of norm < 4: {len(below) + len(belowc)}, "
              f"minimal vectors: {len(mins)}, simple system of type D5 found: {S is not None}, "
              f"its Gram determinant {detS} (= det L: {detS == detL})")

if __name__ == '__main__':
    main(sys.argv[1])
