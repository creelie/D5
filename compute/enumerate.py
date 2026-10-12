"""Enumerate the vertices of R (d = 5) up to Gamma, starting from the D5 vertex J5,1.

For every vertex: its minimal vectors, the extreme rays of its cone, and for every ray either
'unbounded' or the contiguous vertex J + t*Y together with its class.  Results go to a JSON file."""
import sys, json, time
from fractions import Fraction as F
from forms import frac_mat, min_vectors, sym_from_vec, blocks, det, rank, sym_coords
from rays import extreme_rays
from edges import tstar
from unbounded import classify
from equiv import invariant, isometries

N = 6

def fstr(M):
    return [[str(x) for x in row] for row in M]

def main(out):
    J51 = frac_mat([[2 * x for x in r] for r in
                    [[2, -1, -1, 0, 0, 1], [-1, 2, 1, 0, 0, 0], [-1, 1, 2, 0, 0, -1],
                     [0, 0, 0, 2, 0, -1], [0, 0, 0, 0, 2, -1], [1, 0, -1, -1, -1, 2]]])
    verts = []      # dicts: J, mv, inv
    bykey = {}      # invariant -> list of vertex indices

    def add_vertex(J, mv, inv):
        verts.append({'J': J, 'mv': mv, 'inv': inv})
        bykey.setdefault(inv, []).append(len(verts) - 1)
        return len(verts) - 1

    def lookup(J, mv, inv):
        for i in bykey.get(inv, []):
            T = isometries(verts[i]['J'], J, verts[i]['mv'], mv, first=True)
            if T:
                return i, T[0]
        return None, None

    mv0 = min_vectors(J51)
    add_vertex(J51, mv0, invariant(J51, mv0))
    report = []
    i = 0
    t_start = time.time()
    while i < len(verts):
        V = verts[i]
        J, mv = V['J'], V['mv']
        assert rank([sym_coords(k) for k in mv]) == N * (N + 1) // 2, "not a vertex"
        t0 = time.time()
        R = extreme_rays(mv)
        edges = []
        cache = {}
        for y in R:
            Y = sym_from_vec(y, N)
            c, w = classify(Y)
            if c == 'unbounded':
                edges.append({'ray': list(y), 'tstar': None})
                continue
            ts = tstar(J, Y, w)
            M = [[J[a][b] + ts * Y[a][b] for b in range(N)] for a in range(N)]
            mvM = min_vectors(M)
            inv = invariant(M, mvM)
            j, T = lookup(M, mvM, inv)
            if j is None:
                j = add_vertex(M, mvM, inv)
                T = [[1 if a == b else 0 for b in range(N)] for a in range(N)]
            edges.append({'ray': list(y), 'tstar': str(ts), 'to': j, 'T': T})
        Q, r, s = blocks(J)
        rec = {'index': i, 'J': fstr(J), 'detQ': str(det(Q)), 'detJ': str(det(J)),
               'nmin': len(mv), 'nrays': len(R), 'edges': edges}
        report.append(rec)
        nb = sorted({e['to'] for e in edges if e['tstar'] is not None})
        print(f"vertex {i}: detQ={det(Q)} detJ={det(J)} |Min|={len(mv)} rays={len(R)} "
              f"unbounded={sum(1 for e in edges if e['tstar'] is None)} neighbours={nb} "
              f"[{time.time() - t0:.1f}s, total {len(verts)} vertices]", flush=True)
        with open(out, 'w') as fh:
            json.dump(report, fh)
        i += 1
    print(f"done: {len(verts)} vertices in {time.time() - t_start:.0f}s")

if __name__ == '__main__':
    main(sys.argv[1])
