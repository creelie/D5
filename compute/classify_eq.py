"""Classes up to Gamma of the configurations (Q, b) with Q the Q-block of an equality point and b a deep
hole of Q (covering radius exactly 2), i.e. of the optimal sets Lambda u (Lambda + b)."""
import sys, json
from fractions import Fraction as F
from forms import blocks, det, min_vectors, rank, sym_coords
from deepholes import deep_holes
from equiv import invariant, isometries

def config(Q, y):
    d = len(Q)
    Qy = [sum(Q[i][j] * y[j] for j in range(d)) for i in range(d)]
    s = sum(y[i] * Qy[i] for i in range(d))
    J = [row[:] + [Qy[i]] for i, row in enumerate(Q)]
    J.append(Qy + [s])
    return J

def main(path, out):
    pts = json.load(open(path))
    classes = []   # (J, mv, inv, label)
    summary = []
    for P in pts:
        M = [[F(x) for x in row] for row in P['J']]
        Q, r, s = blocks(M)
        mu2, holes, nv = deep_holes(Q)
        assert mu2 == 4
        found = []
        for y in holes:
            J = config(Q, y)
            mv = min_vectors(J)
            inv = invariant(J, mv)
            hit = None
            for c, (J2, mv2, inv2, lab) in enumerate(classes):
                if inv2 == inv and isometries(J2, J, mv2, mv, first=True):
                    hit = c
                    break
            if hit is None:
                classes.append((J, mv, inv, f"P{len(summary)}"))
                hit = len(classes) - 1
            found.append(hit)
        summary.append({'from': P['from'], 't': P['t'], 'classes': sorted(set(found)),
                        'holes': len(holes)})
        print(f"point from vertex {P['from']} t={P['t']}: {len(holes)} deep holes -> classes {sorted(set(found))}", flush=True)
    print("number of classes:", len(classes))
    rec = []
    for J, mv, inv, lab in classes:
        Q, r, s = blocks(J)
        twob = all((2 * x).denominator == 1 for x in __import__('forms').solve(Q, r))
        print(f"  class {lab}: |Min J|={len(mv)} (l=0: {sum(1 for v in mv if v[-1]==0)}, l=1: {sum(1 for v in mv if v[-1]==1)}) "
              f"lattice packing (2b in Lambda): {twob}")
        rec.append({'J': [[str(x) for x in row] for row in J], 'nmin': len(mv), 'lattice': twob})
    json.dump({'points': summary, 'classes': rec}, open(out, 'w'))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
