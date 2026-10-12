"""Second, independent computation of the extreme rays of the vertex cones with lrs (reverse search),
compared as sets of primitive integer vectors with the cddlib (double description) result."""
import sys, json, subprocess, tempfile, os
from fractions import Fraction as F
from math import gcd, lcm
from forms import min_vectors, sym_coords
from rays import extreme_rays

def primitive(v):
    den = 1
    for x in v:
        den = lcm(den, x.denominator)
    iv = [int(x * den) for x in v]
    g = 0
    for x in iv:
        g = gcd(g, abs(x))
    return tuple(x // g for x in iv)

def lrs_rays(mv):
    rows = [[0] + sym_coords(k) for k in mv]
    with tempfile.TemporaryDirectory() as d:
        ine = os.path.join(d, 'c.ine')
        ext = os.path.join(d, 'c.ext')
        with open(ine, 'w') as fh:
            fh.write("cone\nH-representation\nbegin\n%d %d integer\n" % (len(rows), len(rows[0])))
            for r in rows:
                fh.write(" ".join(map(str, r)) + "\n")
            fh.write("end\n")
        subprocess.run(['lrs', ine, ext], check=True, capture_output=True)
        out = []
        inside = False
        for line in open(ext):
            line = line.strip()
            if line == 'begin':
                inside = True
                continue
            if line == 'end':
                break
            if not inside or not line or line.startswith('*') or 'rational' in line:
                continue
            parts = line.split()
            vals = [F(p) for p in parts]
            if vals[0] == 0:
                out.append(primitive(vals[1:]))
            else:
                assert all(x == 0 for x in vals[1:]), "unexpected vertex"
        return out

def main(path):
    data = json.load(open(path))
    total = 0
    for V in data:
        J = [[F(x) for x in row] for row in V['J']]
        mv = min_vectors(J)
        a = set(extreme_rays(mv))
        b = lrs_rays(mv)
        assert len(b) == len(set(b))
        ok = a == set(b)
        total += len(a)
        print(f"vertex {V['index']}: cdd {len(a)} rays, lrs {len(b)} rays, equal: {ok}", flush=True)
        assert ok and len(a) == V['nrays']
    print("all ray sets agree;", total, "rays in all")

if __name__ == '__main__':
    main(sys.argv[1])
