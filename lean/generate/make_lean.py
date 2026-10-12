"""Write the data files of the Lean certificate (lean/D5Lean/Data/*.lean) from cert/certificate.txt
and data/equality_points.json.

For every vertex J the data are: the minimal vectors; bound certificates showing that a lattice
point outside a box has J[k] > 4; the box; 21 minimal vectors whose equations determine J, with a
certificate A B = c I; an elimination certificate for det J and det Q.
For every edge J + tY the data are: the end point (T^T J_to T with T in Gamma) or the rank-one
form of Y; the polynomials det(J + tY), det Q(J + tY) with seven sample points and their
elimination certificates; break points u/wd and the pieces between them, on each of which
det J > 0 or det Q > 512 is shown by a coefficient test.

Usage: python3 make_lean.py certificate.txt equality_points.json outdir   (outdir = lean/)
"""
import sys, os, json, math
from fractions import Fraction as F

E55 = [[1 if (i == 5 and j == 5) else 0 for j in range(6)] for i in range(6)]


def lcm(a, b):
    return a * b // math.gcd(a, b)


def lcm_list(xs):
    r = 1
    for x in xs:
        r = lcm(r, x)
    return r


def qform(J, k):
    return sum(k[i] * J[i][j] * k[j] for i in range(len(k)) for j in range(len(k)))


def symc(k):
    return [(1 if a == b else 2) * k[a] * k[b] for a in range(6) for b in range(a, 6)]


def sym_from_vec(v):
    X = [[0] * 6 for _ in range(6)]
    t = 0
    for a in range(6):
        for b in range(a, 6):
            X[a][b] = X[b][a] = v[t]
            t += 1
    return X


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def transpose(A):
    return [list(r) for r in zip(*A)]


def inverse(A):
    n = len(A)
    M = [[F(x) for x in row] + [F(1 if i == j else 0) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        inv = 1 / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def det(M):
    n = len(M)
    A = [[F(x) for x in row] for row in M]
    d = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return d


def rank_rows(rows):
    """indices of a maximal independent subset of rows (greedy, exact)."""
    basis = []   # reduced rows with pivot columns
    chosen = []
    for idx, r in enumerate(rows):
        v = [F(x) for x in r]
        for piv, b in basis:
            if v[piv] != 0:
                f = v[piv] / b[piv]
                v = [x - f * y for x, y in zip(v, b)]
        nz = next((i for i, x in enumerate(v) if x != 0), None)
        if nz is not None:
            basis.append((nz, v))
            chosen.append(idx)
    return chosen


# ---------- elimination certificates ----------

def elim(A):
    """Integer lower-triangular E with E A upper triangular (leading minors of order <= 5 nonzero)."""
    n = len(A)
    M = [[F(x) for x in row] for row in A]
    L = [[F(1 if i == j else 0) for j in range(n)] for i in range(n)]   # running E (rational)
    for c in range(n - 1):
        assert M[c][c] != 0, "zero pivot"
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
                L[r] = [x - f * y for x, y in zip(L[r], L[c])]
    E = []
    for row in L:
        den = lcm_list([x.denominator for x in row])
        num = [int(x * den) for x in row]
        g = 0
        for x in num:
            g = math.gcd(g, x)
        E.append([x // g for x in num])
    EA = matmul(E, A)
    for i in range(n):
        for j in range(n):
            assert (j <= i or E[i][j] == 0) and (i <= j or EA[i][j] == 0)
    dE = math.prod(E[i][i] for i in range(n))
    assert dE != 0
    return E


def check_det(A, E, d, dq):
    EA = matmul(E, A)
    assert math.prod(EA[i][i] for i in range(6)) == d * math.prod(E[i][i] for i in range(6))
    assert math.prod(EA[i][i] for i in range(5)) == dq * math.prod(E[i][i] for i in range(5))


# ---------- bound certificates ----------

def decompose(M):
    """M rational symmetric PSD: M = sum d_r c_r c_r^T with d_r > 0 (exact)."""
    n = len(M)
    M = [[F(x) for x in row] for row in M]
    terms = []
    while True:
        p = next((i for i in range(n) if M[i][i] != 0), None)
        if p is None:
            assert all(x == 0 for row in M for x in row), "not PSD"
            return terms
        assert M[p][p] > 0, "not PSD"
        col = [M[i][p] for i in range(n)]
        d = 1 / M[p][p]
        terms.append((d, col))
        M = [[M[i][j] - d * col[i] * col[j] for j in range(n)] for i in range(n)]


def bound_certs(J, cd, cn):
    Jp = [[F(cd * J[a][c] + (cn if a == 5 and c == 5 else 0)) for c in range(6)] for a in range(6)]
    Ji = inverse(Jp)
    certs, bs = [], []
    for i in range(5):
        lam = 1 / Ji[i][i]
        R = [[Jp[a][c] - (lam if a == i and c == i else 0) for c in range(6)] for a in range(6)]
        terms = [(lam, [1 if a == i else 0 for a in range(6)])]
        for d, col in decompose(R):
            den = lcm_list([x.denominator for x in col])
            m = [int(x * den) for x in col]
            g = 0
            for x in m:
                g = math.gcd(g, x)
            m = [x // g for x in m]
            # d col col^T = d (g/den)^2 m m^T
            terms.append((d * F(g, den) ** 2, m))
        # smallest b with lam (b+1)^2 > 4 cd + cn
        b = 0
        while not lam * (b + 1) ** 2 > 4 * cd + cn:
            b += 1
        certs.append(terms)
        bs.append(b)
    C = lcm_list([e.denominator for terms in certs for e, _ in terms])
    out = [[(int(e * C), m) for e, m in terms] for terms in certs]
    for i, terms in enumerate(out):
        for a in range(6):
            for c in range(6):
                lhs = C * (cd * J[a][c] + (cn if a == 5 and c == 5 else 0))
                assert lhs == sum(e * m[a] * m[c] for e, m in terms)
        assert C * (4 * cd + cn) < terms[0][0] * (bs[i] + 1) ** 2
    return C, bs, out


def box_size(bs):
    return 2 * math.prod(2 * b + 1 for b in bs)


def best_bounds(J):
    Q = [row[:5] for row in J[:5]]
    r = [J[i][5] for i in range(5)]
    q = [sum(x * y for x, y in zip(row, r)) for row in inverse(Q)]
    sigma = J[5][5] - sum(x * y for x, y in zip(r, q))
    cands = []
    base = max(F(0), -sigma)
    for delta in [F(0), F(1, 16), F(1, 8), F(1, 4), F(1, 2), F(1), F(2), F(4), F(8)]:
        rho = base + delta
        if sigma + rho <= 0:
            continue
        cands.append(rho)
    best = None
    for rho in cands:
        cd, cn = rho.denominator, rho.numerator
        C, bs, certs = bound_certs(J, cd, cn)
        sz = box_size(bs)
        if best is None or sz < best[0]:
            best = (sz, cd, cn, C, bs, certs)
    return best


# ---------- polynomials ----------

def poly_from_values(ts, vals):
    """coefficients (lowest first) of the polynomial of degree < len(ts) through the points."""
    n = len(ts)
    A = [[F(t) ** i for i in range(n)] for t in ts]
    Ai = inverse(A)
    return [sum(Ai[i][j] * vals[j] for j in range(n)) for i in range(n)]


def peval(p, t):
    s = F(0)
    for c in reversed(p):
        s = s * t + c
    return s


def pcomp(p, u, h):
    """coefficients in x of p(u + h x)."""
    out = [F(0)]
    for c in reversed(p):
        # out = out * (u + h x) + c
        new = [F(0)] * (len(out) + 1)
        for i, a in enumerate(out):
            new[i] += a * u
            new[i + 1] += a * h
        new[0] += c
        out = new
    return out


def pos_bound(c):
    i = 0
    while i < len(c) and c[i] == 0:
        i += 1
    if i == len(c):
        return False
    return c[i] + sum(min(x, 0) for x in c[i + 1:]) > 0


def pos_inf(c):
    i = 0
    while i < len(c) and c[i] == 0:
        i += 1
    if i == len(c):
        return False
    return c[i] > 0 and all(x >= 0 for x in c[i + 1:])


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def simple_between(a, b):
    """a rational in [a + (b-a)/3, a + 2(b-a)/3] with small denominator."""
    lo, hi = a + (b - a) / 3, a + 2 * (b - a) / 3
    d = 1
    while True:
        n = math.ceil(lo * d)
        if F(n, d) <= hi:
            return F(n, d)
        d *= 2


class SignFail(Exception):
    pass


def sign_cert(pJ, g, req, tend, eqpts):
    """Break points (sorted rationals from 0) and pieces [(modeJ, dirL)].
    tend: the end t* (bounded) or None (unbounded: the last piece is [u, oo))."""
    def point_ok(t):
        a, b = peval(pJ, t), peval(g, t)
        if not (a > 0 or b >= 0):
            return False
        if a == 0 and b == 0 and t not in eqpts:
            return False
        return True

    def piece(a, b):
        for modeJ in (True, False):
            P = pJ if modeJ else g
            if pos_bound(pcomp(P, a, b - a)):
                return (modeJ, True)
            if pos_bound(pcomp(P, b, a - b)):
                return (modeJ, False)
        return None

    def cover(a, b, depth):
        pc = piece(a, b)
        if pc is not None:
            return [b], [pc]
        if depth > 60:
            raise SignFail(f"no certificate on [{a}, {b}]")
        m = simple_between(a, b)
        if not point_ok(m):
            raise SignFail(f"point {m} fails")
        u1, p1 = cover(a, m, depth + 1)
        u2, p2 = cover(m, b, depth + 1)
        return u1 + u2, p1 + p2

    pts = sorted(set(req))
    tail = None
    if tend is None:
        u = max(pts) if pts else F(0)
        u = max(u, F(1))
        while True:
            if pos_inf(pcomp(pJ, u, 1)):
                tail = True
                break
            if pos_inf(pcomp(g, u, 1)):
                tail = False
                break
            u *= 2
        if u not in pts:
            pts.append(u)
        pts.sort()
    for t in pts:
        if not point_ok(t):
            raise SignFail(f"required point {t} fails")
    us, ps = [pts[0]], []
    for a, b in zip(pts, pts[1:]):
        u1, p1 = cover(a, b, 0)
        us += u1
        ps += p1
    return us, ps, tail


# ---------- parsing ----------

def parse_cert(path):
    lines = open(path).read().split("\n")
    pos = 0

    def nxt():
        nonlocal pos
        l = lines[pos]
        pos += 1
        return l

    w = nxt().split()
    assert w[0] == "VERTICES"
    nv = int(w[1])
    verts = []
    for v in range(nv):
        assert nxt().split() == ["VERTEX", str(v)]
        w = nxt().split()
        assert w[0] == "J"
        Jf = [F(x) for x in w[1:]]
        assert all(x.denominator == 1 for x in Jf)
        J = [[int(Jf[6 * i + j]) for j in range(6)] for i in range(6)]
        w = nxt().split()
        mins = [[int(x) for x in nxt().split()] for _ in range(int(w[1]))]
        w = nxt().split()
        gens = [[int(x) for x in nxt().split()] for _ in range(int(w[1]))]
        w = nxt().split()
        rays = [[int(x) for x in nxt().split()] for _ in range(int(w[1]))]
        w = nxt().split()
        reps = []
        for _ in range(int(w[1])):
            p = nxt().split()
            if p[0] == "U":
                reps.append(dict(ray=rays[int(p[1])], idx=int(p[1]), bounded=False))
            else:
                T = [int(x) for x in p[5:]]
                reps.append(dict(ray=rays[int(p[1])], idx=int(p[1]), bounded=True, tstar=F(p[3]),
                                 to=int(p[4]), T=[T[6 * i:6 * i + 6] for i in range(6)]))
        assert nxt() == "END"
        verts.append(dict(J=J, mins=mins, reps=reps))
    w = nxt().split()
    assert w[0] == "EQUALITY"
    eqs = []
    for _ in range(int(w[1])):
        p = nxt().split()
        eqs.append((int(p[1]), int(p[2]), F(p[3])))
    return verts, eqs


# ---------- vertex and edge data ----------

def vertex_data(V, is_eq):
    J = V['J']
    mins = V['mins']
    for k in mins:
        assert qform(J, k) == 4 and k[5] in (0, 1)
    sz, cd, cn, C, bs, certs = best_bounds(J)
    # box check (exact) mirrors ptOK
    idx = rank_rows([symc(k) for k in mins])
    assert len(idx) == 21
    A = [symc(mins[i]) for i in idx]
    Ai = inverse(A)
    c = lcm_list([x.denominator for row in Ai for x in row])
    Bp = [[int(x * c) for x in row] for row in Ai]
    Bt = transpose(Bp)
    AB = matmul(A, Bp)
    assert all(AB[i][j] == (c if i == j else 0) for i in range(21) for j in range(21))
    E = elim(J)
    dJ = det(J)
    dQ = det([row[:5] for row in J[:5]])
    assert dJ.denominator == 1 and dQ.denominator == 1
    dJ, dQ = int(dJ), int(dQ)
    check_det(J, E, dJ, dQ)
    assert dJ > 0 or dQ >= 512
    assert (dJ == 0 and dQ == 512) == is_eq
    return dict(J=J, mins=mins, cd=cd, cn=cn, C=C, bs=bs, certs=certs, idx=idx, A=A, Bt=Bt, c=c,
                E=E, dJ=dJ, dQ=dQ, eq=is_eq, box=sz)


def detpolys(J, Y):
    ts = list(range(7))
    vJ = [det([[J[a][b] + t * Y[a][b] for b in range(6)] for a in range(6)]) for t in ts]
    vQ = [det([[J[a][b] + t * Y[a][b] for b in range(5)] for a in range(5)]) for t in ts]
    pJ = trim(poly_from_values(ts, vJ))
    pQ = trim(poly_from_values(ts[:6], vQ[:6]))
    assert all(x.denominator == 1 for x in pJ + pQ)
    assert all(peval(pQ, t) == vQ[t] for t in ts)
    return [int(x) for x in pJ], [int(x) for x in pQ]


def rank_one(Y):
    Y5 = [row[:5] for row in Y[:5]]
    y = [Y[a][5] for a in range(5)]
    s = Y[5][5]
    if all(x == 0 for row in Y5 for x in row):
        assert all(x == 0 for x in y) and s >= 0
        return 0, 0, s, 0, [0, 0, 0, 0, 0]
    p = next(a for a in range(5) if Y5[a][a] != 0)
    lam0 = Y5[p][p]
    # w proportional to row p
    w = Y5[p][:]
    g = 0
    for x in w:
        g = math.gcd(g, x)
    w = [x // g for x in w]
    if w[p] < 0:
        w = [-x for x in w]
    lam = F(Y5[p][p], w[p] * w[p])
    assert lam.denominator == 1
    lam = int(lam)
    for a in range(5):
        for b in range(5):
            assert Y5[a][b] == lam * w[a] * w[b]
    mu = F(y[p], w[p])
    assert mu.denominator == 1
    mu = int(mu)
    assert all(y[a] == mu * w[a] for a in range(5))
    n0 = math.floor(F(-mu, lam))
    assert lam * n0 <= -mu < lam * (n0 + 1)
    assert lam * n0 * n0 + 2 * mu * n0 + s >= 0
    assert lam * (n0 + 1) ** 2 + 2 * mu * (n0 + 1) + s >= 0
    return lam, mu, s, n0, w


def edge_data(src, dst, rep, eqpts):
    J = src['J']
    Y = sym_from_vec(rep['ray'])
    pJ, pQ = detpolys(J, Y)
    g = list(pQ) if pQ else [0]
    g = [g[0] - 512] + g[1:]
    d = dict(Y=Y, bounded=rep['bounded'], pJ=pJ, pQ=pQ)
    if rep['bounded']:
        ts = rep['tstar']
        T = rep['T']
        Ti = inverse(T)
        assert all(x.denominator == 1 for row in Ti for x in row)
        Ti = [[int(x) for x in row] for row in Ti]
        Jto = dst['J']
        end = matmul(matmul(transpose(T), Jto), T)
        assert all(J[a][b] + ts * Y[a][b] == end[a][b] for a in range(6) for b in range(6))
        assert all(T[5][j] == 0 for j in range(5)) and T[5][5] in (1, -1)
        kmax = None
        # minimal vectors of the end point are Tinv k for minimal k of J_to
        for k in dst['mins']:
            kk = [sum(Ti[i][j] * k[j] for j in range(6)) for i in range(6)]
            if kk[5] == -1:
                kk = [-x for x in kk]
            if qform(Y, kk) < 0:
                kmax = kk
                break
        assert kmax is not None and kmax[5] in (0, 1)
        assert qform(J, kmax) + ts * qform(Y, kmax) == 4
        d.update(tn=ts.numerator, td=ts.denominator, target=rep['to'], T=T, Tinv=Ti, kmax=kmax)
        tend = ts
        sts = [F(j, 6) * ts for j in range(7)]
    else:
        lam, mu, s, n0, w = rank_one(Y)
        d.update(lam=lam, mu=mu, s=s, n0=n0, w=w)
        tend = None
        sts = [F(j) for j in range(7)]
    sd = lcm_list([x.denominator for x in sts])
    sn = [int(x * sd) for x in sts]
    Es = []
    for x in sn:
        A = [[sd * J[a][b] + x * Y[a][b] for b in range(6)] for a in range(6)]
        E = elim(A)
        hJ = sum(c * x ** i * sd ** (6 - i) for i, c in enumerate(pJ))
        hQ = sum(c * x ** i * sd ** (5 - i) for i, c in enumerate(pQ))
        check_det(A, E, hJ, hQ)
        Es.append(E)
    eqpts = list(eqpts) + ([F(0)] if src['eq'] else []) + \
        ([tend] if rep['bounded'] and dst['eq'] else [])
    req = [F(0)] + ([tend] if tend is not None else []) + list(eqpts)
    us, ps, tail = sign_cert([F(x) for x in pJ], [F(x) for x in g], req, tend,
                             set(F(x) for x in eqpts))
    wd = lcm_list([x.denominator for x in us])
    d.update(sd=sd, sn=sn, Es=Es, wd=wd, us=[int(x * wd) for x in us], ps=ps,
             tailJ=bool(tail) if tail is not None else False,
             eqs=sorted(set(int(x * wd) for x in eqpts)))
    return d


# ---------- Lean output ----------

def L(x):
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, int):
        return str(x) if x >= 0 else f"({x})"
    if isinstance(x, tuple):
        return "(" + ", ".join(L(y) for y in x) + ")"
    if isinstance(x, list):
        return "[" + ", ".join(L(y) for y in x) + "]"
    raise TypeError(x)


def lean_vertex(name, v):
    certs = [[(e, m) for e, m in terms] for terms in v['certs']]
    return f"""def {name} : Vtx where
  J := {L(v['J'])}
  mins := {L(v['mins'])}
  cd := {L(v['cd'])}
  cn := {L(v['cn'])}
  C := {L(v['C'])}
  bs := {L(v['bs'])}
  certs := {L(certs)}
  E := {L(v['E'])}
  dJ := {L(v['dJ'])}
  dQ := {L(v['dQ'])}
  eq := {L(v['eq'])}
"""


def lean_edge(name, e):
    z6 = [[0] * 6 for _ in range(6)]
    f = dict(Y=e['Y'], bounded=e['bounded'], tn=e.get('tn', 0), td=e.get('td', 1),
             target=e.get('target', 0), T=e.get('T', z6), Tinv=e.get('Tinv', z6),
             kmax=e.get('kmax', []), lam=e.get('lam', 0), mu=e.get('mu', 0), s=e.get('s', 0),
             n0=e.get('n0', 0), w=e.get('w', []), pJ=e['pJ'], pQ=e['pQ'], sd=e['sd'], sn=e['sn'],
             Es=e['Es'], wd=e['wd'], us=e['us'], ps=[tuple(p) for p in e['ps']], tailJ=e['tailJ'],
             eqs=e['eqs'])
    body = "\n".join(f"  {k} := {L(v)}" for k, v in f.items())
    return f"def {name} : Edge where\n{body}\n"


HEADER = """/- Generated by lean/generate/make_lean.py from cert/certificate.txt; do not edit. -/
"""


def write_lean(vd, edges, outdir):
    ddir = os.path.join(outdir, "D5Lean", "Data")
    cdir = os.path.join(outdir, "D5Lean", "Cert")
    os.makedirs(ddir, exist_ok=True)
    os.makedirs(cdir, exist_ok=True)
    n = len(vd)
    with open(os.path.join(ddir, "Vertices.lean"), "w") as fh:
        fh.write(HEADER + "import D5Lean.Check\n\nnamespace D5\n\n")
        for i, v in enumerate(vd):
            fh.write(lean_vertex(f"v{i}", v) + "\n")
        fh.write("/-- The 34 vertices. -/\ndef vertices : List Vtx :=\n  [" +
                 ", ".join(f"v{i}" for i in range(n)) + "]\n\n")
        fh.write("/-- Vertex `i` of the list. -/\ndef vtx (i : Nat) : Vtx := vertices.getD i default\n\n")
        fh.write("/-- The check of an edge at vertex `src`: the end point of a bounded edge is one of the\n"
                 f"vertices, number `target`. -/\ndef edgeOK (src : Nat) (e : Edge) : Bool :=\n"
                 f"  (!e.bounded || decide (e.target < {n})) &&\n"
                 "    Edge.ok (vtx src).J (vtx (if e.bounded then e.target else src)).J e\n\n")
        fh.write("end D5\n")
    for i in range(n):
        es = [e for s, e in edges if s == i]
        with open(os.path.join(cdir, f"V{i}.lean"), "w") as fh:
            fh.write(HEADER + "import D5Lean.Data.Vertices\n\nnamespace D5\n\n")
            fh.write(f"theorem v{i}_ok : v{i}.ok = true := by decide +kernel\n\n")
            for j, e in enumerate(es):
                fh.write(lean_edge(f"e{i}_{j}", e) + "\n")
                fh.write(f"theorem e{i}_{j}_ok : edgeOK {i} e{i}_{j} = true := by decide +kernel\n\n")
            fh.write(f"/-- The edges at vertex {i}, one for each orbit of its stabiliser. -/\n")
            fh.write(f"def edges{i} : List Edge :=\n  [" + ", ".join(f"e{i}_{j}" for j in range(len(es))) + "]\n\n")
            fh.write(f"theorem edges{i}_ok : edges{i}.all (edgeOK {i}) = true := by\n")
            fh.write(f"  simp only [edges{i}, List.all_cons, List.all_nil, " +
                     "".join(f"e{i}_{j}_ok, " for j in range(len(es))) + "Bool.and_self, Bool.true_and]\n\n")
            fh.write("end D5\n")
    with open(os.path.join(cdir, "All.lean"), "w") as fh:
        fh.write(HEADER + "".join(f"import D5Lean.Cert.V{i}\n" for i in range(n)) + "\nnamespace D5\n\n")
        fh.write(f"/-- The edges at each vertex. -/\ndef edgeLists : List (List Edge) :=\n  [" +
                 ", ".join(f"edges{i}" for i in range(n)) + "]\n\n")
        fh.write("/-- The edges at vertex `i`. -/\ndef edgesOf (i : Nat) : List Edge := edgeLists.getD i []\n\n")
        fh.write(f"theorem vertices_length : vertices.length = {n} := rfl\n\n")
        fh.write("theorem vertices_ok : vertices.all Vtx.ok = true := by\n")
        fh.write("  simp only [vertices, List.all_cons, List.all_nil, " +
                 "".join(f"v{i}_ok, " for i in range(n)) + "Bool.and_self, Bool.true_and]\n\n")
        fh.write(f"theorem edges_ok : ∀ i < {n}, (edgesOf i).all (edgeOK i) = true\n")
        for i in range(n):
            fh.write(f"  | {i}, _ => edges{i}_ok\n")
        fh.write(f"  | _ + {n}, h => absurd h (by omega)\n")
        fh.write("\nend D5\n")


def main(cpath, epath, outdir):
    verts, eqs = parse_cert(cpath)
    os.makedirs(outdir, exist_ok=True)
    vd = []
    for i, V in enumerate(verts):
        is_eq = any(f == i and r == -1 for f, r, _ in eqs)
        vd.append(vertex_data(V, is_eq))
        print(f"vertex {i}: box {vd[-1]['box']} points", file=sys.stderr)
    total_box = sum(v['box'] for v in vd)
    print(f"box points in all: {total_box}", file=sys.stderr)
    edges = []
    for i, V in enumerate(verts):
        for rep in V['reps']:
            eqp = [t for f, r, t in eqs if f == i and r == rep['idx']]
            dst = vd[rep['to']] if rep['bounded'] else None
            e = edge_data(vd[i], dst, rep, eqp)
            edges.append((i, e))
    print(f"edges: {len(edges)}, pieces: {sum(len(e['ps']) for _, e in edges)}", file=sys.stderr)
    write_lean(vd, edges, outdir)
    return vd, edges


if __name__ == '__main__':
    main(*sys.argv[1:4])
