"""Independent checker (Python, exact rationals) for the certificate of the two-translate problem in R^5.
Usage: python3 check.py certificate.txt"""
import sys
from fractions import Fraction as Fr
from math import gcd, floor

def fail(m):
    print("FAIL:", m, file=sys.stderr)
    sys.exit(1)

def qform(J, k):
    return sum(J[i][j] * k[i] * k[j] for i in range(6) for j in range(6) if k[i] and k[j])

def symc(k):
    return [k[i] * k[j] * (1 if i == j else 2) for i in range(6) for j in range(i, 6)]

def det(A0):
    A = [row[:] for row in A0]
    n = len(A)
    d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[p], A[c] = A[c], A[p]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / A[c][c]
                for j in range(c, n):
                    A[r][j] -= f * A[c][j]
    return d

def ldl(A):
    n = len(A)
    d = [Fr(0)] * n
    m = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        v = A[i][i] - sum(d[k] * m[k][i] ** 2 for k in range(i))
        if v <= 0:
            return None
        d[i] = v
        for j in range(i + 1, n):
            m[i][j] = (A[i][j] - sum(d[k] * m[k][i] * m[k][j] for k in range(i))) / v
    return d, m

def enumerate_short(d, m, c, B, cb):
    """call cb(x) for every integer x with Q[x + c] <= B."""
    n = len(d)
    x = [0] * n

    def rec(i, R):
        if i < 0:
            cb(x)
            return
        u = c[i] + sum(m[i][j] * (x[j] + c[j]) for j in range(i + 1, n))
        x0 = floor(-u)
        v = x0 + 1
        while True:
            val = d[i] * (v + u) ** 2
            if val > R:
                break
            x[i] = v
            rec(i - 1, R - val)
            v += 1
        v = x0
        while True:
            val = d[i] * (v + u) ** 2
            if val > R:
                break
            x[i] = v
            rec(i - 1, R - val)
            v -= 1
        x[i] = 0

    if B >= 0:
        rec(n - 1, B)

def solve(A, b):
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                for j in range(c, n + 1):
                    M[r][j] -= f * M[c][j]
    return [M[i][n] / M[i][i] for i in range(n)]

def in_R(J):
    A = [row[:5] for row in J[:5]]
    f = ldl(A)
    if f is None:
        fail("Q block not positive definite")
    d, m = f
    state = {'ok': True}
    mins = []

    def l0(x):
        if not any(x):
            return
        first = next(v for v in x if v != 0)
        if first < 0:
            return
        k = tuple(x) + (0,)
        v = qform(J, k)
        if v < 4:
            state['ok'] = False
        elif v == 4:
            mins.append(k)

    enumerate_short(d, m, [Fr(0)] * 5, Fr(4), l0)
    r = [J[i][5] for i in range(5)]
    q = solve(A, r)
    sigma = J[5][5] - sum(a * b for a, b in zip(r, q))

    def l1(x):
        k = tuple(x) + (1,)
        v = qform(J, k)
        if v < 4:
            state['ok'] = False
        elif v == 4:
            mins.append(k)

    enumerate_short(d, m, q, Fr(4) - sigma, l1)
    return state['ok'], mins

P61 = (1 << 61) - 1

def rank_mod(rows):
    A = [[x % P61 for x in r] for r in rows]
    if not A:
        return 0
    n = len(A[0])
    rk = 0
    for c in range(n):
        if rk == len(A):
            break
        p = next((r for r in range(rk, len(A)) if A[r][c]), None)
        if p is None:
            continue
        A[p], A[rk] = A[rk], A[p]
        inv = pow(A[rk][c], P61 - 2, P61)
        for r in range(len(A)):
            if r != rk and A[r][c]:
                f = A[r][c] * inv % P61
                A[r] = [(a - f * b) % P61 for a, b in zip(A[r], A[rk])]
        rk += 1
    return rk

# polynomials: lists of Fractions, low degree first
def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p

def evalp(p, t):
    s = Fr(0)
    for a in reversed(p):
        s = s * t + a
    return s

def deriv(p):
    return trim([p[i] * i for i in range(1, len(p))])

def mulp(a, b):
    if not a or not b:
        return []
    c = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)

def divmodp(a, b):
    r = trim(a)
    q = [Fr(0)] * max(len(r) - len(b) + 1, 0)
    while r and len(r) >= len(b):
        f = r[-1] / b[-1]
        sh = len(r) - len(b)
        q[sh] = f
        for i, y in enumerate(b):
            r[sh + i] -= f * y
        r = trim(r)
    return trim(q), r

def gcdp(a, b):
    a, b = trim(a), trim(b)
    while b:
        _, r = divmodp(a, b)
        a, b = b, r
    return [x / a[-1] for x in a] if a else a

def squarefree(p):
    q, r = divmodp(p, gcdp(p, deriv(p)))
    if r:
        fail("sqf")
    return q

def sign(x):
    return (x > 0) - (x < 0)

def sturm(A, B):
    s = [A]
    b = trim(B)
    if not b:
        return s
    s.append(b)
    while True:
        _, r = divmodp(s[-2], s[-1])
        if not r:
            break
        s.append([-x for x in r])
    return s

def changes(s, t):
    prev, ch = 0, 0
    for p in s:
        v = sign(p[-1]) if t is None else sign(evalp(p, t))
        if v == 0:
            continue
        if prev and v != prev:
            ch += 1
        prev = v
    return ch

def taq(P, G, a, b):
    """sum of sign(G(x)) over the roots x of P in (a, b); P squarefree, P(a), P(b) != 0; b None = infinity."""
    if len(P) <= 1:
        return 0
    s = sturm(P, mulp(deriv(P), G))
    return changes(s, a) - changes(s, b)

def interp(v):
    n = len(v)
    c = list(v)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) / j
    p = [c[n - 1]]
    for i in range(n - 2, -1, -1):
        p = mulp(p, [Fr(-i), Fr(1)]) or [Fr(0)]
        p[0] += c[i]
        p = trim(p)
    return trim(p)

def drop_root(P, t):
    if evalp(P, t) != 0:
        return P
    q, _ = divmodp(P, [-t, Fr(1)])
    return q

def act(T, Y):
    """T^T Y T"""
    YT = [[sum(Y[i][k] * T[k][j] for k in range(6) if T[k][j]) for j in range(6)] for i in range(6)]
    return [[sum(T[k][i] * YT[k][j] for k in range(6) if T[k][i]) for j in range(6)] for i in range(6)]

def symm(y):
    Y = [[Fr(0)] * 6 for _ in range(6)]
    t = 0
    for i in range(6):
        for j in range(i, 6):
            Y[i][j] = Y[j][i] = Fr(y[t])
            t += 1
    return Y

def check_gamma(T):
    if any(T[5][j] != 0 for j in range(5)) or abs(T[5][5]) != 1:
        fail("T not in Gamma")
    if abs(det([[Fr(T[i][j]) for j in range(5)] for i in range(5)])) != 1:
        fail("T not in Gamma (det U)")

def main(path):
    toks = open(path).read().split()
    pos = [0]

    def nxt():
        t = toks[pos[0]]
        pos[0] += 1
        return t

    def nint():
        return int(nxt())

    if nxt() != "VERTICES":
        fail("format")
    nv = nint()
    V = []
    for v in range(nv):
        if nxt() != "VERTEX" or nint() != v:
            fail("format vertex")
        nxt()
        flat = [Fr(nxt()) for _ in range(36)]
        J = [flat[6 * i:6 * i + 6] for i in range(6)]
        nxt(); m = nint()
        mins = [tuple(nint() for _ in range(6)) for _ in range(m)]
        nxt(); g = nint()
        gens = []
        for _ in range(g):
            f = [nint() for _ in range(36)]
            gens.append([f[6 * i:6 * i + 6] for i in range(6)])
        nxt(); nr = nint()
        rays = [tuple(nint() for _ in range(21)) for _ in range(nr)]
        nxt(); nreps = nint()
        reps = []
        for _ in range(nreps):
            kind = nxt(); idx = nint(); size = nint()
            if kind == 'B':
                ts = Fr(nxt()); to = nint()
                f = [nint() for _ in range(36)]
                reps.append((kind, idx, size, ts, to, [f[6 * i:6 * i + 6] for i in range(6)]))
            else:
                reps.append((kind, idx, size, None, -1, None))
        if nxt() != "END":
            fail("format end")
        V.append((J, mins, gens, rays, reps))
    nxt(); ne = nint()
    eqpts = []
    for _ in range(ne):
        nxt()
        eqpts.append((nxt(), nxt(), nxt()))

    total_rays = total_edges = eq_found = 0
    eq_mine = []
    for v in range(nv):
        J, lmins, gens, rays, reps = V[v]
        if any(J[i][j] != J[j][i] for i in range(6) for j in range(6)):
            fail("J not symmetric")
        ok, mins = in_R(J)
        if not ok:
            fail("vertex not in R")
        if sorted(lmins) != sorted(mins):
            fail(f"Min(J) differs at vertex {v}")
        if rank_mod([symc(k) for k in mins]) != 21:
            fail("not a vertex")
        for T in gens:
            check_gamma(T)
            if act(T, J) != J:
                fail("generator is not an automorphism")
        posd = {}
        forms = [symc(k) for k in mins]
        for i, y in enumerate(rays):
            if y in posd:
                fail("repeated ray")
            posd[y] = i
            g0 = 0
            for a in y:
                g0 = gcd(g0, abs(a))
            if g0 != 1:
                fail("ray not primitive")
            tight = []
            for a in forms:
                s = sum(p * q for p, q in zip(a, y))
                if s < 0:
                    fail("ray outside the cone")
                if s == 0:
                    tight.append(a)
            if rank_mod(tight) != 20:
                fail("ray not extreme")
        par = list(range(len(rays)))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        for T in gens:
            # y -> vec(T^T sym(y) T) is an integer linear map on the 21 coordinates
            A = []
            for c in range(21):
                e = [0] * 21
                e[c] = 1
                Y = act(T, symm(e))
                A.append([int(Y[a][b]) for a in range(6) for b in range(a, 6)])
            rows = [[A[c][r] for c in range(21)] for r in range(21)]
            for i, y in enumerate(rays):
                yy = [sum(p * q for p, q in zip(row, y) if p) for row in rows]
                j = posd.get(tuple(yy))
                if j is None:
                    fail("image of a ray is not a ray")
                a, b = find(i), find(j)
                if a != b:
                    par[a] = b
        osize = {}
        for i in range(len(rays)):
            r = find(i)
            osize[r] = osize.get(r, 0) + 1
        seen = set()
        for kind, idx, size, ts, to, T in reps:
            o = find(idx)
            if o in seen:
                fail("two representatives in one orbit")
            seen.add(o)
            if osize[o] != size:
                fail("orbit size")
        if len(seen) != len(osize):
            fail("an orbit has no representative")
        total_rays += len(rays)
        dJ = det(J)
        dQ = det([row[:5] for row in J[:5]])
        if dJ <= 0 and dQ < 512:
            fail("vertex with det J <= 0 and det Q < 512")
        if dJ == 0 and dQ == 512:
            eq_mine.append((str(v), "-1", "0"))
        for kind, idx, size, ts, to, T in reps:
            Y = symm(rays[idx])
            Tp = None
            if kind == 'U':
                i0 = next((i for i in range(5) if Y[i][i] != 0), None)
                s = Y[5][5]
                if i0 is None:
                    if any(Y[i][j] != 0 for i in range(5) for j in range(6)):
                        fail("unbounded: Q_Y = 0 but Y != s e e^T")
                    if s < 0:
                        fail("unbounded: s < 0")
                else:
                    if Y[i0][i0] < 0:
                        fail("unbounded: negative diagonal")
                    w = [Y[i][i0] / Y[i0][i0] for i in range(5)]
                    L = 1
                    for x in w:
                        L = L * x.denominator // gcd(L, x.denominator)
                    w = [x * L for x in w]
                    G = 0
                    for x in w:
                        G = gcd(G, abs(x.numerator))
                    w = [x / G for x in w]
                    lam = Y[i0][i0] / w[i0] ** 2
                    if any(Y[i][j] != lam * w[i] * w[j] for i in range(5) for j in range(5)):
                        fail("unbounded: Q_Y not of rank one")
                    mu = Y[i0][5] / w[i0]
                    if any(Y[i][5] != mu * w[i] for i in range(5)):
                        fail("unbounded: r_Y not parallel to w")
                    n0 = floor(-mu / lam)
                    for n in (n0, n0 + 1):
                        if lam * n * n + 2 * mu * n + s < 0:
                            fail("unbounded: negative value")
            else:
                if ts <= 0:
                    fail("t* <= 0")
                M = [[J[i][j] + ts * Y[i][j] for j in range(6)] for i in range(6)]
                ok2, mm = in_R(M)
                if not ok2:
                    fail("end point not in R")
                if not any(qform(Y, k) < 0 for k in mm):
                    fail("t* not maximal")
                if not 0 <= to < nv:
                    fail("target out of range")
                check_gamma(T)
                if act(T, V[to][0]) != M:
                    fail("end point not equivalent to the target vertex")
                Tp = ts
            vJ, vQ = [], []
            for t in range(7):
                A = [[J[i][j] + t * Y[i][j] for j in range(6)] for i in range(6)]
                vJ.append(det(A))
                if t <= 5:
                    vQ.append(det([row[:5] for row in A[:5]]))
            pJ, pQ = interp(vJ), interp(vQ)
            g = list(pQ) if pQ else [Fr(0)]
            g[0] -= 512
            g = trim(g)
            if not pJ:
                if not g:
                    fail("det Q constant 512 on an edge")
                mid = Fr(1) if Tp is None else Tp / 2
                if evalp(g, mid) <= 0:
                    fail("zero edge: det Q < 512 inside")
                h = drop_root(squarefree(g), Fr(0))
                if Tp is not None:
                    h = drop_root(h, Tp)
                if taq(h, [Fr(1)], Fr(0), Tp) != 0:
                    fail("zero edge: det Q - 512 has a root inside")
                total_edges += 1
                continue
            P = drop_root(squarefree(pJ), Fr(0))
            if Tp is not None:
                P = drop_root(P, Tp)
            N = taq(P, [Fr(1)], Fr(0), Tp)
            S = taq(P, g, Fr(0), Tp)
            h = gcdp(P, g) if g else P
            N0 = taq(h, [Fr(1)], Fr(0), Tp)
            if (N - N0 - S) % 2:
                fail("parity")
            if (N - N0 - S) // 2 != 0:
                fail("a point with det J = 0 and det Q < 512")
            if N0 > 0:
                for e in eqpts:
                    if e[0] == str(v) and e[1] == str(idx):
                        t = Fr(e[2])
                        if evalp(pJ, t) != 0 or evalp(g, t) != 0:
                            fail("listed equality point is not one")
                        eq_mine.append(e)
                eq_found += N0
            total_edges += 1
        print(f"vertex {v}: ok", flush=True)
    nvx = sum(1 for e in eq_mine if e[1] == "-1")
    if len(eq_mine) != eq_found + nvx or len(eqpts) != len(eq_mine):
        fail("equality points do not match the certificate")
    print(f"vertices: {nv}, rays: {total_rays}, edge representatives: {total_edges}, equality points: {len(eq_mine)}")
    print("Python checker: all checks passed")

if __name__ == '__main__':
    main(sys.argv[1])
