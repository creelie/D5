// Independent checker (C++, GMP rationals) for the certificate of the two-translate problem in R^5.
// Usage: check certificate.txt
#include <gmpxx.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <functional>
using namespace std;
typedef mpq_class Q;
typedef vector<Q> Poly;   // coefficients, low degree first

static void fail(const string &m) { cerr << "FAIL: " << m << endl; exit(1); }

// ---------- forms ----------
typedef array<array<Q, 6>, 6> M6;
typedef array<long, 6> K6;

static Q qform(const M6 &J, const K6 &k) {
    Q s = 0;
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) s += J[i][j] * k[i] * k[j];
    return s;
}
static vector<long> symc(const K6 &k) {   // linear form Y -> Y[k] in the coordinates y_ij, i <= j
    vector<long> out;
    for (int i = 0; i < 6; i++) for (int j = i; j < 6; j++) out.push_back(k[i] * k[j] * (i == j ? 1 : 2));
    return out;
}
static Q det(vector<vector<Q>> A) {
    int n = A.size(); Q d = 1;
    for (int c = 0; c < n; c++) {
        int p = -1;
        for (int r = c; r < n; r++) if (A[r][c] != 0) { p = r; break; }
        if (p < 0) return 0;
        if (p != c) { swap(A[p], A[c]); d = -d; }
        d *= A[c][c];
        for (int r = c + 1; r < n; r++) if (A[r][c] != 0) {
            Q f = A[r][c] / A[c][c];
            for (int j = c; j < n; j++) A[r][j] -= f * A[c][j];
        }
    }
    return d;
}
static vector<vector<Q>> qblock(const M6 &J) {
    vector<vector<Q>> A(5, vector<Q>(5));
    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) A[i][j] = J[i][j];
    return A;
}
static vector<vector<Q>> full(const M6 &J) {
    vector<vector<Q>> A(6, vector<Q>(6));
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) A[i][j] = J[i][j];
    return A;
}

// Q[x] = sum_i d_i (x_i + sum_{j>i} m_ij x_j)^2 ; returns false if Q is not positive definite
static bool ldl(const vector<vector<Q>> &A, vector<Q> &d, vector<vector<Q>> &m) {
    int n = A.size();
    d.assign(n, 0); m.assign(n, vector<Q>(n, 0));
    for (int i = 0; i < n; i++) {
        Q v = A[i][i];
        for (int k = 0; k < i; k++) v -= d[k] * m[k][i] * m[k][i];
        if (v <= 0) return false;
        d[i] = v;
        for (int j = i + 1; j < n; j++) {
            Q w = A[i][j];
            for (int k = 0; k < i; k++) w -= d[k] * m[k][i] * m[k][j];
            m[i][j] = w / v;
        }
    }
    return true;
}
static mpz_class floorq(const Q &x) { mpz_class r; mpz_fdiv_q(r.get_mpz_t(), x.get_num_mpz_t(), x.get_den_mpz_t()); return r; }

// all integer x with Q[x + c] <= B (Q positive definite, given by its LDL data)
static void enumerate(const vector<Q> &d, const vector<vector<Q>> &m, const vector<Q> &c, const Q &B,
                      const function<void(const vector<long> &)> &cb) {
    int n = d.size();
    vector<long> x(n, 0);
    function<void(int, Q)> rec = [&](int i, Q R) {
        if (i < 0) { cb(x); return; }
        Q u = c[i];
        for (int j = i + 1; j < n; j++) u += m[i][j] * (x[j] + c[j]);
        // d_i (x_i + u)^2 <= R
        // (v + u)^2 increases for v >= x0 + 1 > -u and for v <= x0 <= -u as v moves away
        long x0 = floorq(-u).get_si();
        for (long v = x0 + 1;; v++) { Q t = v + u; Q val = d[i] * t * t; if (val > R) break; x[i] = v; rec(i - 1, R - val); }
        for (long v = x0;; v--) { Q t = v + u; Q val = d[i] * t * t; if (val > R) break; x[i] = v; rec(i - 1, R - val); }
    };
    if (B >= 0) rec(n - 1, B);
}

// minimum of J over M* (and the k attaining 4) ; requires the Q block positive definite.
// Returns false if some value is below 4.
static bool in_R(const M6 &J, vector<K6> &mins) {
    vector<vector<Q>> A = qblock(J);
    vector<Q> d; vector<vector<Q>> m;
    if (!ldl(A, d, m)) fail("Q block not positive definite");
    mins.clear();
    bool ok = true;
    vector<Q> zero(5, 0);
    enumerate(d, m, zero, Q(4), [&](const vector<long> &x) {
        bool nz = false; for (long v : x) if (v) nz = true;
        if (!nz) return;
        int f = 0; while (x[f] == 0) f++;
        if (x[f] < 0) return;                       // one of +-x
        K6 k; for (int i = 0; i < 5; i++) k[i] = x[i]; k[5] = 0;
        Q v = qform(J, k);
        if (v < 4) ok = false; else if (v == 4) mins.push_back(k);
    });
    // l = 1: J[(x,1)] = Q[x + q] + sigma
    vector<Q> r(5); for (int i = 0; i < 5; i++) r[i] = J[i][5];
    // q = Q^{-1} r by Gaussian elimination
    vector<vector<Q>> Aug = A;
    for (int i = 0; i < 5; i++) Aug[i].push_back(r[i]);
    for (int c = 0; c < 5; c++) {
        int p = c; while (Aug[p][c] == 0) p++;
        swap(Aug[p], Aug[c]);
        for (int rr = 0; rr < 5; rr++) if (rr != c && Aug[rr][c] != 0) {
            Q f = Aug[rr][c] / Aug[c][c];
            for (int j = c; j < 6; j++) Aug[rr][j] -= f * Aug[c][j];
        }
    }
    vector<Q> q(5); for (int i = 0; i < 5; i++) q[i] = Aug[i][5] / Aug[i][i];
    Q sigma = J[5][5]; for (int i = 0; i < 5; i++) sigma -= r[i] * q[i];
    enumerate(d, m, q, Q(4) - sigma, [&](const vector<long> &x) {
        K6 k; for (int i = 0; i < 5; i++) k[i] = x[i]; k[5] = 1;
        Q v = qform(J, k);
        if (v < 4) ok = false; else if (v == 4) mins.push_back(k);
    });
    return ok;
}

// ---------- rank modulo a prime ----------
typedef unsigned long long u64;
static const u64 P61 = (1ULL << 61) - 1;
static u64 mulm(u64 a, u64 b) { unsigned __int128 t = (unsigned __int128)a * b; u64 lo = (u64)(t & P61), hi = (u64)(t >> 61); u64 s = lo + hi; if (s >= P61) s -= P61; return s; }
static u64 powm(u64 a, u64 e) { u64 r = 1; while (e) { if (e & 1) r = mulm(r, a); a = mulm(a, a); e >>= 1; } return r; }
static u64 red(long v) { long long t = v % (long long)P61; if (t < 0) t += P61; return (u64)t; }
static int rank_mod(vector<vector<long>> rows) {
    if (rows.empty()) return 0;
    int n = rows[0].size();
    vector<vector<u64>> A;
    for (auto &r : rows) { vector<u64> a; for (long v : r) a.push_back(red(v)); A.push_back(a); }
    int rk = 0;
    for (int c = 0; c < n && rk < (int)A.size(); c++) {
        int p = -1; for (int r = rk; r < (int)A.size(); r++) if (A[r][c]) { p = r; break; }
        if (p < 0) continue;
        swap(A[p], A[rk]);
        u64 inv = powm(A[rk][c], P61 - 2);
        for (int r = 0; r < (int)A.size(); r++) if (r != rk && A[r][c]) {
            u64 f = mulm(A[r][c], inv);
            for (int j = c; j < n; j++) { A[r][j] = (A[r][j] + P61 - mulm(f, A[rk][j])) % P61; }
        }
        rk++;
    }
    return rk;
}

// ---------- polynomials ----------
static void trim(Poly &p) { while (!p.empty() && p.back() == 0) p.pop_back(); }
static Q evalp(const Poly &p, const Q &t) { Q s = 0; for (int i = p.size() - 1; i >= 0; i--) s = s * t + p[i]; return s; }
static Poly deriv(const Poly &p) { Poly d; for (size_t i = 1; i < p.size(); i++) d.push_back(p[i] * (long)i); trim(d); return d; }
static Poly mul(const Poly &a, const Poly &b) { if (a.empty() || b.empty()) return Poly(); Poly c(a.size() + b.size() - 1, 0); for (size_t i = 0; i < a.size(); i++) for (size_t j = 0; j < b.size(); j++) c[i + j] += a[i] * b[j]; trim(c); return c; }
static void divmodp(const Poly &a, const Poly &b, Poly &q, Poly &r) {
    r = a; trim(r); q.assign(r.size() >= b.size() ? r.size() - b.size() + 1 : 0, 0);
    while (!r.empty() && r.size() >= b.size()) {
        Q f = r.back() / b.back(); size_t sh = r.size() - b.size(); q[sh] = f;
        for (size_t i = 0; i < b.size(); i++) r[sh + i] -= f * b[i];
        trim(r);
    }
    trim(q);
}
static Poly gcdp(Poly a, Poly b) { trim(a); trim(b); while (!b.empty()) { Poly q, r; divmodp(a, b, q, r); a = b; b = r; } if (!a.empty()) { Q l = a.back(); for (auto &x : a) x /= l; } return a; }
static Poly squarefree(const Poly &p) { Poly g = gcdp(p, deriv(p)); Poly q, r; divmodp(p, g, q, r); if (!r.empty()) fail("sqf"); return q; }
static int qsgn(const Q &x) { return mpq_sgn(x.get_mpq_t()); }
// number of sign changes of the Sturm-type sequence of (A, B) at t (inf: t = +infinity)
static vector<Poly> sturm(const Poly &A, const Poly &B) {
    vector<Poly> s; s.push_back(A); Poly b = B; trim(b);
    if (b.empty()) return s;
    s.push_back(b);
    while (true) { Poly q, r; divmodp(s[s.size() - 2], s.back(), q, r); if (r.empty()) break; for (auto &x : r) x = -x; s.push_back(r); }
    return s;
}
static int changes(const vector<Poly> &s, const Q *t) {
    int prev = 0, ch = 0;
    for (auto &p : s) {
        int v = t ? qsgn(evalp(p, *t)) : qsgn(p.back());
        if (v == 0) continue;
        if (prev && v != prev) ch++;
        prev = v;
    }
    return ch;
}
// sum over roots x of P in (a, b) of sign(G(x)); P squarefree, P(a), P(b) != 0 ; b == NULL means infinity
static int taq(const Poly &P, const Poly &G, const Q &a, const Q *b) {
    if (P.size() <= 1) return 0;
    Poly R = mul(deriv(P), G);
    vector<Poly> s = sturm(P, R);
    return changes(s, &a) - changes(s, b);
}
// polynomial through (i, v_i), i = 0..n (Newton form expanded)
static Poly interp(const vector<Q> &v) {
    int n = v.size();
    vector<Q> c = v;
    for (int j = 1; j < n; j++) for (int i = n - 1; i >= j; i--) c[i] = (c[i] - c[i - 1]) / Q(j);
    Poly p; p.push_back(c[n - 1]);
    for (int i = n - 2; i >= 0; i--) { Poly f; f.push_back(Q(-i)); f.push_back(Q(1)); p = mul(p, f); if (p.empty()) p.push_back(0); p[0] += c[i]; trim(p); }
    trim(p);
    return p;
}

// ---------- certificate ----------
struct Rep { char kind; int idx; long size; Q tstar; int to; array<long, 36> T; };
struct Vertex { M6 J; vector<K6> mins; vector<array<long, 36>> gens; vector<array<long, 21>> rays; vector<Rep> reps; };

static Q rat(const string &s) { Q x(s); x.canonicalize(); return x; }

static M6 act(const array<long, 36> &T, const M6 &Y) {   // T^T Y T
    M6 R;
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) {
        Q s = 0;
        for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) if (T[a * 6 + i] && T[b * 6 + j]) s += T[a * 6 + i] * Y[a][b] * T[b * 6 + j];
        R[i][j] = s;
    }
    return R;
}
static M6 symm(const array<long, 21> &y) { M6 Y; int t = 0; for (int i = 0; i < 6; i++) for (int j = i; j < 6; j++) { Y[i][j] = y[t]; Y[j][i] = y[t]; t++; } return Y; }
static void check_gamma(const array<long, 36> &T) {
    for (int j = 0; j < 5; j++) if (T[30 + j] != 0) fail("T not in Gamma (last row)");
    if (T[35] != 1 && T[35] != -1) fail("T not in Gamma (eps)");
    vector<vector<Q>> U(5, vector<Q>(5)); for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) U[i][j] = T[i * 6 + j];
    Q d = det(U); if (d != 1 && d != -1) fail("T not in Gamma (det U)");
}

int main(int argc, char **argv) {
    if (argc < 2) { cerr << "usage: check certificate.txt" << endl; return 2; }
    ifstream in(argv[1]);
    string tok;
    vector<Vertex> V;
    vector<array<string, 3>> eqpts;
    in >> tok; if (tok != "VERTICES") fail("format"); int nv; in >> nv;
    for (int v = 0; v < nv; v++) {
        Vertex X; int idx;
        in >> tok >> idx; if (tok != "VERTEX" || idx != v) fail("format vertex");
        in >> tok; for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) { string s; in >> s; X.J[i][j] = rat(s); }
        int m; in >> tok >> m; X.mins.resize(m); for (auto &k : X.mins) for (auto &x : k) in >> x;
        int g; in >> tok >> g; X.gens.resize(g); for (auto &T : X.gens) for (auto &x : T) in >> x;
        int nr; in >> tok >> nr; X.rays.resize(nr); for (auto &y : X.rays) for (auto &x : y) in >> x;
        int nreps; in >> tok >> nreps; X.reps.resize(nreps);
        for (auto &r : X.reps) {
            string kind; in >> kind >> r.idx >> r.size; r.kind = kind[0];
            if (r.kind == 'B') { string s; in >> s; r.tstar = rat(s); in >> r.to; for (auto &x : r.T) in >> x; }
        }
        in >> tok; if (tok != "END") fail("format end");
        V.push_back(X);
    }
    int ne; in >> tok >> ne;
    for (int i = 0; i < ne; i++) { string p, a, b, c; in >> p >> a >> b >> c; eqpts.push_back({a, b, c}); }

    long total_rays = 0, total_edges = 0, eq_found = 0;
    vector<array<string, 3>> eq_mine;
    for (int v = 0; v < nv; v++) {
        Vertex &X = V[v];
        const M6 &J = X.J;
        for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (J[i][j] != J[j][i]) fail("J not symmetric");
        // C1: J in R, Min(J)
        vector<K6> mins;
        if (!in_R(J, mins)) fail("vertex not in R");
        map<K6, int> listed; for (auto &k : X.mins) listed[k]++;
        map<K6, int> mine; for (auto &k : mins) mine[k]++;
        if (listed != mine) fail("Min(J) differs from the certificate at vertex " + to_string(v));
        // C2: vertex
        vector<vector<long>> rows; for (auto &k : mins) rows.push_back(symc(k));
        if (rank_mod(rows) != 21) fail("not a vertex");
        // C3: automorphisms
        for (auto &T : X.gens) { check_gamma(T); M6 R = act(T, J); if (R != J) fail("generator is not an automorphism"); }
        // C4: rays
        map<array<long, 21>, int> pos;
        for (size_t i = 0; i < X.rays.size(); i++) {
            auto &y = X.rays[i];
            if (pos.count(y)) fail("repeated ray");
            pos[y] = i;
            mpz_class g = 0; for (long a : y) g = gcd(g, mpz_class(a)); if (g != 1) fail("ray not primitive");
            vector<vector<long>> tight;
            for (auto &k : mins) {
                vector<long> a = symc(k); long s = 0; for (int t = 0; t < 21; t++) s += a[t] * y[t];
                if (s < 0) fail("ray outside the cone"); if (s == 0) tight.push_back(a);
            }
            if (rank_mod(tight) != 20) fail("ray not extreme");
        }
        // orbits
        vector<int> par(X.rays.size()); for (size_t i = 0; i < par.size(); i++) par[i] = i;
        function<int(int)> find = [&](int a) { while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; } return a; };
        for (auto &T : X.gens) for (size_t i = 0; i < X.rays.size(); i++) {
            M6 Y = act(T, symm(X.rays[i]));
            array<long, 21> y; int t = 0;
            for (int a = 0; a < 6; a++) for (int b = a; b < 6; b++) { if (Y[a][b].get_den() != 1) fail("image not integral"); y[t++] = Y[a][b].get_num().get_si(); }
            auto it = pos.find(y); if (it == pos.end()) fail("image of a ray is not a ray");
            int a = find(i), b = find(it->second); if (a != b) par[a] = b;
        }
        map<int, long> osize; for (size_t i = 0; i < par.size(); i++) osize[find(i)]++;
        map<int, int> seen;
        for (auto &r : X.reps) {
            int o = find(r.idx);
            if (seen.count(o)) fail("two representatives in one orbit");
            seen[o] = 1;
            if (osize[o] != r.size) fail("orbit size");
        }
        if (seen.size() != osize.size()) fail("an orbit has no representative");
        total_rays += X.rays.size();
        // C8
        Q dJ = det(full(J)), dQ = det(qblock(J));
        if (dJ <= 0 && dQ < 512) fail("vertex with det J <= 0 and det Q < 512");
        if (dJ == 0 && dQ == 512) eq_mine.push_back({to_string(v), "-1", "0"});
        // C5 - C7: edges
        for (auto &r : X.reps) {
            M6 Y = symm(X.rays[r.idx]);
            const Q *Tp = nullptr; Q Tv;
            if (r.kind == 'U') {
                // Y[k] >= 0 on M*: Q_Y = lambda w w^T, r_Y = mu w, n = w.x runs over Z
                int i0 = -1; for (int i = 0; i < 5; i++) if (Y[i][i] != 0) { i0 = i; break; }
                Q s = Y[5][5];
                if (i0 < 0) {
                    for (int i = 0; i < 5; i++) for (int j = 0; j < 6; j++) if (Y[i][j] != 0) fail("unbounded: Q_Y = 0 but Y != s e e^T");
                    if (s < 0) fail("unbounded: s < 0");
                } else {
                    if (Y[i0][i0] < 0) fail("unbounded: negative diagonal");
                    // w proportional to column i0
                    vector<Q> w(5); for (int i = 0; i < 5; i++) w[i] = Y[i][i0] / Y[i0][i0];
                    mpz_class L = 1; for (auto &x : w) L = lcm(L, mpz_class(x.get_den()));
                    mpz_class G = 0; for (auto &x : w) { x *= L; G = gcd(G, x.get_num()); } for (auto &x : w) x /= G;
                    Q lam = Y[i0][i0] / (w[i0] * w[i0]);
                    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) if (Y[i][j] != lam * w[i] * w[j]) fail("unbounded: Q_Y not of rank one");
                    Q mu = Y[i0][5] / w[i0];
                    for (int i = 0; i < 5; i++) if (Y[i][5] != mu * w[i]) fail("unbounded: r_Y not parallel to w");
                    mpz_class n0 = floorq(-mu / lam);
                    for (int dn = 0; dn < 2; dn++) { Q n = Q(n0 + dn); if (lam * n * n + 2 * mu * n + s < 0) fail("unbounded: negative value"); }
                }
            } else {
                if (r.tstar <= 0) fail("t* <= 0");
                M6 M; for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) M[i][j] = J[i][j] + r.tstar * Y[i][j];
                vector<K6> mm;
                if (!in_R(M, mm)) fail("end point not in R");
                bool cut = false;
                for (auto &k : mm) if (qform(Y, k) < 0) cut = true;
                if (!cut) { cerr << "vertex " << v << " rep idx " << r.idx << " t* " << r.tstar << " |Min M| " << mm.size() << endl; for (auto &k : mm) { for (long a : k) cerr << a << " "; cerr << " Y[k]=" << qform(Y, k) << " M[k]=" << qform(M, k) << endl; } fail("t* not maximal"); }
                if (r.to < 0 || r.to >= nv) fail("target out of range");
                check_gamma(r.T);
                if (act(r.T, V[r.to].J) != M) fail("end point not equivalent to the target vertex");
                Tv = r.tstar; Tp = &Tv;
            }
            // C7: det Q >= 512 where det J = 0 on the open edge
            vector<Q> vJ, vQ;
            for (int t = 0; t <= 6; t++) {
                M6 A; for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) A[i][j] = J[i][j] + Q(t) * Y[i][j];
                vJ.push_back(det(full(A)));
                if (t <= 5) vQ.push_back(det(qblock(A)));
            }
            Poly pJ = interp(vJ), pQ = interp(vQ);
            Poly g = pQ; if (g.empty()) g.push_back(0); g[0] -= 512; trim(g);
            Q zero = 0;
            if (pJ.empty()) {
                // the whole edge lies in det J = 0
                if (g.empty()) fail("det Q constant 512 on an edge");
                Poly gs = squarefree(g);
                Q a = 0; Poly one; one.push_back(1);
                // roots of g in [0, T]
                Q mid = Tp ? *Tp / 2 : Q(1);
                if (qsgn(evalp(g, mid)) <= 0) fail("zero edge: det Q < 512 inside");
                // remove endpoint roots of gs
                Poly q, rem; Poly lin; lin.push_back(0); lin.push_back(1);
                Poly h = gs;
                if (evalp(h, zero) == 0) { divmodp(h, lin, q, rem); h = q; }
                if (Tp && evalp(h, *Tp) == 0) { Poly l2; l2.push_back(-*Tp); l2.push_back(1); divmodp(h, l2, q, rem); h = q; }
                if (taq(h, one, a, Tp) != 0) fail("zero edge: det Q - 512 has a root inside");
                total_edges++;
                continue;
            }
            Poly P = squarefree(pJ);
            Poly q, rem, lin; lin.push_back(0); lin.push_back(1);
            if (evalp(P, zero) == 0) { divmodp(P, lin, q, rem); P = q; }
            if (Tp && evalp(P, *Tp) == 0) { Poly l2; l2.push_back(-*Tp); l2.push_back(1); divmodp(P, l2, q, rem); P = q; }
            Poly one; one.push_back(1);
            int N = taq(P, one, zero, Tp);
            int S = taq(P, g, zero, Tp);
            Poly h = g.empty() ? P : gcdp(P, g);
            int N0 = taq(h, one, zero, Tp);
            if ((N - N0 - S) % 2 != 0) fail("parity");
            int neg = (N - N0 - S) / 2;
            if (neg != 0) fail("a point with det J = 0 and det Q < 512");
            if (N0 > 0) {
                // equality points: the roots of h in (0, T) (rational here)
                for (int c = 0; c < (int)eqpts.size(); c++) {
                    if (eqpts[c][0] == to_string(v) && eqpts[c][1] == to_string(r.idx)) {
                        Q t = rat(eqpts[c][2]);
                        if (evalp(pJ, t) != 0 || evalp(g, t) != 0) fail("listed equality point is not one");
                        eq_mine.push_back(eqpts[c]);
                    }
                }
                eq_found += N0;
            }
            total_edges++;
        }
    }
    eq_found += 0;
    long eq_vertices = 0; for (auto &e : eq_mine) if (e[1] == "-1") eq_vertices++;
    if ((long)eq_mine.size() != eq_found + eq_vertices || (long)eqpts.size() != (long)eq_mine.size())
        fail("equality points do not match the certificate");
    cout << "vertices: " << nv << ", rays: " << total_rays << ", edge representatives: " << total_edges
         << ", equality points: " << eq_mine.size() << endl;
    cout << "C++ checker: all checks passed" << endl;
    return 0;
}
