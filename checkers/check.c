/* Independent checker (C, GMP rationals) for the certificate of the two-translate problem in R^5.
   Usage: check_c certificate.txt */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void fail(const char *m) { fprintf(stderr, "FAIL: %s\n", m); exit(1); }

/* ---------- small helpers on mpq matrices ---------- */
typedef mpq_t M6[6][6];

static void m6_init(M6 A) { for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) mpq_init(A[i][j]); }
static void m6_clear(M6 A) { for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) mpq_clear(A[i][j]); }

static void qform(mpq_t out, M6 J, const long *k) {
    mpq_t t, c; mpq_init(t); mpq_init(c); mpq_set_ui(out, 0, 1);
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (k[i] && k[j]) {
        mpq_set_si(c, k[i] * k[j], 1); mpq_mul(t, J[i][j], c); mpq_add(out, out, t);
    }
    mpq_clear(t); mpq_clear(c);
}
static void symc(const long *k, long *out) {
    int t = 0;
    for (int i = 0; i < 6; i++) for (int j = i; j < 6; j++) out[t++] = k[i] * k[j] * (i == j ? 1 : 2);
}
/* determinant of the leading n x n block of a 6 x 6 matrix */
static void det_n(mpq_t out, M6 A0, int n) {
    mpq_t A[6][6], f, t; mpq_init(f); mpq_init(t);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { mpq_init(A[i][j]); mpq_set(A[i][j], A0[i][j]); }
    mpq_set_ui(out, 1, 1);
    for (int c = 0; c < n; c++) {
        int p = -1;
        for (int r = c; r < n; r++) if (mpq_sgn(A[r][c]) != 0) { p = r; break; }
        if (p < 0) { mpq_set_ui(out, 0, 1); break; }
        if (p != c) { for (int j = 0; j < n; j++) mpq_swap(A[p][j], A[c][j]); mpq_neg(out, out); }
        mpq_mul(out, out, A[c][c]);
        for (int r = c + 1; r < n; r++) if (mpq_sgn(A[r][c]) != 0) {
            mpq_div(f, A[r][c], A[c][c]);
            for (int j = c; j < n; j++) { mpq_mul(t, f, A[c][j]); mpq_sub(A[r][j], A[r][j], t); }
        }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) mpq_clear(A[i][j]);
    mpq_clear(f); mpq_clear(t);
}

/* ---------- short vectors ---------- */
typedef struct { mpq_t d[5]; mpq_t m[5][5]; } LDL;
static int ldl5(LDL *L, M6 A) {
    mpq_t v, t; mpq_init(v); mpq_init(t);
    for (int i = 0; i < 5; i++) { mpq_init(L->d[i]); for (int j = 0; j < 5; j++) mpq_init(L->m[i][j]); }
    for (int i = 0; i < 5; i++) {
        mpq_set(v, A[i][i]);
        for (int k = 0; k < i; k++) { mpq_mul(t, L->m[k][i], L->m[k][i]); mpq_mul(t, t, L->d[k]); mpq_sub(v, v, t); }
        if (mpq_sgn(v) <= 0) { mpq_clear(v); mpq_clear(t); return 0; }
        mpq_set(L->d[i], v);
        for (int j = i + 1; j < 5; j++) {
            mpq_set(L->m[i][j], A[i][j]);
            for (int k = 0; k < i; k++) { mpq_mul(t, L->m[k][i], L->m[k][j]); mpq_mul(t, t, L->d[k]); mpq_sub(L->m[i][j], L->m[i][j], t); }
            mpq_div(L->m[i][j], L->m[i][j], v);
        }
    }
    mpq_clear(v); mpq_clear(t);
    return 1;
}
static void ldl_clear(LDL *L) { for (int i = 0; i < 5; i++) { mpq_clear(L->d[i]); for (int j = 0; j < 5; j++) mpq_clear(L->m[i][j]); } }

typedef void (*cbfun)(const long *x, void *ctx);
static void rec_enum(LDL *L, mpq_t *c, long *x, int i, mpq_t R, cbfun cb, void *ctx) {
    if (i < 0) { cb(x, ctx); return; }
    mpq_t u, t, w, val, R2; mpq_init(u); mpq_init(t); mpq_init(w); mpq_init(val); mpq_init(R2);
    mpq_set(u, c[i]);
    for (int j = i + 1; j < 5; j++) { mpq_set_si(t, x[j], 1); mpq_add(t, t, c[j]); mpq_mul(t, t, L->m[i][j]); mpq_add(u, u, t); }
    mpz_t fl; mpz_init(fl);
    mpq_neg(t, u); mpz_fdiv_q(fl, mpq_numref(t), mpq_denref(t));
    long x0 = mpz_get_si(fl);
    for (int dir = 0; dir < 2; dir++) {
        long v = dir == 0 ? x0 + 1 : x0;
        while (1) {
            mpq_set_si(w, v, 1); mpq_add(w, w, u); mpq_mul(val, w, w); mpq_mul(val, val, L->d[i]);
            if (mpq_cmp(val, R) > 0) break;
            x[i] = v; mpq_sub(R2, R, val); rec_enum(L, c, x, i - 1, R2, cb, ctx);
            v += dir == 0 ? 1 : -1;
        }
    }
    x[i] = 0;
    mpz_clear(fl); mpq_clear(u); mpq_clear(t); mpq_clear(w); mpq_clear(val); mpq_clear(R2);
}

typedef struct { M6 *J; int ok; long (*mins)[6]; int nmin, cap; int layer; } RCtx;
static void collect(const long *x, void *vctx) {
    RCtx *c = (RCtx *)vctx;
    long k[6];
    if (c->layer == 0) {
        int nz = -1; for (int i = 0; i < 5; i++) if (x[i]) { nz = i; break; }
        if (nz < 0 || x[nz] < 0) return;
        for (int i = 0; i < 5; i++) k[i] = x[i];
        k[5] = 0;
    } else { for (int i = 0; i < 5; i++) k[i] = x[i]; k[5] = 1; }
    mpq_t v; mpq_init(v); qform(v, *c->J, k);
    int cmp = mpq_cmp_si(v, 4, 1);
    if (cmp < 0) c->ok = 0;
    else if (cmp == 0) {
        if (c->nmin == c->cap) { c->cap = c->cap ? 2 * c->cap : 64; c->mins = realloc(c->mins, sizeof(long[6]) * c->cap); }
        memcpy(c->mins[c->nmin++], k, sizeof(long[6]));
    }
    mpq_clear(v);
}
/* returns 1 if all values on M* are >= 4; fills the minimal vectors */
static int in_R(M6 J, long (**mins)[6], int *nmin) {
    LDL L;
    if (!ldl5(&L, J)) fail("Q block not positive definite");
    RCtx ctx = { (M6 *)J, 1, NULL, 0, 0, 0 };
    mpq_t c[5], B; for (int i = 0; i < 5; i++) mpq_init(c[i]); mpq_init(B);
    long x[5] = {0, 0, 0, 0, 0};
    mpq_set_ui(B, 4, 1);
    rec_enum(&L, c, x, 4, B, collect, &ctx);
    /* q = Q^{-1} r, sigma = s - r.q */
    mpq_t A[5][6], f, t; mpq_init(f); mpq_init(t);
    for (int i = 0; i < 5; i++) for (int j = 0; j < 6; j++) { mpq_init(A[i][j]); mpq_set(A[i][j], J[i][j < 5 ? j : 5]); }
    for (int col = 0; col < 5; col++) {
        int p = col; while (mpq_sgn(A[p][col]) == 0) p++;
        for (int j = 0; j < 6; j++) mpq_swap(A[p][j], A[col][j]);
        for (int r = 0; r < 5; r++) if (r != col && mpq_sgn(A[r][col]) != 0) {
            mpq_div(f, A[r][col], A[col][col]);
            for (int j = col; j < 6; j++) { mpq_mul(t, f, A[col][j]); mpq_sub(A[r][j], A[r][j], t); }
        }
    }
    mpq_t sigma; mpq_init(sigma); mpq_set(sigma, J[5][5]);
    for (int i = 0; i < 5; i++) { mpq_div(c[i], A[i][5], A[i][i]); mpq_mul(t, J[i][5], c[i]); mpq_sub(sigma, sigma, t); }
    mpq_set_ui(B, 4, 1); mpq_sub(B, B, sigma);
    ctx.layer = 1;
    if (mpq_sgn(B) >= 0) rec_enum(&L, c, x, 4, B, collect, &ctx);
    for (int i = 0; i < 5; i++) for (int j = 0; j < 6; j++) mpq_clear(A[i][j]);
    for (int i = 0; i < 5; i++) mpq_clear(c[i]);
    mpq_clear(B); mpq_clear(f); mpq_clear(t); mpq_clear(sigma);
    ldl_clear(&L);
    *mins = ctx.mins; *nmin = ctx.nmin;
    return ctx.ok;
}

/* ---------- rank modulo 2^61 - 1 ---------- */
typedef unsigned long long u64;
static const u64 P61 = (1ULL << 61) - 1;
static u64 mulm(u64 a, u64 b) { unsigned __int128 t = (unsigned __int128)a * b; u64 lo = (u64)(t & P61), hi = (u64)(t >> 61); u64 s = lo + hi; if (s >= P61) s -= P61; return s; }
static u64 powm(u64 a, u64 e) { u64 r = 1; while (e) { if (e & 1) r = mulm(r, a); a = mulm(a, a); e >>= 1; } return r; }
static int rank_mod(long (*rows)[21], int nr) {
    u64 (*A)[21] = malloc(sizeof(u64[21]) * (nr ? nr : 1));
    for (int i = 0; i < nr; i++) for (int j = 0; j < 21; j++) { long long v = rows[i][j] % (long long)P61; if (v < 0) v += P61; A[i][j] = v; }
    int rk = 0;
    for (int c = 0; c < 21 && rk < nr; c++) {
        int p = -1; for (int r = rk; r < nr; r++) if (A[r][c]) { p = r; break; }
        if (p < 0) continue;
        for (int j = 0; j < 21; j++) { u64 s = A[p][j]; A[p][j] = A[rk][j]; A[rk][j] = s; }
        u64 inv = powm(A[rk][c], P61 - 2);
        for (int r = 0; r < nr; r++) if (r != rk && A[r][c]) {
            u64 f = mulm(A[r][c], inv);
            for (int j = c; j < 21; j++) A[r][j] = (A[r][j] + P61 - mulm(f, A[rk][j])) % P61;
        }
        rk++;
    }
    free(A);
    return rk;
}

/* ---------- polynomials (degree <= 12 suffices) ---------- */
#define PMAX 16
typedef struct { int n; mpq_t c[PMAX]; } Poly;   /* n = number of coefficients after trimming */
static void p_init(Poly *p) { p->n = 0; for (int i = 0; i < PMAX; i++) mpq_init(p->c[i]); }
static void p_clear(Poly *p) { for (int i = 0; i < PMAX; i++) mpq_clear(p->c[i]); }
static void p_trim(Poly *p) { while (p->n > 0 && mpq_sgn(p->c[p->n - 1]) == 0) p->n--; }
static void p_set(Poly *a, const Poly *b) { a->n = b->n; for (int i = 0; i < PMAX; i++) mpq_set(a->c[i], b->c[i]); }
static void p_eval(mpq_t out, const Poly *p, const mpq_t t) {
    mpq_set_ui(out, 0, 1);
    for (int i = p->n - 1; i >= 0; i--) { mpq_mul(out, out, t); mpq_add(out, out, p->c[i]); }
}
static void p_deriv(Poly *d, const Poly *p) {
    Poly r; p_init(&r); mpq_t k; mpq_init(k);
    for (int i = 1; i < p->n; i++) { mpq_set_si(k, i, 1); mpq_mul(r.c[i - 1], p->c[i], k); }
    r.n = p->n > 0 ? p->n - 1 : 0; p_trim(&r); p_set(d, &r); p_clear(&r); mpq_clear(k);
}
static void p_mul(Poly *out, const Poly *a, const Poly *b) {
    Poly r; p_init(&r); mpq_t t; mpq_init(t);
    if (a->n && b->n) {
        if (a->n + b->n - 1 > PMAX) fail("degree");
        for (int i = 0; i < a->n; i++) for (int j = 0; j < b->n; j++) { mpq_mul(t, a->c[i], b->c[j]); mpq_add(r.c[i + j], r.c[i + j], t); }
        r.n = a->n + b->n - 1; p_trim(&r);
    }
    p_set(out, &r); p_clear(&r); mpq_clear(t);
}
static void p_divmod(Poly *q, Poly *r, const Poly *a, const Poly *b) {
    Poly Q, R; p_init(&Q); p_init(&R); p_set(&R, a); p_trim(&R);
    mpq_t f, t; mpq_init(f); mpq_init(t);
    Q.n = R.n >= b->n ? R.n - b->n + 1 : 0;
    while (R.n > 0 && R.n >= b->n) {
        mpq_div(f, R.c[R.n - 1], b->c[b->n - 1]); int sh = R.n - b->n;
        mpq_set(Q.c[sh], f);
        for (int i = 0; i < b->n; i++) { mpq_mul(t, f, b->c[i]); mpq_sub(R.c[sh + i], R.c[sh + i], t); }
        p_trim(&R);
    }
    p_trim(&Q);
    if (q) p_set(q, &Q);
    if (r) p_set(r, &R);
    p_clear(&Q); p_clear(&R); mpq_clear(f); mpq_clear(t);
}
static void p_gcd(Poly *g, const Poly *a0, const Poly *b0) {
    Poly a, b, r; p_init(&a); p_init(&b); p_init(&r);
    p_set(&a, a0); p_set(&b, b0); p_trim(&a); p_trim(&b);
    while (b.n > 0) { p_divmod(NULL, &r, &a, &b); p_set(&a, &b); p_set(&b, &r); }
    if (a.n > 0) { mpq_t l; mpq_init(l); mpq_set(l, a.c[a.n - 1]); for (int i = 0; i < a.n; i++) mpq_div(a.c[i], a.c[i], l); mpq_clear(l); }
    p_set(g, &a); p_clear(&a); p_clear(&b); p_clear(&r);
}
static void p_squarefree(Poly *out, const Poly *p) {
    Poly d, g, r; p_init(&d); p_init(&g); p_init(&r);
    p_deriv(&d, p); p_gcd(&g, p, &d); p_divmod(out, &r, p, &g);
    if (r.n) fail("sqf");
    p_clear(&d); p_clear(&g); p_clear(&r);
}
static void p_drop_root(Poly *P, const mpq_t t) {
    mpq_t v; mpq_init(v); p_eval(v, P, t);
    if (mpq_sgn(v) == 0) {
        Poly lin; p_init(&lin); mpq_neg(lin.c[0], t); mpq_set_ui(lin.c[1], 1, 1); lin.n = 2;
        p_divmod(P, NULL, P, &lin); p_clear(&lin);
    }
    mpq_clear(v);
}
/* sign changes of the Sturm sequence of (A, B) at t (t == NULL: +infinity) */
static int changes_at(Poly *s, int ns, const mpq_t *t) {
    int prev = 0, ch = 0; mpq_t v; mpq_init(v);
    for (int i = 0; i < ns; i++) {
        int sg;
        if (t) { p_eval(v, &s[i], *t); sg = mpq_sgn(v); } else sg = mpq_sgn(s[i].c[s[i].n - 1]);
        if (sg == 0) continue;
        if (prev && sg != prev) ch++;
        prev = sg;
    }
    mpq_clear(v);
    return ch;
}
/* sum over roots x of P in (a, b) of sign(G(x)); b == NULL means +infinity */
static int taq(const Poly *P, const Poly *G, const mpq_t a, const mpq_t *b) {
    if (P->n <= 1) return 0;
    Poly s[PMAX + 2]; int ns = 0;
    for (int i = 0; i < PMAX + 2; i++) p_init(&s[i]);
    Poly d; p_init(&d); p_deriv(&d, P);
    p_set(&s[0], P); p_mul(&s[1], &d, G); ns = s[1].n ? 2 : 1;
    while (ns >= 2) {
        Poly r; p_init(&r); p_divmod(NULL, &r, &s[ns - 2], &s[ns - 1]);
        if (r.n == 0) { p_clear(&r); break; }
        for (int i = 0; i < r.n; i++) mpq_neg(r.c[i], r.c[i]);
        p_set(&s[ns++], &r); p_clear(&r);
    }
    mpq_t at; mpq_init(at); mpq_set(at, a);
    int res = changes_at(s, ns, (const mpq_t *)&at) - changes_at(s, ns, b);
    mpq_clear(at);
    for (int i = 0; i < PMAX + 2; i++) p_clear(&s[i]);
    p_clear(&d);
    return res;
}
static void p_interp(Poly *out, mpq_t *v, int n) {
    mpq_t c[8], t, k; mpq_init(t); mpq_init(k);
    for (int i = 0; i < n; i++) { mpq_init(c[i]); mpq_set(c[i], v[i]); }
    for (int j = 1; j < n; j++) for (int i = n - 1; i >= j; i--) { mpq_sub(c[i], c[i], c[i - 1]); mpq_set_si(k, j, 1); mpq_div(c[i], c[i], k); }
    Poly p, f; p_init(&p); p_init(&f);
    mpq_set(p.c[0], c[n - 1]); p.n = 1; p_trim(&p);
    for (int i = n - 2; i >= 0; i--) {
        mpq_set_si(f.c[0], -i, 1); mpq_set_ui(f.c[1], 1, 1); f.n = 2;
        p_mul(&p, &p, &f); if (p.n == 0) { mpq_set_ui(p.c[0], 0, 1); p.n = 1; }
        mpq_add(p.c[0], p.c[0], c[i]); p_trim(&p);
    }
    p_set(out, &p); p_clear(&p); p_clear(&f);
    for (int i = 0; i < n; i++) mpq_clear(c[i]);
    mpq_clear(t); mpq_clear(k);
}

/* ---------- certificate ---------- */
typedef struct { char kind; int idx; long size; mpq_t tstar; int to; long T[36]; } Rep;
typedef struct { M6 J; int nmin; long (*mins)[6]; int ngen; long (*gens)[36]; int nray; long (*rays)[21]; int nrep; Rep *reps; } Vertex;

static FILE *fin;
static char tokbuf[4096];
static const char *tok(void) { if (fscanf(fin, "%4095s", tokbuf) != 1) fail("unexpected end of file"); return tokbuf; }
static long ltok(void) { return atol(tok()); }
static void qtok(mpq_t q) { if (mpq_set_str(q, tok(), 10) != 0) fail("bad rational"); mpq_canonicalize(q); }

static void act(M6 out, const long *T, M6 Y) {   /* T^T Y T */
    mpq_t s, t, c; mpq_init(s); mpq_init(t); mpq_init(c);
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) {
        mpq_set_ui(s, 0, 1);
        for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) if (T[a * 6 + i] && T[b * 6 + j]) {
            mpq_set_si(c, T[a * 6 + i] * T[b * 6 + j], 1); mpq_mul(t, Y[a][b], c); mpq_add(s, s, t);
        }
        mpq_set(out[i][j], s);
    }
    mpq_clear(s); mpq_clear(t); mpq_clear(c);
}
static int m6_eq(M6 A, M6 B) { for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (!mpq_equal(A[i][j], B[i][j])) return 0; return 1; }
static void symm(M6 Y, const long *y) { int t = 0; for (int i = 0; i < 6; i++) for (int j = i; j < 6; j++) { mpq_set_si(Y[i][j], y[t], 1); mpq_set_si(Y[j][i], y[t], 1); t++; } }
static void check_gamma(const long *T) {
    for (int j = 0; j < 5; j++) if (T[30 + j]) fail("T not in Gamma (last row)");
    if (T[35] != 1 && T[35] != -1) fail("T not in Gamma (eps)");
    M6 U; m6_init(U);
    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) mpq_set_si(U[i][j], T[i * 6 + j], 1);
    mpq_t d; mpq_init(d); det_n(d, U, 5);
    if (mpq_cmp_si(d, 1, 1) != 0 && mpq_cmp_si(d, -1, 1) != 0) fail("T not in Gamma (det U)");
    mpq_clear(d); m6_clear(U);
}
static int cmp6(const void *a, const void *b) { return memcmp(a, b, sizeof(long[6])); }

/* hash table for rays */
typedef struct { long *keys; int *vals; int cap; } Hash;
static u64 hray(const long *y) { u64 h = 1469598103934665603ULL; for (int i = 0; i < 21; i++) { h ^= (u64)(y[i] + 1000003); h *= 1099511628211ULL; } return h; }
static void h_init(Hash *H, int n) { H->cap = 1; while (H->cap < 4 * n + 4) H->cap <<= 1; H->keys = malloc(sizeof(long) * 21 * H->cap); H->vals = malloc(sizeof(int) * H->cap); for (int i = 0; i < H->cap; i++) H->vals[i] = -1; }
static int h_find(Hash *H, const long *y) { u64 i = hray(y) & (H->cap - 1); while (H->vals[i] >= 0) { if (!memcmp(H->keys + 21 * i, y, sizeof(long) * 21)) return H->vals[i]; i = (i + 1) & (H->cap - 1); } return -1; }
static void h_add(Hash *H, const long *y, int v) { u64 i = hray(y) & (H->cap - 1); while (H->vals[i] >= 0) { if (!memcmp(H->keys + 21 * i, y, sizeof(long) * 21)) fail("repeated ray"); i = (i + 1) & (H->cap - 1); } memcpy(H->keys + 21 * i, y, sizeof(long) * 21); H->vals[i] = v; }
static int *par;
static int findr(int a) { while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; } return a; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: check_c certificate.txt\n"); return 2; }
    fin = fopen(argv[1], "r"); if (!fin) fail("cannot open certificate");
    if (strcmp(tok(), "VERTICES")) fail("format");
    int nv = ltok();
    Vertex *V = calloc(nv, sizeof(Vertex));
    for (int v = 0; v < nv; v++) {
        Vertex *X = &V[v];
        if (strcmp(tok(), "VERTEX") || ltok() != v) fail("format vertex");
        tok(); m6_init(X->J);
        for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) qtok(X->J[i][j]);
        tok(); X->nmin = ltok(); X->mins = malloc(sizeof(long[6]) * X->nmin);
        for (int i = 0; i < X->nmin; i++) for (int j = 0; j < 6; j++) X->mins[i][j] = ltok();
        tok(); X->ngen = ltok(); X->gens = malloc(sizeof(long[36]) * (X->ngen ? X->ngen : 1));
        for (int i = 0; i < X->ngen; i++) for (int j = 0; j < 36; j++) X->gens[i][j] = ltok();
        tok(); X->nray = ltok(); X->rays = malloc(sizeof(long[21]) * X->nray);
        for (int i = 0; i < X->nray; i++) for (int j = 0; j < 21; j++) X->rays[i][j] = ltok();
        tok(); X->nrep = ltok(); X->reps = calloc(X->nrep, sizeof(Rep));
        for (int i = 0; i < X->nrep; i++) {
            Rep *r = &X->reps[i]; r->kind = tok()[0]; r->idx = ltok(); r->size = ltok(); mpq_init(r->tstar);
            if (r->kind == 'B') { qtok(r->tstar); r->to = ltok(); for (int j = 0; j < 36; j++) r->T[j] = ltok(); }
        }
        if (strcmp(tok(), "END")) fail("format end");
    }
    tok(); int ne = ltok();
    long (*eqv)[2] = malloc(sizeof(long[2]) * (ne ? ne : 1)); mpq_t *eqt = malloc(sizeof(mpq_t) * (ne ? ne : 1));
    for (int i = 0; i < ne; i++) { tok(); eqv[i][0] = ltok(); eqv[i][1] = ltok(); mpq_init(eqt[i]); qtok(eqt[i]); }
    int *eq_used = calloc(ne ? ne : 1, sizeof(int));

    long total_rays = 0, total_edges = 0, eq_found = 0, eq_vert = 0;
    mpq_t a, b, zero, val; mpq_init(a); mpq_init(b); mpq_init(zero); mpq_init(val);
    for (int v = 0; v < nv; v++) {
        Vertex *X = &V[v];
        for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) if (!mpq_equal(X->J[i][j], X->J[j][i])) fail("J not symmetric");
        long (*mins)[6]; int nmin;
        if (!in_R(X->J, &mins, &nmin)) fail("vertex not in R");
        if (nmin != X->nmin) fail("Min(J) differs (count)");
        long (*lm)[6] = malloc(sizeof(long[6]) * nmin); memcpy(lm, X->mins, sizeof(long[6]) * nmin);
        qsort(lm, nmin, sizeof(long[6]), cmp6); qsort(mins, nmin, sizeof(long[6]), cmp6);
        if (memcmp(lm, mins, sizeof(long[6]) * nmin)) fail("Min(J) differs");
        long (*forms)[21] = malloc(sizeof(long[21]) * nmin);
        for (int i = 0; i < nmin; i++) symc(mins[i], forms[i]);
        if (rank_mod(forms, nmin) != 21) fail("not a vertex");
        M6 R; m6_init(R);
        for (int g = 0; g < X->ngen; g++) { check_gamma(X->gens[g]); act(R, X->gens[g], X->J); if (!m6_eq(R, X->J)) fail("generator is not an automorphism"); }
        Hash H; h_init(&H, X->nray);
        long (*tight)[21] = malloc(sizeof(long[21]) * nmin);
        for (int i = 0; i < X->nray; i++) {
            long *y = X->rays[i];
            h_add(&H, y, i);
            mpz_t g; mpz_init_set_ui(g, 0); for (int j = 0; j < 21; j++) { mpz_t t; mpz_init_set_si(t, y[j]); mpz_gcd(g, g, t); mpz_clear(t); }
            if (mpz_cmp_ui(g, 1)) fail("ray not primitive");
            mpz_clear(g);
            int nt = 0;
            for (int k = 0; k < nmin; k++) {
                long s = 0; for (int j = 0; j < 21; j++) s += forms[k][j] * y[j];
                if (s < 0) fail("ray outside the cone");
                if (s == 0) memcpy(tight[nt++], forms[k], sizeof(long[21]));
            }
            if (rank_mod(tight, nt) != 20) fail("ray not extreme");
        }
        par = malloc(sizeof(int) * X->nray); for (int i = 0; i < X->nray; i++) par[i] = i;
        M6 Y, Z; m6_init(Y); m6_init(Z);
        for (int g = 0; g < X->ngen; g++) for (int i = 0; i < X->nray; i++) {
            symm(Y, X->rays[i]); act(Z, X->gens[g], Y);
            long yy[21]; int t = 0;
            for (int p = 0; p < 6; p++) for (int q = p; q < 6; q++) { if (mpz_cmp_ui(mpq_denref(Z[p][q]), 1)) fail("image not integral"); yy[t++] = mpz_get_si(mpq_numref(Z[p][q])); }
            int j = h_find(&H, yy); if (j < 0) fail("image of a ray is not a ray");
            int ra = findr(i), rb = findr(j); if (ra != rb) par[ra] = rb;
        }
        long *osize = calloc(X->nray, sizeof(long)); int norb = 0;
        for (int i = 0; i < X->nray; i++) { int r = findr(i); if (osize[r]++ == 0) norb++; }
        int *seen = calloc(X->nray, sizeof(int));
        for (int i = 0; i < X->nrep; i++) {
            int o = findr(X->reps[i].idx);
            if (seen[o]) fail("two representatives in one orbit");
            seen[o] = 1;
            if (osize[o] != X->reps[i].size) fail("orbit size");
        }
        if (norb != X->nrep) fail("an orbit has no representative");
        total_rays += X->nray;
        mpq_t dJ, dQ; mpq_init(dJ); mpq_init(dQ); det_n(dJ, X->J, 6); det_n(dQ, X->J, 5);
        if (mpq_sgn(dJ) <= 0 && mpq_cmp_si(dQ, 512, 1) < 0) fail("vertex with det J <= 0 and det Q < 512");
        if (mpq_sgn(dJ) == 0 && mpq_cmp_si(dQ, 512, 1) == 0) {
            int hit = 0; for (int e = 0; e < ne; e++) if (eqv[e][0] == v && eqv[e][1] == -1) { eq_used[e]++; hit = 1; }
            if (!hit) fail("equality vertex not listed");
            eq_vert++;
        }
        for (int ri = 0; ri < X->nrep; ri++) {
            Rep *r = &X->reps[ri];
            symm(Y, X->rays[r->idx]);
            int bounded = 0;
            if (r->kind == 'U') {
                int i0 = -1; for (int i = 0; i < 5; i++) if (mpq_sgn(Y[i][i])) { i0 = i; break; }
                if (i0 < 0) {
                    for (int i = 0; i < 5; i++) for (int j = 0; j < 6; j++) if (mpq_sgn(Y[i][j])) fail("unbounded: Q_Y = 0 but Y != s e e^T");
                    if (mpq_sgn(Y[5][5]) < 0) fail("unbounded: s < 0");
                } else {
                    if (mpq_sgn(Y[i0][i0]) < 0) fail("unbounded: negative diagonal");
                    /* w: column i0 made primitive integral */
                    mpz_t w[5], L, G; mpz_init_set_ui(L, 1); mpz_init_set_ui(G, 0);
                    mpq_t wq[5], lam, mu, t, n, f;
                    for (int i = 0; i < 5; i++) { mpq_init(wq[i]); mpq_div(wq[i], Y[i][i0], Y[i0][i0]); mpz_lcm(L, L, mpq_denref(wq[i])); }
                    for (int i = 0; i < 5; i++) { mpz_init(w[i]); mpz_mul(w[i], mpq_numref(wq[i]), L); mpz_divexact(w[i], w[i], mpq_denref(wq[i])); mpz_gcd(G, G, w[i]); }
                    for (int i = 0; i < 5; i++) { mpz_divexact(w[i], w[i], G); mpq_set_z(wq[i], w[i]); }
                    mpq_init(lam); mpq_init(mu); mpq_init(t); mpq_init(n); mpq_init(f);
                    mpq_mul(t, wq[i0], wq[i0]); mpq_div(lam, Y[i0][i0], t);
                    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) { mpq_mul(t, wq[i], wq[j]); mpq_mul(t, t, lam); if (!mpq_equal(t, Y[i][j])) fail("unbounded: Q_Y not of rank one"); }
                    mpq_div(mu, Y[i0][5], wq[i0]);
                    for (int i = 0; i < 5; i++) { mpq_mul(t, mu, wq[i]); if (!mpq_equal(t, Y[i][5])) fail("unbounded: r_Y not parallel to w"); }
                    mpq_div(t, mu, lam); mpq_neg(t, t);
                    mpz_t fl; mpz_init(fl); mpz_fdiv_q(fl, mpq_numref(t), mpq_denref(t));
                    for (int dn = 0; dn < 2; dn++) {
                        mpq_set_z(n, fl); if (dn) { mpq_set_ui(t, 1, 1); mpq_add(n, n, t); }
                        mpq_mul(f, n, n); mpq_mul(f, f, lam); mpq_mul(t, mu, n); mpq_add(f, f, t); mpq_add(f, f, t); mpq_add(f, f, Y[5][5]);
                        if (mpq_sgn(f) < 0) fail("unbounded: negative value");
                    }
                    mpz_clear(fl);
                    for (int i = 0; i < 5; i++) { mpz_clear(w[i]); mpq_clear(wq[i]); }
                    mpz_clear(L); mpz_clear(G); mpq_clear(lam); mpq_clear(mu); mpq_clear(t); mpq_clear(n); mpq_clear(f);
                }
            } else {
                if (mpq_sgn(r->tstar) <= 0) fail("t* <= 0");
                M6 M; m6_init(M); mpq_t t; mpq_init(t);
                for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) { mpq_mul(t, r->tstar, Y[i][j]); mpq_add(M[i][j], X->J[i][j], t); }
                long (*mm)[6]; int nm;
                if (!in_R(M, &mm, &nm)) fail("end point not in R");
                int cut = 0; for (int i = 0; i < nm; i++) { qform(t, Y, mm[i]); if (mpq_sgn(t) < 0) cut = 1; }
                if (!cut) fail("t* not maximal");
                free(mm);
                if (r->to < 0 || r->to >= nv) fail("target out of range");
                check_gamma(r->T);
                act(Z, r->T, V[r->to].J);
                if (!m6_eq(Z, M)) fail("end point not equivalent to the target vertex");
                m6_clear(M); mpq_clear(t);
                bounded = 1;
            }
            /* polynomials det J(t), det Q(t) */
            mpq_t vJ[7], vQ[6], tt; mpq_init(tt);
            M6 A; m6_init(A);
            for (int s = 0; s <= 6; s++) {
                mpq_init(vJ[s]); if (s <= 5) mpq_init(vQ[s]);
                for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) { mpq_set_si(tt, s, 1); mpq_mul(tt, tt, Y[i][j]); mpq_add(A[i][j], X->J[i][j], tt); }
                det_n(vJ[s], A, 6); if (s <= 5) det_n(vQ[s], A, 5);
            }
            Poly pJ, pQ, g, P, h, one, tmp; p_init(&pJ); p_init(&pQ); p_init(&g); p_init(&P); p_init(&h); p_init(&one); p_init(&tmp);
            p_interp(&pJ, vJ, 7); p_interp(&pQ, vQ, 6);
            p_set(&g, &pQ); if (g.n == 0) g.n = 1; mpq_set_si(tt, 512, 1); mpq_sub(g.c[0], g.c[0], tt); p_trim(&g);
            mpq_set_ui(one.c[0], 1, 1); one.n = 1;
            const mpq_t *Tp = bounded ? (const mpq_t *)&r->tstar : NULL;
            if (pJ.n == 0) {
                if (g.n == 0) fail("det Q constant 512 on an edge");
                if (bounded) { mpq_set_ui(tt, 2, 1); mpq_div(tt, r->tstar, tt); } else mpq_set_ui(tt, 1, 1);
                p_eval(val, &g, tt); if (mpq_sgn(val) <= 0) fail("zero edge: det Q < 512 inside");
                p_squarefree(&h, &g); p_drop_root(&h, zero); if (bounded) p_drop_root(&h, r->tstar);
                if (taq(&h, &one, zero, Tp) != 0) fail("zero edge: det Q - 512 has a root inside");
            } else {
                p_squarefree(&P, &pJ); p_drop_root(&P, zero); if (bounded) p_drop_root(&P, r->tstar);
                int N = taq(&P, &one, zero, Tp);
                int S = taq(&P, &g, zero, Tp);
                if (g.n) p_gcd(&h, &P, &g); else p_set(&h, &P);
                int N0 = taq(&h, &one, zero, Tp);
                if ((N - N0 - S) % 2) fail("parity");
                if ((N - N0 - S) / 2 != 0) fail("a point with det J = 0 and det Q < 512");
                if (N0 > 0) {
                    for (int e = 0; e < ne; e++) if (eqv[e][0] == v && eqv[e][1] == r->idx) {
                        p_eval(val, &pJ, eqt[e]); if (mpq_sgn(val)) fail("listed equality point is not one");
                        p_eval(val, &g, eqt[e]); if (mpq_sgn(val)) fail("listed equality point is not one");
                        eq_used[e]++; eq_found++;
                    }
                    long cnt = 0; for (int e = 0; e < ne; e++) if (eqv[e][0] == v && eqv[e][1] == r->idx) cnt++;
                    if (cnt != N0) fail("equality points on an edge do not match the certificate");
                }
            }
            for (int s = 0; s <= 6; s++) { mpq_clear(vJ[s]); if (s <= 5) mpq_clear(vQ[s]); }
            mpq_clear(tt); m6_clear(A);
            p_clear(&pJ); p_clear(&pQ); p_clear(&g); p_clear(&P); p_clear(&h); p_clear(&one); p_clear(&tmp);
            total_edges++;
        }
        mpq_clear(dJ); mpq_clear(dQ);
        free(osize); free(seen); free(par); free(H.keys); free(H.vals); free(tight); free(forms); free(lm); free(mins);
        m6_clear(R); m6_clear(Y); m6_clear(Z);
    }
    for (int e = 0; e < ne; e++) if (eq_used[e] != 1) fail("equality points do not match the certificate");
    printf("vertices: %d, rays: %ld, edge representatives: %ld, equality points: %ld\n", nv, total_rays, total_edges, eq_found + eq_vert);
    printf("C checker: all checks passed\n");
    return 0;
}
