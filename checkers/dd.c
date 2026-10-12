/*
 * Independent enumeration of the extreme rays of the vertex cones, by the double
 * description method, written without cdd or lrs.
 *
 * For a vertex J with minimal vectors k_1, ..., k_m the cone is
 *     C_J = { Y symmetric : Y[k_i] >= 0, i = 1..m },
 * in the 21 coordinates y_ab (a <= b) of Y, where Y[k] = sum_a y_aa k_a^2 + 2 sum_{a<b} y_ab k_a k_b.
 * The program starts from the simplicial cone of 21 independent constraints, adds the
 * others one at a time, and keeps exactly the extreme rays: after adding a, a pair (p, n)
 * with a.p > 0 > a.n gives the new ray (a.p) n - (a.n) p when p and n are adjacent, that
 * is when no other ray is tight on every constraint on which both are tight (the
 * combinatorial adjacency test of Fukuda and Prodon). All arithmetic is exact
 * (64-bit with checked 128-bit products; the initial rays by GMP).
 * The rays found are compared, as sets of primitive integer vectors, with the RAYS of the
 * certificate.
 *
 * Usage: dd certificate.txt
 * Build: gcc -O2 -std=c11 -o dd dd.c -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

#define D 21
#define MAXK 64

typedef struct { int64_t c[D]; uint64_t tight; } ray;

static void fail(const char *m) { fprintf(stderr, "FAIL: %s\n", m); exit(1); }

static int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

static void normalize(int64_t *c) {
    int64_t g = 0;
    for (int i = 0; i < D; i++) g = gcd64(g, c[i]);
    if (g == 0) fail("zero ray");
    for (int i = 0; i < D; i++) c[i] /= g;
}

static int64_t dot(const int64_t *a, const int64_t *b) {
    __int128 s = 0;
    for (int i = 0; i < D; i++) s += (__int128)a[i] * b[i];
    if (s > ((__int128)1 << 62) || s < -((__int128)1 << 62)) fail("overflow in a dot product");
    return (int64_t)s;
}

static void symc(const int *k, int64_t *out) {
    int t = 0;
    for (int a = 0; a < 6; a++)
        for (int b = a; b < 6; b++) out[t++] = (int64_t)(a == b ? 1 : 2) * k[a] * k[b];
}

/* rank of a set of integer rows modulo the prime 2^61 - 1 */
static const uint64_t P = (1ULL << 61) - 1;
static uint64_t mulmod(uint64_t a, uint64_t b) { return (uint64_t)((unsigned __int128)a * b % P); }
static uint64_t powmod(uint64_t a, uint64_t e) {
    uint64_t r = 1;
    while (e) { if (e & 1) r = mulmod(r, a); a = mulmod(a, a); e >>= 1; }
    return r;
}
static int rank_mod(int64_t rows[][D], int n, int *pivrows) {
    static uint64_t A[MAXK][D];
    int idx[MAXK];
    for (int i = 0; i < n; i++) {
        idx[i] = i;
        for (int j = 0; j < D; j++) {
            int64_t v = rows[i][j] % (int64_t)P;
            A[i][j] = (uint64_t)(v < 0 ? v + (int64_t)P : v);
        }
    }
    int rk = 0;
    for (int c = 0; c < D && rk < n; c++) {
        int p = -1;
        for (int r = rk; r < n; r++) if (A[r][c]) { p = r; break; }
        if (p < 0) continue;
        for (int j = 0; j < D; j++) { uint64_t t = A[p][j]; A[p][j] = A[rk][j]; A[rk][j] = t; }
        int ti = idx[p]; idx[p] = idx[rk]; idx[rk] = ti;
        uint64_t inv = powmod(A[rk][c], P - 2);
        for (int r = 0; r < n; r++) {
            if (r != rk && A[r][c]) {
                uint64_t f = mulmod(A[r][c], inv);
                for (int j = 0; j < D; j++) A[r][j] = (A[r][j] + P - mulmod(f, A[rk][j])) % P;
            }
        }
        if (pivrows) pivrows[rk] = idx[rk];
        rk++;
    }
    return rk;
}

/* a nonzero integer vector x with B x = 0 for the 20 x 21 integer matrix B of rank 20 */
static void kernel_vector(int64_t B[20][D], int64_t *x) {
    mpz_t M[20][D], t, u;
    mpz_init(t); mpz_init(u);
    for (int i = 0; i < 20; i++) for (int j = 0; j < D; j++) mpz_init_set_si(M[i][j], B[i][j]);
    int piv[20], r = 0;
    /* fraction-free reduction to reduced echelon form (every pivot column cleared above and below) */
    for (int c = 0; c < D && r < 20; c++) {
        int p = -1;
        for (int i = r; i < 20; i++) if (mpz_sgn(M[i][c])) { p = i; break; }
        if (p < 0) continue;
        if (p != r) for (int j = 0; j < D; j++) mpz_swap(M[p][j], M[r][j]);
        for (int i = 0; i < 20; i++) {
            if (i == r || !mpz_sgn(M[i][c])) continue;
            mpz_set(u, M[i][c]);
            for (int j = 0; j < D; j++) {
                mpz_mul(M[i][j], M[i][j], M[r][c]);
                mpz_mul(t, u, M[r][j]);
                mpz_sub(M[i][j], M[i][j], t);
            }
            mpz_set_ui(t, 0);
            for (int j = 0; j < D; j++) mpz_gcd(t, t, M[i][j]);
            if (mpz_sgn(t)) for (int j = 0; j < D; j++) mpz_divexact(M[i][j], M[i][j], t);
        }
        piv[r++] = c;
    }
    if (r != 20) fail("initial constraints are dependent");
    int fr = -1;
    for (int c = 0, k = 0; c < D; c++) { if (k < 20 && piv[k] == c) k++; else { fr = c; break; } }
    /* x_fr = L, x_piv[i] = -M[i][fr] * L / M[i][piv[i]], L = lcm of the pivots */
    mpz_t L, xv[D];
    mpz_init_set_ui(L, 1);
    for (int i = 0; i < 20; i++) mpz_lcm(L, L, M[i][piv[i]]);
    for (int j = 0; j < D; j++) mpz_init_set_ui(xv[j], 0);
    mpz_set(xv[fr], L);
    for (int i = 0; i < 20; i++) {
        mpz_divexact(t, L, M[i][piv[i]]);
        mpz_mul(t, t, M[i][fr]);
        mpz_neg(xv[piv[i]], t);
    }
    mpz_set_ui(t, 0);
    for (int j = 0; j < D; j++) mpz_gcd(t, t, xv[j]);
    for (int j = 0; j < D; j++) {
        mpz_divexact(xv[j], xv[j], t);
        if (!mpz_fits_slong_p(xv[j])) fail("initial ray too large");
        x[j] = mpz_get_si(xv[j]);
    }
    for (int i = 0; i < 20; i++) for (int j = 0; j < D; j++) mpz_clear(M[i][j]);
    for (int j = 0; j < D; j++) mpz_clear(xv[j]);
    mpz_clear(L); mpz_clear(t); mpz_clear(u);
}

static int cmp_ray(const void *a, const void *b) {
    const int64_t *x = ((const ray *)a)->c, *y = ((const ray *)b)->c;
    for (int i = 0; i < D; i++) if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
    return 0;
}

static int popc(uint64_t x) { return __builtin_popcountll(x); }

/* the extreme rays of { y : a_i . y >= 0 }, i < m */
static ray *enumerate(int64_t a[][D], int m, int *nout) {
    int basis[D];
    if (rank_mod(a, m, basis) != D) fail("constraints do not have rank 21");
    int inb[MAXK] = {0}, order[MAXK], no = 0;
    for (int i = 0; i < D; i++) inb[basis[i]] = 1;
    for (int i = 0; i < m; i++) if (!inb[i]) order[no++] = i;
    size_t cap = 1024, n = 0;
    ray *R = malloc(cap * sizeof(ray));
    for (int i = 0; i < D; i++) {
        int64_t B[20][D];
        for (int j = 0, r = 0; j < D; j++) if (j != i) memcpy(B[r++], a[basis[j]], sizeof(B[0]));
        ray x; kernel_vector(B, x.c);
        int64_t s = dot(a[basis[i]], x.c);
        if (s == 0) fail("initial ray tight on its own constraint");
        if (s < 0) for (int j = 0; j < D; j++) x.c[j] = -x.c[j];
        x.tight = 0;
        for (int j = 0; j < D; j++) if (j != i) x.tight |= 1ULL << basis[j];
        R[n++] = x;
    }
    int64_t *sv = NULL;
    for (int step = 0; step < no; step++) {
        int64_t *ak = a[order[step]];
        uint64_t bit = 1ULL << order[step];
        sv = realloc(sv, n * sizeof(int64_t));
        size_t *pos = malloc(n * sizeof(size_t)), *neg = malloc(n * sizeof(size_t)), np = 0, nn = 0;
        for (size_t i = 0; i < n; i++) {
            sv[i] = dot(ak, R[i].c);
            if (sv[i] > 0) pos[np++] = i; else if (sv[i] < 0) neg[nn++] = i;
        }
        size_t ncap = n + 1024, nn2 = 0;
        ray *N = malloc(ncap * sizeof(ray));
        for (size_t i = 0; i < n; i++) {
            if (sv[i] < 0) continue;
            N[nn2] = R[i];
            if (sv[i] == 0) N[nn2].tight |= bit;
            nn2++;
        }
        for (size_t ip = 0; ip < np; ip++) {
            const ray *p = &R[pos[ip]];
            for (size_t in = 0; in < nn; in++) {
                const ray *q = &R[neg[in]];
                uint64_t S = p->tight & q->tight;
                if (popc(S) < D - 2) continue;
                int adj = 1;
                for (size_t r = 0; r < n; r++) {
                    if (r == pos[ip] || r == neg[in]) continue;
                    if ((S & ~R[r].tight) == 0) { adj = 0; break; }
                }
                if (!adj) continue;
                ray x;
                int64_t sp = sv[pos[ip]], sq = sv[neg[in]];
                for (int j = 0; j < D; j++) {
                    __int128 v = (__int128)sp * q->c[j] - (__int128)sq * p->c[j];
                    if (v > ((__int128)1 << 62) || v < -((__int128)1 << 62)) fail("overflow in a new ray");
                    x.c[j] = (int64_t)v;
                }
                normalize(x.c);
                x.tight = S | bit;
                if (nn2 == ncap) { ncap *= 2; N = realloc(N, ncap * sizeof(ray)); }
                N[nn2++] = x;
            }
        }
        free(pos); free(neg); free(R);
        R = N; n = nn2;
    }
    free(sv);
    *nout = (int)n;
    return R;
}

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: dd certificate.txt\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) fail("cannot open the certificate");
    static char line[1 << 16];
    int nv = -1, total = 0;
    if (!fgets(line, sizeof line, f) || sscanf(line, "VERTICES %d", &nv) != 1) fail("format");
    for (int v = 0; v < nv; v++) {
        int idx, m, g, nr, nreps;
        if (!fgets(line, sizeof line, f) || sscanf(line, "VERTEX %d", &idx) != 1 || idx != v) fail("format vertex");
        if (!fgets(line, sizeof line, f)) fail("format J");
        if (!fgets(line, sizeof line, f) || sscanf(line, "MIN %d", &m) != 1 || m > MAXK) fail("format MIN");
        int64_t a[MAXK][D];
        for (int i = 0; i < m; i++) {
            int k[6];
            if (!fgets(line, sizeof line, f) || sscanf(line, "%d %d %d %d %d %d", &k[0], &k[1], &k[2], &k[3], &k[4], &k[5]) != 6)
                fail("format minimal vector");
            symc(k, a[i]);
        }
        if (!fgets(line, sizeof line, f) || sscanf(line, "GENS %d", &g) != 1) fail("format GENS");
        for (int i = 0; i < g; i++) if (!fgets(line, sizeof line, f)) fail("format generator");
        if (!fgets(line, sizeof line, f) || sscanf(line, "RAYS %d", &nr) != 1) fail("format RAYS");
        ray *C = malloc((size_t)nr * sizeof(ray));
        for (int i = 0; i < nr; i++) {
            if (!fgets(line, sizeof line, f)) fail("format ray");
            char *s = line, *e;
            for (int j = 0; j < D; j++) { C[i].c[j] = strtoll(s, &e, 10); if (e == s) fail("format ray entry"); s = e; }
        }
        if (!fgets(line, sizeof line, f) || sscanf(line, "REPS %d", &nreps) != 1) fail("format REPS");
        for (int i = 0; i < nreps; i++) if (!fgets(line, sizeof line, f)) fail("format rep");
        if (!fgets(line, sizeof line, f) || strncmp(line, "END", 3) != 0) fail("format END");
        int n;
        ray *R = enumerate(a, m, &n);
        qsort(R, (size_t)n, sizeof(ray), cmp_ray);
        qsort(C, (size_t)nr, sizeof(ray), cmp_ray);
        int same = (n == nr);
        for (int i = 0; same && i < n; i++) if (cmp_ray(&R[i], &C[i]) != 0) same = 0;
        printf("vertex %d: %d constraints, %d extreme rays, certificate %d, equal: %s\n", v, m, n, nr, same ? "yes" : "no");
        fflush(stdout);
        if (!same) fail("the extreme rays differ from the certificate");
        total += n;
        free(R); free(C);
    }
    fclose(f);
    printf("dd: all %d ray sets agree with the certificate (%d rays)\n", nv, total);
    return 0;
}
