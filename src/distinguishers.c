#include "f2ld.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int64_t ipow(int p, int k) { int64_t T = 1; for (int i = 0; i < k; i++) T *= p; return T; }

static int64_t lin_index(const int32_t *a, int k, int p) {
    int64_t idx = 0, mul = 1;
    for (int d = 0; d < k; d++) { idx += (int64_t)a[d] * mul; mul *= p; }
    return idx;
}

static void decode_index(int64_t idx, int k, int p, int32_t *out) {
    for (int d = 0; d < k; d++) { out[d] = (int32_t)(idx % p); idx /= p; }
}

void result_free(result_t *r) { if (r->guess) { free(r->guess); r->guess = NULL; } }

result_t distinguisher_enum(const samples_t *s) {
    result_t r; memset(&r, 0, sizeof(r));
    int p = s->p, k = s->k;
    int64_t T = ipow(p, k);
    double *cph = malloc(sizeof(double) * p), *sph = malloc(sizeof(double) * p);
    for (int m = 0; m < p; m++) { cph[m] = cos(-2.0 * M_PI * m / p); sph[m] = sin(-2.0 * M_PI * m / p); }
    double t0 = now_seconds();
    double best = -1e300; int64_t bestc = 0;
    int32_t sd[256];
    for (int64_t c = 0; c < T; c++) {
        decode_index(c, k, p, sd);
        double acc = 0.0;
        for (int64_t j = 0; j < s->N; j++) {
            int64_t inner = 0;
            const int32_t *a = &s->addr[j * k];
            for (int d = 0; d < k; d++) inner += (int64_t)a[d] * sd[d];
            int m = (int)(((inner % p) + p) % p);
            acc += creal(s->w[j]) * cph[m] - cimag(s->w[j]) * sph[m];
        }
        if (acc > best) { best = acc; bestc = c; }
    }
    r.runtime_s = now_seconds() - t0;
    r.guess = malloc(sizeof(int32_t) * k);
    decode_index(bestc, k, p, r.guess);
    r.peak = best; r.found = 1; r.Teff = T; r.Neff = s->N;
    r.peak_bytes = (double)(2 * p) * sizeof(double);
    free(cph); free(sph);
    return r;
}

result_t distinguisher_dense_fft(const samples_t *s) {
    result_t r; memset(&r, 0, sizeof(r));
    int p = s->p, k = s->k;
    int64_t T = ipow(p, k);
    double _Complex *g = calloc((size_t)T, sizeof(double _Complex));
    double t0 = now_seconds();
    for (int64_t j = 0; j < s->N; j++) g[lin_index(&s->addr[j * k], k, p)] += s->w[j];
    ndft_double(g, p, k, 0);
    double best = -1e300; int64_t bestc = 0;
    for (int64_t c = 0; c < T; c++) if (creal(g[c]) > best) { best = creal(g[c]); bestc = c; }
    r.runtime_s = now_seconds() - t0;
    r.guess = malloc(sizeof(int32_t) * k);
    decode_index(bestc, k, p, r.guess);
    r.peak = best; r.found = 1; r.Teff = T; r.Neff = s->N;
    r.peak_bytes = (double)T * sizeof(double _Complex);
    free(g);
    return r;
}

result_t distinguisher_matzov(const samples_t *s, int k_enum) {
    result_t r; memset(&r, 0, sizeof(r));
    int p = s->p, k = s->k;
    if (k_enum < 0) k_enum = 0;
    if (k_enum > k - 1) k_enum = k - 1;
    int k_fft = k - k_enum;
    int64_t Tf = ipow(p, k_fft);
    int64_t Te = ipow(p, k_enum);
    double _Complex *g = malloc(sizeof(double _Complex) * (size_t)Tf);
    double *cph = malloc(sizeof(double) * p), *sph = malloc(sizeof(double) * p);
    for (int m = 0; m < p; m++) { cph[m] = cos(-2.0 * M_PI * m / p); sph[m] = sin(-2.0 * M_PI * m / p); }
    double t0 = now_seconds();
    double best = -1e300; int64_t best_ge = 0, best_cf = 0;
    int32_t ge[256];
    for (int64_t e = 0; e < Te; e++) {
        decode_index(e, k_enum, p, ge);
        for (int64_t t = 0; t < Tf; t++) g[t] = 0.0;
        for (int64_t j = 0; j < s->N; j++) {
            const int32_t *a = &s->addr[j * k];
            int64_t inner_e = 0;
            for (int d = 0; d < k_enum; d++) inner_e += (int64_t)a[d + k_fft] * ge[d];
            int m = (int)(((inner_e % p) + p) % p);
            double wr = creal(s->w[j]) * cph[m] - cimag(s->w[j]) * sph[m];
            double wi = creal(s->w[j]) * sph[m] + cimag(s->w[j]) * cph[m];
            g[lin_index(a, k_fft, p)] += wr + I * wi;
        }
        ndft_double(g, p, k_fft, 0);
        for (int64_t c = 0; c < Tf; c++) if (creal(g[c]) > best) { best = creal(g[c]); best_cf = c; best_ge = e; }
    }
    r.runtime_s = now_seconds() - t0;
    r.guess = malloc(sizeof(int32_t) * k);
    decode_index(best_cf, k_fft, p, r.guess);
    decode_index(best_ge, k_enum, p, r.guess + k_fft);
    r.peak = best; r.found = 1; r.Teff = Tf; r.Neff = s->N;
    r.peak_bytes = (double)Tf * sizeof(double _Complex);
    free(g); free(cph); free(sph);
    return r;
}

static int cmp_i64(const void *x, const void *y) {
    int64_t a = *(const int64_t *)x, b = *(const int64_t *)y;
    return (a > b) - (a < b);
}

result_t distinguisher_sparse(const samples_t *s) {
    result_t r; memset(&r, 0, sizeof(r));
    int p = s->p, k = s->k;
    int64_t T = ipow(p, k);
    int64_t *keys = malloc(sizeof(int64_t) * (size_t)s->N);
    for (int64_t j = 0; j < s->N; j++) keys[j] = lin_index(&s->addr[j * k], k, p);
    int64_t *order = malloc(sizeof(int64_t) * (size_t)s->N);
    for (int64_t j = 0; j < s->N; j++) order[j] = j;
    int64_t *sk = malloc(sizeof(int64_t) * (size_t)s->N);
    memcpy(sk, keys, sizeof(int64_t) * (size_t)s->N);
    qsort(sk, s->N, sizeof(int64_t), cmp_i64);
    int64_t m = 0;
    for (int64_t j = 0; j < s->N; j++) if (j == 0 || sk[j] != sk[j - 1]) m++;
    int64_t *ukey = malloc(sizeof(int64_t) * (size_t)m);
    double _Complex *uw = calloc((size_t)m, sizeof(double _Complex));
    int64_t idx = -1;
    for (int64_t j = 0; j < s->N; j++) {
        if (j == 0 || sk[j] != sk[j - 1]) { idx++; ukey[idx] = sk[j]; }
    }
    for (int64_t j = 0; j < s->N; j++) {
        int64_t lo = 0, hi = m - 1, pos = 0;
        while (lo <= hi) { int64_t mid = (lo + hi) / 2; if (ukey[mid] < keys[j]) lo = mid + 1; else { pos = mid; hi = mid - 1; } }
        uw[pos] += s->w[j];
    }
    int32_t (*udig)[256] = malloc(sizeof(int32_t) * 256 * (size_t)m);
    for (int64_t t = 0; t < m; t++) decode_index(ukey[t], k, p, udig[t]);
    double *cph = malloc(sizeof(double) * p), *sph = malloc(sizeof(double) * p);
    for (int mm = 0; mm < p; mm++) { cph[mm] = cos(-2.0 * M_PI * mm / p); sph[mm] = sin(-2.0 * M_PI * mm / p); }
    double t0 = now_seconds();
    double best = -1e300; int64_t bestc = 0;
    int32_t sd[256];
    for (int64_t c = 0; c < T; c++) {
        decode_index(c, k, p, sd);
        double acc = 0.0;
        for (int64_t t = 0; t < m; t++) {
            int64_t inner = 0;
            for (int d = 0; d < k; d++) inner += (int64_t)udig[t][d] * sd[d];
            int mm = (int)(((inner % p) + p) % p);
            acc += creal(uw[t]) * cph[mm] - cimag(uw[t]) * sph[mm];
        }
        if (acc > best) { best = acc; bestc = c; }
    }
    r.runtime_s = now_seconds() - t0;
    r.guess = malloc(sizeof(int32_t) * k);
    decode_index(bestc, k, p, r.guess);
    r.peak = best; r.found = 1; r.Teff = T; r.Neff = m;
    r.peak_bytes = (double)m * (sizeof(int64_t) + sizeof(double _Complex) + k * sizeof(int32_t));
    free(keys); free(order); free(sk); free(ukey); free(uw); free(udig); free(cph); free(sph);
    return r;
}

result_t distinguisher_ours_split(const samples_t *s, const hintset_t *hs, int k_enum) {
    result_t r; memset(&r, 0, sizeof(r));
    int64_t Teff = 0, Neff = 0;
    int folded = (hs && hs->count > 0);
    samples_t *f = folded ? fold_hints(s, hs, &Teff, &Neff) : (samples_t *)s;
    int p = f->p, k = f->k;
    if (k_enum < 0) k_enum = 0;
    if (k_enum > k - 1) k_enum = k - 1;
    int k_fft = k - k_enum;
    int64_t Tf = ipow(p, k_fft);
    int64_t Te = ipow(p, k_enum);
    const int Q = 1 << 14;
    int32_t *re = malloc(sizeof(int32_t) * (size_t)Tf);
    int32_t *im = malloc(sizeof(int32_t) * (size_t)Tf);
    double *cph = malloc(sizeof(double) * p), *sph = malloc(sizeof(double) * p);
    for (int m = 0; m < p; m++) { cph[m] = cos(-2.0 * M_PI * m / p); sph[m] = sin(-2.0 * M_PI * m / p); }
    double t0 = now_seconds();
    int64_t best = INT64_MIN; int64_t best_cf = 0, best_ge = 0;
    int32_t ge[256];
    for (int64_t e = 0; e < Te; e++) {
        decode_index(e, k_enum, p, ge);
        memset(re, 0, sizeof(int32_t) * (size_t)Tf);
        memset(im, 0, sizeof(int32_t) * (size_t)Tf);
        for (int64_t j = 0; j < f->N; j++) {
            const int32_t *a = &f->addr[j * k];
            int64_t inner_e = 0;
            for (int d = 0; d < k_enum; d++) inner_e += (int64_t)a[d + k_fft] * ge[d];
            int m = (int)(((inner_e % p) + p) % p);
            double wr = creal(f->w[j]) * cph[m] - cimag(f->w[j]) * sph[m];
            double wi = creal(f->w[j]) * sph[m] + cimag(f->w[j]) * cph[m];
            int64_t c = lin_index(a, k_fft, p);
            re[c] += (int32_t)lround(wr * Q);
            im[c] += (int32_t)lround(wi * Q);
        }
        int shift = 0;
        ndft_fixed(re, im, p, k_fft, &shift);
        for (int64_t c = 0; c < Tf; c++) if ((int64_t)re[c] > best) { best = re[c]; best_cf = c; best_ge = e; }
    }
    r.runtime_s = now_seconds() - t0;
    r.guess = malloc(sizeof(int32_t) * k);
    decode_index(best_cf, k_fft, p, r.guess);
    decode_index(best_ge, k_enum, p, r.guess + k_fft);
    r.peak = (double)best; r.found = 1; r.Teff = Tf; r.Neff = f->N;
    r.peak_bytes = (double)Tf * 2 * sizeof(int32_t);
    free(re); free(im); free(cph); free(sph);
    if (folded) samples_free(f);
    return r;
}

result_t distinguisher_ours(const samples_t *s, const hintset_t *hs) {
    return distinguisher_ours_split(s, hs, 0);
}

result_t distinguisher_prior(const samples_t *s, const hintset_t *hs, double lambda) {
    result_t r; memset(&r, 0, sizeof(r));
    int p = s->p, k = s->k;
    int64_t T = ipow(p, k);
    double _Complex *g = calloc((size_t)T, sizeof(double _Complex));
    double *bias = NULL;
    double t0 = now_seconds();
    for (int64_t j = 0; j < s->N; j++) g[lin_index(&s->addr[j * k], k, p)] += s->w[j];
    ndft_double(g, p, k, 0);
    if (hs && hs->count > 0 && lambda > 0.0) {
        double **lp = calloc(k, sizeof(double *));
        for (int i = 0; i < hs->count; i++) {
            hint_t *h = &hs->h[i];
            if (h->kind == HINT_APPROX && h->prior && h->coord >= 0 && h->coord < k) {
                double *v = malloc(sizeof(double) * p);
                for (int u = 0; u < p; u++) { double pr = h->prior[u] < 1e-9 ? 1e-9 : h->prior[u]; v[u] = lambda * log(pr); }
                if (lp[h->coord]) free(lp[h->coord]);
                lp[h->coord] = v;
            }
        }
        int32_t sd[256];
        double sc = 0.0;
        for (int64_t j = 0; j < s->N; j++) sc += cabs(s->w[j]);
        sc = sc > 0 ? sc : 1.0;
        for (int64_t c = 0; c < T; c++) {
            decode_index(c, k, p, sd);
            double add = 0.0;
            for (int d = 0; d < k; d++) if (lp[d]) add += lp[d][sd[d]];
            g[c] += sc * add / (double)s->N;
        }
        for (int d = 0; d < k; d++) free(lp[d]);
        free(lp);
    }
    double best = -1e300; int64_t bestc = 0;
    for (int64_t c = 0; c < T; c++) if (creal(g[c]) > best) { best = creal(g[c]); bestc = c; }
    r.runtime_s = now_seconds() - t0;
    r.guess = malloc(sizeof(int32_t) * k);
    decode_index(bestc, k, p, r.guess);
    r.peak = best; r.found = 1; r.Teff = T; r.Neff = s->N;
    r.peak_bytes = (double)T * sizeof(double _Complex);
    free(g); free(bias);
    return r;
}
