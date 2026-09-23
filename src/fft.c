#include "f2ld.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void fft_radix2(double _Complex *a, int64_t n, int inv) {
    for (int64_t i = 1, j = 0; i < n; i++) {
        int64_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { double _Complex t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int64_t len = 2; len <= n; len <<= 1) {
        double ang = (inv ? 2.0 : -2.0) * M_PI / (double)len;
        double _Complex wl = cos(ang) + I * sin(ang);
        for (int64_t i = 0; i < n; i += len) {
            double _Complex w = 1.0;
            for (int64_t k = 0; k < len / 2; k++) {
                double _Complex u = a[i + k];
                double _Complex v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wl;
            }
        }
    }
}

static void dft_naive_line(double _Complex *buf, int p, int inv, const double _Complex *root) {
    double _Complex tmp[64];
    for (int u = 0; u < p; u++) {
        double _Complex acc = 0.0;
        for (int t = 0; t < p; t++) acc += buf[t] * root[(u * t) % p];
        tmp[u] = acc;
    }
    (void)inv;
    for (int u = 0; u < p; u++) buf[u] = tmp[u];
}

static int is_pow2(int x) { return x && ((x & (x - 1)) == 0); }

void ndft_double(double _Complex *a, int p, int k, int inv) {
    int64_t T = 1;
    for (int i = 0; i < k; i++) T *= p;
    double _Complex *root = NULL;
    if (!is_pow2(p)) {
        root = malloc(sizeof(double _Complex) * p);
        for (int r = 0; r < p; r++) {
            double ang = (inv ? 2.0 : -2.0) * M_PI * r / p;
            root[r] = cos(ang) + I * sin(ang);
        }
    }
    double _Complex *buf = malloc(sizeof(double _Complex) * p);
    int64_t stride = 1;
    for (int d = 0; d < k; d++) {
        int64_t block = stride * p;
        for (int64_t base = 0; base < T; base += block) {
            for (int64_t off = 0; off < stride; off++) {
                int64_t start = base + off;
                for (int j = 0; j < p; j++) buf[j] = a[start + (int64_t)j * stride];
                if (is_pow2(p)) fft_radix2(buf, p, inv);
                else dft_naive_line(buf, p, inv, root);
                for (int j = 0; j < p; j++) a[start + (int64_t)j * stride] = buf[j];
            }
        }
        stride = block;
    }
    free(buf);
    if (root) free(root);
}

#define FXBITS 15
#define FXONE (1 << FXBITS)

static void fft_radix2_fixed(int64_t *re, int64_t *im, int n, const int32_t *tw_c, const int32_t *tw_s) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            int64_t t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1;
        int step = n / len;
        for (int i = 0; i < n; i += len) {
            for (int k = 0; k < half; k++) {
                int idx = k * step;
                int32_t c = tw_c[idx], s = tw_s[idx];
                int64_t vr = re[i + k + half], vi = im[i + k + half];
                int64_t tr = (vr * c - vi * s) >> FXBITS;
                int64_t ti = (vr * s + vi * c) >> FXBITS;
                int64_t ur = re[i + k], ui = im[i + k];
                re[i + k] = ur + tr; im[i + k] = ui + ti;
                re[i + k + half] = ur - tr; im[i + k + half] = ui - ti;
            }
        }
    }
}

static void dft_naive_fixed(int64_t *re, int64_t *im, int p, const int32_t *tc, const int32_t *ts) {
    int64_t or_[64], oi_[64];
    for (int u = 0; u < p; u++) {
        int64_t ar = 0, ai = 0;
        for (int t = 0; t < p; t++) {
            int idx = (u * t) % p;
            int64_t c = tc[idx], s = ts[idx];
            ar += (re[t] * c - im[t] * s) >> FXBITS;
            ai += (re[t] * s + im[t] * c) >> FXBITS;
        }
        or_[u] = ar; oi_[u] = ai;
    }
    for (int u = 0; u < p; u++) { re[u] = or_[u]; im[u] = oi_[u]; }
}

void ndft_fixed(int32_t *re32, int32_t *im32, int p, int k, int *shift_out) {
    int64_t T = 1;
    for (int i = 0; i < k; i++) T *= p;
    int pw2 = p && ((p & (p - 1)) == 0);
    int32_t *tw_c = malloc(sizeof(int32_t) * p);
    int32_t *tw_s = malloc(sizeof(int32_t) * p);
    for (int r = 0; r < p; r++) {
        double ang = -2.0 * M_PI * r / p;
        tw_c[r] = (int32_t)lround(cos(ang) * FXONE);
        tw_s[r] = (int32_t)lround(sin(ang) * FXONE);
    }
    int64_t *br = malloc(sizeof(int64_t) * p);
    int64_t *bi = malloc(sizeof(int64_t) * p);
    int total_shift = 0;
    int64_t stride = 1;
    for (int d = 0; d < k; d++) {
        int64_t block = stride * p;
        int64_t maxmag = 0;
        for (int64_t base = 0; base < T; base += block) {
            for (int64_t off = 0; off < stride; off++) {
                int64_t start = base + off;
                for (int j = 0; j < p; j++) { br[j] = re32[start + j * stride]; bi[j] = im32[start + j * stride]; }
                if (pw2) fft_radix2_fixed(br, bi, p, tw_c, tw_s);
                else dft_naive_fixed(br, bi, p, tw_c, tw_s);
                for (int j = 0; j < p; j++) {
                    re32[start + j * stride] = (int32_t)br[j];
                    im32[start + j * stride] = (int32_t)bi[j];
                    int64_t m = llabs(br[j]); if (m > maxmag) maxmag = m;
                    m = llabs(bi[j]); if (m > maxmag) maxmag = m;
                }
            }
        }
        int sh = 0;
        while (maxmag > (1 << 29)) { maxmag >>= 1; sh++; }
        if (sh) {
            for (int64_t t = 0; t < T; t++) { re32[t] >>= sh; im32[t] >>= sh; }
            total_shift += sh;
        }
        stride = block;
    }
    free(br); free(bi); free(tw_c); free(tw_s);
    if (shift_out) *shift_out = total_shift;
}
