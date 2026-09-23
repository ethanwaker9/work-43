#include "f2ld.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

hintset_t *hintset_new(void) {
    hintset_t *hs = calloc(1, sizeof(hintset_t));
    hs->cap = 8;
    hs->h = calloc(hs->cap, sizeof(hint_t));
    return hs;
}

void hintset_free(hintset_t *hs) {
    if (!hs) return;
    for (int i = 0; i < hs->count; i++) free(hs->h[i].prior);
    free(hs->h); free(hs);
}

static void hintset_push(hintset_t *hs, hint_t h) {
    if (hs->count == hs->cap) {
        hs->cap *= 2;
        hs->h = realloc(hs->h, hs->cap * sizeof(hint_t));
    }
    hs->h[hs->count++] = h;
}

void leak_hamming(const samples_t *s, hintset_t *hs, double snr, int nhints, uint64_t seed) {
    uint64_t st = seed ? seed : 0xdeadbeefULL;
    int p = s->p;
    int used[256];
    memset(used, 0, sizeof(used));
    for (int t = 0; t < nhints && t < s->k; t++) {
        int d;
        do { d = (int)(xorshift64(&st) % (uint64_t)s->k); } while (used[d]);
        used[d] = 1;
        int truev = s->secret[d];
        double *post = malloc(sizeof(double) * p);
        double Z = 0.0;
        for (int v = 0; v < p; v++) {
            int hw_t = __builtin_popcount((unsigned)truev);
            int hw_v = __builtin_popcount((unsigned)v);
            double d2 = (double)(hw_t - hw_v);
            double leak = d2 + gaussrand(&st) / (snr > 1e-6 ? snr : 1e-6);
            post[v] = exp(-0.5 * leak * leak * snr * snr);
            Z += post[v];
        }
        for (int v = 0; v < p; v++) post[v] /= Z;
        hint_t h;
        h.coord = d; h.modulus = p; h.kind = HINT_APPROX; h.prior = post;
        double best = -1; int bestv = 0;
        for (int v = 0; v < p; v++) if (post[v] > best) { best = post[v]; bestv = v; }
        h.value = bestv;
        hintset_push(hs, h);
    }
}

samples_t *fold_hints(const samples_t *in, const hintset_t *hs, int64_t *Teff_out, int64_t *Neff_out) {
    int p = in->p, k = in->k;
    int elim[256];
    int elimval[256];
    memset(elim, 0, sizeof(elim));
    double conf_thr = 0.90;
    int ne = 0;
    if (hs) {
        for (int i = 0; i < hs->count; i++) {
            hint_t *h = &hs->h[i];
            int d = h->coord;
            if (d < 0 || d >= k || elim[d]) continue;
            int pin = 0, val = h->value;
            if (h->kind == HINT_PERFECT) pin = 1;
            else if (h->kind == HINT_MODULAR && h->modulus == p) pin = 1;
            else if (h->kind == HINT_APPROX && h->prior) {
                double best = 0; for (int v = 0; v < p; v++) if (h->prior[v] > best) best = h->prior[v];
                if (best >= conf_thr) pin = 1;
            }
            if (pin) { elim[d] = 1; elimval[d] = val; ne++; }
        }
    }
    int kp = k - ne;
    if (kp < 1) kp = 1;
    int map[256]; int mc = 0;
    for (int d = 0; d < k; d++) if (!elim[d]) map[mc++] = d;
    if (mc == 0) { map[0] = 0; mc = 1; elim[0] = 0; }
    samples_t *out = samples_alloc(in->N, mc, p);
    for (int j = 0; j < mc; j++) out->secret[j] = in->secret[map[j]];
    for (int64_t j = 0; j < in->N; j++) {
        int64_t elim_inner = 0;
        for (int d = 0; d < k; d++) if (elim[d]) elim_inner += (int64_t)in->addr[j * k + d] * elimval[d];
        double ang = -2.0 * M_PI * (double)(((elim_inner % p) + p) % p) / (double)p;
        double _Complex rot = cos(ang) + I * sin(ang);
        out->w[j] = in->w[j] * rot;
        for (int c = 0; c < mc; c++) out->addr[j * mc + c] = in->addr[j * k + map[c]];
    }
    int64_t T = 1; for (int c = 0; c < mc; c++) T *= p;
    if (Teff_out) *Teff_out = T;
    if (Neff_out) *Neff_out = in->N;
    return out;
}
