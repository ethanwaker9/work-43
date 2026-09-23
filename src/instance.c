#include "f2ld.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

samples_t *samples_alloc(int64_t N, int k, int p) {
    samples_t *s = calloc(1, sizeof(samples_t));
    s->N = N; s->k = k; s->p = p;
    s->addr = malloc(sizeof(int32_t) * (size_t)N * k);
    s->w = malloc(sizeof(double _Complex) * (size_t)N);
    s->secret = calloc(k, sizeof(int32_t));
    return s;
}

void samples_free(samples_t *s) {
    if (!s) return;
    free(s->addr); free(s->w); free(s->secret); free(s);
}

void plant_instance(samples_t *s, double contrast, uint64_t seed) {
    uint64_t st = seed ? seed : 0x9e3779b97f4a7c15ULL;
    if (contrast > 0.999999) contrast = 0.999999;
    if (contrast < 1e-9) contrast = 1e-9;
    double se = sqrt(-log(contrast) / (2.0 * M_PI * M_PI));
    for (int d = 0; d < s->k; d++) s->secret[d] = (int32_t)(xorshift64(&st) % (uint64_t)s->p);
    for (int64_t j = 0; j < s->N; j++) {
        int64_t inner = 0;
        for (int d = 0; d < s->k; d++) {
            int32_t a = (int32_t)(xorshift64(&st) % (uint64_t)s->p);
            s->addr[j * s->k + d] = a;
            inner += (int64_t)a * s->secret[d];
        }
        double phase = 2.0 * M_PI * (double)(inner % s->p) / (double)s->p;
        double eps = gaussrand(&st) * se;
        double ang = phase + 2.0 * M_PI * eps;
        s->w[j] = cos(ang) + I * sin(ang);
    }
}

double instance_contrast(const params_t *pp, double blocksize_delta) {
    double ell = blocksize_delta;
    double lat = ell * pp->sigma_e / (double)pp->q;
    double round = 1.0 / (double)(pp->p) / sqrt(12.0);
    double tau = lat * lat + round * round;
    return exp(-2.0 * M_PI * M_PI * tau);
}
