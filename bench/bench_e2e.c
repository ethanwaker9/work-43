#include "../src/f2ld.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static int eq(const int32_t *a, const int32_t *b, int k) { return memcmp(a, b, k * sizeof(int32_t)) == 0; }
static double medd(double *v, int n) {
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (v[j] < v[i]) { double t = v[i]; v[i] = v[j]; v[j] = t; }
    return v[n / 2];
}

int main(int argc, char **argv) {
    int p = argc > 1 ? atoi(argv[1]) : 2;
    int k = argc > 2 ? atoi(argv[2]) : 22;
    int64_t N = argc > 3 ? atoll(argv[3]) : 1048576;
    double contrast = argc > 4 ? atof(argv[4]) : 0.14;
    int hmax = argc > 5 ? atoi(argv[5]) : 16;
    int hstep = argc > 6 ? atoi(argv[6]) : 2;
    int trials = argc > 7 ? atoi(argv[7]) : 5;
    double snr = argc > 8 ? atof(argv[8]) : 6.0;
    int mat_kenum = argc > 9 ? atoi(argv[9]) : 6;
    uint64_t seed = argc > 10 ? strtoull(argv[10], NULL, 10) : 4242;

    printf("scenario,method,h,kprime,T,table_bytes,time_s,success\n");

    double *ts = malloc(sizeof(double) * trials);
    int64_t Tfull = 1; for (int i = 0; i < k; i++) Tfull *= p;

    int do_sparse = (Tfull <= 100000);
    const char *bm[3] = {"dense", "matzov", "sparse"};
    for (int b = 0; b < 3; b++) {
        if (b == 2 && !do_sparse) continue;
        int ok = 1; double tb = 0;
        for (int t = 0; t < trials; t++) {
            samples_t *s = samples_alloc(N, k, p);
            plant_instance(s, contrast, seed + 100 + t);
            result_t r;
            if (b == 0) r = distinguisher_dense_fft(s);
            else if (b == 1) r = distinguisher_matzov(s, mat_kenum);
            else r = distinguisher_sparse(s);
            ts[t] = r.runtime_s; tb = r.peak_bytes; ok &= eq(r.guess, s->secret, k);
            result_free(&r); samples_free(s);
        }
        printf("leaky,%s,0,%d,%lld,%.0f,%.6f,%d\n", bm[b], k, (long long)Tfull, tb, medd(ts, trials), ok);
    }

    for (int h = 0; h <= hmax; h += hstep) {
        int ok = 1; double tb = 0; int kp = k;
        int64_t Tp = Tfull;
        for (int t = 0; t < trials; t++) {
            samples_t *s = samples_alloc(N, k, p);
            plant_instance(s, contrast, seed + 200 + t);
            hintset_t *hs = hintset_new();
            if (h > 0) leak_hamming(s, hs, snr, h, seed + 300 + t);
            int64_t Te, Ne; samples_t *fref = fold_hints(s, hs, &Te, &Ne);
            kp = (int)llround(log((double)Te) / log((double)p));
            Tp = Te;
            result_t r = distinguisher_ours(s, hs);
            ts[t] = r.runtime_s; tb = r.peak_bytes;
            ok &= eq(r.guess, fref->secret, kp);
            result_free(&r); samples_free(s); samples_free(fref); hintset_free(hs);
        }
        printf("leaky,ours,%d,%d,%lld,%.0f,%.6f,%d\n", h, kp, (long long)Tp, tb, medd(ts, trials), ok);
    }
    free(ts);
    return 0;
}
