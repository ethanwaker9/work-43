#include "../src/f2ld.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static int eq(const int32_t *a, const int32_t *b, int k) { return memcmp(a, b, k * sizeof(int32_t)) == 0; }

static double median(double *v, int n) {
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (v[j] < v[i]) { double t = v[i]; v[i] = v[j]; v[j] = t; }
    return v[n / 2];
}

int main(int argc, char **argv) {
    if (argc < 8) {
        fprintf(stderr, "usage: %s method p k N contrast k_enum trials [seed]\n", argv[0]);
        fprintf(stderr, "method in {enum,dense,matzov,sparse,ours}\n");
        return 2;
    }
    const char *method = argv[1];
    int p = atoi(argv[2]), k = atoi(argv[3]);
    int64_t N = atoll(argv[4]);
    double contrast = atof(argv[5]);
    int k_enum = atoi(argv[6]);
    int trials = atoi(argv[7]);
    uint64_t seed = argc > 8 ? strtoull(argv[8], NULL, 10) : 777;

    double *ts = malloc(sizeof(double) * trials);
    double table_bytes = 0; int ok_all = 1;
    for (int t = 0; t < trials; t++) {
        samples_t *s = samples_alloc(N, k, p);
        plant_instance(s, contrast, seed + t);
        result_t r;
        if (!strcmp(method, "enum")) r = distinguisher_enum(s);
        else if (!strcmp(method, "dense")) r = distinguisher_dense_fft(s);
        else if (!strcmp(method, "matzov")) r = distinguisher_matzov(s, k_enum);
        else if (!strcmp(method, "sparse")) r = distinguisher_sparse(s);
        else if (!strcmp(method, "ours")) r = distinguisher_ours(s, NULL);
        else if (!strcmp(method, "ourss")) r = distinguisher_ours_split(s, NULL, k_enum);
        else { fprintf(stderr, "bad method\n"); return 2; }
        ts[t] = r.runtime_s;
        table_bytes = r.peak_bytes;
        int ok = eq(r.guess, s->secret, k);
        ok_all &= ok;
        result_free(&r);
        samples_free(s);
    }
    double tmed = median(ts, trials);
    double rss = peak_rss_mb();
    int64_t Tfull = 1; for (int i = 0; i < k; i++) Tfull *= p;
    printf("%s,%d,%d,%lld,%lld,%.4f,%.8f,%.0f,%.3f,%d\n",
           method, p, k, (long long)N, (long long)Tfull, contrast, tmed, table_bytes, rss, ok_all);
    free(ts);
    return 0;
}
