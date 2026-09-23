#include "../src/f2ld.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static int eq(const int32_t *a, const int32_t *b, int k) { return memcmp(a, b, k * sizeof(int32_t)) == 0; }

int main(int argc, char **argv) {
    int p = argc > 1 ? atoi(argv[1]) : 2;
    int k = argc > 2 ? atoi(argv[2]) : 14;
    double contrast = argc > 3 ? atof(argv[3]) : 0.10;
    int trials = argc > 4 ? atoi(argv[4]) : 200;
    int nhints = argc > 5 ? atoi(argv[5]) : 14;
    double snr = argc > 6 ? atof(argv[6]) : 1.2;
    double lambda = argc > 7 ? atof(argv[7]) : 3.0;
    uint64_t seed = argc > 8 ? strtoull(argv[8], NULL, 10) : 20260922;
    int64_t Nmin = argc > 9 ? atoll(argv[9]) : 2000;
    int64_t Nmax = argc > 10 ? atoll(argv[10]) : 200000;
    int points = argc > 11 ? atoi(argv[11]) : 12;

    printf("N,sr_plain,sr_prior\n");
    for (int pt = 0; pt < points; pt++) {
        double frac = points > 1 ? (double)pt / (points - 1) : 0.0;
        int64_t N = (int64_t)llround(Nmin * pow((double)Nmax / Nmin, frac));
        int okp = 0, okr = 0;
        for (int t = 0; t < trials; t++) {
            samples_t *s = samples_alloc(N, k, p);
            plant_instance(s, contrast, seed + (uint64_t)pt * 100003 + t);
            hintset_t *hs = hintset_new();
            leak_hamming(s, hs, snr, nhints, seed + (uint64_t)pt * 777 + t + 1);
            result_t rp = distinguisher_dense_fft(s);
            result_t rr = distinguisher_prior(s, hs, lambda);
            okp += eq(rp.guess, s->secret, k);
            okr += eq(rr.guess, s->secret, k);
            result_free(&rp); result_free(&rr);
            samples_free(s); hintset_free(hs);
        }
        printf("%lld,%.4f,%.4f\n", (long long)N, (double)okp / trials, (double)okr / trials);
        fflush(stdout);
    }
    return 0;
}
