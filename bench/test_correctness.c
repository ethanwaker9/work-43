#include "../src/f2ld.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

static int eq(const int32_t *a, const int32_t *b, int k) { return memcmp(a, b, k * sizeof(int32_t)) == 0; }

int main(void) {
    int fails = 0;
    int p = 4, k = 5;
    int64_t N = 40000;
    double contrast = 0.25;
    samples_t *s = samples_alloc(N, k, p);
    plant_instance(s, contrast, 12345);

    result_t re = distinguisher_enum(s);
    result_t rf = distinguisher_dense_fft(s);
    result_t rm = distinguisher_matzov(s, 2);
    result_t rs = distinguisher_sparse(s);
    result_t ro = distinguisher_ours(s, NULL);

    printf("secret:"); for (int d = 0; d < k; d++) printf(" %d", s->secret[d]); printf("\n");
    printf("enum   :"); for (int d = 0; d < k; d++) printf(" %d", re.guess[d]); printf("  ok=%d t=%.4f\n", eq(re.guess, s->secret, k), re.runtime_s);
    printf("dense  :"); for (int d = 0; d < k; d++) printf(" %d", rf.guess[d]); printf("  ok=%d t=%.4f\n", eq(rf.guess, s->secret, k), rf.runtime_s);
    printf("matzov :"); for (int d = 0; d < k; d++) printf(" %d", rm.guess[d]); printf("  ok=%d t=%.4f\n", eq(rm.guess, s->secret, k), rm.runtime_s);
    printf("sparse :"); for (int d = 0; d < k; d++) printf(" %d", rs.guess[d]); printf("  ok=%d t=%.4f (m=%lld)\n", eq(rs.guess, s->secret, k), rs.runtime_s, (long long)rs.Neff);
    printf("ours   :"); for (int d = 0; d < k; d++) printf(" %d", ro.guess[d]); printf("  ok=%d t=%.4f\n", eq(ro.guess, s->secret, k), ro.runtime_s);

    fails += !eq(re.guess, s->secret, k);
    fails += !eq(rf.guess, s->secret, k);
    fails += !eq(rm.guess, s->secret, k);
    fails += !eq(rs.guess, s->secret, k);
    fails += !eq(ro.guess, s->secret, k);

    printf("\n[hint-folding test] pin 2 coords via high-SNR approx hints\n");
    hintset_t *hs = hintset_new();
    leak_hamming(s, hs, 6.0, 2, 999);
    int64_t Te, Ne; samples_t *fref = fold_hints(s, hs, &Te, &Ne);
    int kp = (int)llround(log((double)Te) / log((double)p));
    result_t rh = distinguisher_ours(s, hs);
    int ok = memcmp(rh.guess, fref->secret, kp * sizeof(int32_t)) == 0;
    printf("folded k'=%d  T'=%lld (full T=%lld)  ok=%d t=%.5f bytes=%.0f\n",
           kp, (long long)Te, (long long)llround(pow(p, k)), ok, rh.runtime_s, rh.peak_bytes);
    fails += !ok;

    printf("\nRESULT: %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    result_free(&re); result_free(&rf); result_free(&rm); result_free(&rs); result_free(&ro); result_free(&rh);
    samples_free(s); samples_free(fref); hintset_free(hs);
    return fails ? 1 : 0;
}
