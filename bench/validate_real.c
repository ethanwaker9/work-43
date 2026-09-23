#include "../src/f2ld.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int q, n, m, n2, k1;
static int A[64][64];
static int secret[64];
static int errv[64];
static int bvec[64];

static int rq(uint64_t *st, int mod) { return (int)(xorshift64(st) % (uint64_t)mod); }
static int small(uint64_t *st) { return (int)(xorshift64(st) % 3) - 1; }

int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 20;
    q = argc > 2 ? atoi(argv[2]) : 97;
    m = argc > 3 ? atoi(argv[3]) : 46;
    int bs = argc > 4 ? atoi(argv[4]) : 20;
    uint64_t seed = argc > 5 ? strtoull(argv[5], NULL, 10) : 12345;
    k1 = 1;
    n2 = n - k1;
    uint64_t st = seed;

    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) A[i][j] = rq(&st, q);
    for (int j = 0; j < n; j++) secret[j] = small(&st);
    for (int i = 0; i < m; i++) {
        long acc = 0;
        for (int j = 0; j < n; j++) acc += (long)A[i][j] * secret[j];
        errv[i] = small(&st);
        bvec[i] = (int)(((acc + errv[i]) % q + q) % q);
    }

    long C = (long)q * 4096;
    int dim = m + n2;
    char path[512];
    snprintf(path, sizeof(path), "/tmp/f2ld_lat_%llu.txt", (unsigned long long)seed);
    FILE *f = fopen(path, "w");
    fprintf(f, "[");
    for (int i = 0; i < m; i++) {
        fprintf(f, "[");
        for (int col = 0; col < m; col++) fprintf(f, "%d ", col == i ? 1 : 0);
        for (int t = 0; t < n2; t++) fprintf(f, "%ld ", C * A[i][k1 + t]);
        fprintf(f, "]\n");
    }
    for (int t = 0; t < n2; t++) {
        fprintf(f, "[");
        for (int col = 0; col < m; col++) fprintf(f, "0 ");
        for (int u = 0; u < n2; u++) fprintf(f, "%ld ", u == t ? C * q : 0);
        fprintf(f, "]\n");
    }
    fprintf(f, "]\n");
    fclose(f);

    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "fplll -a bkz -b %d < %s", bs, path);
    FILE *pp = popen(cmd, "r");
    if (!pp) { fprintf(stderr, "fplll failed\n"); return 1; }

    int cap = 4096;
    samples_t *S = samples_alloc(cap, k1, q);
    int64_t cnt = 0;
    double biassum = 0; int biasn = 0;
    int row[256];
    int ridx = 0, val, sign, reading = 0;
    int c;
    long *rowbuf = malloc(sizeof(long) * dim);
    (void)row;
    char line[65536];
    while (fgets(line, sizeof(line), pp)) {
        char *pt = line;
        while (*pt == '[' || *pt == ' ') pt++;
        int nc = 0; long vals[512];
        char *tok = strtok(pt, " ]\n[");
        while (tok && nc < dim) { vals[nc++] = atol(tok); tok = strtok(NULL, " ]\n["); }
        if (nc < dim) continue;
        int ok = 1;
        for (int t = 0; t < n2; t++) if (vals[m + t] != 0) { ok = 0; break; }
        if (!ok) continue;
        long norm = 0; for (int col = 0; col < m; col++) norm += vals[col] * vals[col];
        if (norm == 0) continue;
        long y1 = 0, cc = 0;
        for (int i = 0; i < m; i++) {
            y1 += vals[i] * A[i][0];
            cc += vals[i] * bvec[i];
        }
        y1 = ((y1 % q) + q) % q;
        cc = ((cc % q) + q) % q;
        long ers = 0; for (int i = 0; i < m; i++) ers += vals[i] * errv[i];
        biassum += cos(2.0 * M_PI * (double)(((ers % q) + q) % q) / q); biasn++;
        if (cnt < cap) {
            S->addr[cnt] = (int32_t)y1;
            double ang = 2.0 * M_PI * (double)cc / q;
            S->w[cnt] = cos(ang) + I * sin(ang);
            cnt++;
        }
        (void)ridx; (void)val; (void)sign; (void)reading; (void)c; (void)rowbuf;
    }
    pclose(pp);
    S->N = cnt;
    S->secret[0] = secret[0];

    double best = -1e300; int arg = 0;
    for (int s1 = 0; s1 < q; s1++) {
        double acc = 0;
        for (int64_t j = 0; j < cnt; j++) {
            int mm = (int)((((long)S->addr[j] * s1) % q + q) % q);
            double ang = -2.0 * M_PI * mm / q;
            acc += creal(S->w[j]) * cos(ang) - cimag(S->w[j]) * sin(ang);
        }
        if (acc > best) { best = acc; arg = s1; }
    }
    double meas_bias = biasn ? biassum / biasn : 0;
    int true_s = ((secret[0] % q) + q) % q;
    printf("n=%d q=%d m=%d bs=%d samples=%lld measured_bias=%.4f recovered=%d true=%d success=%d\n",
           n, q, m, bs, (long long)cnt, meas_bias, arg, true_s, arg == true_s);
    samples_free(S);
    free(rowbuf);
    remove(path);
    return arg == true_s ? 0 : 2;
}
