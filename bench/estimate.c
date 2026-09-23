#include "../src/f2ld.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct { const char *name; int n; int q; double sigma; int beta; int klat; int p; } kyber_t;

static double delta_of_beta(int beta) {
    double b = beta;
    return pow((pow(M_PI * b, 1.0 / b) * b / (2.0 * M_PI * M_E)), 1.0 / (2.0 * (b - 1.0)));
}

int main(void) {
    kyber_t sets[3] = {
        {"Kyber-512", 512, 3329, 1.22, 406, 40, 2},
        {"Kyber-768", 768, 3329, 1.00, 623, 46, 2},
        {"Kyber-1024", 1024, 3329, 1.00, 873, 52, 2},
    };
    double m4_sram = 512.0 * 1024.0;
    double m0_sram = 32.0 * 1024.0;

    printf("set,n,q,beta,coreSVP_classical_bits,coreSVP_quantum_bits,delta\n");
    for (int i = 0; i < 3; i++) {
        kyber_t *K = &sets[i];
        double cls = 0.292 * K->beta + 16.4;
        double qnt = 0.265 * K->beta + 16.4;
        printf("%s,%d,%d,%d,%.1f,%.1f,%.5f\n", K->name, K->n, K->q, K->beta, cls, qnt, delta_of_beta(K->beta));
    }

    printf("\nset,h,klat_eff,mem_dense_bytes,mem_ours_bytes,guess_bits_removed,fits_M4,fits_M0\n");
    for (int i = 0; i < 3; i++) {
        kyber_t *K = &sets[i];
        for (int h = 0; h <= K->klat - 4; h += 4) {
            int ke = K->klat - h;
            double log2Td = ke * log2((double)K->p);
            double mem_dense = pow(2.0, log2Td) * (double)sizeof(double _Complex);
            double mem_ours = pow(2.0, log2Td) * 2.0 * (double)sizeof(int32_t);
            double removed = h * log2((double)K->p);
            int fm4 = mem_ours <= m4_sram;
            int fm0 = mem_ours <= m0_sram;
            printf("%s,%d,%d,%.3e,%.3e,%.1f,%d,%d\n", K->name, h, ke, mem_dense, mem_ours, removed, fm4, fm0);
        }
    }
    return 0;
}
