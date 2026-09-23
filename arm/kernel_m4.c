#include <stdint.h>

#ifndef KBITS
#define KBITS 14
#endif
#define T (1u << KBITS)

static int32_t re[T];
static int32_t im[T];

void f2ld_accumulate(const uint32_t *addr, const int32_t *wre, const int32_t *wim, uint32_t N) {
    for (uint32_t j = 0; j < N; j++) {
        uint32_t c = addr[j] & (T - 1u);
        re[c] += wre[j];
        im[c] += wim[j];
    }
}

int32_t f2ld_wht_argmax(void) {
    for (uint32_t len = 1; len < T; len <<= 1) {
        for (uint32_t i = 0; i < T; i += (len << 1)) {
            for (uint32_t k = 0; k < len; k++) {
                int32_t ur = re[i + k], ui = im[i + k];
                int32_t vr = re[i + k + len], vi = im[i + k + len];
                re[i + k] = ur + vr; im[i + k] = ui + vi;
                re[i + k + len] = ur - vr; im[i + k + len] = ui - vi;
            }
        }
        int32_t mx = 0;
        for (uint32_t t = 0; t < T; t++) { int32_t a = re[t] < 0 ? -re[t] : re[t]; if (a > mx) mx = a; }
        if (mx > (1 << 29)) for (uint32_t t = 0; t < T; t++) { re[t] >>= 1; im[t] >>= 1; }
    }
    int32_t best = re[0]; int32_t arg = 0;
    for (uint32_t t = 1; t < T; t++) if (re[t] > best) { best = re[t]; arg = (int32_t)t; }
    return arg;
}

int main(void) {
    for (uint32_t t = 0; t < T; t++) { re[t] = (int32_t)((t * 2654435761u) >> 20); im[t] = 0; }
    volatile int32_t s = f2ld_wht_argmax();
    return (int)(s & 0x7f);
}
