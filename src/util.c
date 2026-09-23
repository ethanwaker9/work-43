#include "f2ld.h"
#include <time.h>
#include <sys/resource.h>
#include <math.h>

double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

double peak_rss_mb(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);
}

uint64_t xorshift64(uint64_t *st) {
    uint64_t x = *st;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *st = x;
    return x;
}

double gaussrand(uint64_t *st) {
    double u1 = (double)(xorshift64(st) >> 11) / 9007199254740992.0;
    double u2 = (double)(xorshift64(st) >> 11) / 9007199254740992.0;
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}
