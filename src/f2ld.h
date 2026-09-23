#ifndef F2LD_H
#define F2LD_H

#include <stdint.h>
#include <stddef.h>
#include <complex.h>

typedef struct {
    int k;
    int p;
    int logp;
    int64_t T;
    int n;
    int q;
    double sigma_s;
    double sigma_e;
} params_t;

typedef struct {
    int32_t *addr;
    double _Complex *w;
    int64_t N;
    int k;
    int p;
    int32_t *secret;
} samples_t;

typedef struct {
    int coord;
    int modulus;
    int value;
    int kind;
    double *prior;
} hint_t;

typedef struct {
    hint_t *h;
    int count;
    int cap;
} hintset_t;

enum { HINT_PERFECT = 0, HINT_MODULAR = 1, HINT_APPROX = 2 };

typedef struct {
    int found;
    int32_t *guess;
    double peak;
    double runtime_s;
    double peak_bytes;
    int64_t Teff;
    int64_t Neff;
} result_t;

void fft_radix2(double _Complex *a, int64_t n, int inv);
void ndft_double(double _Complex *a, int p, int k, int inv);
void ndft_fixed(int32_t *re, int32_t *im, int p, int k, int *shift_out);

samples_t *samples_alloc(int64_t N, int k, int p);
void samples_free(samples_t *s);
void plant_instance(samples_t *s, double contrast, uint64_t seed);
double instance_contrast(const params_t *pp, double blocksize_delta);

hintset_t *hintset_new(void);
void hintset_free(hintset_t *hs);
void leak_hamming(const samples_t *s, hintset_t *hs, double snr, int nhints, uint64_t seed);
samples_t *fold_hints(const samples_t *in, const hintset_t *hs, int64_t *Teff_out, int64_t *Neff_out);

result_t distinguisher_enum(const samples_t *s);
result_t distinguisher_dense_fft(const samples_t *s);
result_t distinguisher_matzov(const samples_t *s, int k_enum);
result_t distinguisher_sparse(const samples_t *s);
result_t distinguisher_ours(const samples_t *s, const hintset_t *hs);
result_t distinguisher_ours_split(const samples_t *s, const hintset_t *hs, int k_enum);
result_t distinguisher_prior(const samples_t *s, const hintset_t *hs, double lambda);
void result_free(result_t *r);

double now_seconds(void);
double peak_rss_mb(void);
uint64_t xorshift64(uint64_t *st);
double gaussrand(uint64_t *st);

#endif
