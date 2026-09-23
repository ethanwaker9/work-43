# F2LD: Hint-Aware Dual Lattice Attacks on Embedded Systems

This repository contains the implementation and benchmarks for efficient side-channel dual lattice key recovery on embedded systems. Our code implements the F2LD distinguisher together with the four existing dual-attack distinguishers it is compared against, and
validates end-to-end key recovery on real reduced lattices.

## Files and Contents

```
src/                core library (C11)
  f2ld.h            public types and API
  fft.c             double and fixed point k-dimensional transforms over Z_p
  instance.c        LWE sample generation, planting, contrast model
  hints.c           side-channel leakage simulation and hint folding
  distinguishers.c  enumeration, dense FFT, MATZOV hybrid, sparse, F2LD
  util.c            timing, peak memory, RNG
bench/
  test_correctness.c  correctness self test for all distinguishers
  bench_kernel.c      each method time and space microbenchmark
  bench_e2e.c         hint-assisted key recovery sweep (embedded scenario)
  bench_success.c     success rate vs samples, with and without the prior
  estimate.c          security and memory estimator for Kyber parameter sets
  validate_real.c     E2E dual attack on fplll-reduced lattices
arm/
  kernel_m4.c         Cortex-M4 fixed point kernel for footprint measurement
scripts/
  run_all.sh          regenerates every CSV in results/
results/              generated CSV data
```

## Requirements

- A C11 compiler (`cc`/`clang`/`gcc`).
- `python3` with `numpy` and `matplotlib` for the plots.
- `fplll` on the `PATH` for `validate_real` (`brew install fplll` or
  `apt install fplll`).
- `arm-none-eabi-gcc` for the Cortex-M4 footprint (optional).

## Build
```
make
```
This builds all benchmarks into `bench/`. To check correctness:
```
./bench/test_correctness
```

## Running the Experiments
```
bash scripts/run_all.sh        # writes results/*.csv
python3 scripts/plots.py       # writes ../final_paper/figs/*.eps
```
`run_all.sh` runs the kernel microbenchmark, the hint sweep, the
success versus samples experiment, the Kyber estimator, the Cortex-M4 footprint,
and the real lattice validation. 
```
./bench/bench_kernel <method> <p> <k> <N> <contrast> <k_enum> <trials>
    method in {enum, dense, matzov, sparse, ours, ourss}

./bench/bench_e2e <p> <k> <N> <contrast> <hmax> <hstep> <trials> <snr> <k_enum> <seed>

./bench/bench_success <p> <k> <contrast> <trials> <nhints> <snr> <lambda> <seed> <Nmin> <Nmax> <points>

./bench/estimate

./bench/validate_real <n> <q> <m> <blocksize> <seed>
```
Each `bench_kernel` invocation prints one CSV line
`method,p,k,N,T,contrast,time_s,table_bytes,rss_mb,ok`. Running one method per
process keeps the reported peak memory clean.

## Cortex-M4 footprint
```
arm-none-eabi-gcc -mcpu=cortex-m4 -mthumb -mfloat-abi=soft -O3 \
    -ffreestanding -c -DKBITS=14 arm/kernel_m4.c -o arm/k.o
arm-none-eabi-size arm/k.o
```
The `p = 2` transform is a Walsh Hadamard transform, so the kernel uses only
additions and subtractions; the code stays near 424 bytes and the table sets
the static memory requirement.
