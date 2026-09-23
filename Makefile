CC ?= cc
CFLAGS ?= -O3 -march=native -std=c11 -Wall -Wextra -Isrc
LDFLAGS ?= -lm

SRC = src/fft.c src/util.c src/instance.c src/hints.c src/distinguishers.c
HDR = src/f2ld.h

BIN = bench/test_correctness bench/bench_kernel bench/bench_e2e bench/bench_success bench/estimate bench/validate_real

all: $(BIN)

bench/%: bench/%.c $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ $< $(SRC) $(LDFLAGS)

clean:
	rm -f $(BIN)

.PHONY: all clean
