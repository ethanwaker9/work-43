#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/.."
make all >/dev/null 2>&1
mkdir -p results

echo "[1/6] kernel micro-benchmark (5 methods vs T)"
HDR="method,p,k,N,T,contrast,time_s,table_bytes,rss_mb,ok"
echo "$HDR" > results/kernel.csv
N=262144; CON=0.15; TR=3
for k in 8 10 12 14 16 18 20 22 24; do
  ./bench/bench_kernel dense  2 $k $N $CON 0 $TR >> results/kernel.csv
  ./bench/bench_kernel ours   2 $k $N $CON 0 $TR >> results/kernel.csv
  ke=$(( k/3 )); [ $ke -lt 1 ] && ke=1
  ./bench/bench_kernel matzov 2 $k $N $CON $ke $TR >> results/kernel.csv
  if [ $k -le 12 ]; then ./bench/bench_kernel enum 2 $k $N $CON 0 $TR >> results/kernel.csv; fi
  if [ $k -le 14 ]; then ./bench/bench_kernel sparse 2 $k $N $CON 0 $TR >> results/kernel.csv; fi
done

echo "[2/6] end-to-end hint sweep (embedded scenario)"
./bench/bench_e2e 2 22 1048576 0.14 18 2 5 6.0 6 4242 > results/e2e.csv

echo "[3/6] success vs samples with/without Bayesian prior"
./bench/bench_success 2 16 0.05 300 16 1.7 9.0 20260922 1200 40000 14 > results/success.csv

echo "[4/6] security / memory estimator for Kyber"
./bench/estimate > results/estimate.csv

echo "[5/6] Cortex-M4 footprint"
echo "kbits,text_bytes,bss_bytes,table_bytes,wht_ops" > results/arm.csv
for kb in 8 10 12 14 16; do
  arm-none-eabi-gcc -mcpu=cortex-m4 -mthumb -mfloat-abi=soft -O3 -ffreestanding -c -DKBITS=$kb -o arm/k_$kb.o arm/kernel_m4.c
  sz=$(arm-none-eabi-size arm/k_$kb.o | tail -1)
  txt=$(echo $sz | awk '{print $1}'); bss=$(echo $sz | awk '{print $3}')
  echo "$kb,$txt,$bss,$((8*(1<<kb))),$(( (1<<kb)*kb ))" >> results/arm.csv
done
rm -f arm/k_*.o

echo "[6/6] real fplll-reduced dual attack validation"
echo "seed,samples,measured_bias,success" > results/validate.csv
for sd in $(seq 1 24); do
  out=$(./bench/validate_real 20 97 46 20 $sd || true)
  smp=$(echo "$out" | grep -oE "samples=[0-9]+" | cut -d= -f2)
  b=$(echo "$out" | grep -oE "measured_bias=[0-9.]+" | cut -d= -f2)
  s=$(echo "$out" | grep -oE "success=[01]" | cut -d= -f2)
  echo "$sd,$smp,$b,$s" >> results/validate.csv
done

echo "done. CSVs in results/"
