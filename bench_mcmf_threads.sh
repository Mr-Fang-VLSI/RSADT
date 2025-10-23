#!/bin/bash
# --------------------------------------------------------
# Weighted-MCMF (OC-shortest) ???????
# ???: 1 2 4 8 24,??: 4x4 ~ 15x15
# ??: ??????,???????????RSS
# --------------------------------------------------------

exe=./build/oc_shortest
if [ ! -x "$exe" ]; then
  echo "? ???? oc_shortest (make)"
  exit 1
fi

THREADS=(1 2 4 8 24)
SIZES=$(seq 4 15)

echo "=== Weighted-MCMF ????? (real RSS) ==="
printf "%8s" "Thread"
for m in $SIZES; do printf "%10s" "${m}x${m}"; done
echo
echo "---------------------------------------------------------------"

for thr in "${THREADS[@]}"; do
  printf "%8d" "$thr"
  for m in $SIZES; do
    h=$m
    t0=$(date +%s%N)
    /usr/bin/time -v $exe $m $h 1 0 $thr > /tmp/mcmf_out.log 2> /tmp/mcmf_time.log
    t1=$(date +%s%N)
    sec=$(printf "%.3f" "$(echo "scale=6;($t1-$t0)/1000000000" | bc -l)")
    mem_kb=$(grep "Maximum resident set size" /tmp/mcmf_time.log | awk '{print $(NF)}')
    mem_mb=$(echo "scale=1; $mem_kb/1024" | bc)
    printf "%10s" "${sec}s"
  done
  echo
done

echo "---------------------------------------------------------------"
echo "? ???????? wall time,????"
echo "   ???????? RSS ??????"
