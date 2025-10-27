#!/bin/bash
# ==========================================
#  OC-Shortest benchmark (4x4 ~ 16x16)
#  Threads: <=12 ? 4 threads, >12 ? 8 threads
#  Record HPWL, pure runtime, max RAM (MB)
# ==========================================

BIN=./build/oc_shortest
WH=wH_unit.txt
WV=wV_unit.txt
OUT=bench_result.txt

# ???????
echo "# m  h  threads  HPWL  runtime_ms  max_RAM_MB" > $OUT

# ???1????(??16x16??)
yes 1 | head -n $((16*15)) > $WH
yes 1 | head -n $((15*16)) > $WV

# ??????
for size in 4 5 6 7 8 9 10 11 12 13 14 15 16; do
  m=$size
  h=$size
  if [ $size -le 12 ]; then
    thr=24
  else
    thr=24
  fi

  echo "? Running ${m}x${h} (threads=$thr)..."

  # ?? /usr/bin/time ???????
  # “Solve-only time”??????,??????
  /usr/bin/time -v $BIN $m $h 1 -1 $thr \
      --pow2=1 --check-oc=0 --verify=1 \
      --rsad=0 --wH $WH --wV $WV > tmp_${m}x${h}.log 2> tmp_${m}x${h}_time.log

  # ?? HPWL ? runtime
  hpwl=$(grep "\[Verify\]" tmp_${m}x${h}.log | awk '{print $2}' | cut -d'=' -f2)
  runtime=$(grep "Solve-only time" tmp_${m}x${h}.log | awk '{print $3}')
  ram=$(grep "Maximum resident set size" tmp_${m}x${h}_time.log | awk '{print $6}')

  echo "$m $h $thr $hpwl $runtime $ram" >> $OUT
done

echo -e "\n? Benchmark done. Results saved to $OUT"
column -t $OUT
