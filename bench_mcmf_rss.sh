#!/bin/bash
# --------------------------------------------------------
# Weighted-MCMF Benchmark (real RSS) ???
# --------------------------------------------------------

exe=./build/oc_shortest
if [ ! -x "$exe" ]; then
  echo "? ???? oc_shortest (make)"
  exit 1
fi

echo "=== Weighted-MCMF Benchmark (real RSS) 4x4~15x15 ==="
printf "%5s %5s %10s %12s %12s\n" "m" "h" "HPWL" "Time(s)" "Peak(MB)"
echo "-----------------------------------------------------------"

for m in $(seq 4 15); do
  h=$m
  /usr/bin/time -v $exe $m $h > /tmp/mcmf_out.log 2> /tmp/mcmf_time.log

  # ?? HPWL ???
  hpwl=$(grep -Eo "HPWL=[0-9]+" /tmp/mcmf_out.log | sed 's/HPWL=//' | head -1)

  # ?? wall time(?)
  t0=$(date +%s%N)
  /usr/bin/time -v $exe $m $h > /tmp/mcmf_out.log 2> /tmp/mcmf_time.log
  t1=$(date +%s%N)
  sec=$(printf "%.3f" "$(echo "scale=6;($t1-$t0)/1000000000" | bc -l)")

  # ???? RSS
  mem_kb=$(grep "Maximum resident set size" /tmp/mcmf_time.log | awk '{print $(NF)}')
  mem_mb=$(echo "scale=2;$mem_kb/1024" | bc)

  printf "%5d %5d %10s %12s %12s\n" "$m" "$h" "$hpwl" "$sec" "$mem_mb"
done

echo "-----------------------------------------------------------"
echo "? Weighted-MCMF Benchmark finished."
