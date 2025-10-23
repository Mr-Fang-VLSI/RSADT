#!/bin/bash
# -------------------------------
# Benchmark DP (no Max-T, unweighted)
# Usage: bash bench_dp.sh > bench.log
# -------------------------------

exe=./dp_model
if [ ! -x "$exe" ]; then
  echo "? please compile dp_model first (make)"
  exit 1
fi

echo "=== DP Benchmark: from 4x4 to 15x15 ==="
printf "%5s %5s %12s %12s %12s\n" "m" "h" "HPWL" "Time(s)" "Peak(MB)"
echo "-----------------------------------------------------------"

for m in $(seq 4 15); do
  h=$m
  out=$($exe $m $h | grep -E "Best HPWL|Time =" | tr '\n' ' ')
  # parse
  hpwl=$(echo "$out" | grep -oE "Best HPWL = [0-9]+" | awk '{print $4}')
  time=$(echo "$out" | grep -oE "Time = [0-9.]+ s" | awk '{print $3}')
  mem=$(echo "$out" | grep -oE "PeakRSS = [0-9.]+ MB" | awk '{print $3}')
  printf "%5d %5d %12s %12s %12s\n" "$m" "$h" "$hpwl" "$time" "$mem"
done

echo "-----------------------------------------------------------"
echo "? Benchmark finished. Results above."
