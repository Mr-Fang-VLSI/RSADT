#!/usr/bin/env bash
set -euo pipefail

BIN="${BIN:-build/oc_observe}"
ROUNDS="${ROUNDS:-200}"
ALPHA="${ALPHA:-0.5}"
ETA="${ETA:-0.5}"
TOP="${TOP:-0.10}"
CAP="${CAP:-32}"
POW2="${POW2:-1}"
PROG="${PROG:-1}"

pairs=("8 8" "12 8" "16 8")
for p in "${pairs[@]}"; do
  set -- $p; m=$1; h=$2
  T=$(python3 - <<PY
m=${m}; h=${h}
print(int(round(1.5*min(m,h))))
PY
)
  OUT="metrics_m${m}_h${h}_T${T}.txt"
  echo "== Run m=${m} h=${h} T=${T} rounds=${ROUNDS} =="

  ( set -x; "${BIN}" ${m} ${h} --rounds="${ROUNDS}" --T="${T}" \
      --alpha="${ALPHA}" --eta="${ETA}" --top="${TOP}" --cap="${CAP}" \
      --pow2="${POW2}" --progress="${PROG}" --outfile "${OUT}" )
done

echo "[DONE] metrics written to metrics_m*_h*_T*.txt"
