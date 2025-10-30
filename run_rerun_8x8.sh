#!/usr/bin/env bash
set -euo pipefail

# ======= 基本配置 =======
M=8; H=8
ROUNDS="${ROUNDS:-200}"
DV="${DV:-1}"
POW2="${POW2:-1}"
PROGRESS="${PROGRESS:-1}"
BIN="${BIN:-build/oc_singlecol}"

# ======= 惯性加权参数（可自由调节） =======
ALPHA="${ALPHA:-0.5}"        # 动量衰减系数 α
ETA="${ETA:-0.5}"            # 步长 η
TOP_RATIO="${TOP_RATIO:-0.15}"  # Top-p% 更新比例
CAP_MAX="${CAP_MAX:-8}"     # 权重上限 cap_max
CAP_MIN="${CAP_MIN:-0}"      # 权重下限 cap_min
SCALE="${SCALE:-1000}"       # 输出整数缩放比例

timestamp=$(date +%Y%m%d_%H%M%S)
RESULT_DIR="results_8x8_${timestamp}"
LOG_DIR="${RESULT_DIR}/logs"
SUMMARY="${RESULT_DIR}/summary_8x8.txt"
TARBALL="${RESULT_DIR}.tar.gz"

echo "[INFO] Creating new result folder: ${RESULT_DIR}"
mkdir -p "${LOG_DIR}"

# ======= 准备单位权重 =======
cat > "${RESULT_DIR}/wH_unit.txt" <<'EOF'
1 1 1 1 1 1 1
1 1 1 1 1 1 1
1 1 1 1 1 1 1
1 1 1 1 1 1 1
1 1 1 1 1 1 1
1 1 1 1 1 1 1
1 1 1 1 1 1 1
1 1 1 1 1 1 1
EOF

cat > "${RESULT_DIR}/wV_unit.txt" <<'EOF'
1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1
EOF

tr -d '\r' < "${RESULT_DIR}/wH_unit.txt" > "${RESULT_DIR}/wH.tmp" && mv "${RESULT_DIR}/wH.tmp" "${RESULT_DIR}/wH_unit.txt"
tr -d '\r' < "${RESULT_DIR}/wV_unit.txt" > "${RESULT_DIR}/wV.tmp" && mv "${RESULT_DIR}/wV.tmp" "${RESULT_DIR}/wV_unit.txt"

# ======= summary 表头 =======
echo "T  HPWL_eq(last)  maxL(last)  Meet_MaxT" > "${SUMMARY}"

# ======= T 从 8 到 20 =======
for T in $(seq 8 20); do
  echo "========== Run T=${T} (m=${M},h=${H},rounds=${ROUNDS}) =========="
  LOG="${LOG_DIR}/singlecol_8x8_T_${T}.log"
  OUT_NAME="placement${M}_Col_1_T_${T}.txt"

  # === 执行 oc_singlecol ===
  ( set -x
    time "${BIN}" "${M}" "${H}" \
      --rounds="${ROUNDS}" --dV="${DV}" --T="${T}" \
      --pow2="${POW2}" --progress="${PROGRESS}" \
      --alpha="${ALPHA}" --eta="${ETA}" --top="${TOP_RATIO}" \
      --cap="${CAP_MAX}" --scale="${SCALE}" \
      --wH "${RESULT_DIR}/wH_unit.txt" --wV "${RESULT_DIR}/wV_unit.txt"
  ) 2>&1 | tee "${LOG}"

  # === 移动 placement 文件 ===
  if [[ -f "${OUT_NAME}" ]]; then
    mv "${OUT_NAME}" "${RESULT_DIR}/"
    echo "[OK] Moved ${OUT_NAME} -> ${RESULT_DIR}/"
  else
    echo "[WARN] Expected output not found for T=${T} (see ${LOG})"
  fi

  # === 解析日志最后一轮数据 ===
  HPWL=$(grep "HPWL_eq=" "${LOG}" | tail -n1 | sed -E 's/.*HPWL_eq=([0-9]+).*/\1/')
  MAXL=$(grep "maxL="    "${LOG}" | tail -n1 | sed -E 's/.*maxL=([0-9]+).*/\1/')
  [[ -z "${HPWL}" ]] && HPWL="NA"
  [[ -z "${MAXL}" ]] && MAXL="NA"
  MEET="NA"
  if [[ "${MAXL}" != "NA" ]]; then
    (( MAXL <= T )) && MEET="YES" || MEET="NO"
  fi

  printf "%-3s %-14s %-10s %-3s\n" "${T}" "${HPWL}" "${MAXL}" "${MEET}" >> "${SUMMARY}"
done

# ======= 压缩结果目录 =======
echo "[INFO] Compressing ${RESULT_DIR} -> ${TARBALL}"
tar -czf "${TARBALL}" "${RESULT_DIR}"
echo "[DONE] All runs complete. Archive: ${TARBALL}"
