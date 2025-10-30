#!/usr/bin/env bash
set -euo pipefail

# ======================================================
#                基本配置参数（可被覆盖）
# ======================================================
M=16; H=8; T=12
ROUNDS="${ROUNDS:-120}"
BIN="${BIN:-build/oc_policy_singlecol}"
FREEZE="${FREEZE:-1}"      # 可行即冻结
CAP="${CAP:-64}"
SCALE="${SCALE:-1000}"
DV="${DV:-1}"
POW2="${POW2:-1}"
PROGRESS="${PROGRESS:-1}"

# --- Lagrange 参数 ---
RHO_LAG="${RHO_LAG:-0.3}"

# --- IRL1 参数 ---
LAM_IRL="${LAM_IRL:-2.0}"
DELTA_IRL="${DELTA_IRL:-2.0}"

# --- Soft 参数 (DreamPlace-style + 条件收紧 + 冷启动) ---
ALPHA_SOFT="${ALPHA_SOFT:-0.9}"
ETA_SOFT="${ETA_SOFT:-1.0}"
GAM_SOFT="${GAM_SOFT:-0.3}"
BET_SOFT="${BET_SOFT:-3.0}"
SIGMA_SOFT="${SIGMA_SOFT:-0.8}"
SHRINK_SOFT="${SHRINK_SOFT:-0.9}"

# ======================================================
#              结果目录与输出文件
# ======================================================
timestamp=$(date +%Y%m%d_%H%M%S)
RESULT_DIR="results_16x8_T12_cmp_${timestamp}"
LOG_DIR="${RESULT_DIR}/logs"
SUMMARY="${RESULT_DIR}/compare_16x8_T12.txt"
TARBALL="${RESULT_DIR}.tar.gz"

mkdir -p "${LOG_DIR}"
echo "[INFO] 创建结果目录: ${RESULT_DIR}"

# ======================================================
#        输出表头 (maxL 最大/最小/轮次 + HPWL_eq/W)
# ======================================================
echo "policy,maxL_max,round_of_maxL_max,maxL_min,round_of_maxL_min,HPWL_eq_final,HPWL_w_final" > "${SUMMARY}"

# ======================================================
#        函数：从日志中解析指标
# ======================================================
parse_log() {
  local log="$1"
  local Tval="$2"
  awk -v T="${Tval}" '
  BEGIN {
    nr=0; maxLmax=-1e9; maxLmin=1e9;
  }
  # 匹配 Round 行，兼容 busybox/mawk
  /\[Round[[:space:]]+[0-9]+\]/ {
    # round 号
    r=-1;
    if ($0 ~ /\[Round[[:space:]]+[0-9]+\]/) {
      sub(/^.*\[Round[[:space:]]+/, "", $0);
      r=$1+0;
    }
    # WNS
    w=0;
    if (index($0, "WNS=") > 0) {
      split($0, arr, "WNS=");
      split(arr[2], tmp, " ");
      w=tmp[1]+0;
    }
    # HPWL_eq
    he=0;
    if (index($0, "HPWL_eq=") > 0) {
      split($0, arr, "HPWL_eq=");
      split(arr[2], tmp, " ");
      he=tmp[1]+0;
    }
    # HPWL_w
    hw=0;
    if (index($0, "HPWL_w=") > 0) {
      split($0, arr, "HPWL_w=");
      split(arr[2], tmp, " ");
      hw=tmp[1]+0;
    }
    if (r >= 0) {
      maxL = T - w;
      if (nr == 0 || maxL > maxLmax) { maxLmax = maxL; rmax = r; }
      if (nr == 0 || maxL < maxLmin) { maxLmin = maxL; rmin = r; }
      nr++;
      last_he = he;
      last_hw = hw;
    }
  }
  END {
    if (nr > 0)
      printf "%d %d %d %d %d %d\n", maxLmax, rmax, maxLmin, rmin, last_he, last_hw;
    else exit 1;
  }' "${log}"
}


# ======================================================
#        函数：执行某个策略并写入结果
# ======================================================
run_policy() {
  local policy="$1"
  local log="${LOG_DIR}/policy_${policy}.log"
  local out_name="placement${M}_Col_1_T_${T}.txt"
  local out_renamed="${RESULT_DIR}/placement${M}_Col_1_T_${T}_${policy}.txt"

  echo "===== 运行策略 ${policy} (m=${M},h=${H},T=${T},rounds=${ROUNDS}) ====="

  if [[ "${policy}" == "lagrange" ]]; then
    ( set -x
      time "${BIN}" "${M}" "${H}" "${ROUNDS}" \
        --T="${T}" --policy=lagrange --freeze="${FREEZE}" \
        --rho="${RHO_LAG}" \
        --cap="${CAP}" --scale="${SCALE}" \
        --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}"
    ) 2>&1 | tee "${log}"

  elif [[ "${policy}" == "irl1" ]]; then
    ( set -x
      time "${BIN}" "${M}" "${H}" "${ROUNDS}" \
        --T="${T}" --policy=irl1 --freeze="${FREEZE}" \
        --lam="${LAM_IRL}" --delta="${DELTA_IRL}" \
        --cap="${CAP}" --scale="${SCALE}" \
        --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}"
    ) 2>&1 | tee "${log}"

  elif [[ "${policy}" == "soft" ]]; then
    ( set -x
      time "${BIN}" "${M}" "${H}" "${ROUNDS}" \
        --T="${T}" --policy=soft --freeze="${FREEZE}" \
        --alpha="${ALPHA_SOFT}" --eta="${ETA_SOFT}" \
        --gamma="${GAM_SOFT}" --beta="${BET_SOFT}" \
        --cap="${CAP}" --scale="${SCALE}" \
        --homotopy=1 --sigma="${SIGMA_SOFT}" --shrink="${SHRINK_SOFT}" \
        --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}"
    ) 2>&1 | tee "${log}"

  else
    echo "[ERR] 未知策略 ${policy}"; exit 3
  fi

  # 移动 placement 文件
  if [[ -f "${out_name}" ]]; then
    mv "${out_name}" "${out_renamed}"
    echo "[OK] 已移动 ${out_renamed}"
  else
    echo "[WARN] 未找到输出文件 (${out_name})，请检查日志。"
  fi

  # 解析日志，写入对比表
  if RES=$(parse_log "${log}" "${T}"); then
    set -- ${RES}
    MAXL_MAX="$1"; RMAX="$2"; MAXL_MIN="$3"; RMIN="$4"; HE_LAST="$5"; HW_LAST="$6"
    echo "${policy},${MAXL_MAX},${RMAX},${MAXL_MIN},${RMIN},${HE_LAST},${HW_LAST}" >> "${SUMMARY}"
  else
    echo "${policy},NA,NA,NA,NA,NA,NA" >> "${SUMMARY}"
  fi
}

# ======================================================
#                 运行三种策略
# ======================================================
run_policy "lagrange"
run_policy "irl1"
run_policy "soft"

# ======================================================
#                 打包结果
# ======================================================
echo "[INFO] 压缩结果目录 -> ${TARBALL}"
tar -czf "${TARBALL}" "${RESULT_DIR}"
echo "[DONE] 对比表已生成: ${SUMMARY}"
