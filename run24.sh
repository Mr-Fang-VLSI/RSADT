#!/usr/bin/env bash
set -euo pipefail

# ============================================
# Usage:
#   ./run24.sh <Tmin> <Tmax> [ROUNDS]
#
# Example:
#   ./run24.sh 8 20 150
#
# 固定跑 24x24, K=4，多列封装器 build/oc_multicol
# T 从 Tmin 到 Tmax（步长 1），你可以用环境变量调 soft / irl1 / lagrange 的参数。
# ============================================

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <Tmin> <Tmax> [ROUNDS]"
  echo "Example: $0 8 20 150"
  exit 1
fi

TMIN="$1"
TMAX="$2"
ROUNDS="${3:-150}"

M=24
H=24
K="${K:-4}"                      # 默认 4 列，你也可以临时 K=5 ./run24.sh ...

BIN="${BIN:-build/oc_multicol}"  # 多列封装器
POLICY="${POLICY:-soft}"         # soft | irl1 | lagrange
PROGRESS="${PROGRESS:-1}"        # 透传给 oc_policy_singlecol

# ============================================
# soft 参数 —— 默认用你 8x8 上那组
#   alpha=0.5 eta=0.5 gamma=0.4 beta=3.0
#   sigma=0.85 shrink=0.92 cap=8 scale=1000
# ============================================
ALPHA_SOFT="${ALPHA_SOFT:-0.5}"
ETA_SOFT="${ETA_SOFT:-0.5}"
GAM_SOFT="${GAM_SOFT:-0.4}"
BET_SOFT="${BET_SOFT:-3.0}"
SIGMA_SOFT="${SIGMA_SOFT:-0.85}"
SHRINK_SOFT="${SHRINK_SOFT:-0.92}"
HOMOTOPY_SOFT="${HOMOTOPY_SOFT:-1}"

# irl1 参数（用不到也保持兼容）
LAM_IRL="${LAM_IRL:-2.0}"
DELTA_IRL="${DELTA_IRL:-2.0}"
HOMOTOPY_IRL1="${HOMOTOPY_IRL1:-1}"
SIGMA_IRL1="${SIGMA_IRL1:-0.8}"
SHRINK_IRL1="${SHRINK_IRL1:-0.9}"

# lagrange 参数（用不到也保持兼容）
RHO_LAG="${RHO_LAG:-0.3}"
HOMOTOPY_LAG="${HOMOTOPY_LAG:-0}"
SIGMA_LAG="${SIGMA_LAG:-0.8}"
SHRINK_LAG="${SHRINK_LAG:-0.9}"

# 通用参数
FREEZE="${FREEZE:-1}"
CAP="${CAP:-8}"                  # 注意：默认 cap=8
SCALE="${SCALE:-1000}"
DV="${DV:-1}"
POW2="${POW2:-1}"

# ============================================
# 结果目录 & 全局 summary
# ============================================
timestamp=$(date +%Y%m%d_%H%M%S)
ROOT="results_24x24_K${K}_${timestamp}"
mkdir -p "${ROOT}"

GLOBAL="${ROOT}/summary_24x24_K${K}_all.txt"
# 和 run_multicol_K2K4 的列基本对齐
echo "M,H,K,h_prime,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,tighten_steps,final_Tprime,avg_solve_ms,last_solve_ms,last_picked,policy,alpha,eta,gamma,beta,lam,delta,rho,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path" > "${GLOBAL}"

# ============================================
# 解析日志：抓 maxL、HPWL、solve_ms 等
# （和 run_multicol_K2K4 保持一致）
# ============================================
parse_log_metrics() {
  local log="$1"
  local Tval="$2"
  awk -v T="${Tval}" '
  BEGIN {
    tighten=0; finalTp="NA";
    nr=0; first_meet=-1;
    last_r=-1; last_ml="NA"; last_he="NA"; last_hw="NA"; last_ms="NA"; last_wt="NA";
    sum_ms=0; cnt=0;
  }
  /\[Tighten\]/ { tighten++ }

  /\[Round[[:space:]]+[0-9]+\]/ {
    line=$0
    tmp=line
    sub(/^.*\[Round[[:space:]]+/, "", tmp)
    split(tmp, ff, /[^0-9]+/)
    r=(ff[1]==""? -1: ff[1]+0)
    if (r<0) next

    Tp=""; if (index(line,"T_prime=")>0){ split(line,a,"T_prime="); split(a[2],b," "); Tp=b[1]+0 }
    ml=""; if (index(line,"maxL=")>0){ split(line,a2,"maxL="); split(a2[2],b2," "); ml=b2[1]+0 }
    he=""; if (index(line,"HPWL_eq=")>0){ split(line,a3,"HPWL_eq="); split(a3[2],b3," "); he=b3[1]+0 }
    hw=""; if (index(line,"HPWL_w=")>0){ split(line,a4,"HPWL_w="); split(a4[2],b4," "); hw=b4[1]+0 }
    ms=""; if (index(line,"solve_ms=")>0){ split(line,a5,"solve_ms="); split(a5[2],b5," "); ms=b5[1]+0 }
    wt=""; if (index(line,"WNS*(T)=")>0){ split(line,a6,"WNS*(T)="); split(a6[2],b6," "); wt=b6[1]+0 }

    last_r=r; if (Tp!="") finalTp=Tp;
    if (ml!=""){
      last_ml=ml;
      if (first_meet<0 && ml<=T) first_meet=r;
    }
    if (he!="") last_he=he;
    if (hw!="") last_hw=hw;
    if (ms!=""){ last_ms=ms; sum_ms+=ms; cnt++; }
    if (wt!="") last_wt=wt;

    nr++;
  }
  END {
    avg_ms=(cnt>0 ? sum_ms/cnt : 0);
    meet=(last_ml!="NA" && last_ml<=T) ? "YES" : "NO";
    printf "%s %s %s %s %s %s %s %s %s %d\n",
      last_r, last_ml, last_he, last_hw, last_wt, meet, (first_meet<0?"NA":first_meet), tighten, finalTp, avg_ms;
  }' "${log}"
}

# ============================================
# 主流程：固定 24x24，只扫 T
# ============================================
if [[ ! -x "${BIN}" ]]; then
  echo "[ERR] Binary not found: ${BIN}. Please: make multicol"
  exit 2
fi

# 检查 K 是否可行
if (( H % K != 0 )); then
  echo "[ERR] H=${H} not divisible by K=${K}"
  exit 3
fi
HP=$(( H / K ))    # strip 高度（也是 h_prime）

LOG_DIR="${ROOT}/logs"
mkdir -p "${LOG_DIR}"

SUMMARY="${ROOT}/summary_24x24_K${K}.txt"
echo "M,H,K,h_prime,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,tighten_steps,final_Tprime,avg_solve_ms,last_solve_ms,last_picked,policy,alpha,eta,gamma,beta,lam,delta,rho,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path" > "${SUMMARY}"

echo "[INFO] 24x24, K=${K}, strip=${M}x${HP}, T∈[${TMIN},${TMAX}], rounds=${ROUNDS}, policy=${POLICY}"

for (( T="${TMIN}"; T<="${TMAX}"; T+=1 )); do
  LOG="${LOG_DIR}/multicol_24x24_K${K}_T_${T}.log"
  OUT="placement${M}_Col_${K}_T_${T}.txt"
  OUT_MOVED="${ROOT}/placement24x24_Col_${K}_T_${T}.txt"

  echo "----- [Run] 24x24  K=${K}  T=${T}  rounds=${ROUNDS} -----"
  SECONDS=0

  # 组装命令
  if [[ "${POLICY}" == "soft" ]]; then
    CMD=( "${BIN}" "${M}" "${H}" "${ROUNDS}"
          "--K=${K}" "--T=${T}" "--policy=soft"
          "--alpha=${ALPHA_SOFT}" "--eta=${ETA_SOFT}"
          "--gamma=${GAM_SOFT}" "--beta=${BET_SOFT}"
          "--homotopy=${HOMOTOPY_SOFT}" "--sigma=${SIGMA_SOFT}" "--shrink=${SHRINK_SOFT}"
          "--freeze=${FREEZE}" "--cap=${CAP}" "--scale=${SCALE}"
          "--dV=${DV}" "--pow2=${POW2}" "--progress=${PROGRESS}" "--snake=1" )
  elif [[ "${POLICY}" == "irl1" ]]; then
    CMD=( "${BIN}" "${M}" "${H}" "${ROUNDS}"
          "--K=${K}" "--T=${T}" "--policy=irl1"
          "--lam=${LAM_IRL}" "--delta=${DELTA_IRL}"
          "--homotopy=${HOMOTOPY_IRL1}" "--sigma=${SIGMA_IRL1}" "--shrink=${SHRINK_IRL1}"
          "--freeze=${FREEZE}" "--cap=${CAP}" "--scale=${SCALE}"
          "--dV=${DV}" "--pow2=${POW2}" "--progress=${PROGRESS}" "--snake=1" )
  elif [[ "${POLICY}" == "lagrange" ]]; then
    CMD=( "${BIN}" "${M}" "${H}" "${ROUNDS}"
          "--K=${K}" "--T=${T}" "--policy=lagrange"
          "--rho=${RHO_LAG}"
          "--homotopy=${HOMOTOPY_LAG}" "--sigma=${SIGMA_LAG}" "--shrink=${SHRINK_LAG}"
          "--freeze=${FREEZE}" "--cap=${CAP}" "--scale=${SCALE}"
          "--dV=${DV}" "--pow2=${POW2}" "--progress=${PROGRESS}" "--snake=1" )
  else
    echo "[ERR] Unknown POLICY=${POLICY}"
    exit 4
  fi

  # 运行
  "${CMD[@]}" 2>&1 | tee "${LOG}"
  RUNTIME=${SECONDS}

  # 移动 placement
  if [[ -f "${OUT}" ]]; then
    mv "${OUT}" "${OUT_MOVED}"
  else
    echo "[WARN] Not found placement: ${OUT} (see ${LOG})"
    OUT_MOVED="NA"
  fi

  # 解析日志
  if RES=$(parse_log_metrics "${LOG}" "${T}"); then
    # fields: last_r last_ml last_he last_hw last_wt meet first_meet tighten finalTp avg_ms
    set -- ${RES}
    LAST_R="$1"; LAST_ML="$2"; LAST_HE="$3"; LAST_HW="$4"; LAST_WT="$5"; MEET="$6"; FIRST_MEET="$7"; TIGHTEN="$8"; FINAL_TP="$9"; AVG_MS="${10}"
    echo "${M},${H},${K},${HP},${T},${ROUNDS},${RUNTIME},${LAST_R},${LAST_ML},${LAST_HE},${LAST_HW},${LAST_WT},${MEET},${FIRST_MEET},${TIGHTEN},${FINAL_TP},${AVG_MS},${LAST_MS:-NA},${LAST_PK:-NA},${POLICY},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${LAM_IRL},${DELTA_IRL},${RHO_LAG},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED}" >> "${SUMMARY}"
    echo "${M},${H},${K},${HP},${T},${ROUNDS},${RUNTIME},${LAST_R},${LAST_ML},${LAST_HE},${LAST_HW},${LAST_WT},${MEET},${FIRST_MEET},${TIGHTEN},${FINAL_TP},${AVG_MS},${LAST_MS:-NA},${LAST_PK:-NA},${POLICY},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${LAM_IRL},${DELTA_IRL},${RHO_LAG},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED}" >> "${GLOBAL}"
  else
    echo "[WARN] Parse failed for ${LOG}"
    echo "${M},${H},${K},${HP},${T},${ROUNDS},${RUNTIME},NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,${POLICY},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${LAM_IRL},${DELTA_IRL},${RHO_LAG},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED}" >> "${SUMMARY}"
    echo "${M},${H},${K},${HP},${T},${ROUNDS},${RUNTIME},NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,${POLICY},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${LAM_IRL},${DELTA_IRL},${RHO_LAG},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED}" >> "${GLOBAL}"
  fi
done

# 打包
tar -czf "${ROOT}.tar.gz" "${ROOT}"
echo "[DONE] All finished. Global summary:"
echo "  ${GLOBAL}"
