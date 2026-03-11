#!/usr/bin/env bash
set -euo pipefail

# ===============================
# 参数解析
# ===============================
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <M> <H> [ROUNDS]"
  echo "Example: $0 16 16 150"
  exit 1
fi

M="$1"; H="$2"; ROUNDS="${3:-150}"

BIN="${BIN:-build/oc_multicol}"          # 多列封装器
POLICY="${POLICY:-soft}"                 # soft | lagrange | irl1
PROGRESS="${PROGRESS:-1}"                # oc_policy_singlecol 的 --progress 透传（通过 oc_multicol 执行）

# soft 参数（独立环境变量）
ALPHA_SOFT="${ALPHA_SOFT:-0.9}"
ETA_SOFT="${ETA_SOFT:-1.0}"
GAM_SOFT="${GAM_SOFT:-0.3}"
BET_SOFT="${BET_SOFT:-3.0}"
SIGMA_SOFT="${SIGMA_SOFT:-0.8}"
SHRINK_SOFT="${SHRINK_SOFT:-0.9}"
HOMOTOPY_SOFT="${HOMOTOPY_SOFT:-1}"

# irl1 参数（独立）
LAM_IRL="${LAM_IRL:-2.0}"
DELTA_IRL="${DELTA_IRL:-2.0}"
HOMOTOPY_IRL1="${HOMOTOPY_IRL1:-1}"
SIGMA_IRL1="${SIGMA_IRL1:-0.8}"
SHRINK_IRL1="${SHRINK_IRL1:-0.9}"

# lagrange 参数（独立）
RHO_LAG="${RHO_LAG:-0.3}"
HOMOTOPY_LAG="${HOMOTOPY_LAG:-0}"   # 默认关
SIGMA_LAG="${SIGMA_LAG:-0.8}"
SHRINK_LAG="${SHRINK_LAG:-0.9}"

# 通用
FREEZE="${FREEZE:-1}"
CAP="${CAP:-64}"
SCALE="${SCALE:-1000}"
DV="${DV:-1}"
POW2="${POW2:-1}"

# ===============================
# 结果目录
# ===============================
timestamp=$(date +%Y%m%d_%H%M%S)
ROOT="results_multicol_${M}x${H}_${timestamp}"
mkdir -p "${ROOT}"

# 全局汇总
GLOBAL="${ROOT}/summary_multicol_all.txt"
echo "M,H,K,h_prime,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,tighten_steps,final_Tprime,avg_solve_ms,last_solve_ms,last_picked,policy,alpha,eta,gamma,beta,lam,delta,rho,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path" > "${GLOBAL}"

# ===============================
# 解析函数（POSIX awk）
# ===============================
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

# ===============================
# 跑某个 K
# ===============================
run_for_K() {
  local K="$1"

  if (( H % K != 0 )); then
    echo "[INFO] Skip K=${K}: h=${H} not divisible by ${K}"
    return
  fi

  local HP=$(( H / K ))
  local MINV=$(( M < HP ? M : HP ))
  local TMIN=${MINV}
  local TMAX=$(( (5*MINV) / 2 ))   # floor(2.5*MINV)

  local RUN_DIR="${ROOT}/K${K}_${M}x${H}_T${TMIN}-${TMAX}"
  local LOG_DIR="${RUN_DIR}/logs"
  mkdir -p "${LOG_DIR}"

  local SUMMARY="${RUN_DIR}/summary_multicol_${M}x${H}_K${K}.txt"
  echo "M,H,K,h_prime,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,tighten_steps,final_Tprime,avg_solve_ms,last_solve_ms,last_picked,policy,alpha,eta,gamma,beta,lam,delta,rho,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path" > "${SUMMARY}"

  echo "[INFO] >>> Run K=${K}: base strip ${M}x${HP}, T∈[${TMIN},${TMAX}] rounds=${ROUNDS} policy=${POLICY}"

  for (( T="${TMIN}"; T<="${TMAX}"; T+=1 )); do
    local LOG="${LOG_DIR}/multicol_${M}x${H}_K${K}_T_${T}.log"
    local OUT="placement${M}_Col_${K}_T_${T}.txt"
    local OUT_MOVED="${RUN_DIR}/placement${M}_Col_${K}_T_${T}.txt"

    echo "----- [Run] ${M}x${H}  K=${K}  T=${T}  rounds=${ROUNDS} -----"
    SECONDS=0

    if [[ "${POLICY}" == "soft" ]]; then
      CMD=( "${BIN}" "${M}" "${H}" "${ROUNDS}"
            "--K=${K}" "--T=${T}" "--policy=soft"
            "--alpha=${ALPHA_SOFT}" "--eta=${ETA_SOFT}" "--gamma=${GAM_SOFT}" "--beta=${BET_SOFT}"
            "--homotopy=${HOMOTOPY_SOFT}" "--sigma=${SIGMA_SOFT}" "--shrink=${SHRINK_SOFT}"
            "--freeze=${FREEZE}" "--cap=${CAP}" "--scale=${SCALE}"
            "--dV=${DV}" "--pow2=${POW2}" "--progress=${PROGRESS}" "--snake=1" )
    elif [[ "${POLICY}" == "irl1" ]]; then
      CMD=( "${BIN}" "${M}" "${H}" "${ROUNDS}"
            "--K=${K}" "--T=${T}" "--policy=irl1"
            "--lam=${LAM_IRL}" "--delta=${DELTA_IRL}"
            "--homotopy=${HOMOTOPY_IRL1}" "--sigma=${SIGMA_IRL1}" "--shrink=${SHRINK_IRL1}"
            "--freeze=${FREEZE}" "--cap=${CAP}" "--scale=${SCALE}"
            "--dV=${DV}" "--pow2=${POW2}" "--progress=${PROGRESS}" "--snake=1")
    elif [[ "${POLICY}" == "lagrange" ]]; then
      CMD=( "${BIN}" "${M}" "${H}" "${ROUNDS}"
            "--K=${K}" "--T=${T}" "--policy=lagrange"
            "--rho=${RHO_LAG}"
            "--homotopy=${HOMOTOPY_LAG}" "--sigma=${SIGMA_LAG}" "--shrink=${SHRINK_LAG}"
            "--freeze=${FREEZE}" "--cap=${CAP}" "--scale=${SCALE}"
            "--dV=${DV}" "--pow2=${POW2}" "--progress=${PROGRESS}" "--snake=1")
    else
      echo "[ERR] Unknown POLICY=${POLICY}"; exit 3
    fi

    "${CMD[@]}" 2>&1 | tee "${LOG}"
    local RUNTIME=${SECONDS}

    if [[ -f "${OUT}" ]]; then
      mv "${OUT}" "${OUT_MOVED}"
    else
      echo "[WARN] Not found: ${OUT} (see ${LOG})"
      OUT_MOVED="NA"
    fi

    if RES=$(parse_log_metrics "${LOG}" "${T}"); then
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

  tar -czf "${RUN_DIR}.tar.gz" -C "${ROOT}" "$(basename "${RUN_DIR}")"
  echo "[INFO] Packed: ${RUN_DIR}.tar.gz"
}

# ===============================
# 主流程
# ===============================
if [[ ! -x "${BIN}" ]]; then
  echo "[ERR] Binary not found: ${BIN}. Please: make policy multicol"
  exit 2
fi

echo "[INFO] Output root: ${ROOT}"

# 视需要解开 2 列，这里默认跑 4 列和 5 列
# run_for_K 2
run_for_K 5
# run_for_K 5

tar -czf "${ROOT}.tar.gz" "${ROOT}"
echo "[DONE] All finished. Global summary:"
echo "  ${GLOBAL}"
