#!/usr/bin/env bash
set -euo pipefail

# ======================================================
#                 命令行参数解析
# ======================================================
usage() {
  cat <<EOF
Usage:
  $0 --pairs "8x8 16x8" --Tmin 8 --Tmax 20 [--Tstep 1] [--rounds 150]

Options:
  --pairs    "<MxH> <MxH> ..."   要跑的尺寸列表，例如 "8x8 16x8"
  --Tmin     <int>               T 起始值（含）
  --Tmax     <int>               T 结束值（含）
  --Tstep    <int>               T 步长，默认 1
  --rounds   <int>               每个 T 的最大迭代轮数，默认 150
  --out-prefix <str>             输出目录前缀，默认 results_soft
  --bin      <path>              可执行文件，默认 build/oc_policy_singlecol
  --progress <0|1>               打印 DP 进度，默认 1

示例：
  $0 --pairs "8x8 16x8" --Tmin 8 --Tmax 20 --rounds 150
EOF
  exit 1
}

PAIRS=""
TMIN=""
TMAX=""
TSTEP=1
ROUNDS=150
BIN="build/oc_policy_singlecol"
OUT_PREFIX="results_soft"
PROGRESS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pairs)   PAIRS="$2"; shift 2;;
    --Tmin)    TMIN="$2"; shift 2;;
    --Tmax)    TMAX="$2"; shift 2;;
    --Tstep)   TSTEP="$2"; shift 2;;
    --rounds)  ROUNDS="$2"; shift 2;;
    --bin)     BIN="$2"; shift 2;;
    --out-prefix) OUT_PREFIX="$2"; shift 2;;
    --progress) PROGRESS="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "[ERR] Unknown arg: $1"; usage;;
  esac
done

[[ -z "${PAIRS}" || -z "${TMIN}" || -z "${TMAX}" ]] && usage

# ======================================================
#          ⭐ 你指定的 Soft Policy 参数（已默认化）⭐
# ======================================================
ALPHA_SOFT="${ALPHA_SOFT:-0.50}"
ETA_SOFT="${ETA_SOFT:-0.50}"
GAM_SOFT="${GAM_SOFT:-0.40}"
BET_SOFT="${BET_SOFT:-3.0}"

SIGMA_SOFT="${SIGMA_SOFT:-0.85}"
SHRINK_SOFT="${SHRINK_SOFT:-0.92}"

CAP="${CAP:-8}"          # 你给的 cap=8
SCALE="${SCALE:-1000}"   # 你给的 scale=1000

DV="${DV:-1}"
POW2="${POW2:-1}"
FREEZE="${FREEZE:-1}"

echo "[INFO] Soft policy parameters:"
echo "       alpha=${ALPHA_SOFT}, eta=${ETA_SOFT}, gamma=${GAM_SOFT}, beta=${BET_SOFT}"
echo "       sigma=${SIGMA_SOFT}, shrink=${SHRINK_SOFT}, cap=${CAP}, scale=${SCALE}"

# ======================================================
#                结果目录与汇总文件
# ======================================================
timestamp=$(date +%Y%m%d_%H%M%S)
RESULT_ROOT="${OUT_PREFIX}_${timestamp}"
mkdir -p "${RESULT_ROOT}"
echo "[INFO] 输出根目录: ${RESULT_ROOT}"

GLOBAL_SUMMARY="${RESULT_ROOT}/summary_soft_all.txt"
echo "M,H,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,init_Tprime,tighten_steps,final_Tprime,min_maxL,round_of_min,max_maxL,round_of_max,avg_solve_ms,last_solve_ms,last_picked,alpha,eta,gamma,beta,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path,tprime_path" > "${GLOBAL_SUMMARY}"

# ======================================================
#                解析日志
# ======================================================
parse_log_metrics() {
  local log="$1"
  local Tval="$2"
  awk -v T="${Tval}" '
  BEGIN {
    init_tp="NA"; tighten=0;
    nr=0; minL=1e18; rmin=-1; maxL=-1e18; rmax=-1;
    last_r=-1; first_meet=-1;
    sum_ms=0; cnt=0;
    last_ml="NA"; last_he="NA"; last_hw="NA"; last_ms="NA"; last_pk="NA"; last_WNS_T="NA"; last_Tp="NA";
  }
  /\[Init-Homotopy\]/ {
    if (index($0,"T_prime0=")>0) {
      split($0,a,"T_prime0="); split(a[2],b," "); init_tp=b[1];
    }
  }
  /\[Tighten\]/ { tighten++ }

  /\[Round[[:space:]]+[0-9]+\]/ {
    line=$0
    tmp=line
    sub(/^.*\[Round[[:space:]]+/, "", tmp)
    split(tmp, ff, /[^0-9]+/)
    r=(ff[1]==""? -1: ff[1]+0)
    if (r<0) next

    Tp=""; if (index(line,"T_prime=")>0){ split(line,a,"T_prime="); split(a[2],b," "); Tp=b[1] }
    ml=""; if (index(line,"maxL=")>0){ split(line,a2,"maxL="); split(a2[2],b2," "); ml=b2[1] }
    he=""; if (index(line,"HPWL_eq=")>0){ split(line,a3,"HPWL_eq="); split(a3[2],b3," "); he=b3[1] }
    hw=""; if (index(line,"HPWL_w=")>0){ split(line,a4,"HPWL_w="); split(a4[2],b4," "); hw=b4[1] }
    ms=""; if (index(line,"solve_ms=")>0){ split(line,a5,"solve_ms="); split(a5[2],b5," "); ms=b5[1] }
    pk=""; if (index(line,"picked=")>0){ split(line,a6,"picked="); split(a6[2],b6," "); pk=b6[1] }

    if (ml != ""){
      if (ml < minL){ minL=ml; rmin=r; }
      if (ml > maxL){ maxL=ml; rmax=r; }
      if (first_meet<0 && ml <= T) first_meet=r;
    }

    last_r=r; last_ml=ml; last_he=he; last_hw=hw;
    last_ms=ms; last_pk=pk; last_Tp=Tp;

    if (ml!="") last_WNS_T = (T - ml);

    if (ms!="") { sum_ms+=ms; cnt++; }
  }
  END {
    avg_ms=(cnt>0?sum_ms/cnt:0);
    printf "%s %d %s %s %s %s %s %s %s %s %s %s %s %.3f %s %s\n",
      init_tp, tighten, last_Tp, last_r, last_ml, last_he, last_hw, last_WNS_T,
      (first_meet<0?"NA":first_meet), minL, rmin, maxL, rmax, avg_ms, last_ms, last_pk;
  }' "${log}"
}

# ======================================================
#          dump T_prime evolution（每个 T）
# ======================================================
dump_tprime_evolution() {
  local log="$1"
  local outfile="$2"
  local T="$3"
  {
    echo "# round,T_prime,maxL,WNS(T'),WNS(T),picked,HPWL_eq,HPWL_w,solve_ms"
    awk -v TT="${T}" '
      /\[Round[[:space:]]+[0-9]+\]/ {
        line=$0
        tmp=line
        sub(/^.*\[Round[[:space:]]+/, "", tmp)
        split(tmp, ff, /[^0-9]+/)
        r=ff[1]

        Tp=""; ml=""; he=""; hw=""; ms=""; pk="";
        if(index(line,"T_prime=")>0){ split(line,a,"T_prime="); split(a[2],b," "); Tp=b[1] }
        if(index(line,"maxL=")>0){ split(line,a2,"maxL="); split(a2[2],b2," "); ml=b2[1] }
        if(index(line,"HPWL_eq=")>0){ split(line,a3,"HPWL_eq="); split(a3[2],b3," "); he=b3[1] }
        if(index(line,"HPWL_w=")>0){ split(line,a4,"HPWL_w="); split(a4[2],b4," "); hw=b4[1] }
        if(index(line,"solve_ms=")>0){ split(line,a5,"solve_ms="); split(a5[2],b5," "); ms=b5[1] }
        if(index(line,"picked=")>0){ split(line,a6,"picked="); split(a6[2],b6," "); pk=b6[1] }

        wnsp = (Tp=="" || ml=="" ? "" : (Tp - ml))
        wnst = (ml=="" ? "" : (TT - ml))

        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n", r,Tp,ml,wnsp,wnst,pk,he,hw,ms
      }
    ' "${log}"
  } > "${outfile}"
}

# ======================================================
#                 主循环：逐 case × T
# ======================================================
if [[ ! -x "${BIN}" ]]; then
  echo "[ERR] 找不到可执行文件: ${BIN}"
  exit 2
fi

for pair in ${PAIRS}; do
  M="${pair%x*}"; H="${pair#*x}"

  RUN_DIR="${RESULT_ROOT}/soft_${M}x${H}_T${TMIN}-${TMAX}"
  LOG_DIR="${RUN_DIR}/logs"
  mkdir -p "${LOG_DIR}"

  SUMMARY="${RUN_DIR}/summary_soft_${M}x${H}.txt"
  echo "M,H,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,init_Tprime,tighten_steps,final_Tprime,min_maxL,round_of_min,max_maxL,round_of_max,avg_solve_ms,last_solve_ms,last_picked,alpha,eta,gamma,beta,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path,tprime_path" > "${SUMMARY}"

  echo "[INFO] >>> Case ${M}x${H}, T=[${TMIN}, ${TMAX}] step=${TSTEP}, rounds=${ROUNDS}"

  for (( T="${TMIN}"; T<="${TMAX}"; T+=TSTEP )); do
    LOG="${LOG_DIR}/soft_${M}x${H}_T_${T}.log"
    OUT_NAME="placement${M}_Col_1_T_${T}.txt"
    OUT_MOVED="${RUN_DIR}/placement${M}x${H}_Col_1_T_${T}.txt"
    TP_EV="${RUN_DIR}/tprime_soft_${M}x${H}_T_${T}.txt"

    echo "----- [Run] ${M}x${H}, T=${T}, rounds=${ROUNDS} -----"
    SECONDS=0

    "${BIN}" "${M}" "${H}" "${ROUNDS}" \
      --T="${T}" --policy=soft \
      --alpha="${ALPHA_SOFT}" --eta="${ETA_SOFT}" \
      --gamma="${GAM_SOFT}" --beta="${BET_SOFT}" \
      --homotopy=1 --sigma="${SIGMA_SOFT}" --shrink="${SHRINK_SOFT}" \
      --cap="${CAP}" --scale="${SCALE}" --freeze="${FREEZE}" \
      --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}" \
      2>&1 | tee "${LOG}"

    runtime_sec=${SECONDS}

    [[ -f "${OUT_NAME}" ]] && mv "${OUT_NAME}" "${OUT_MOVED}" || OUT_MOVED="NA"

    dump_tprime_evolution "${LOG}" "${TP_EV}" "${T}"

    if RES=$(parse_log_metrics "${LOG}" "${T}"); then
      set -- ${RES}
      INIT_TP="$1"; TIGHT="$2"; FINAL_TP="$3"; LAST_R="$4"; LAST_MAXL="$5"; LAST_HE="$6"; LAST_HW="$7"; LAST_WT="$8"
      FIRST_MEET="$9"; MINL="${10}"; RMIN="${11}"; MAXL="${12}"; RMAX="${13}"; AVG_MS="${14}"; LAST_MS="${15}"; LAST_PK="${16}"

      MEET="NO"; [[ "${LAST_MAXL}" != "NA" && "${LAST_MAXL}" -le "${T}" ]] && MEET="YES"

      echo "${M},${H},${T},${ROUNDS},${runtime_sec},${LAST_R},${LAST_MAXL},${LAST_HE},${LAST_HW},${LAST_WT},${MEET},${FIRST_MEET},${INIT_TP},${TIGHT},${FINAL_TP},${MINL},${RMIN},${MAXL},${RMAX},${AVG_MS},${LAST_MS},${LAST_PK},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED},${TP_EV}" >> "${SUMMARY}"

      echo "${M},${H},${T},${ROUNDS},${runtime_sec},${LAST_R},${LAST_MAXL},${LAST_HE},${LAST_HW},${LAST_WT},${MEET},${FIRST_MEET},${INIT_TP},${TIGHT},${FINAL_TP},${MINL},${RMIN},${MAXL},${RMAX},${AVG_MS},${LAST_MS},${LAST_PK},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED},${TP_EV}" >> "${GLOBAL_SUMMARY}"
    fi
  done

  tar -czf "${RUN_DIR}.tar.gz" -C "${RESULT_ROOT}" "$(basename "${RUN_DIR}")"
done

tar -czf "${RESULT_ROOT}.tar.gz" "${RESULT_ROOT}"

echo "[DONE] 全部完成。全局汇总：${GLOBAL_SUMMARY}"
