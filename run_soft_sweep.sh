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

你也可以通过环境变量覆盖 soft 的参数：
  ALPHA_SOFT=0.9 ETA_SOFT=1.0 GAM_SOFT=0.3 BET_SOFT=3.0 \\
  SIGMA_SOFT=0.8 SHRINK_SOFT=0.9 CAP=64 SCALE=1000 \\
  DV=1 POW2=1 FREEZE=1  $0 --pairs "8x8 16x8" --Tmin 8 --Tmax 20

示例（你当前的首要测试）：
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
#                   软策略的默认参数
#   可用环境变量覆盖：ALPHA_SOFT、ETA_SOFT、GAM_SOFT、BET_SOFT、
#                    SIGMA_SOFT、SHRINK_SOFT、CAP、SCALE、DV、POW2、FREEZE
# ======================================================
ALPHA_SOFT="${ALPHA_SOFT:-0.9}"
ETA_SOFT="${ETA_SOFT:-1.0}"
GAM_SOFT="${GAM_SOFT:-0.3}"
BET_SOFT="${BET_SOFT:-3.0}"
SIGMA_SOFT="${SIGMA_SOFT:-0.8}"
SHRINK_SOFT="${SHRINK_SOFT:-0.9}"

CAP="${CAP:-64}"
SCALE="${SCALE:-1000}"

DV="${DV:-1}"
POW2="${POW2:-1}"
FREEZE="${FREEZE:-1}"

# ======================================================
#                结果目录与汇总文件
# ======================================================
timestamp=$(date +%Y%m%d_%H%M%S)
RESULT_ROOT="${OUT_PREFIX}_${timestamp}"
mkdir -p "${RESULT_ROOT}"
echo "[INFO] 输出根目录: ${RESULT_ROOT}"

# 汇总表头（记录指标尽量全面）
# 注：参数也记录，便于复现实验
GLOBAL_SUMMARY="${RESULT_ROOT}/summary_soft_all.txt"
echo "M,H,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,init_Tprime,tighten_steps,final_Tprime,min_maxL,round_of_min,max_maxL,round_of_max,avg_solve_ms,last_solve_ms,last_picked,alpha,eta,gamma,beta,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path,tprime_path" > "${GLOBAL_SUMMARY}"

# ======================================================
#       解析日志：抽取每轮与最终指标（POSIX awk）
# ======================================================
# 输出 16 个字段（空格分隔）：
# 1:init_Tprime 2:tighten_steps 3:final_Tprime 4:last_round 5:last_maxL 6:last_HPWL_eq
# 7:last_HPWL_w 8:last_WNS_T 9:first_meet_T_round 10:min_maxL 11:round_of_min
# 12:max_maxL 13:round_of_max 14:avg_solve_ms 15:last_solve_ms 16:last_picked
parse_log_metrics() {
  local log="$1"
  local Tval="$2"
  awk -v T="${Tval}" '
  BEGIN {
    init_tp="NA"; tighten=0;
    nr=0; minL=1e18; rmin=-1; maxL=-1e18; rmax=-1;
    last_r=-1; first_meet=-1;
    sum_ms=0; cnt=0;
    last_ml="NA"; last_he="NA"; last_hw="NA"; last_ms="NA"; last_pk="NA"; last_Tp="NA"; last_Tt="NA"; last_WNS_T="NA";
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

    # Extract tokens
    Tp=""; if (index(line,"T_prime=")>0){ split(line,a,"T_prime="); split(a[2],b," "); Tp=b[1]+0 }
    Tt=""; if (index(line,"T_target=")>0){ split(line,a,"T_target="); split(a[2],b," "); Tt=b[1]+0 }
    ml=""; if (index(line,"maxL=")>0){ split(line,a,"maxL="); split(a[2],b," "); ml=b[1]+0 }
    he=""; if (index(line,"HPWL_eq=")>0){ split(line,a,"HPWL_eq="); split(a[2],b," "); he=b[1]+0 }
    hw=""; if (index(line,"HPWL_w=")>0){ split(line,a,"HPWL_w="); split(a[2],b," "); hw=b[1]+0 }
    ms=""; if (index(line,"solve_ms=")>0){ split(line,a,"solve_ms="); split(a[2],b," "); ms=b[1]+0 }
    pk=""; if (index(line,"picked=")>0){ split(line,a,"picked="); split(a[2],b," "); pk=b[1]+0 }
    # WNS*(T)= 可选字段；没有就用 T - maxL
    wT=""; if (index(line,"WNS*(T)=")>0){ split(line,a,"WNS*(T)="); split(a[2],b," "); wT=b[1]+0 }

    # update min/max & first meet
    if (ml != "") {
      if (ml < minL){ minL=ml; rmin=r; }
      if (ml > maxL){ maxL=ml; rmax=r; }
      if (first_meet<0 && ml <= T){ first_meet=r; }
    }

    # last values
    last_r=r; last_ml=(ml==""?"NA":ml); last_he=(he==""?"NA":he); last_hw=(hw==""?"NA":hw);
    last_ms=(ms==""?"NA":ms); last_pk=(pk==""?"NA":pk); last_Tp=(Tp==""?"NA":Tp); last_Tt=(Tt==""?"NA":Tt);
    if (wT!="") last_WNS_T=wT; else if (ml!="") last_WNS_T=(T - ml);

    if (ms!="") { sum_ms += ms; cnt++; }
  }
  END {
    avg_ms = (cnt>0 ? sum_ms/cnt : 0);
    printf "%s %d %s %s %s %s %s %s %s %s %s %s %s %.3f %s %s\n",
      init_tp, tighten, last_Tp, last_r, last_ml, last_he, last_hw, last_WNS_T,
      (first_meet<0 ? "NA" : first_meet), minL, rmin, maxL, rmax, avg_ms, last_ms, last_pk;
  }' "${log}"
}

# ======================================================
#     生成 T_prime 演化 TXT：round + tighten（CSV）
# ======================================================
dump_tprime_evolution() {
  local log="$1"
  local outfile="$2"
  local tgtT="$3"
  {
    echo "# T_prime evolution (target T=${tgtT})"
    echo "type,round,T_prime,T_target,maxL,WNS_Tprime,WNS_T,picked,HPWL_eq,HPWL_w,solve_ms,old_T_prime,new_T_prime,sigma_maxL,shrink_val"
    awk -v TT="${tgtT}" '
      BEGIN { last_round=""; }
      /\[Round[[:space:]]+[0-9]+\]/ {
        line=$0
        tmp=line
        sub(/^.*\[Round[[:space:]]+/, "", tmp)
        split(tmp, ff, /[^0-9]+/)
        r=(ff[1]==""? -1: ff[1]+0); if (r<0) next
        last_round=r

        Tp=""; if (index(line,"T_prime=")>0){ split(line,a,"T_prime="); split(a[2],b," "); Tp=b[1]+0 }
        Tt=""; if (index(line,"T_target=")>0){ split(line,a,"T_target="); split(a[2],b," "); Tt=b[1]+0 }
        ml=""; if (index(line,"maxL=")>0){ split(line,a,"maxL="); split(a[2],b," "); ml=b[1]+0 }
        pick=""; if (index(line,"picked=")>0){ split(line,a,"picked="); split(a[2],b," "); pick=b[1]+0 }
        he=""; if (index(line,"HPWL_eq=")>0){ split(line,a,"HPWL_eq="); split(a[2],b," "); he=b[1]+0 }
        hw=""; if (index(line,"HPWL_w=")>0){ split(line,a,"HPWL_w="); split(a[2],b," "); hw=b[1]+0 }
        ms=""; if (index(line,"solve_ms=")>0){ split(line,a,"solve_ms="); split(a[2],b," "); ms=b[1]+0 }


        # 计算两种 WNS
        wnsT  = (ml=="" ? "" : (TT - ml))
        wnsp  = (Tp=="" || ml=="" ? "" : (Tp - ml))

        printf "round,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,,,,\n", r,Tp,Tt,ml,wnsp,wnsT,pick,he,hw,ms
      }
      /\[Tighten\]/ {
        line=$0
        old=""; newv=""; smax=""; shr="";
        if (index(line,"T_prime:")>0){
          split(line, a, "T_prime:");
          rest=a[2]; sub(/^[[:space:]]+/,"",rest)
          split(rest, b, /[[:space:]]*->[[:space:]]*/)
          old=b[1]
          newpart=b[2]; sub(/^[[:space:]]+/,"",newpart)
          split(newpart, c, /[^0-9]+/); newv=c[1]
          if(index(line,"sigma*maxL=")>0){ split(line,d,"sigma*maxL="); split(d[2],e,/[), ]/); smax=e[1]; }
          if(index(line,"shrink=")>0){ split(line,f,"shrink="); split(f[2],g,/[), ]/); shr=g[1]; }
        }
        r=(last_round==""? "": last_round)
        printf "tighten,%s,,,,,,,,,,%s,%s,%s,%s\n", r, old, newv, smax, shr
      }
    ' "${log}"
  } > "${outfile}"
}

# ======================================================
#                 主循环：逐尺寸 x 逐 T
# ======================================================
if [[ ! -x "${BIN}" ]]; then
  echo "[ERR] 可执行文件不存在: ${BIN}"
  exit 2
fi

for pair in ${PAIRS}; do
  M="${pair%x*}"; H="${pair#*x}"
  if ! [[ "${M}" =~ ^[0-9]+$ && "${H}" =~ ^[0-9]+$ ]]; then
    echo "[WARN] 跳过非法尺寸：${pair}"; continue
  fi

  RUN_DIR="${RESULT_ROOT}/soft_${M}x${H}_T${TMIN}-${TMAX}"
  LOG_DIR="${RUN_DIR}/logs"
  mkdir -p "${LOG_DIR}"

  SUMMARY="${RUN_DIR}/summary_soft_${M}x${H}.txt"
  echo "M,H,T,rounds_total,total_runtime_sec,last_round,last_maxL,last_HPWL_eq,last_HPWL_w,last_WNS_T,meet_T,first_meet_T_round,init_Tprime,tighten_steps,final_Tprime,min_maxL,round_of_min,max_maxL,round_of_max,avg_solve_ms,last_solve_ms,last_picked,alpha,eta,gamma,beta,sigma,shrink,cap,scale,dV,pow2,freeze,log_path,placement_path,tprime_path" > "${SUMMARY}"

  echo "[INFO] >>> Size ${M}x${H}, T in [${TMIN}, ${TMAX}] step ${TSTEP}, rounds=${ROUNDS}"

  for (( T="${TMIN}"; T<="${TMAX}"; T+=TSTEP )); do
    LOG="${LOG_DIR}/soft_${M}x${H}_T_${T}.log"
    OUT_NAME="placement${M}_Col_1_T_${T}.txt"
    OUT_MOVED="${RUN_DIR}/placement${M}x${H}_Col_1_T_${T}.txt"
    TPEV="${RUN_DIR}/tprime_evolution_soft_${M}x${H}_T_${T}.txt"

    echo "----- [Run] ${M}x${H}  T=${T}  rounds=${ROUNDS} -----"
    SECONDS=0
    "${BIN}" "${M}" "${H}" "${ROUNDS}" \
      --T="${T}" --policy=soft --freeze="${FREEZE}" \
      --alpha="${ALPHA_SOFT}" --eta="${ETA_SOFT}" \
      --gamma="${GAM_SOFT}" --beta="${BET_SOFT}" \
      --cap="${CAP}" --scale="${SCALE}" \
      --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}" \
      --homotopy=1 --sigma="${SIGMA_SOFT}" --shrink="${SHRINK_SOFT}" \
      2>&1 | tee "${LOG}"
    runtime_sec=${SECONDS}

    # 移动放置文件
    if [[ -f "${OUT_NAME}" ]]; then
      mv "${OUT_NAME}" "${OUT_MOVED}"
    else
      OUT_MOVED="NA"
      echo "[WARN] 未发现放置文件：${OUT_NAME}"
    fi

    # 生成 T_prime 演化
    dump_tprime_evolution "${LOG}" "${TPEV}" "${T}"

    # 解析日志得到指标
    if RES=$(parse_log_metrics "${LOG}" "${T}"); then
      # 16个字段：
      # init_Tprime tighten_steps final_Tprime last_round last_maxL last_HPWL_eq last_HPWL_w last_WNS_T first_meet_T_round min_maxL rmin maxL rmax avg_ms last_ms last_pk
      set -- ${RES}
      INIT_TP="$1"; TIGHT="$2"; FINAL_TP="$3"; LAST_R="$4"; LAST_MAXL="$5"; LAST_HE="$6"; LAST_HW="$7"; LAST_WT="$8"; FIRST_MEET="$9"
      MIN_L="${10}"; RMIN="${11}"; MAX_L="${12}"; RMAX="${13}"; AVG_MS="${14}"; LAST_MS="${15}"; LAST_PK="${16}"

      MEET_T="NO"; [[ "${LAST_MAXL}" != "NA" && "${LAST_MAXL}" -le "${T}" ]] && MEET_T="YES"

      # 写入当前尺寸的 summary
      echo "${M},${H},${T},${ROUNDS},${runtime_sec},${LAST_R},${LAST_MAXL},${LAST_HE},${LAST_HW},${LAST_WT},${MEET_T},${FIRST_MEET},${INIT_TP},${TIGHT},${FINAL_TP},${MIN_L},${RMIN},${MAX_L},${RMAX},${AVG_MS},${LAST_MS},${LAST_PK},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED},${TPEV}" >> "${SUMMARY}"

      # 写入全局汇总
      echo "${M},${H},${T},${ROUNDS},${runtime_sec},${LAST_R},${LAST_MAXL},${LAST_HE},${LAST_HW},${LAST_WT},${MEET_T},${FIRST_MEET},${INIT_TP},${TIGHT},${FINAL_TP},${MIN_L},${RMIN},${MAX_L},${RMAX},${AVG_MS},${LAST_MS},${LAST_PK},${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED},${TPEV}" >> "${GLOBAL_SUMMARY}"
    else
      echo "[WARN] 解析失败：${LOG}"
      echo "${M},${H},${T},${ROUNDS},${runtime_sec},NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED},${TPEV}" >> "${SUMMARY}"
      echo "${M},${H},${T},${ROUNDS},${runtime_sec},NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,${ALPHA_SOFT},${ETA_SOFT},${GAM_SOFT},${BET_SOFT},${SIGMA_SOFT},${SHRINK_SOFT},${CAP},${SCALE},${DV},${POW2},${FREEZE},${LOG},${OUT_MOVED},${TPEV}" >> "${GLOBAL_SUMMARY}"
    fi
  done

  # 每个尺寸打包
  tar -czf "${RUN_DIR}.tar.gz" -C "${RESULT_ROOT}" "$(basename "${RUN_DIR}")"
  echo "[INFO] 打包尺寸 ${M}x${H}: ${RUN_DIR}.tar.gz"
done

# 根目录也打包一份
tar -czf "${RESULT_ROOT}.tar.gz" "${RESULT_ROOT}"
echo "[DONE] 全部完成。全局汇总：${GLOBAL_SUMMARY}"
