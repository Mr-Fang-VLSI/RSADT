#!/usr/bin/env bash
set -euo pipefail

# ======================================================
#                参数解析与默认值
# ======================================================
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <M> <H> <T> [ROUNDS]"
  echo "Example: $0 12 8 12 120"
  exit 1
fi

M="$1"; H="$2"; T="$3"
ROUNDS="${4:-120}"

BIN="${BIN:-build/oc_policy_singlecol}"   # 可执行文件
FREEZE="${FREEZE:-1}"                      # 可行即冻结
CAP="${CAP:-64}"
SCALE="${SCALE:-1000}"
DV="${DV:-1}"
POW2="${POW2:-1}"
PROGRESS="${PROGRESS:-1}"

# --- Lagrange 参数（可被环境覆盖） ---
RHO_LAG="${RHO_LAG:-0.3}"

# --- IRL1（改进 IRL1+）参数（可被环境覆盖；单独的 homotopy，不影响 soft） ---
LAM_IRL="${LAM_IRL:-2.0}"
DELTA_IRL="${DELTA_IRL:-2.0}"
IRL1_HOMOTOPY="${IRL1_HOMOTOPY:-1}"
IRL1_SIGMA="${IRL1_SIGMA:-0.8}"
IRL1_SHRINK="${IRL1_SHRINK:-0.9}"

# --- Soft (DreamPlace + 条件收紧 + 冷启动) 参数（可被环境覆盖） ---
ALPHA_SOFT="${ALPHA_SOFT:-0.9}"
ETA_SOFT="${ETA_SOFT:-1.0}"
GAM_SOFT="${GAM_SOFT:-0.3}"
BET_SOFT="${BET_SOFT:-3.0}"
SOFT_HOMOTOPY="${SOFT_HOMOTOPY:-1}"
SIGMA_SOFT="${SIGMA_SOFT:-0.8}"
SHRINK_SOFT="${SHRINK_SOFT:-0.9}"

# ======================================================
#              结果目录与输出文件
# ======================================================
timestamp=$(date +%Y%m%d_%H%M%S)
RESULT_DIR="results_${M}x${H}_T${T}_cmp_${timestamp}"
LOG_DIR="${RESULT_DIR}/logs"
SUMMARY="${RESULT_DIR}/compare_${M}x${H}_T${T}.txt"
TARBALL="${RESULT_DIR}.tar.gz"

mkdir -p "${LOG_DIR}"
echo "[INFO] 创建结果目录: ${RESULT_DIR}"

# ======================================================
#        输出表头 (maxL 最大/最小/轮次 + HPWL_eq/W + runtime)
# ======================================================
echo "policy,maxL_max,round_of_maxL_max,maxL_min,round_of_maxL_min,HPWL_eq_final,HPWL_w_final,runtime_sec" > "${SUMMARY}"

# ======================================================
#        解析函数（POSIX兼容：直接解析 maxL= / HPWL_*）
# ======================================================
parse_log() {
  local log="$1"
  awk '
  BEGIN {
    nr=0; maxLmax=-1e9; maxLmin=1e9;
    rmax=0; rmin=0; last_he="NA"; last_hw="NA";
  }
  /\[Round[[:space:]]+[0-9]+\]/ {
    line=$0

    # round 号
    r=-1
    tmp=line
    sub(/^.*\[Round[[:space:]]+/, "", tmp)
    split(tmp, ff, /[^0-9]+/)
    if (ff[1] != "") { r = ff[1] + 0 }

    # maxL
    ml_found=0; ml=0
    if (index(line, "maxL=") > 0) {
      split(line, a, "maxL=")
      split(a[2], b, " ")
      ml = b[1] + 0
      ml_found=1
    }

    # HPWL_eq
    he=0
    if (index(line, "HPWL_eq=") > 0) {
      split(line, c, "HPWL_eq=")
      split(c[2], d, " ")
      he = d[1] + 0
      last_he = he
    }

    # HPWL_w
    hw=0
    if (index(line, "HPWL_w=") > 0) {
      split(line, e, "HPWL_w=")
      split(e[2], f, " ")
      hw = f[1] + 0
      last_hw = hw
    }

    if (r >= 0 && ml_found==1) {
      if (nr==0 || ml > maxLmax) { maxLmax = ml; rmax = r; }
      if (nr==0 || ml < maxLmin) { maxLmin = ml; rmin = r; }
      nr++
    }
  }
  END {
    if (nr > 0) {
      printf "%d %d %d %d %s %s\n", maxLmax, rmax, maxLmin, rmin, last_he, last_hw;
    } else {
      exit 1
    }
  }' "${log}"
}

# ======================================================
#   生成 T_prime 演化 TXT（round 行 + tighten 事件，单文件）
#   文件格式：CSV
#   columns:
#   type,round,T_prime,T_target,maxL,WNS_Tprime,WNS_T,picked,HPWL_eq,HPWL_w,solve_ms,old_T_prime,new_T_prime,sigma_maxL,shrink_val
#   - round 行：填充前 11 列，后 4 列留空
#   - tighten 行：填充 type/round 与最后 4 列，中间 9 列留空
# ======================================================
dump_tprime_evolution() {
  local log="$1"
  local outfile="$2"
  {
    echo "# T_prime evolution (M=${M}, H=${H}, T=${T}, rounds=${ROUNDS})"
    echo "type,round,T_prime,T_target,maxL,WNS_Tprime,WNS_T,picked,HPWL_eq,HPWL_w,solve_ms,old_T_prime,new_T_prime,sigma_maxL,shrink_val"
    awk '
      BEGIN {
        last_round="";
        empties=""; for (i=0;i<9;i++) empties = empties ",";
      }
      /\[Round[[:space:]]+[0-9]+\]/ {
        line=$0
        tmp=line
        sub(/^.*\[Round[[:space:]]+/, "", tmp)
        split(tmp, ff, /[^0-9]+/)
        r = ff[1] + 0
        last_round = r

        Tp=""; if (index(line,"T_prime=")>0){ split(line,a,"T_prime="); split(a[2],b," "); Tp=b[1] }
        Tt=""; if (index(line,"T_target=")>0){ split(line,a,"T_target="); split(a[2],b," "); Tt=b[1] }
        ml=""; if (index(line,"maxL=")>0){ split(line,a,"maxL="); split(a[2],b," "); ml=b[1] }
        # WNS(T\047)= 这里绕过单引号困难，先不抓；只抓 WNS*(T)=
        wnsT=""; if (index(line,"WNS*(T)=")>0){ split(line,a,"WNS*(T)="); split(a[2],b," "); wnsT=b[1] }
        pick=""; if (index(line,"picked=")>0){ split(line,a,"picked="); split(a[2],b," "); pick=b[1] }
        he="";   if (index(line,"HPWL_eq=")>0){ split(line,a,"HPWL_eq="); split(a[2],b," "); he=b[1] }
        hw="";   if (index(line,"HPWL_w=")>0){  split(line,a,"HPWL_w=");  split(a[2],b," "); hw=b[1] }
        ms="";   if (index(line,"solve_ms=")>0){ split(line,a,"solve_ms="); split(a[2],b," "); ms=b[1] }

        # 输出 round 行（共 15 列；后 4 列空）
        printf "round,%s,%s,%s,%s,,%s,%s,%s,%s,%s,,,,\n", r,Tp,Tt,ml,wnsT,pick,he,hw,ms
      }
      /\[Tighten\]/ {
        line=$0
        old=""; newv=""; smax=""; shr="";
        if (index(line,"T_prime:")>0){
          split(line, a, "T_prime:");
          rest=a[2]; sub(/^[[:space:]]+/, "", rest);
          # 切 old -> new
          split(rest, b, /[[:space:]]*->[[:space:]]*/)
          old=b[1]
          newpart=b[2]; sub(/^[[:space:]]+/, "", newpart)
          split(newpart, c, /[^0-9]+/); newv=c[1]
          if(index(line,"sigma*maxL=")>0){ split(line,d,"sigma*maxL="); split(d[2],e,/[), ]/); smax=e[1]; }
          if(index(line,"shrink=")>0){ split(line,f,"shrink="); split(f[2],g,/[), ]/); shr=g[1]; }
        }
        r=(last_round==""? "": last_round)
        # 输出 tighten 行：type,round + 9 个空列 + old,new,sigma,shrink
        # empties = ",,,,,,,,,"
        printf "tighten,%s%s%s,%s,%s,%s\n", r, empties, old, newv, smax, shr
      }
    ' "${log}"
  } > "${outfile}"
  echo "[OK] T_prime 演化写入: ${outfile}"
}

# ======================================================
#        运行策略并写入结果（记录 runtime）+ 生成 T_prime 演化
# ======================================================
run_policy() {
  local policy="$1"
  local log="${LOG_DIR}/policy_${policy}.log"
  local out_name="placement${M}_Col_1_T_${T}.txt"
  local out_renamed="${RESULT_DIR}/placement${M}_Col_1_T_${T}_${policy}.txt"
  local tprime_txt="${RESULT_DIR}/tprime_evolution_${policy}.txt"

  echo "===== 运行策略 ${policy} (m=${M},h=${H},T=${T},rounds=${ROUNDS}) ====="

  # runtime 计时（秒）
  local runtime_sec
  SECONDS=0

  if [[ "${policy}" == "lagrange" ]]; then
    "${BIN}" "${M}" "${H}" "${ROUNDS}" \
      --T="${T}" --policy=lagrange --freeze="${FREEZE}" \
      --rho="${RHO_LAG}" \
      --cap="${CAP}" --scale="${SCALE}" \
      --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}" \
      2>&1 | tee "${log}"
    runtime_sec=${SECONDS}

  elif [[ "${policy}" == "irl1" ]]; then
    # IRL1 推荐：开启独立 homotopy，不影响 soft
    "${BIN}" "${M}" "${H}" "${ROUNDS}" \
      --T="${T}" --policy=irl1 --freeze="${FREEZE}" \
      --lam="${LAM_IRL}" --delta="${DELTA_IRL}" \
      --cap="${CAP}" --scale="${SCALE}" \
      --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}" \
      $( [[ "${IRL1_HOMOTOPY}" == "1" ]] && echo --homotopy=1 --sigma="${IRL1_SIGMA}" --shrink="${IRL1_SHRINK}" ) \
      2>&1 | tee "${log}"
    runtime_sec=${SECONDS}

  elif [[ "${policy}" == "soft" ]]; then
    # Soft：保持独立 homotopy 参数
    "${BIN}" "${M}" "${H}" "${ROUNDS}" \
      --T="${T}" --policy=soft --freeze="${FREEZE}" \
      --alpha="${ALPHA_SOFT}" --eta="${ETA_SOFT}" \
      --gamma="${GAM_SOFT}" --beta="${BET_SOFT}" \
      --cap="${CAP}" --scale="${SCALE}" \
      --dV="${DV}" --pow2="${POW2}" --progress="${PROGRESS}" \
      $( [[ "${SOFT_HOMOTOPY}" == "1" ]] && echo --homotopy=1 --sigma="${SIGMA_SOFT}" --shrink="${SHRINK_SOFT}" ) \
      2>&1 | tee "${log}"
    runtime_sec=${SECONDS}

  else
    echo "[ERR] 未知策略 ${policy}"; exit 3
  fi

  # 移动 placement 文件（若存在）
  if [[ -f "${out_name}" ]]; then
    mv "${out_name}" "${out_renamed}"
    echo "[OK] 已移动 ${out_renamed}"
  else
    echo "[WARN] 未找到输出文件 (${out_name})，请检查日志：${log}"
  fi

  # 生成 T_prime 演化文件
  dump_tprime_evolution "${log}" "${tprime_txt}"

  # 解析日志，写入对比表
  if RES=$(parse_log "${log}"); then
    set -- ${RES}
    local MAXL_MAX="$1" RMAX="$2" MAXL_MIN="$3" RMIN="$4" HE_LAST="$5" HW_LAST="$6"
    echo "${policy},${MAXL_MAX},${RMAX},${MAXL_MIN},${RMIN},${HE_LAST},${HW_LAST},${runtime_sec}" >> "${SUMMARY}"
  else
    echo "${policy},NA,NA,NA,NA,NA,NA,${runtime_sec}" >> "${SUMMARY}"
  fi
}

# ======================================================
#                 运行三种策略
# ======================================================
if [[ ! -x "${BIN}" ]]; then
  echo "[ERR] 找不到可执行文件: ${BIN} ；请先 make 构建。"
  exit 2
fi

run_policy "lagrange"
run_policy "irl1"
run_policy "soft"

# ======================================================
#                 打包结果
# ======================================================
echo "[INFO] 压缩结果目录 -> ${TARBALL}"
tar -czf "${TARBALL}" "${RESULT_DIR}"
echo "[DONE] 对比表已生成: ${SUMMARY}"
echo "[HINT] 每个策略的 T_prime 演化文件位于: ${RESULT_DIR}/tprime_evolution_<policy>.txt"
