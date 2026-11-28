#!/usr/bin/env bash
# Memory and CPU monitor for a UCloud job, with a safety stop. Start once per job, in the background:
#   nohup bash "/work/YoavLevy#3923/riseup_preproccesing/python/Yoav/scripts/monitor_job.sh" >/dev/null 2>&1 &
# Watch it:   tail -f "/work/YoavLevy#3923/logs/job_monitor.log"
#
# Reads the job's own counters (Linux cgroup v2, the same numbers as UCloud's usage chart):
#   memory.max      the job's memory limit; above it the kernel kills processes (or the whole job);
#   memory.stat     anon = memory the programs really hold (page cache is left out: the kernel frees it when needed);
#   cpu.stat/max    CPU time used and the number of CPUs of the job.
# Every INTERVAL s one line: time, anon memory (GB), limit, share of the limit, CPU use (in CPUs), OOM kills so far.
# Above WARN_GB a warning; above STOP_GB the largest python/R process is stopped (SIGTERM, then SIGKILL after 10 s),
# so that one runaway worker does not bring the whole job down. Settings via the environment, e.g. STOP_GB=40.
set -u
INTERVAL=${INTERVAL:-10}
WARN_GB=${WARN_GB:-38}
STOP_GB=${STOP_GB:-44}
LOG=${LOG:-"/work/YoavLevy#3923/logs/job_monitor.log"}
CG=/sys/fs/cgroup
mkdir -p "$(dirname "$LOG")"

limit=$(cat $CG/memory.max)
read -r quota period < $CG/cpu.max
cpus=$(( quota == 0 ? $(nproc) : quota / period ))
echo "$(date '+%F %T') start: limit $((limit / 1000000000)) GB, $cpus CPUs, warn $WARN_GB GB, stop $STOP_GB GB, pid $$" >> "$LOG"

prev_usage=$(awk '/^usage_usec/ {print $2}' $CG/cpu.stat)
prev_time=$(date +%s%N)
while true; do
  sleep "$INTERVAL"
  anon=$(awk '/^anon / {print $2}' $CG/memory.stat)
  usage=$(awk '/^usage_usec/ {print $2}' $CG/cpu.stat)
  now=$(date +%s%N)
  ooms=$(awk '/^oom_kill / {print $2}' $CG/memory.events)
  cpu_used=$(awk -v u=$((usage - prev_usage)) -v t=$(((now - prev_time) / 1000)) 'BEGIN {printf "%.1f", u / t}')
  prev_usage=$usage; prev_time=$now
  gb=$(awk -v a="$anon" 'BEGIN {printf "%.1f", a / 1e9}')
  pct=$(awk -v a="$anon" -v l="$limit" 'BEGIN {printf "%.0f", 100 * a / l}')
  line="$(date '+%F %T') mem ${gb}/$((limit / 1000000000)) GB (${pct}%) cpu ${cpu_used}/${cpus} oom_kills ${ooms}"
  if awk -v g="$gb" -v s="$STOP_GB" 'BEGIN {exit !(g >= s)}'; then
    # largest process by resident memory among python and R
    read -r pid rss cmd < <(ps -eo pid=,rss=,comm= --sort=-rss | awk '$3 ~ /^(python|python3|R|Rscript)/ {print; exit}')
    if [ -n "${pid:-}" ]; then
      line="$line STOP: pid $pid ($cmd, $((rss / 1000000)) GB) stopped"
      kill -TERM "$pid" 2>/dev/null
      ( sleep 10; kill -KILL "$pid" 2>/dev/null ) &
    fi
  elif awk -v g="$gb" -v w="$WARN_GB" 'BEGIN {exit !(g >= w)}'; then
    line="$line WARNING above ${WARN_GB} GB"
  fi
  echo "$line" >> "$LOG"
done
