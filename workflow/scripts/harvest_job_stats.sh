#!/bin/bash
# Harvest per-rule Slurm resource stats for a pipeline run.
#
# WHY. Two rules (abc_call_macs_peaks, abc_create_neighborhoods) killed run
# 43001627 by timing out at the profile's flat runtime="6h" on the eight
# largest catlas-human clusters -- and because retries do NOT scale walltime
# (verified: all 3 attempts got exactly 06:00:00), each doomed cluster burned
# ~18 h. Nothing in the pipeline recorded how long rules actually take or how
# much memory they actually use, so the right limits could only be guessed.
# This closes that gap: every completed run leaves a per-rule profile behind.
#
# Reads ONLY sacct plus the slurm executor's own log tree
#   .snakemake/slurm_logs/rule_<rule>/<wildcards>/<jobid>.log
# so it needs no pipeline change and can run during or after a run.
#
# USAGE
#   workflow/scripts/harvest_job_stats.sh [outfile.tsv]
#
# Emits one row per Slurm job (rule, wildcards, jobid, state, elapsed, ...)
# and prints a per-rule summary to stderr.

set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOGROOT="$REPO/.snakemake/slurm_logs"
OUT="${1:-$REPO/resources/job_stats_$(date +%Y%m%d_%H%M%S).tsv}"

[ -d "$LOGROOT" ] || { echo "no $LOGROOT -- has a slurm run happened here?" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"

# rule + wildcards + jobid from the log tree
map=$(mktemp); trap 'rm -f "$map"' EXIT
find "$LOGROOT" -name "*.log" -printf "%p\n" 2>/dev/null | while read -r p; do
    rel="${p#$LOGROOT/}"                 # rule_<rule>/<wildcards>/<jobid>.log
    rule="${rel%%/*}"; rule="${rule#rule_}"
    rest="${rel#*/}"; jid="$(basename "$rest" .log)"; wc_="${rest%/*}"
    [[ "$jid" =~ ^[0-9]+$ ]] && printf "%s\t%s\t%s\n" "$jid" "$rule" "$wc_"
done | sort -u > "$map"

n=$(wc -l < "$map")
echo "[harvest] $n slurm jobs found in the log tree" >&2
[ "$n" -eq 0 ] && { echo "[harvest] nothing to do" >&2; exit 0; }

# sacct in chunks -- the id list can be long
printf "rule\twildcards\tjobid\tstate\telapsed_s\ttimelimit_s\tmaxrss_mb\treqmem\tcpus\tpartition\tnode\n" > "$OUT"
cut -f1 "$map" | paste -sd, - | tr ',' '\n' | split -l 400 - "${map}.chunk."
for c in "${map}".chunk.*; do
    ids=$(paste -sd, "$c")
    # -X would hide MaxRSS (it lives on the .batch/step rows), so take all rows
    # and fold the per-step max back onto the parent job id.
    sacct -j "$ids" -n -P -o JobID,State,ElapsedRaw,TimelimitRaw,MaxRSS,ReqMem,AllocCPUS,Partition,NodeList 2>/dev/null
done | awk -F'|' -v mapf="$map" '
  BEGIN { while ((getline line < mapf) > 0) { split(line,a,"\t"); rule[a[1]]=a[2]; wc[a[1]]=a[3] } }
  {
    id=$1; sub(/\..*$/,"",id)                 # 123.batch -> 123
    if (!(id in rule)) next
    if ($1 !~ /\./) { st[id]=$2; el[id]=$3+0; tl[id]=$4*60; rm[id]=$6; cp[id]=$7; pa[id]=$8; nd[id]=$9 }
    m=$5; if (m ~ /K$/) { sub(/K$/,"",m); m=m/1024 } else if (m ~ /M$/) { sub(/M$/,"",m); m=m+0 }
      else if (m ~ /G$/) { sub(/G$/,"",m); m=m*1024 } else m=0
    if (m > mx[id]) mx[id]=m
  }
  END { for (id in st) printf "%s\t%s\t%s\t%s\t%d\t%d\t%.1f\t%s\t%s\t%s\t%s\n",
          rule[id], wc[id], id, st[id], el[id], tl[id], mx[id], rm[id], cp[id], pa[id], nd[id] }
' | sort -k1,1 -k5,5nr >> "$OUT"
rm -f "${map}".chunk.*

echo "[harvest] wrote $OUT" >&2
echo >&2
awk -F'\t' 'NR>1 {
    n[$1]++; if($5+0>mxel[$1]) mxel[$1]=$5+0; tot[$1]+=$5+0
    if($7+0>mxrss[$1]) mxrss[$1]=$7+0
    if($4 ~ /TIMEOUT/) to[$1]++
    if($4 ~ /FAILED|CANCELLED|NODE_FAIL|OUT_OF_ME/) bad[$1]++
    if($6+0>lim[$1]) lim[$1]=$6+0
  }
  END {
    printf "%-46s %5s %10s %10s %10s %8s %7s\n","rule","jobs","max_h","mean_h","limit_h","maxRSS_G","TMOUT"
    for (r in n) printf "%-46s %5d %10.2f %10.2f %10.2f %8.1f %7d\n",
        r, n[r], mxel[r]/3600, tot[r]/n[r]/3600, lim[r]/3600, mxrss[r]/1024, to[r]+0
  }' "$OUT" | { read -r h; echo "$h"; sort -k3,3nr; } >&2
