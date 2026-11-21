#!/bin/bash
# Which nodes can the 16-rank jobs run on?  (Slurm only; run on the login node.)
#   bash hw.sh [partition]                 node types of the partition: nodes, nodes with free CPUs, CPUs per node, memory (MB), features
#   bash hw.sh model <feature> [partition] CPU model of a node with that feature (runs `lscpu` as a 1-task job)
#   bash hw.sh test <feature> [partition]  estimated start of a 16-rank job on such nodes (sbatch --test-only, nothing is submitted)
# then:  HWCONSTRAINT=<feature> ARRAY=3-5 UPS="150 250 400" bash submit_ladder.sh      (HWEXCL=1: whole node; HWEXCLUDE=<node,...>)
set -u
mode=list; [ "${1:-}" = model ] || [ "${1:-}" = test ] && { mode=$1; shift; }
if [ "$mode" = list ]; then
  P=${1:-general}
  echo "partition $P:  nodes  free  | cpus | memory MB | features"
  sinfo -h -N -p "$P" -o "%f|%c|%m|%t" | awk -F'|' '{ k = $2 " | " $3 " | " $1; tot[k]++; if ($4 ~ /^(idle|mix)/) free[k]++ }
    END { for (k in tot) printf "  %5d %5d  | %s\n", tot[k], free[k] + 0, k }' | sort -k5,5 -k7
  echo "feature names in the last column; a node type is picked with HWCONSTRAINT=<one feature>. 'bash hw.sh model <feature>' shows its CPU."
  exit 0
fi
F=${1:?feature}; P=${2:-general}
if [ "$mode" = model ]; then
  srun -p "$P" -n 1 --constraint="$F" --time=2 bash -c 'echo "$(hostname -s): $(lscpu | sed -n "s/^Model name: *//p" | head -1), $(lscpu | sed -n "s/^Thread(s) per core: *//p") thread(s) per core"'
else
  sbatch --test-only -p "$P" -N 1 -n "${NP:-16}" --mem="${MEM:-48G}" --hint=nomultithread --constraint="$F" --time=24:00:00 --wrap=true
fi
