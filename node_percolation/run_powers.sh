#!/bin/bash
set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: run_powers.sh [options] MIN_EXP MAX_EXP [extra node_percolation args]

Options (override defaults used for every run):
  --strategy STR       Ranking strategy (default: ci)
  --avg-degree C       ER average degree c (default: 4)
  --L-min L            Minimum collective influence radius L (default: 0)
  --L-max L            Maximum collective influence radius (default: same as --L-min)
  --m-graphs NUM       Number of graphs per size (default: auto: 10k for 2^10<=N<2^17, 1k for 2^17<=2^22, else 1k)
  --m-realizations NUM Number of realizations per graph (default: 1)
  --threads NUM        OpenMP threads (default: auto)
  --output PREFIX      Output prefix passed via -o (default: NR)
  --extra "ARGS..."    Additional arguments appended to every run

All remaining arguments after MAX_EXP are forwarded verbatim.
EOF
}

STRATEGY="ci"
AVG_DEGREE=4
CI_RADIUS_MIN=0
CI_RADIUS_MAX=
DEFAULT_M_GRAPHS=1000
M_GRAPHS=$DEFAULT_M_GRAPHS
M_GRAPHS_OVERRIDE=0
M_REALIZATIONS=1
THREADS=0
OUTPUT_PREFIX="NR"
declare -a GLOBAL_EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --strategy)
      STRATEGY="$2"; shift 2 ;;
    --avg-degree)
      AVG_DEGREE="$2"; shift 2 ;;
    --L-min|--ci-radius|--L)
      CI_RADIUS_MIN="$2"
      CI_RADIUS_MAX="$2"
      shift 2 ;;
    --L-max|--ci-max)
      CI_RADIUS_MAX="$2"; shift 2 ;;
    --m-graphs)
      M_GRAPHS="$2"
      M_GRAPHS_OVERRIDE=1
      shift 2 ;;
    --m-realizations)
      M_REALIZATIONS="$2"; shift 2 ;;
    --threads)
      THREADS="$2"; shift 2 ;;
    --output)
      OUTPUT_PREFIX="$2"; shift 2 ;;
    --extra)
      GLOBAL_EXTRA+=("$2"); shift 2 ;;
    --help)
      show_usage; exit 0 ;;
    --*)
      echo "Unknown option: $1" >&2
      exit 1 ;;
    *)
      break ;;
  esac
done

if [[ $# -lt 2 ]]; then
  show_usage >&2
  exit 1
fi

MIN_EXP=$1
MAX_EXP=$2
shift 2

if [[ "$MIN_EXP" -gt "$MAX_EXP" ]]; then
  echo "Error: MIN_EXP must be <= MAX_EXP" >&2
  exit 1
fi

if [[ -z "${CI_RADIUS_MAX:-}" ]]; then
  CI_RADIUS_MAX="$CI_RADIUS_MIN"
fi

if (( CI_RADIUS_MIN > CI_RADIUS_MAX )); then
  echo "Error: --L-min must be <= --L-max" >&2
  exit 1
fi

BIN="$(dirname "$0")/node_percolation"

for ((pow=MIN_EXP; pow<=MAX_EXP; pow++)); do
  N=$((1 << pow))
  if (( M_GRAPHS_OVERRIDE )); then
    RUN_M_GRAPHS="$M_GRAPHS"
  else
    if (( pow >= 10 && pow < 25 )); then
      RUN_M_GRAPHS=10000
    elif (( pow >= 25 )); then
      RUN_M_GRAPHS=1000
    else
      RUN_M_GRAPHS="$DEFAULT_M_GRAPHS"
    fi
  fi
  for ((L=CI_RADIUS_MIN; L<=CI_RADIUS_MAX; L++)); do
    echo "== Running N=$N (2^$pow), L=$L =="
    CMD=("$BIN" -g ER -N "$N" -c "$AVG_DEGREE" -s "$STRATEGY" -L "$L" \
         -M_graphs "$RUN_M_GRAPHS" -M_realizations "$M_REALIZATIONS" -o "$OUTPUT_PREFIX")
    if [[ "$THREADS" -gt 0 ]]; then
      CMD+=(-n_threads "$THREADS")
    fi
    if ((${#GLOBAL_EXTRA[@]})); then
      CMD+=("${GLOBAL_EXTRA[@]}")
    fi
    CMD+=("$@")
    "${CMD[@]}"
  done
done
