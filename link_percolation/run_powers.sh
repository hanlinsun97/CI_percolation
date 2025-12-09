#!/bin/bash
set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: run_powers.sh [options] MIN_EXP MAX_EXP [extra link_percolation args]

Options (override defaults used for every run):
  --avg-degree C       ER average degree c (default: 4)
  --m-graphs NUM       Number of graphs per size (default: auto: 10k for 2^10<=N<2^17, 1k otherwise)
  --m-realizations NUM Number of realizations per graph (default: 1)
  --threads NUM        OpenMP threads (default: auto)
  --output PREFIX      Output prefix passed via -o (default: LR)
  --static-core        Enable static 2-core removal (adds -static_core and suffixes output prefix)
  --random             Remove edges uniformly at random from all edges (adds -random and suffixes output prefix)
  --seed, -r SEED      Random seed passed to binary (default: time)
  --extra "ARGS..."    Additional arguments appended to every run

All remaining arguments after MAX_EXP are forwarded verbatim.
EOF
}

AVG_DEGREE=4
DEFAULT_M_GRAPHS=10
M_GRAPHS=$DEFAULT_M_GRAPHS
M_GRAPHS_OVERRIDE=0
M_REALIZATIONS=1000
THREADS=0
OUTPUT_PREFIX="LR"
SEED=""
STATIC_CORE=0
RANDOM_MODE=0
declare -a GLOBAL_EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --avg-degree)
      AVG_DEGREE="$2"; shift 2 ;;
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
    --static-core)
      STATIC_CORE=1; shift 1 ;;
    --random)
      RANDOM_MODE=1; shift 1 ;;
    --seed|-r)
      SEED="$2"; shift 2 ;;
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

# Adjust output prefix when static-core mode is on to disambiguate files.
if (( STATIC_CORE )) && [[ "$OUTPUT_PREFIX" != *"_stat" ]]; then
  OUTPUT_PREFIX="${OUTPUT_PREFIX}_stat"
fi
if (( RANDOM_MODE )) && [[ "$OUTPUT_PREFIX" != *"_rand" ]]; then
  OUTPUT_PREFIX="${OUTPUT_PREFIX}_rand"
fi

BIN="$(dirname "$0")/link_percolation"

re_num='^[0-9]+$'
if ! [[ "$MIN_EXP" =~ $re_num && "$MAX_EXP" =~ $re_num ]]; then
  echo "Error: MIN_EXP/MAX_EXP must be integers" >&2
  exit 1
fi

for ((pow=MIN_EXP; pow<=MAX_EXP; pow++)); do
  N=$((1 << pow))
  if (( M_GRAPHS_OVERRIDE )); then
    RUN_M_GRAPHS="$M_GRAPHS"
  else
    if (( pow >= 10 && pow < 17 )); then
      RUN_M_GRAPHS=1000
    else
      RUN_M_GRAPHS="$DEFAULT_M_GRAPHS"
    fi
  fi

  echo "== Running N=$N (2^$pow) =="
  CMD=("$BIN" -g ER -N "$N" -c "$AVG_DEGREE" \
       -M_graphs "$RUN_M_GRAPHS" -M_realizations "$M_REALIZATIONS" -o "$OUTPUT_PREFIX")
  if (( STATIC_CORE )); then
    CMD+=(-static_core)
  fi
  if (( RANDOM_MODE )); then
    CMD+=(-random)
  fi
  if [[ -n "$SEED" ]]; then
    CMD+=(-r "$SEED")
  fi
  if [[ "$THREADS" -gt 0 ]]; then
    CMD+=(-n_threads "$THREADS")
  fi
  if ((${#GLOBAL_EXTRA[@]})); then
    CMD+=("${GLOBAL_EXTRA[@]}")
  fi
  # Append user trailing args (after MAX_EXP)
  CMD+=("$@")
  "${CMD[@]}"
done
