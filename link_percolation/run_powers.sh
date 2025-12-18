#!/bin/bash
set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: run_powers.sh [options] MIN_EXP MAX_EXP [extra link_percolation args]

Options (override defaults used for every run):
  --avg-degree C       ER average degree c (default: 4)
  --m-graphs NUM       Number of graphs per size (default: 1000)
  --m-realizations NUM Number of realizations per graph (default: 10)
  --threads NUM        OpenMP threads (default: auto)
  --samples NUM        Output samples in [0,1] (0=full grid, default: binary default)
  --output PREFIX      Output prefix basename (default: LR)
  --adaptive           Adaptive removal (default; random edge in CURRENT 2-core; writes into ./adaptive/)
  --static-core        Static removal (random edge in ORIGINAL 2-core; writes into ./static/; adds -static_core)
  --random             Random removal (random edge among ALL edges; writes into ./random/; adds -random)
  --adaptive-biased    Adaptive biased (max k_i+k_j in CURRENT 2-core; writes into ./adaptive_biased/; adds -adaptive_biased)
  --static-biased      Static biased (max k_i+k_j in ORIGINAL 2-core; writes into ./static_biased/; adds -static_biased)
  --random-biased      Random biased (max deg0_i+deg0_j in ORIGINAL graph; writes into ./random_biased/; adds -random_biased)
  --random-biased-adaptive Random biased adaptive (max k_i+k_j in CURRENT graph; writes into ./random_biased_adaptive/; adds -random_biased_adaptive)
  --seed, -r SEED      Random seed passed to binary (default: time)
  --extra "ARGS..."    Additional arguments appended to every run

All remaining arguments after MAX_EXP are forwarded verbatim.
EOF
}

AVG_DEGREE=4
M_GRAPHS=1000
M_REALIZATIONS=10
THREADS=0
SAMPLES=""
OUTPUT_PREFIX="LR"
SEED=""
STRATEGY="adaptive"
declare -a GLOBAL_EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --avg-degree)
      AVG_DEGREE="$2"; shift 2 ;;
    --m-graphs)
      M_GRAPHS="$2"
      shift 2 ;;
    --m-realizations)
      M_REALIZATIONS="$2"; shift 2 ;;
    --threads)
      THREADS="$2"; shift 2 ;;
    --samples)
      SAMPLES="$2"; shift 2 ;;
    --output)
      OUTPUT_PREFIX="$2"; shift 2 ;;
    --adaptive)
      STRATEGY="adaptive"; shift 1 ;;
    --static-core)
      STRATEGY="static"; shift 1 ;;
    --random)
      STRATEGY="random"; shift 1 ;;
    --adaptive-biased)
      STRATEGY="adaptive_biased"; shift 1 ;;
    --static-biased)
      STRATEGY="static_biased"; shift 1 ;;
    --random-biased)
      STRATEGY="random_biased"; shift 1 ;;
    --random-biased-adaptive)
      STRATEGY="random_biased_adaptive"; shift 1 ;;
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

case "$STRATEGY" in
  adaptive|static|random|adaptive_biased|static_biased|random_biased|random_biased_adaptive) ;;
  *)
    echo "Error: unknown strategy '$STRATEGY'" >&2
    exit 1 ;;
esac

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN="$BASE_DIR/link_percolation"

OUT_DIR="$BASE_DIR/$STRATEGY"
mkdir -p "$OUT_DIR"
OUT_PREFIX="$OUT_DIR/$OUTPUT_PREFIX"

re_num='^[0-9]+$'
if ! [[ "$MIN_EXP" =~ $re_num && "$MAX_EXP" =~ $re_num ]]; then
  echo "Error: MIN_EXP/MAX_EXP must be integers" >&2
  exit 1
fi

for ((pow=MIN_EXP; pow<=MAX_EXP; pow++)); do
  N=$((1 << pow))
  RUN_M_GRAPHS="$M_GRAPHS"

  echo "== Running N=$N (2^$pow) =="
  CMD=("$BIN" -g ER -N "$N" -c "$AVG_DEGREE" \
       -M_graphs "$RUN_M_GRAPHS" -M_realizations "$M_REALIZATIONS" -o "$OUT_PREFIX")
  case "$STRATEGY" in
    adaptive) ;;
    static) CMD+=(-static_core) ;;
    random) CMD+=(-random) ;;
    adaptive_biased) CMD+=(-adaptive_biased) ;;
    static_biased) CMD+=(-static_biased) ;;
    random_biased) CMD+=(-random_biased) ;;
    random_biased_adaptive) CMD+=(-random_biased_adaptive) ;;
  esac
  if [[ -n "$SEED" ]]; then
    CMD+=(-r "$SEED")
  fi
  if [[ "$THREADS" -gt 0 ]]; then
    CMD+=(-n_threads "$THREADS")
  fi
  if [[ -n "$SAMPLES" ]]; then
    if ! [[ "$SAMPLES" =~ ^[0-9]+$ ]]; then
      echo "Error: --samples must be an integer >= 0" >&2
      exit 1
    fi
    CMD+=(-samples "$SAMPLES")
  fi
  if ((${#GLOBAL_EXTRA[@]})); then
    CMD+=("${GLOBAL_EXTRA[@]}")
  fi
  # Append user trailing args (after MAX_EXP)
  CMD+=("$@")
  "${CMD[@]}"
done
