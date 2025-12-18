#!/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
RUN_POWERS="$BASE_DIR/run_powers.sh"

# All available removal strategies supported by run_powers.sh
STRATEGIES=(
  --adaptive
  --static-core
  --random
  --adaptive-biased
  --static-biased
  --random-biased
)

show_usage() {
  cat <<EOF
Usage: $(basename "$0") [common options] MIN_EXP MAX_EXP [extra link_percolation args]

Runs \`run_powers.sh\` sequentially for all strategies:
  ${STRATEGIES[*]}

Common options (forwarded to run_powers.sh for every strategy):
  --avg-degree C
  --m-graphs NUM        (default: 1000)
  --m-realizations NUM  (default: 10)
  --threads NUM
  --output PREFIX
  --seed, -r SEED
  --extra "ARGS..."

Example (matches your commands):
  $(basename "$0") --m-graphs 1000 --m-realizations 10 10 24
EOF
}

COMMON_ARGS=()
M_GRAPHS=1000
M_REALIZATIONS=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --avg-degree|--threads|--output|--seed|-r|--extra)
      COMMON_ARGS+=("$1" "$2"); shift 2 ;;
    --m-graphs)
      M_GRAPHS="$2"; shift 2 ;;
    --m-realizations)
      M_REALIZATIONS="$2"; shift 2 ;;
    --help)
      show_usage; exit 0 ;;
    --adaptive|--static-core|--random|--adaptive-biased|--static-biased|--random-biased)
      echo "Error: do not pass a strategy flag; this script runs all strategies." >&2
      exit 1 ;;
    --*)
      echo "Unknown option: $1" >&2
      show_usage >&2
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

# Always set these explicitly (default or user-provided) to keep runs consistent.
COMMON_ARGS+=(--m-graphs "$M_GRAPHS" --m-realizations "$M_REALIZATIONS")

for strategy in "${STRATEGIES[@]}"; do
  echo "== ${strategy} =="
  "$RUN_POWERS" "$strategy" "${COMMON_ARGS[@]}" "$MIN_EXP" "$MAX_EXP" "$@"
done
