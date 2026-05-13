#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/cpp/build-release/hyperelliptic_cpp}"
OUT_DIR="${OUT_DIR:-$ROOT/results/p3_random}"
MAX_SPARSITY="${MAX_SPARSITY:-1}"
PROGRESS_INTERVAL="${PROGRESS_INTERVAL:-10000}"
IRREDUCIBLE_MEMORY_BUDGET_MB="${IRREDUCIBLE_MEMORY_BUDGET_MB:-1024}"
LIMIT="${LIMIT:-170000000}"
GENUS_START="${GENUS_START:-11}"
GENUS_END="${GENUS_END:-16}"
SKIP_COMPLETE="${SKIP_COMPLETE:-1}"

if [[ ! -x "$BIN" ]]; then
  echo "missing executable: $BIN" >&2
  echo "build it with: cmake -S cpp -B cpp/build-release -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build-release -j" >&2
  exit 1
fi

if [[ "$MAX_SPARSITY" -lt 0 ]]; then
  SPARSITY_LABEL="all"
else
  SPARSITY_LABEL="s_${MAX_SPARSITY}"
fi

is_complete_sqlite() {
  local db="$1"
  [[ -s "$db" ]] || return 1
  command -v sqlite3 >/dev/null 2>&1 || return 1
  local ok
  ok="$(sqlite3 "file:$db?mode=ro" \
    "SELECT CASE WHEN processed >= CAST(total_coefficient_vectors AS INTEGER) THEN 1 ELSE 0 END FROM enumeration_summary WHERE id = 1;" \
    2>/dev/null || true)"
  [[ "$ok" == "1" ]]
}

mkdir -p "$OUT_DIR"

for ((g = GENUS_START; g <= GENUS_END; ++g)); do
  out="$OUT_DIR/p3_g${g}_${SPARSITY_LABEL}_random.sqlite"

  if [[ "$SKIP_COMPLETE" == "1" ]] && is_complete_sqlite "$out"; then
    echo "skip complete: $out"
    continue
  fi

  if [[ -e "$out" ]] && ! is_complete_sqlite "$out"; then
    stamp="$(date +%Y%m%d_%H%M%S)"
    echo "existing incomplete output moved to: ${out}.incomplete.${stamp}"
    mv "$out" "${out}.incomplete.${stamp}"
    rm -f "$out-wal" "$out-shm"
  fi

  echo "run: p=3 genus=$g limit=$LIMIT out=$out"
  "$BIN" \
    --p 3 \
    --genus "$g" \
    --enumeration-mode random \
    --max-sparsity "$MAX_SPARSITY" \
    --limit "$LIMIT" \
    --progress-interval "$PROGRESS_INTERVAL" \
    --random-seed "$g" \
    --irreducible-memory-budget-mb "$IRREDUCIBLE_MEMORY_BUDGET_MB" \
    --out "$out"
done
