#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/cpp/build-release/hyperelliptic_cpp}"
OUT_ROOT="${OUT_ROOT:-$ROOT/results}"
MAX_SPARSITY="${MAX_SPARSITY:-1}"
PROGRESS_INTERVAL="${PROGRESS_INTERVAL:-100000}"
IRREDUCIBLE_MEMORY_BUDGET_MB="${IRREDUCIBLE_MEMORY_BUDGET_MB:-1024}"
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
  ok="$(sqlite3 "file:$db?mode=ro&immutable=1" \
    "SELECT CASE WHEN processed = CAST(total_coefficient_vectors AS INTEGER) THEN 1 ELSE 0 END FROM enumeration_summary WHERE id = 1;" \
    2>/dev/null || true)"
  [[ "$ok" == "1" ]]
}

run_one() {
  local g="$1"
  local dir="$OUT_ROOT/p3_enumerate"
  local out="$dir/p3_g${g}_${SPARSITY_LABEL}.sqlite"

  mkdir -p "$dir"

  if [[ "$SKIP_COMPLETE" == "1" ]] && is_complete_sqlite "$out"; then
    echo "skip complete: $out"
    return 0
  fi

  if [[ -e "$out" ]] && ! is_complete_sqlite "$out"; then
    local stamp
    stamp="$(date +%Y%m%d_%H%M%S)"
    echo "existing incomplete output moved to: ${out}.incomplete.${stamp}"
    mv "$out" "${out}.incomplete.${stamp}"
    rm -f "$out-wal" "$out-shm"
  fi

  echo "run: p=3 genus=$g out=$out"
  "$BIN" \
    --p 3 \
    --genus "$g" \
    --enumeration-mode enumerate \
    --max-sparsity "$MAX_SPARSITY" \
    --progress-interval "$PROGRESS_INTERVAL" \
    --irreducible-memory-budget-mb "$IRREDUCIBLE_MEMORY_BUDGET_MB" \
    --out "$out"
}

for g in 9 10 11 12 13; do
  run_one "$g"
done
