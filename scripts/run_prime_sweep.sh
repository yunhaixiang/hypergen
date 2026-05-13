#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/cpp/build-release/hyperelliptic_cpp}"
OUT_ROOT="${OUT_ROOT:-$ROOT/results/prime_sweep}"
MAX_SPARSITY="${MAX_SPARSITY:-1}"
PROGRESS_INTERVAL="${PROGRESS_INTERVAL:-10000}"
IRREDUCIBLE_MEMORY_BUDGET_MB="${IRREDUCIBLE_MEMORY_BUDGET_MB:-1024}"
EXTRA_MAX_PRIME="${EXTRA_MAX_PRIME:-43}"
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

is_prime() {
  local n="$1"
  if (( n < 2 )); then return 1; fi
  if (( n == 2 )); then return 0; fi
  if (( n % 2 == 0 )); then return 1; fi
  local d=3
  while (( d * d <= n )); do
    if (( n % d == 0 )); then return 1; fi
    d=$((d + 2))
  done
  return 0
}

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
  local p="$1"
  local g="$2"
  local dir="$OUT_ROOT/p${p}_enumerate"
  local out="$dir/p${p}_g${g}_${SPARSITY_LABEL}.sqlite"
  local log_dir="$OUT_ROOT/logs"
  local log="$log_dir/p${p}_g${g}_${SPARSITY_LABEL}.log"

  mkdir -p "$dir" "$log_dir"

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

  echo "run: p=$p genus=$g out=$out"
  "$BIN" \
    --p "$p" \
    --genus "$g" \
    --enumeration-mode enumerate \
    --max-sparsity "$MAX_SPARSITY" \
    --progress-interval "$PROGRESS_INTERVAL" \
    --irreducible-memory-budget-mb "$IRREDUCIBLE_MEMORY_BUDGET_MB" \
    --out "$out" \
    2>&1 | tee "$log"
}

run_range() {
  local p="$1"
  local g_start="$2"
  local g_end="$3"
  local g
  for ((g = g_start; g <= g_end; ++g)); do
    run_one "$p" "$g"
  done
}

# Fixed requested ranges. Genus 1 is skipped.
run_range 5 2 5
run_range 7 2 4
run_range 11 2 3

# Extra prime sweep: p=13,17,19,... for genus 2 only.
# Set EXTRA_MAX_PRIME=0 to keep going until interrupted.
p=13
while true; do
  if (( EXTRA_MAX_PRIME > 0 && p > EXTRA_MAX_PRIME )); then
    break
  fi
  if is_prime "$p"; then
    run_one "$p" 2
  fi
  p=$((p + 2))
done
