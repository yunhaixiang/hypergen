# Hypergen

`hypergen` is the C++ branch-divisor search runner for generating sparse hyperelliptic-curve data. The current repository layout is flat:

```text
cpp/                              C++ source and CMake file
legacy/python_prototype/           frozen old Python prototype, reference only
results/p3_enumerate_batch/        generated SQLite data through genus 8
```

The Python prototype is kept only as a legacy snapshot and is not updated. The
active implementation is the C++ runner under `cpp/`.

## Build

Requirements:

- CMake
- `pkg-config`
- FLINT
- SQLite3

From the repository root:

Use a release build for actual runs:

```bash
cmake -S cpp -B cpp/build-release -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build-release -j
```

The release executable is:

```text
cpp/build-release/hyperelliptic_cpp
```

## Run

Deterministic branch-divisor enumeration:

```bash
cpp/build-release/hyperelliptic_cpp \
  --p 3 \
  --genus 7 \
  --enumeration-mode enumerate \
  --max-sparsity 1 \
  --out results/p3_enumerate_batch/p3_g7_s_1.sqlite
```

Batch run over a genus range:

```bash
cpp/build-release/hyperelliptic_cpp \
  --p 3 \
  --genus-start 7 \
  --genus-end 10 \
  --enumeration-mode enumerate \
  --max-sparsity 1 \
  --out-dir results/p3_enumerate_batch
```

Random sparse search:

```bash
cpp/build-release/hyperelliptic_cpp \
  --p 5 \
  --genus 11 \
  --enumeration-mode random \
  --max-sparsity 1 \
  --limit 10000 \
  --out results/p5_g11_s_1_random.sqlite
```

If `--out` is omitted, single-genus runs write under `results/`. If `--out-dir` is omitted, batch runs write under `results/`.

## Modes

`--enumeration-mode enumerate` deterministically enumerates squarefree branch divisors by factorization pattern and irreducible-factor choices. It includes both normalized degree `2g+1` and degree `2g+2` models. Degree `2g+2` branch patterns with an `F_p`-linear factor are skipped because they are represented by odd-degree models.

`--enumeration-mode random` samples random factorized branch divisors for higher-genus sparse search. It first builds all feasible factorization-pattern strata, then adaptively gives more samples to patterns whose previous samples produced sparse curves. The optional `--random-max-factors N` restricts this to patterns with at most `N` irreducible factors; if omitted, there is no hard factor-count cap. Even without a cap, the sampling weights strongly favor patterns with fewer irreducible factors. If `--limit` is omitted, random mode runs until interrupted.

`--limit N` has mode-dependent meaning:

- in `enumerate` mode, stop after processing `N` deterministic candidates
- in `random` mode, sample `N` random candidates

Both modes use Hasse-Witt filtering before exact point counts when `--max-sparsity` is supplied. Exact L-polynomial coefficients are computed by point counting over extension fields and Newton identities.

Random-mode progress includes a `random_patterns` line showing the number of active strata and the currently best observed pattern by sparse-hit rate.

## Memory Budget

Irreducible-factor tables are materialized only while the estimated memory use stays within:

```bash
--irreducible-memory-budget-mb 1024
```

The default is `1024` MB. Higher degrees beyond the budget are streamed in `enumerate` mode or sampled directly in `random` mode.

## Sparsity

For a genus `g` curve, the runner computes the stored part of the L-polynomial

```text
1 + a_1 T + a_2 T^2 + ... + a_g T^g + ...
```

and records only `a_1, ..., a_g`; the remaining coefficients are determined by
Poincare duality. The `sparsity` of a curve is the number of nonzero entries
among `a_1, ..., a_g`.

`--max-sparsity N` keeps only curves with sparsity at most `N`. The value used
for a run is recorded in `enumeration_summary.max_sparsity`, while each stored
curve has its actual value in `sparse_curves.sparsity`. If `--max-sparsity` is
omitted, no sparsity filter is applied and `max_sparsity` is stored as `NULL`.

## Progress Output

The runner prints progress blocks of the form:

```text
prime: p
genus: g
progress: processed/total
progress_percent: P%
elapsed: Dd HHh MMm SSs
estimated_remaining: Dd HHh MMm SSs
sparse_presentations: N
sparse_isomorphism_classes: M
canonicalized_isomorphism_classes: K
-
```

`total` is the actual deterministic branch-divisor presentation count in `enumerate` mode, even if `--limit` is used. In unbounded `random` mode, `total`, `progress_percent`, and `estimated_remaining` are printed as `?`.

`sparse_presentations` is the sum of normalized enumerated orbit sizes for sparse isomorphism classes.

## SQLite Output

Each generated SQLite file is intentionally lean and keeps only sparse output plus run metadata. The tables are:

- `sparse_curves`: sparse survivors, including readable expanded coefficients, readable factorized branch data, exact L-polynomial coefficients, and sparsity.
- `enumeration_summary`: one-row run summary with `prime`, `genus`, `max_sparsity`, `enumeration_mode`, `total_coefficient_vectors`, `processed`, `sparse_presentations`, `sparse_isomorphism_classes`, and timing/status fields.
- `enumeration_progress`: progress snapshots at the print interval.

The output database does not store all candidates, all duplicates, all Hasse-Witt failures, or all orbit-cache rows.

On normal completion or graceful interrupt, the runner writes the summary row,
checkpoints SQLite, switches the database back to `DELETE` journal mode, and
removes transient `*.sqlite-wal` and `*.sqlite-shm` sidecar files.

## Included Data

Current included data is under:

```text
results/p3_enumerate_batch/
```

Summary:

```text
g=2  total=560       sparse_presentations=538   sparse_classes=69
g=3  total=5124      sparse_presentations=1652  sparse_classes=140
g=4  total=46296     sparse_presentations=3789  sparse_classes=340
g=5  total=416972    sparse_presentations=2858  sparse_classes=255
g=6  total=3753216   sparse_presentations=5911  sparse_classes=763
g=7  total=33779604  sparse_presentations=7878  sparse_classes=1022
g=8  total=304017320 sparse_presentations=12554 sparse_classes=1542
```

## Useful Inspection Commands

List sparse rows:

```bash
sqlite3 results/p3_enumerate_batch/p3_g7_s_1.sqlite \
  "SELECT count(*) FROM sparse_curves;"
```

View the run summary:

```bash
sqlite3 results/p3_enumerate_batch/p3_g7_s_1.sqlite \
  "SELECT * FROM enumeration_summary;"
```
