# Hypergen

`hypergen` is the C++ branch-divisor search runner for generating sparse hyperelliptic-curve data. The current repository layout is flat:

```text
cpp/                              C++ source and CMake file
results/p3_enumerate_batch/        generated SQLite data through genus 7
```

The Python prototype is not part of this repository. This copy is intended to be the C++ runner plus generated data.

## Build

Requirements:

- CMake
- `pkg-config`
- FLINT
- SQLite3

From the repository root:

```bash
cmake -S cpp -B cpp/build
cmake --build cpp/build -j
```

The executable is:

```text
cpp/build/hyperelliptic_cpp
```

## Run

Deterministic branch-divisor enumeration:

```bash
cpp/build/hyperelliptic_cpp \
  --p 3 \
  --genus 7 \
  --enumeration-mode enumerate \
  --max-sparsity 1 \
  --out results/p3_enumerate_batch/p3_g7_s_1.sqlite
```

Batch run over a genus range:

```bash
cpp/build/hyperelliptic_cpp \
  --p 3 \
  --genus-start 7 \
  --genus-end 10 \
  --enumeration-mode enumerate \
  --max-sparsity 1 \
  --out-dir results/p3_enumerate_batch
```

Random sparse search:

```bash
cpp/build/hyperelliptic_cpp \
  --p 5 \
  --genus 11 \
  --enumeration-mode random \
  --max-sparsity 1 \
  --random-steps 10000 \
  --out results/p5_g11_s_1_random.sqlite
```

If `--out` is omitted, single-genus runs write under `results/`. If `--out-dir` is omitted, batch runs write under `results/`.

## Modes

`--enumeration-mode enumerate` deterministically enumerates squarefree branch divisors by factorization pattern and irreducible-factor choices. It includes both normalized degree `2g+1` and degree `2g+2` models. Degree `2g+2` branch patterns with an `F_p`-linear factor are skipped because they are represented by odd-degree models.

`--enumeration-mode random` samples random factorized branch divisors for higher-genus sparse search. If neither `--random-steps` nor `--limit` is supplied, it runs until interrupted.

Both modes use Hasse-Witt filtering before exact point counts when `--max-sparsity` is supplied. Exact L-polynomial coefficients are computed by point counting over extension fields and Newton identities.

## Memory Budget

Irreducible-factor tables are materialized only while the estimated memory use stays within:

```bash
--irreducible-memory-budget-mb 1024
```

The default is `1024` MB. Higher degrees beyond the budget are streamed in `enumerate` mode or sampled directly in `random` mode.

## Progress Output

The runner prints progress blocks of the form:

```text
prime: p
genus: g
progress: processed/total
sparse_presentations: N
sparse_isomorphism_classes: M
canonicalized_isomorphism_classes: K
-
```

`total` is the actual deterministic branch-divisor presentation count in `enumerate` mode, even if `--limit` is used.

`sparse_presentations` is the sum of normalized enumerated orbit sizes for sparse isomorphism classes.

## SQLite Output

Each generated SQLite file is intentionally lean and keeps only sparse output plus run metadata. The tables are:

- `sparse_curves`: sparse survivors, including readable expanded coefficients, readable factorized branch data, exact L-polynomial coefficients, and sparsity.
- `enumeration_summary`: one-row run summary with `prime`, `genus`, `max_sparsity`, `enumeration_mode`, `total_coefficient_vectors`, `processed`, `sparse_presentations`, `sparse_isomorphism_classes`, and timing/status fields.
- `enumeration_progress`: progress snapshots at the print interval.

The output database does not store all candidates, all duplicates, all Hasse-Witt failures, or all orbit-cache rows.

## Included Data

Current included data is under:

```text
results/p3_enumerate_batch/
```

Summary:

```text
g=1  total=60        sparse_presentations=60    sparse_classes=12
g=2  total=560       sparse_presentations=538   sparse_classes=69
g=3  total=5124      sparse_presentations=1652  sparse_classes=140
g=4  total=46296     sparse_presentations=3789  sparse_classes=340
g=5  total=416972    sparse_presentations=2858  sparse_classes=255
g=6  total=3753216   sparse_presentations=5911  sparse_classes=763
g=7  total=33779604  sparse_presentations=7878  sparse_classes=1022
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
