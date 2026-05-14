# Hypergen

`hypergen` is the C++ branch-divisor search runner for generating sparse hyperelliptic-curve data. The current repository layout is flat:

```text
cpp/                              C++ source and CMake file
legacy/python_prototype/           frozen old Python prototype, reference only
results/p*_enumerate/              deterministic branch-divisor data by prime
results/p*_random/                 random sparse-search data by prime
scripts/                           batch-run helper scripts
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
  --out results/p3_enumerate/p3_g7_s_1.sqlite
```

Batch run over a genus range:

```bash
cpp/build-release/hyperelliptic_cpp \
  --p 3 \
  --genus-start 7 \
  --genus-end 10 \
  --enumeration-mode enumerate \
  --max-sparsity 1 \
  --out-dir results/p3_enumerate
```

Random sparse search:

```bash
cpp/build-release/hyperelliptic_cpp \
  --p 5 \
  --genus 11 \
  --enumeration-mode random \
  --max-sparsity 1 \
  --limit 10000 \
  --out results/p5_random/p5_g11_s_1_random.sqlite
```

If `--out` is omitted, single-genus runs write under `results/`. If `--out-dir` is omitted, batch runs write under `results/`.

## Modes

`--enumeration-mode enumerate` deterministically enumerates squarefree branch divisors by branch divisor type and irreducible-factor choices. It includes both normalized degree `2g+1` and degree `2g+2` models. Degree `2g+2` branch types with an `F_p`-linear factor are skipped because they are represented by odd-degree models.

`--enumeration-mode random` samples random factorized branch divisors for higher-genus sparse search. It first builds all feasible branch-divisor-type strata, then adaptively gives more samples to types whose previous samples produced sparse curves. The optional `--random-max-factors N` restricts this to types with at most `N` irreducible finite factors; if omitted, there is no hard factor-count cap. Even without a cap, the sampling weights strongly favor types with fewer irreducible factors. If `--limit` is omitted, random mode runs until interrupted.

Before either mode starts choosing actual irreducible factors, it precomputes
the normalized branch divisor types that can occur over `F_p`. A type is kept
only if `F_p` has enough irreducible polynomials of each requested degree. When
`--max-sparsity` is supplied, the type list is also filtered using the mod-2
branch divisor theorem:

```text
L_C(T) = product_i (1 + T^{d_i}) / (1 + T)^2 mod 2
```

Here `d_i` are the Frobenius orbit degrees of the branch points; the stored
infinity marker `0` counts as degree `1` for this formula. If this mod-2
polynomial already has more than `max_sparsity` nonzero coefficients among
`a_1,...,a_{g-1}`, then no integral L-polynomial with that branch divisor type
can satisfy the sparsity bound, so the entire type is skipped.

When no exact sparse curves have been found yet, random mode also uses
Hasse-Witt pass rates as a proxy signal for choosing branch divisor types.
This matters in higher genus, where exact sparse hits can be rare enough that
the sampler otherwise has no early feedback.

`--limit N` has mode-dependent meaning:

- in `enumerate` mode, stop after processing `N` deterministic candidates
- in `random` mode, sample `N` random candidates

Both modes use Hasse-Witt filtering before exact point counts when `--max-sparsity` is supplied. Exact L-polynomial coefficients are computed by point counting over extension fields and Newton identities.

Random-mode progress includes a `random_patterns` line showing the number of active strata and the currently best observed pattern by sparse-hit rate, breaking ties by Hasse-Witt pass rate.

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
among `a_1, ..., a_{g-1}`. The middle coefficient `a_g` is still stored in the
L-polynomial data, but it is not counted toward sparsity.

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

`total` is the actual deterministic branch-divisor presentation count in `enumerate` mode after feasibility and mod-2 branch-type filtering, even if `--limit` is used. In unbounded `random` mode, `total`, `progress_percent`, and `estimated_remaining` are printed as `?`.

`sparse_presentations` is the sum of normalized enumerated orbit sizes for sparse isomorphism classes.

## SQLite Output

Each generated SQLite file is intentionally lean and keeps only sparse output plus run metadata. The tables are:

- `sparse_curves`: sparse survivors, including readable expanded coefficients, readable factorized branch data, exact L-polynomial coefficients, and sparsity.
- `enumeration_summary`: one-row run summary with `prime`, `genus`, `max_sparsity`, `enumeration_mode`, `total_coefficient_vectors`, `processed`, `sparse_presentations`, `sparse_isomorphism_classes`, and timing/status fields.
- `enumeration_progress`: progress snapshots at the print interval.

`sparse_curves.branch_factorization_pattern` is stored as a partition of the
branch polynomial degree. For example, `[1,11,11]` means one linear factor and
two irreducible degree-11 factors.

`sparse_curves.branch_divisor_type` stores the full branch divisor type, using
`0` for the branch point at infinity. For example, `[0,1,11,11]` means infinity,
one linear factor, and two irreducible degree-11 factors.

The output database does not store all candidates, all duplicates, all Hasse-Witt failures, or all orbit-cache rows.

On normal completion or graceful interrupt, the runner writes the summary row,
checkpoints SQLite, switches the database back to `DELETE` journal mode, and
removes transient `*.sqlite-wal` and `*.sqlite-shm` sidecar files.

## Included Data

Current generated data is under:

```text
results/p3_enumerate/
results/p3_random/
results/p5_enumerate/
results/p5_random/
results/p7_enumerate/
results/p11_enumerate/
results/p13_enumerate/
results/p17_enumerate/
```

All rows below use `max_sparsity = 1`.

Complete deterministic `enumerate` outputs:

| prime | genus | file | total processed | sparse presentations | sparse classes |
| --- | ---: | --- | ---: | ---: | ---: |
| 3 | 2 | `results/p3_enumerate/p3_g2_s_1.sqlite` | 560 | 538 | 69 |
| 3 | 3 | `results/p3_enumerate/p3_g3_s_1.sqlite` | 5,124 | 1,652 | 140 |
| 3 | 4 | `results/p3_enumerate/p3_g4_s_1.sqlite` | 46,296 | 3,789 | 340 |
| 3 | 5 | `results/p3_enumerate/p3_g5_s_1.sqlite` | 416,972 | 2,858 | 255 |
| 3 | 7 | `results/p3_enumerate/p3_g7_s_1.sqlite` | 33,779,604 | 7,878 | 1,022 |
| 3 | 8 | `results/p3_enumerate/p3_g8_s_1.sqlite` | 304,017,320 | 12,554 | 1,542 |
| 5 | 2 | `results/p5_enumerate/p5_g2_s_1.sqlite` | 12,460 | 11,855 | 285 |
| 5 | 3 | `results/p5_enumerate/p5_g3_s_1.sqlite` | 313,390 | 88,678 | 1,378 |
| 5 | 4 | `results/p5_enumerate/p5_g4_s_1.sqlite` | 7,841,152 | 363,490 | 6,355 |
| 5 | 5 | `results/p5_enumerate/p5_g5_s_1.sqlite` | 196,044,530 | 491,152 | 11,470 |
| 7 | 2 | `results/p7_enumerate/p7_g2_s_1.sqlite` | 93,282 | 89,418 | 749 |
| 7 | 3 | `results/p7_enumerate/p7_g3_s_1.sqlite` | 4,585,140 | 1,165,453 | 6,208 |
| 7 | 4 | `results/p7_enumerate/p7_g4_s_1.sqlite` | 224,743,932 | 6,482,570 | 39,431 |
| 11 | 2 | `results/p11_enumerate/p11_g2_s_1.sqlite` | 1,381,380 | 1,337,216 | 2,813 |
| 11 | 3 | `results/p11_enumerate/p11_g3_s_1.sqlite` | 167,357,190 | 36,163,204 | 47,078 |
| 13 | 2 | `results/p13_enumerate/p13_g2_s_1.sqlite` | 3,739,580 | 3,631,447 | 4,589 |

Interrupted or legacy partial outputs:

| prime | genus | mode | file | processed | sparse presentations | sparse classes | note |
| --- | ---: | --- | --- | ---: | ---: | ---: | --- |
| 5 | 6 | enumerate | `results/p5_enumerate/p5_g6_s_1.sqlite` | 164,519,787 / 4,901,145,620 | 9,549 | 635 | interrupted |
| 5 | 7 | random | `results/p5_random/p5_g7_s_1_random.sqlite` | 144,698,904 | 121,496 | 2,764 | interrupted |
| 3 | 9 | random | `results/p3_random/p3_g9_s_1_random.sqlite` | 174,220,000 | 2,742 | 478 | old file with progress rows but no summary row |
| 3 | 11 | random | `results/p3_random/p3_g11_s_1_random.sqlite` | 141,934,429 / 170,000,000 | 0 | 0 | interrupted |
| 17 | 2 | enumerate | `results/p17_enumerate/p17_g2_s_1.sqlite` | 100,000 | 102,815 | 690 | old file with progress rows but no summary row |

Files that should be regenerated before being treated as clean:

```text
results/p3_enumerate/p3_g6_s_1.sqlite
results/p3_random/p3_g10_s_1_random.sqlite
```

Those two files contain sparse rows, but their `enumeration_summary` metadata
does not match the filename/content.

## Useful Inspection Commands

List sparse rows:

```bash
sqlite3 results/p3_enumerate/p3_g7_s_1.sqlite \
  "SELECT count(*) FROM sparse_curves;"
```

View the run summary:

```bash
sqlite3 results/p3_enumerate/p3_g7_s_1.sqlite \
  "SELECT * FROM enumeration_summary;"
```
