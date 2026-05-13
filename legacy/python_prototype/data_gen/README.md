# Hyperelliptic Data Generation

Legacy note: this is a frozen snapshot of the old Python prototype. It is kept
for reference only and is not the active implementation. New development should
use the C++ runner in `hypergen/cpp/`.

This directory has a Python data-generation implementation.

Current Python status: basic prime fields, finite-field polynomials, hyperelliptic model validation, point counting over extensions, Hasse-Witt sparsity filtering, L-polynomial coefficient computation, and sparsity-limited early stopping are implemented.

The Python implementation is in `hyperelliptic.py`.

There is also an experimental C++ implementation in `data_gen/cpp/`. It targets the current branch-divisor pipeline and links against FLINT and SQLite:

```bash
cmake -S data_gen/cpp -B data_gen/cpp/build
cmake --build data_gen/cpp/build -j
```

Example C++ runs:

```bash
data_gen/cpp/build/hyperelliptic_cpp --p 3 --genus 2 --enumeration-mode enumerate --max-sparsity 1 --limit 1000
data_gen/cpp/build/hyperelliptic_cpp --p 5 --genus 10 --enumeration-mode random --random-steps 10000 --max-sparsity 1
```

The C++ runner currently implements the two active modes, `enumerate` and `random`, branch-divisor generation, Hasse-Witt filtering, factorized `PGL_2(F_p)` branch canonicalization with expanded binary-form fallback/storage, exact L-polynomial point counting with integer-encoded extension fields, and SQLite sparse-output tables. It uses FLINT for modular polynomial irreducibility, polynomial products/powers, and Hasse-Witt characteristic polynomials. `--irreducible-memory-budget-mb` controls how many irreducible-factor tables are materialized in memory; the default is `1024`, and degrees beyond the estimated budget are streamed or randomly sampled instead of loaded. The Python implementation remains the reference for the fuller cache/resume schema.

Implemented basic Python structures:

- `PrimeField` stores an odd prime characteristic and normalizes integer representatives.
- `Polynomial` stores coefficients over a `PrimeField` in low-to-high degree order and provides derivative, remainder, gcd, monic normalization, and squarefreeness checks.
- `FiniteExtension` builds a simple polynomial-basis model of `F_{p^r}` using a monic irreducible modulus, with tuple-based compatibility methods and integer-encoded arithmetic for hot point-counting paths.
- `PointCountingContext` caches finite extensions, field elements, quadratic residues, and powers of `x` for reuse across many curves.
- `EnumerationContext` owns enumeration-level caches, uses rational-branch-count plus SQLite orbit lookup for isomorphism matching, canonicalizes binary forms under `PGL_2(F_p)` up to square scalar on orbit-cache misses, tracks seen isomorphism classes, and caches mod-`p` and exact L-polynomial results by canonical key.
- `IrreduciblePolynomialCache` stores complete tables of monic irreducible polynomials by `(prime, degree)` in SQLite, so small-prime factor/enumeration helpers can reuse them across runs.
- `HyperellipticCurve` stores a model `y^2 = f(x)`, validates squarefreeness, computes Hasse-Witt data, counts points over extensions, and computes `a_1, ..., a_g`.

The sparsity-limited method returns `None` as soon as the sparsity among `a_1, ..., a_{g-1}` exceeds the requested limit:

```python
curve.l_polynomial_coefficients_with_sparsity_limit(max_sparsity=1)
```

The Hasse-Witt filter gives a fast safe rejection test modulo `p`:

```python
curve.hasse_witt_matrix()
curve.l_polynomial_coefficients_mod_p()
curve.passes_hasse_witt_sparsity_filter(max_sparsity=1)
```

When a sparsity bound is supplied, enumeration always runs this Hasse-Witt mod-`p` test before canonicalization. Mod-`p` sparsity failures are counted as `rejected_hasse_witt_uncanonicalized`, are not inserted into the canonical-class tables, and never reach exact extension point counts.

For enumeration over a fixed field and degree, share one point-counting context:

```python
field = PrimeField(5)
context = PointCountingContext(field, polynomial_degree=5)
curve = HyperellipticCurve(Polynomial(field, [1, 1, 0, 0, 0, 1]), point_counting_context=context)
```

For full enumeration, use `EnumerationContext` instead:

```python
context = EnumerationContext(prime=5, genus=2)
polynomial = context.polynomial([1, 1, 0, 0, 0, 1])

if context.is_new_isomorphism_class(polynomial):
    mod_p = context.l_polynomial_coefficients_mod_p(polynomial)
    exact = context.l_polynomial_coefficients_with_sparsity_limit(polynomial, max_sparsity=1)
```

To write enumeration output to SQLite while running, pass `sqlite_path` and stream coefficient vectors through the output helper:

```python
context = EnumerationContext(prime=5, genus=2, sqlite_path="curves.sqlite")
stats = context.process_polynomials_for_output(coefficient_vectors, max_sparsity=1)
timing = context.timing_summary()
context.close_sqlite()
```

Reusing the same `sqlite_path` resumes from previous canonical-class rows by loading `curve_cache` into the in-memory indexes at startup:

```python
context = EnumerationContext(prime=5, genus=2, sqlite_path="curves.sqlite")
```

There is also a command-line runner:

```bash
python3 -m data_gen.hyperelliptic --p 5 --genus 2 --max-sparsity 2
```

For a batch run over increasing genus at a fixed prime, use `--genus-start` and `--genus-end`. This writes one SQLite database per genus under `--out-dir`:

```bash
python3 -m data_gen.hyperelliptic --p 5 --genus-start 2 --genus-end 8 --max-sparsity 1 --enumeration-mode random --random-steps 10000 --out-dir data_gen/results/p5_batch
```

Omit `--max-sparsity` to compute without a sparsity restriction:

```bash
python3 -m data_gen.hyperelliptic --p 5 --genus 2
```

For high-genus sparse search, use a sparsity bound; Hasse-Witt filtering is automatic:

```bash
python3 -m data_gen.hyperelliptic --p 5 --genus 20 --max-sparsity 1
```

First-time irreducible generation uses aperiodic necklace representatives: each Lyndon representative is interpreted in a normal basis of `F_{p^d}`, then converted to its Frobenius-orbit minimal polynomial. To warm the persistent irreducible-polynomial cache for a small prime and degrees `1..N`, run:

```bash
python3 -m data_gen.hyperelliptic --p 5 --precompute-irreducibles-up-to-degree 8
```

The default cache is prime-specific, `data_gen/irreducibles/irreducibles_p{p}.sqlite`; override it with `--irreducible-cache PATH`. Enumeration contexts also use this cache lazily when they need monic irreducible polynomials.

Irreducible table materialization is budgeted. The default memory budget is `1024` MB and can be changed with `--irreducible-memory-budget-mb`. During branch enumeration, degrees already present in SQLite are loaded into memory when the budget permits; missing degrees are streamed or sampled instead of being generated into the cache.

Precomputed small-prime cache files are kept under `data_gen/irreducibles/`, for example `data_gen/irreducibles/irreducibles_p5.sqlite`. This folder is separate from the runtime `data_gen/cache/` directory and is ignored as generated data.

By default, deterministic `enumerate` mode uses SQLite BLOB orbit lookup and stores orbit keys for complete canonical-class enumeration. If `--max-sparsity` is set, Hasse-Witt runs before canonicalization; failures are counted as `rejected_hasse_witt_uncanonicalized` and are not inserted into the canonical-class tables. Hasse-Witt survivors in `enumerate` mode use SQLite orbit lookup before full canonicalization. In `random` mode, Hasse-Witt survivors compute the exact L-polynomial first and canonicalize only if the curve passes the sparsity limit; exact sparse failures are counted as `rejected_exact_uncanonicalized`.

The current enumeration modes always enumerate both degree families: degree `2g+1` models first, then degree `2g+2` models. The runner uses a fixed leading-coefficient normalization: degree `2g+1` models are enumerated monic, and degree `2g+2` models are enumerated with leading coefficient `1` and the smallest nonsquare in `F_p`. Since odd models are included, squarefree degree `2g+2` models with an `F_p`-rational branch point are skipped as `covered_by_odd_model`; a `PGL_2(F_p)` transform can move that branch point to infinity, giving a degree `2g+1` representative.

The default enumeration mode is `--enumeration-mode enumerate`, which enumerates squarefree branch divisors by factorization pattern and irreducible-factor choices. The other active mode is `--enumeration-mode random`, which samples random factorized branch divisors for high-genus sparse search. Before branch enumeration starts, it loads irreducibles already present in the SQLite cache; missing degrees are streamed by the necklace path for deterministic enumeration and sampled with Sage for random search. Even degree `2g+2` branch patterns with linear factors are skipped because those curves are represented by odd degree models. In bounded sparse branch runs, the Hasse-Witt prefilter is computed directly from the irreducible factors, so mod-`p` sparsity failures are rejected without first expanding `f(x)`. In `random` mode, exact point counts also run before canonicalization, so nonsparse exact failures are skipped without expanding the polynomial or running PGL2 canonicalization. Sparse branch hits store readable factorized branch data in `sparse_curves.branch_factors` plus the readable expanded coefficients and L-polynomial. Branch candidates now use the same factorized PGL2 canonicalization path in all degrees; if a factorized transform cannot represent a branch divisor in the normalized odd/even model family, the code falls back to expanded binary-form canonicalization. Exact L-polynomial computation for branch candidates uses integer-encoded extension-field elements and cached quadratic-character arrays `(k, q) -> [chi(q(a))]`, then multiplies those arrays to compute point counts. Exact repeated branch keys are memoized before Hasse-Witt, and branch orbit lookups also keep a factorization-pattern-indexed in-memory cache. The factorized path caches small-degree PGL2 action matrices and transformed small factors in memory, and also persists small-factor transforms in `branch_factor_transform_cache` when SQLite output is enabled. Random mode uses `--random-seed` for reproducibility and `--random-max-factors` to bias toward a small number of Frobenius orbits. In random mode, SQLite stores only reduced orbit-cache rows for exact sparse hits instead of every transformed orbit representative; deterministic `enumerate` mode still stores full orbit rows for exhaustive duplicate lookup. If `--random-steps` or `--limit` is provided, random mode stops after that many samples; if both are omitted, it runs until interrupted and records the unknown total as `-1` in SQLite. Random branch sampling uses SQLite random access for cached irreducible degrees; for missing degrees it uses Sage's random irreducible generator through an in-process Sage import when available, otherwise through a persistent `sage -python` helper process.

The runner prints progress as:

```text
prime: p
genus: g
progress: processed/total
sparse_presentations: N
sparse_isomorphism_classes: M
canonicalized_isomorphism_classes: K
-
```

The progress line always reports `canonicalized_isomorphism_classes`, since Hasse-Witt failures are not canonicalized in bounded sparse runs.

The printed timing summary includes separate counters for Hasse-Witt filtering, factorized PGL2 work, polynomial expansion, exact L-polynomial computation, ground invariants, SQLite loading, and SQLite writing. CLI runs batch internal SQLite writes, while progress and final summary rows are still flushed explicitly.

The SQLite output uses:

- `orbit_cache`: BLOB lookup table for `(rational_branch_count, ground_point_count, hasse_witt_lpoly_mod_p, orbit_key) -> canonical_key` in bounded sparse runs. Older databases without the Hasse-Witt column are still readable.
- `curve_cache`: per-canonical-key computation results for the output file's sparsity bound, including rational branch count, L-polynomial data, and rejection status.
- `sparse_curves`: sparse survivors with their expanded coefficient form, readable factorized branch form, and exact `a_1, ..., a_g` coefficients.
- `enumeration_summary`: one-row run summary with `prime`, `genus`, `max_sparsity`, enumeration settings, `total_coefficient_vectors`, progress counts, and timing fields. In deterministic `enumerate` mode, `total_coefficient_vectors` is the full branch-divisor presentation count even if `--limit` stops the run early; bounded random runs store the requested sample count, and unbounded random runs store `-1`.
- `enumeration_progress`: progress snapshots at the same cadence as printed progress, including cumulative counts and delta counts since the previous snapshot. Use `delta_canonicalized_isomorphism_classes` and `delta_sparse_isomorphism_classes` to locate spike and flat regions.
- `branch_orbit_cache` and `branch_curve_cache`: branch mode caches for factorized PGL2 orbit matching before polynomial expansion.
- `branch_factor_transform_cache`: persistent cache for small irreducible-factor PGL2 transforms.

In `curve_cache` and `sparse_curves`, `coefficients` are readable JSON integer lists. In `sparse_curves`, `branch_factors`, `branch_factorization_pattern`, and `lpoly` are also readable JSON. Internal cache keys and intermediate L-polynomial fields are stored as compact BLOBs.

With `--max-sparsity`, Hasse-Witt survivors use SQLite orbit lookup to skip repeated canonicalization when a presentation is already in a stored `PGL_2` orbit. The difference is that Hasse-Witt failures are not canonicalized or stored.

## Python Tests

From the repository root:

```bash
python3 -m unittest data_gen.tests.test_hyperelliptic
```
