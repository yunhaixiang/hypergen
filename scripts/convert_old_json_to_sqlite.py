#!/usr/bin/env python3
"""Convert legacy old JSON curve data to the current SQLite layout.

Run with Sage's Python:

    sage -python scripts/convert_old_json_to_sqlite.py
"""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
from pathlib import Path

from sage.all import GF, PolynomialRing


SCHEMA = """
CREATE TABLE sparse_curves (
    canonical_key TEXT PRIMARY KEY,
    coefficients TEXT NOT NULL,
    branch_factors TEXT,
    branch_infinity_branch INTEGER,
    branch_leading_coefficient INTEGER,
    branch_factorization_pattern TEXT,
    branch_divisor_type TEXT,
    lpoly TEXT NOT NULL,
    sparsity INTEGER NOT NULL,
    rational_branch_count INTEGER NOT NULL
);

CREATE TABLE enumeration_summary (
    id INTEGER PRIMARY KEY CHECK (id = 1),
    prime INTEGER NOT NULL,
    genus INTEGER NOT NULL,
    max_sparsity INTEGER,
    enumeration_mode TEXT NOT NULL,
    irreducible_memory_budget_mb INTEGER,
    limit_count INTEGER,
    total_coefficient_vectors TEXT,
    processed INTEGER NOT NULL,
    sparse_presentations TEXT,
    sparse_isomorphism_classes INTEGER NOT NULL,
    canonicalized_isomorphism_classes INTEGER NOT NULL,
    elapsed_seconds REAL NOT NULL,
    status_counts TEXT NOT NULL
);

CREATE TABLE enumeration_progress (
    position INTEGER PRIMARY KEY,
    processed INTEGER NOT NULL,
    sparse_presentations TEXT,
    sparse_isomorphism_classes INTEGER NOT NULL,
    canonicalized_isomorphism_classes INTEGER NOT NULL,
    elapsed_seconds REAL NOT NULL
);
"""


def json_compact(value) -> str:
    return json.dumps(value, separators=(",", ":"))


def trim_coefficients(coeffs: list[int]) -> list[int]:
    out = list(coeffs)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def coefficient_support_degree(coeffs: list[int], p: int) -> int:
    for index in range(len(coeffs) - 1, -1, -1):
        if coeffs[index] % p:
            return index
    return 0


def factor_data(coeffs: list[int], p: int, genus: int):
    field = GF(p)
    ring = PolynomialRing(field, "x")
    x = ring.gen()
    poly = sum(field(value) * (x**degree) for degree, value in enumerate(coeffs))
    factors = []
    counts = [0] * (2 * genus + 3)
    partition = []
    rational_finite = 0
    for factor, exponent in poly.factor():
        factor_coeffs = [int(c) for c in factor.list()]
        factor_degree = factor.degree()
        for _ in range(int(exponent)):
            factors.append(factor_coeffs)
            if factor_degree >= len(counts):
                counts.extend([0] * (factor_degree - len(counts) + 1))
            counts[factor_degree] += 1
            partition.append(factor_degree)
            if factor_degree == 1:
                rational_finite += 1
    factors.sort(key=lambda item: (len(item), item))
    partition.sort()
    return factors, partition, rational_finite


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)


def convert_file(path: Path, out_root: Path) -> Path:
    match = re.match(r"p(\d+)_g(\d+)_pgl2save_c\.json$", path.name)
    if not match:
        raise ValueError(f"unexpected old-data filename: {path.name}")
    p = int(match.group(1))
    genus = int(match.group(2))
    data = json.loads(path.read_text())
    curves = data.get("curves") or []

    out_dir = out_root / f"p{p}_old"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"p{p}_g{genus}_s_0_old.sqlite"
    if out_path.exists():
        out_path.unlink()

    conn = sqlite3.connect(out_path)
    try:
        create_schema(conn)
        seen_keys: set[str] = set()
        for curve in curves:
            coeffs = trim_coefficients([int(c) % p for c in curve.get("f_coeffs", [])])
            canonical_coeffs = trim_coefficients([int(c) % p for c in curve.get("canonical_f_coeffs", coeffs)])
            degree = coefficient_support_degree(coeffs, p)
            leading = coeffs[degree] % p if coeffs else 0
            infinity_branch = 1 if degree <= 2 * genus + 1 else 0
            factors, pattern, rational_finite = factor_data(coeffs, p, genus)
            rational_branch_count = rational_finite + infinity_branch
            branch_divisor_type = ([0] if infinity_branch else []) + pattern

            reduction_key = curve.get("reduction_key")
            canonical_key = ",".join(str(int(v)) for v in reduction_key) if reduction_key else ",".join(str(v) for v in canonical_coeffs)
            if canonical_key in seen_keys:
                canonical_key = f"{canonical_key}#old_index={curve.get('index')}"
            seen_keys.add(canonical_key)

            middle = int(curve.get("middle_coefficient", 0))
            lpoly = [0] * max(genus - 1, 0) + [middle]

            conn.execute(
                """
                INSERT INTO sparse_curves (
                    canonical_key,
                    coefficients,
                    branch_factors,
                    branch_infinity_branch,
                    branch_leading_coefficient,
                    branch_factorization_pattern,
                    branch_divisor_type,
                    lpoly,
                    sparsity,
                    rational_branch_count
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    canonical_key,
                    json_compact(coeffs),
                    json_compact(factors),
                    infinity_branch,
                    leading,
                    json_compact(pattern),
                    json_compact(branch_divisor_type),
                    json_compact(lpoly),
                    0,
                    rational_branch_count,
                ),
            )

        sparse_classes = int(data.get("isomorphism_class_count") or len(curves))
        canonical_classes = int(data.get("reduction_class_count") or sparse_classes)
        presentations = int(data.get("total_presentations_found") or len(curves))
        status = {
            "source": str(path),
            "source_format": "old_pgl2save_c_json",
            "search_status": data.get("search_status"),
            "complete_list": data.get("complete_list"),
            "presentation_reduction": data.get("presentation_reduction"),
            "monic_enforced": data.get("monic_enforced"),
            "allow_nonmonic": data.get("allow_nonmonic"),
            "hasse_witt_prefilter": data.get("hasse_witt_prefilter"),
            "workers": data.get("workers"),
            "worker_index": data.get("worker_index"),
            "sparse": sparse_classes,
        }
        conn.execute(
            """
            INSERT INTO enumeration_summary (
                id,
                prime,
                genus,
                max_sparsity,
                enumeration_mode,
                irreducible_memory_budget_mb,
                limit_count,
                total_coefficient_vectors,
                processed,
                sparse_presentations,
                sparse_isomorphism_classes,
                canonicalized_isomorphism_classes,
                elapsed_seconds,
                status_counts
            ) VALUES (1, ?, ?, 0, 'old-support', NULL, NULL, ?, ?, ?, ?, ?, 0.0, ?)
            """,
            (
                p,
                genus,
                str(presentations),
                presentations,
                str(presentations),
                sparse_classes,
                canonical_classes,
                json_compact(status),
            ),
        )
        conn.execute(
            """
            INSERT INTO enumeration_progress (
                position,
                processed,
                sparse_presentations,
                sparse_isomorphism_classes,
                canonicalized_isomorphism_classes,
                elapsed_seconds
            ) VALUES (?, ?, ?, ?, ?, 0.0)
            """,
            (presentations, presentations, str(presentations), sparse_classes, canonical_classes),
        )
        conn.commit()
        result = conn.execute("PRAGMA integrity_check").fetchone()[0]
        if result != "ok":
            raise RuntimeError(f"{out_path}: integrity_check failed: {result}")
    finally:
        conn.close()
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-dir", type=Path, default=Path("results/old/batch_results_c"))
    parser.add_argument("--out-root", type=Path, default=Path("results"))
    args = parser.parse_args()

    files = sorted(args.old_dir.glob("p*_g*_pgl2save_c.json"))
    if not files:
        raise SystemExit(f"no old JSON files found in {args.old_dir}")

    for index, path in enumerate(files, start=1):
        out_path = convert_file(path, args.out_root)
        print(f"converted {index}/{len(files)} {path.name} -> {out_path}")


if __name__ == "__main__":
    main()
