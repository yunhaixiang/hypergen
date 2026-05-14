#!/usr/bin/env python3
"""Add and populate sparse_curves.branch_divisor_type.

The value is a JSON partition-like list. A leading 0 denotes the branch point at
infinity. For example, [0,1,11,11] means infinity plus one linear factor and two
degree-11 irreducible factors.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path


def finite_pattern_to_partition(raw: str | None) -> list[int]:
    if raw is None:
        return []
    values = json.loads(raw)
    if not isinstance(values, list):
        return []
    if values and values[0] == 0:
        partition: list[int] = []
        for degree, multiplicity in enumerate(values):
            if degree == 0:
                continue
            for _ in range(int(multiplicity)):
                partition.append(degree)
        return partition
    return [int(value) for value in values]


def branch_divisor_type(infinity_branch: int | None, finite_raw: str | None) -> str:
    values: list[int] = []
    if infinity_branch:
        values.append(0)
    values.extend(finite_pattern_to_partition(finite_raw))
    return json.dumps(values, separators=(",", ":"))


def column_exists(conn: sqlite3.Connection, table: str, column: str) -> bool:
    return any(row[1] == column for row in conn.execute(f"PRAGMA table_info({table})"))


def migrate_database(path: Path) -> int:
    conn = sqlite3.connect(path)
    try:
        has_table = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sparse_curves'"
        ).fetchone()
        if not has_table:
            return 0
        if not column_exists(conn, "sparse_curves", "branch_divisor_type"):
            conn.execute("ALTER TABLE sparse_curves ADD COLUMN branch_divisor_type TEXT")
        updates: list[tuple[str, int]] = []
        for rowid, infinity_branch, finite_raw, current in conn.execute(
            """
            SELECT rowid,
                   branch_infinity_branch,
                   branch_factorization_pattern,
                   branch_divisor_type
            FROM sparse_curves
            """
        ):
            desired = branch_divisor_type(infinity_branch, finite_raw)
            if current != desired:
                updates.append((desired, rowid))
        if updates:
            conn.executemany(
                "UPDATE sparse_curves SET branch_divisor_type = ? WHERE rowid = ?",
                updates,
            )
        conn.commit()
        result = conn.execute("PRAGMA integrity_check").fetchone()[0]
        if result != "ok":
            raise RuntimeError(f"{path}: integrity_check failed: {result}")
        return len(updates)
    finally:
        conn.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--root", type=Path, default=Path("results"))
    parser.add_argument("--skip", action="append", default=[])
    args = parser.parse_args()

    paths = args.paths or sorted(args.root.glob("*/*.sqlite"))
    skip = {str(Path(item)) for item in args.skip}

    total = 0
    for path in paths:
        if str(path) in skip:
            print(f"skipped {path}")
            continue
        changed = migrate_database(path)
        total += changed
        if changed:
            print(f"updated {changed:6d} rows in {path}")
    print(f"updated {total} rows total")


if __name__ == "__main__":
    main()
