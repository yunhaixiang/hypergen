#!/usr/bin/env python3
"""Rewrite branch_factorization_pattern from degree-count arrays to partitions."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path


def counts_to_partition(values: list[int]) -> list[int]:
    partition: list[int] = []
    for degree, multiplicity in enumerate(values):
        if degree == 0:
            continue
        for _ in range(int(multiplicity)):
            partition.append(degree)
    return partition


def convert_value(raw: str | None) -> tuple[str | None, bool]:
    if raw is None:
        return None, False
    try:
        values = json.loads(raw)
    except json.JSONDecodeError:
        return raw, False
    if not isinstance(values, list) or not values:
        return raw, False
    if not all(isinstance(value, int) for value in values):
        return raw, False
    if values[0] != 0:
        return raw, False
    converted = counts_to_partition(values)
    return json.dumps(converted, separators=(",", ":")), converted != values


def migrate_database(path: Path) -> int:
    conn = sqlite3.connect(path)
    try:
        has_table = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sparse_curves'"
        ).fetchone()
        if not has_table:
            return 0
        columns = [row[1] for row in conn.execute("PRAGMA table_info(sparse_curves)")]
        if "branch_factorization_pattern" not in columns:
            return 0

        updates: list[tuple[str | None, int]] = []
        for rowid, raw in conn.execute(
            "SELECT rowid, branch_factorization_pattern FROM sparse_curves"
        ):
            converted, changed = convert_value(raw)
            if changed:
                updates.append((converted, rowid))
        if not updates:
            return 0
        conn.executemany(
            "UPDATE sparse_curves SET branch_factorization_pattern = ? WHERE rowid = ?",
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
        normalized = str(path)
        if normalized in skip:
            print(f"skipped {path}")
            continue
        changed = migrate_database(path)
        total += changed
        if changed:
            print(f"updated {changed:6d} rows in {path}")
    print(f"updated {total} rows total")


if __name__ == "__main__":
    main()
