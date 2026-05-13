from __future__ import annotations

import argparse
from collections import OrderedDict
from dataclasses import dataclass, field
from itertools import combinations, product
import json
from math import comb
import os
from pathlib import Path
import random
import shutil
import sqlite3
import subprocess
from time import perf_counter
from typing import Iterable, Optional
import zlib


DEFAULT_IRREDUCIBLE_MEMORY_BUDGET_MB = 1024
BYTES_PER_MEGABYTE = 1024 * 1024
BRANCH_FACTOR_ACTION_MATRIX_CACHE_MAX_DEGREE = 8
BRANCH_FACTOR_TRANSFORM_CACHE_MAX_ENTRIES = 200_000
BRANCH_PRODUCT_CHARACTER_CACHE_MAX_ENTRIES = 100_000
BRANCH_HASSE_WITT_PRODUCT_CACHE_MAX_ENTRIES = 100_000
SQLITE_WRITE_BATCH_SIZE = 100
SAGE_RANDOM_IRREDUCIBLE_SERVER_CODE = r"""
import json
import sys
from sage.all import GF, PolynomialRing, set_random_seed

rings = {}
while True:
    line = sys.stdin.readline()
    if not line:
        break
    try:
        request = json.loads(line)
        prime = int(request["prime"])
        degree = int(request["degree"])
        count = int(request["count"])
        seed = int(request["seed"])
        set_random_seed(seed)
        ring = rings.get(prime)
        if ring is None:
            ring = PolynomialRing(GF(prime), "x")
            rings[prime] = ring
        polynomials = []
        seen = set()
        while len(polynomials) < count:
            polynomial = ring.irreducible_element(degree, algorithm="random")
            leading = int(polynomial[degree]) % prime
            leading_inverse = pow(leading, -1, prime)
            coefficients = tuple((int(polynomial[index]) * leading_inverse) % prime for index in range(degree + 1))
            if coefficients not in seen:
                seen.add(coefficients)
                polynomials.append(list(coefficients))
        print(json.dumps({"ok": True, "polynomials": polynomials}, separators=(",", ":")), flush=True)
    except Exception as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, separators=(",", ":")), flush=True)
"""


def _is_odd_prime(value: int) -> bool:
    if value < 3 or value % 2 == 0:
        return False
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


@dataclass
class CanonicalRecord:
    canonical_key: tuple[int, ...]
    coefficients: tuple[int, ...]
    rational_branch_count: int
    ground_point_count: Optional[int] = None
    hasse_witt_lpoly_mod_p: Optional[tuple[int, ...]] = None
    orbit_size: int = 1
    status_by_max_sparsity: dict[Optional[int], str] = field(default_factory=dict)
    sparsity_by_max_sparsity: dict[Optional[int], Optional[int]] = field(default_factory=dict)
    exact_lpoly_by_max_sparsity: dict[Optional[int], Optional[tuple[int, ...]]] = field(default_factory=dict)


@dataclass(frozen=True)
class BranchDivisorCandidate:
    leading_coefficient: int
    factors: tuple[tuple[int, ...], ...]
    infinity_branch: bool
    factorization_pattern: tuple[int, ...]


@dataclass(frozen=True)
class PrimeField:
    prime: int

    def __post_init__(self) -> None:
        if not _is_odd_prime(self.prime):
            raise ValueError("prime field characteristic must be an odd prime")

    def normalize(self, value: int) -> int:
        return value % self.prime

    def inverse(self, value: int) -> int:
        value %= self.prime
        if value == 0:
            raise ValueError("zero has no inverse")
        return pow(value, -1, self.prime)


def _trim(coefficients: list[int]) -> list[int]:
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def _poly_trim_mod(coefficients: list[int], p: int) -> tuple[int, ...]:
    reduced = [coefficient % p for coefficient in coefficients]
    return tuple(_trim(reduced))


def _poly_add_mod(lhs: tuple[int, ...], rhs: tuple[int, ...], p: int) -> tuple[int, ...]:
    length = max(len(lhs), len(rhs))
    return _poly_trim_mod(
        [(lhs[i] if i < len(lhs) else 0) + (rhs[i] if i < len(rhs) else 0) for i in range(length)],
        p,
    )


def _poly_sub_mod(lhs: tuple[int, ...], rhs: tuple[int, ...], p: int) -> tuple[int, ...]:
    length = max(len(lhs), len(rhs))
    return _poly_trim_mod(
        [(lhs[i] if i < len(lhs) else 0) - (rhs[i] if i < len(rhs) else 0) for i in range(length)],
        p,
    )


def _poly_mul_mod(lhs: tuple[int, ...], rhs: tuple[int, ...], p: int) -> tuple[int, ...]:
    if lhs == (0,) or rhs == (0,):
        return (0,)
    result = [0] * (len(lhs) + len(rhs) - 1)
    for i, lhs_coefficient in enumerate(lhs):
        if lhs_coefficient == 0:
            continue
        for j, rhs_coefficient in enumerate(rhs):
            if rhs_coefficient == 0:
                continue
            result[i + j] += lhs_coefficient * rhs_coefficient
    return _poly_trim_mod(result, p)


def _poly_pow_mod(base: tuple[int, ...], exponent: int, p: int) -> tuple[int, ...]:
    result = (1,)
    while exponent > 0:
        if exponent % 2 == 1:
            result = _poly_mul_mod(result, base, p)
        base = _poly_mul_mod(base, base, p)
        exponent //= 2
    return result


def _poly_product_tree_mod(polynomials: Iterable[tuple[int, ...]], p: int) -> tuple[int, ...]:
    level = list(polynomials)
    if not level:
        return (1,)
    while len(level) > 1:
        next_level = []
        for index in range(0, len(level), 2):
            if index + 1 == len(level):
                next_level.append(level[index])
            else:
                next_level.append(_poly_mul_mod(level[index], level[index + 1], p))
        level = next_level
    return level[0]


def _poly_exact_div_mod(numerator: tuple[int, ...], denominator: tuple[int, ...], p: int) -> tuple[int, ...]:
    if denominator == (0,):
        raise ZeroDivisionError("division by zero polynomial")
    if numerator == (0,):
        return (0,)

    rem = list(numerator)
    quotient = [0] * max(1, len(numerator) - len(denominator) + 1)
    denominator_degree = len(denominator) - 1
    denominator_leading_inverse = pow(denominator[-1], -1, p)

    while len(rem) >= len(denominator) and rem != [0]:
        shift = len(rem) - len(denominator)
        scale = rem[-1] * denominator_leading_inverse % p
        quotient[shift] = scale
        for i in range(denominator_degree + 1):
            rem[shift + i] = (rem[shift + i] - scale * denominator[i]) % p
        _trim(rem)

    if rem != [0]:
        raise ValueError("polynomial division was not exact")
    return _poly_trim_mod(quotient, p)


def _poly_scalar_mul_mod(polynomial: tuple[int, ...], scalar: int, p: int) -> tuple[int, ...]:
    return _poly_trim_mod([scalar * coefficient for coefficient in polynomial], p)


def _poly_exact_quotient_mod(numerator: tuple[int, ...], denominator: tuple[int, ...], p: int) -> tuple[int, ...]:
    return _poly_exact_div_mod(_poly_trim_mod(list(numerator), p), _poly_trim_mod(list(denominator), p), p)


def _poly_remainder_mod(numerator: tuple[int, ...], denominator: tuple[int, ...], p: int) -> tuple[int, ...]:
    if denominator == (0,):
        raise ZeroDivisionError("division by zero polynomial")
    rem = list(_poly_trim_mod(list(numerator), p))
    denominator = _poly_trim_mod(list(denominator), p)
    denominator_degree = len(denominator) - 1
    denominator_leading_inverse = pow(denominator[-1], -1, p)

    while len(rem) >= len(denominator) and rem != [0]:
        shift = len(rem) - len(denominator)
        scale = rem[-1] * denominator_leading_inverse % p
        for i in range(denominator_degree + 1):
            rem[shift + i] = (rem[shift + i] - scale * denominator[i]) % p
        _trim(rem)

    return _poly_trim_mod(rem, p)


def _poly_monic_mod(polynomial: tuple[int, ...], p: int) -> tuple[int, ...]:
    polynomial = _poly_trim_mod(list(polynomial), p)
    if polynomial == (0,):
        return polynomial
    leading_inverse = pow(polynomial[-1], -1, p)
    return _poly_scalar_mul_mod(polynomial, leading_inverse, p)


def _poly_gcd_mod(lhs: tuple[int, ...], rhs: tuple[int, ...], p: int) -> tuple[int, ...]:
    lhs = _poly_trim_mod(list(lhs), p)
    rhs = _poly_trim_mod(list(rhs), p)
    while rhs != (0,):
        lhs, rhs = rhs, _poly_remainder_mod(lhs, rhs, p)
    return _poly_monic_mod(lhs, p)


def _mobius(value: int) -> int:
    remaining = value
    prime_factor_count = 0
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            remaining //= divisor
            prime_factor_count += 1
            if remaining % divisor == 0:
                return 0
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1 if divisor == 2 else 2
    if remaining > 1:
        prime_factor_count += 1
    return -1 if prime_factor_count % 2 else 1


def monic_irreducible_count(prime: int, degree: int) -> int:
    if degree < 1:
        raise ValueError("degree must be positive")
    return sum(
        _mobius(divisor) * prime ** (degree // divisor)
        for divisor in range(1, degree + 1)
        if degree % divisor == 0
    ) // degree


def estimated_irreducible_tuple_memory_bytes(prime: int, degree: int) -> int:
    count = monic_irreducible_count(prime, degree)
    return count * (40 + 36 * (degree + 1))


def _memory_budget_bytes(memory_budget_mb: int) -> int:
    if memory_budget_mb < 1:
        raise ValueError("irreducible memory budget must be positive")
    return memory_budget_mb * BYTES_PER_MEGABYTE


def _poly_mul_remainder_mod(lhs: tuple[int, ...], rhs: tuple[int, ...], modulus: tuple[int, ...], p: int) -> tuple[int, ...]:
    return _poly_remainder_mod(_poly_mul_mod(lhs, rhs, p), modulus, p)


def _poly_pow_remainder_mod(base: tuple[int, ...], exponent: int, modulus: tuple[int, ...], p: int) -> tuple[int, ...]:
    result = (1,)
    base = _poly_remainder_mod(base, modulus, p)
    while exponent > 0:
        if exponent % 2 == 1:
            result = _poly_mul_remainder_mod(result, base, modulus, p)
        base = _poly_mul_remainder_mod(base, base, modulus, p)
        exponent //= 2
    return result


def _determinant_polynomial_mod(matrix: list[list[tuple[int, ...]]], p: int) -> tuple[int, ...]:
    n = len(matrix)
    if n == 0:
        return (1,)
    if n == 1:
        return matrix[0][0]

    working = [[entry for entry in row] for row in matrix]
    previous_pivot = (1,)
    sign = 1

    for k in range(n - 1):
        pivot_row = None
        pivot_column = None
        for i in range(k, n):
            for j in range(k, n):
                if working[i][j] != (0,):
                    pivot_row = i
                    pivot_column = j
                    break
            if pivot_row is not None:
                break

        if pivot_row is None or pivot_column is None:
            return (0,)

        if pivot_row != k:
            working[k], working[pivot_row] = working[pivot_row], working[k]
            sign = -sign
        if pivot_column != k:
            for row in working:
                row[k], row[pivot_column] = row[pivot_column], row[k]
            sign = -sign

        pivot = working[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = _poly_sub_mod(
                    _poly_mul_mod(working[i][j], pivot, p),
                    _poly_mul_mod(working[i][k], working[k][j], p),
                    p,
                )
                working[i][j] = _poly_exact_div_mod(numerator, previous_pivot, p)

        previous_pivot = pivot

    determinant = working[n - 1][n - 1]
    if sign == -1:
        determinant = _poly_scalar_mul_mod(determinant, -1, p)
    return determinant


def _l_polynomial_coefficients_mod_p_from_hasse_witt_matrix(
    matrix: tuple[tuple[int, ...], ...],
    prime: int,
) -> list[int]:
    genus = len(matrix)
    polynomial_matrix = []
    for i, row in enumerate(matrix):
        polynomial_row = []
        for j, entry in enumerate(row):
            if i == j:
                polynomial_row.append(_poly_trim_mod([1, -entry], prime))
            else:
                polynomial_row.append(_poly_trim_mod([0, -entry], prime))
        polynomial_matrix.append(polynomial_row)

    determinant = _determinant_polynomial_mod(polynomial_matrix, prime)
    return [determinant[i] if i < len(determinant) else 0 for i in range(1, genus + 1)]


def hasse_witt_matrix_from_branch_factors(
    factors: Iterable[tuple[int, ...]],
    leading_coefficient: int,
    genus: int,
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    exponent = (prime - 1) // 2
    h = (pow(leading_coefficient, exponent, prime),)
    for factor in factors:
        h = _poly_mul_mod(h, _poly_pow_mod(factor, exponent, prime), prime)

    rows = []
    for i in range(1, genus + 1):
        row = []
        for j in range(1, genus + 1):
            coefficient_index = prime * i - j
            row.append(h[coefficient_index] if coefficient_index < len(h) else 0)
        rows.append(tuple(row))
    return tuple(rows)


def l_polynomial_coefficients_mod_p_from_branch_factors(
    factors: Iterable[tuple[int, ...]],
    leading_coefficient: int,
    genus: int,
    prime: int,
) -> list[int]:
    return _l_polynomial_coefficients_mod_p_from_hasse_witt_matrix(
        hasse_witt_matrix_from_branch_factors(factors, leading_coefficient, genus, prime),
        prime,
    )


def _affine_polynomial_to_binary_form(polynomial: Polynomial, binary_degree: int) -> tuple[int, ...]:
    if polynomial.degree > binary_degree:
        raise ValueError("polynomial degree exceeds binary form degree")
    return tuple(polynomial.coefficient(i) for i in range(binary_degree + 1))


def _precompute_pgl2(prime: int) -> tuple[tuple[int, int, int, int], ...]:
    representatives = {}
    for a, b, c, d in product(range(prime), repeat=4):
        determinant = (a * d - b * c) % prime
        if determinant == 0:
            continue

        entries = (a, b, c, d)
        first_nonzero = next(entry for entry in entries if entry != 0)
        scale = pow(first_nonzero, -1, prime)
        representative = tuple(entry * scale % prime for entry in entries)
        representatives[representative] = representative

    return tuple(sorted(representatives))


def _linear_power_coefficients(alpha: int, beta: int, exponent: int, prime: int) -> tuple[int, ...]:
    return tuple(
        comb(exponent, k) * pow(alpha, k, prime) * pow(beta, exponent - k, prime) % prime
        for k in range(exponent + 1)
    )


def _transform_binary_form(binary_form: tuple[int, ...], matrix: tuple[int, int, int, int], prime: int) -> tuple[int, ...]:
    a, b, c, d = matrix
    degree = len(binary_form) - 1
    result = [0] * (degree + 1)

    powers_ax_bz = [_linear_power_coefficients(a, b, exponent, prime) for exponent in range(degree + 1)]
    powers_cx_dz = [_linear_power_coefficients(c, d, exponent, prime) for exponent in range(degree + 1)]

    for i, coefficient in enumerate(binary_form):
        if coefficient == 0:
            continue
        left = powers_ax_bz[i]
        right = powers_cx_dz[degree - i]
        for left_x_degree, left_coefficient in enumerate(left):
            if left_coefficient == 0:
                continue
            for right_x_degree, right_coefficient in enumerate(right):
                if right_coefficient == 0:
                    continue
                x_degree = left_x_degree + right_x_degree
                result[x_degree] = (result[x_degree] + coefficient * left_coefficient * right_coefficient) % prime

    return tuple(result)


def _pgl2_action_matrix(matrix: tuple[int, int, int, int], degree: int, prime: int) -> tuple[tuple[int, ...], ...]:
    a, b, c, d = matrix
    columns = []

    powers_ax_bz = [_linear_power_coefficients(a, b, exponent, prime) for exponent in range(degree + 1)]
    powers_cx_dz = [_linear_power_coefficients(c, d, exponent, prime) for exponent in range(degree + 1)]

    for i in range(degree + 1):
        column = [0] * (degree + 1)
        left = powers_ax_bz[i]
        right = powers_cx_dz[degree - i]
        for left_x_degree, left_coefficient in enumerate(left):
            if left_coefficient == 0:
                continue
            for right_x_degree, right_coefficient in enumerate(right):
                if right_coefficient == 0:
                    continue
                x_degree = left_x_degree + right_x_degree
                column[x_degree] = (column[x_degree] + left_coefficient * right_coefficient) % prime
        columns.append(tuple(column))

    return tuple(columns)


def _apply_binary_form_action_matrix(binary_form: tuple[int, ...], action_matrix: tuple[tuple[int, ...], ...], prime: int) -> tuple[int, ...]:
    result = [0] * len(binary_form)
    for coefficient, column in zip(binary_form, action_matrix):
        if coefficient == 0:
            continue
        for j, matrix_entry in enumerate(column):
            if matrix_entry != 0:
                result[j] = (result[j] + coefficient * matrix_entry) % prime
    return tuple(result)


def _normalize_binary_form_up_to_square_scalar(binary_form: tuple[int, ...], prime: int) -> tuple[int, ...]:
    if all(coefficient == 0 for coefficient in binary_form):
        raise ValueError("zero binary form cannot be normalized")

    square_scalars = sorted({value * value % prime for value in range(1, prime)})
    return min(tuple(scalar * coefficient % prime for coefficient in binary_form) for scalar in square_scalars)


def _smallest_nonsquare(prime: int) -> int:
    squares = {value * value % prime for value in range(1, prime)}
    for value in range(2, prime):
        if value not in squares:
            return value
    raise ValueError("prime field has no nonsquare representative")


def _pack_int_tuple(values: tuple[int, ...] | list[int]) -> str:
    return json.dumps(list(values), separators=(",", ":"))


def _unpack_int_tuple(text: str | bytes) -> tuple[int, ...]:
    if isinstance(text, bytes):
        text = text.decode("ascii")
    return tuple(json.loads(text))


def _pack_general_int_tuple_blob(values: tuple[int, ...] | list[int]) -> bytes:
    return zlib.compress(json.dumps(list(values), separators=(",", ":")).encode("ascii"))


def _unpack_general_int_tuple_blob(value: str | bytes | None) -> Optional[tuple[int, ...]]:
    if value is None:
        return None
    if isinstance(value, str):
        return _unpack_int_tuple(value)
    try:
        return tuple(json.loads(zlib.decompress(value).decode("ascii")))
    except zlib.error:
        return _unpack_int_tuple(value)


def _pack_branch_factors_text(factors: Iterable[tuple[int, ...]]) -> str:
    return json.dumps([list(factor) for factor in factors], separators=(",", ":"))


def _pack_branch_pattern_text(pattern: tuple[int, ...] | list[int]) -> str:
    return json.dumps(list(pattern), separators=(",", ":"))


def _parse_max_sparsity_from_sqlite_path(path: Optional[Path]) -> Optional[int]:
    if path is None:
        return None
    stem = path.stem
    marker = "_s_"
    if marker not in stem:
        marker = "_sparsity"
        if marker not in stem:
            return None
    suffix = stem.rsplit(marker, 1)[1]
    digits = []
    for character in suffix:
        if not character.isdigit():
            break
        digits.append(character)
    return int("".join(digits)) if digits else None


def default_irreducible_cache_path(prime: int) -> Path:
    return Path(__file__).resolve().parent / "irreducibles" / f"irreducibles_p{prime}.sqlite"


class IrreduciblePolynomialCache:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.path)
        self.connection.execute("PRAGMA journal_mode=DELETE")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self._initialize_schema()

    def _initialize_schema(self) -> None:
        self.connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS irreducible_cache_metadata (
                prime INTEGER NOT NULL,
                degree INTEGER NOT NULL,
                complete INTEGER NOT NULL,
                generated_at REAL NOT NULL,
                count INTEGER NOT NULL,
                PRIMARY KEY (prime, degree)
            );

            CREATE TABLE IF NOT EXISTS irreducible_polynomials (
                prime INTEGER NOT NULL,
                degree INTEGER NOT NULL,
                position INTEGER NOT NULL,
                coefficients BLOB NOT NULL,
                PRIMARY KEY (prime, degree, position),
                UNIQUE (prime, degree, coefficients)
            );

            CREATE INDEX IF NOT EXISTS idx_irreducible_polynomials_prime_degree
                ON irreducible_polynomials (prime, degree);
            """
        )
        self.connection.commit()

    def get_complete(self, prime: int, degree: int) -> Optional[tuple[tuple[int, ...], ...]]:
        row = self.connection.execute(
            """
            SELECT complete, count
            FROM irreducible_cache_metadata
            WHERE prime = ? AND degree = ?
            """,
            (prime, degree),
        ).fetchone()
        if row is None or not row[0]:
            return None

        rows = self.connection.execute(
            """
            SELECT coefficients
            FROM irreducible_polynomials
            WHERE prime = ? AND degree = ?
            ORDER BY position
            """,
            (prime, degree),
        ).fetchall()
        if len(rows) != row[1]:
            return None

        return tuple(_unpack_general_int_tuple_blob(coefficients) or () for (coefficients,) in rows)

    def has_complete(self, prime: int, degree: int) -> bool:
        row = self.connection.execute(
            """
            SELECT complete
            FROM irreducible_cache_metadata
            WHERE prime = ? AND degree = ?
            """,
            (prime, degree),
        ).fetchone()
        return row is not None and bool(row[0])

    def complete_count(self, prime: int, degree: int) -> Optional[int]:
        row = self.connection.execute(
            """
            SELECT complete, count
            FROM irreducible_cache_metadata
            WHERE prime = ? AND degree = ?
            """,
            (prime, degree),
        ).fetchone()
        if row is None or not row[0]:
            return None
        return int(row[1])

    def get_at_position(self, prime: int, degree: int, position: int) -> tuple[int, ...]:
        row = self.connection.execute(
            """
            SELECT coefficients
            FROM irreducible_polynomials
            WHERE prime = ? AND degree = ? AND position = ?
            """,
            (prime, degree, position),
        ).fetchone()
        if row is None:
            raise IndexError("irreducible polynomial position is outside the cached degree block")
        unpacked = _unpack_general_int_tuple_blob(row[0])
        if unpacked is None:
            raise ValueError("cached irreducible polynomial entry is empty")
        return unpacked

    def iter_complete(self, prime: int, degree: int) -> Iterable[tuple[int, ...]]:
        if not self.has_complete(prime, degree):
            return
        rows = self.connection.execute(
            """
            SELECT coefficients
            FROM irreducible_polynomials
            WHERE prime = ? AND degree = ?
            ORDER BY position
            """,
            (prime, degree),
        )
        for (coefficients,) in rows:
            unpacked = _unpack_general_int_tuple_blob(coefficients)
            if unpacked is not None:
                yield unpacked

    def store_complete(self, prime: int, degree: int, polynomials: tuple[tuple[int, ...], ...]) -> None:
        with self.connection:
            self.connection.execute(
                "DELETE FROM irreducible_polynomials WHERE prime = ? AND degree = ?",
                (prime, degree),
            )
            self.connection.executemany(
                """
                INSERT INTO irreducible_polynomials (prime, degree, position, coefficients)
                VALUES (?, ?, ?, ?)
                """,
                (
                    (prime, degree, position, _pack_general_int_tuple_blob(coefficients))
                    for position, coefficients in enumerate(polynomials)
                ),
            )
            self.connection.execute(
                """
                INSERT OR REPLACE INTO irreducible_cache_metadata (
                    prime,
                    degree,
                    complete,
                    generated_at,
                    count
                )
                VALUES (?, ?, 1, ?, ?)
                """,
                (prime, degree, perf_counter(), len(polynomials)),
            )

    def close(self) -> None:
        self.connection.close()


@dataclass(frozen=True)
class Polynomial:
    field: PrimeField
    coefficients: tuple[int, ...]

    def __init__(self, field: PrimeField, coefficients: list[int] | tuple[int, ...]) -> None:
        if not coefficients:
            coefficients = [0]
        if any(coefficient < 0 or coefficient >= field.prime for coefficient in coefficients):
            raise ValueError(f"polynomial coefficients must be in 0..{field.prime - 1}")
        object.__setattr__(self, "field", field)
        object.__setattr__(self, "coefficients", tuple(_trim(list(coefficients))))

    def coefficient(self, index: int) -> int:
        return self.coefficients[index] if index < len(self.coefficients) else 0

    @property
    def degree(self) -> int:
        return len(self.coefficients) - 1

    def is_zero(self) -> bool:
        return self.coefficients == (0,)

    def is_monic(self) -> bool:
        return not self.is_zero() and self.coefficients[-1] == 1

    def derivative(self) -> Polynomial:
        if self.degree == 0:
            return Polynomial(self.field, [0])
        return Polynomial(
            self.field,
            [self.field.normalize(i * self.coefficients[i]) for i in range(1, len(self.coefficients))],
        )

    def remainder(self, divisor: Polynomial) -> Polynomial:
        self._require_same_field(divisor)
        if divisor.is_zero():
            raise ValueError("division by zero polynomial")

        rem = list(self.coefficients)
        divisor_degree = divisor.degree
        divisor_leading_inverse = self.field.inverse(divisor.coefficients[-1])

        while len(rem) - 1 >= divisor_degree and rem != [0]:
            shift = len(rem) - 1 - divisor_degree
            scale = self.field.normalize(rem[-1] * divisor_leading_inverse)
            for i in range(divisor_degree + 1):
                rem[shift + i] = self.field.normalize(rem[shift + i] - scale * divisor.coefficient(i))
            _trim(rem)

        return Polynomial(self.field, rem)

    def monic(self) -> Polynomial:
        if self.is_zero():
            return self
        leading_inverse = self.field.inverse(self.coefficients[-1])
        return Polynomial(self.field, [self.field.normalize(coefficient * leading_inverse) for coefficient in self.coefficients])

    def gcd(self, other: Polynomial) -> Polynomial:
        self._require_same_field(other)
        current = self
        while not other.is_zero():
            current, other = other, current.remainder(other)
        return current.monic()

    def is_squarefree(self) -> bool:
        return not self.is_zero() and self.gcd(self.derivative()).degree == 0

    def _require_same_field(self, other: Polynomial) -> None:
        if self.field.prime != other.field.prime:
            raise ValueError("polynomials are over different fields")


class FiniteExtension:
    def __init__(self, prime: int, degree: int) -> None:
        if degree < 1:
            raise ValueError("extension degree must be positive")
        self.prime = prime
        self.degree = degree
        self.modulus = self._find_irreducible_polynomial(prime, degree)
        self._elements: Optional[tuple[tuple[int, ...], ...]] = None
        self._int_elements: Optional[range] = None
        self._squares: Optional[set[tuple[int, ...]]] = None
        self._int_quadratic_characters: Optional[tuple[int, ...]] = None
        self._point_contributions: Optional[dict[tuple[int, ...], int]] = None
        self._frobenius_matrix: Optional[tuple[tuple[int, ...], ...]] = None

    @property
    def size(self) -> int:
        return self.prime**self.degree

    def zero(self) -> tuple[int, ...]:
        return (0,) * self.degree

    def one(self) -> tuple[int, ...]:
        return (1,) + (0,) * (self.degree - 1)

    def constant(self, value: int) -> tuple[int, ...]:
        return (value % self.prime,) + (0,) * (self.degree - 1)

    def elements(self) -> tuple[tuple[int, ...], ...]:
        if self._elements is None:
            self._elements = tuple(product(range(self.prime), repeat=self.degree))
        return self._elements

    def int_elements(self) -> range:
        if self._int_elements is None:
            self._int_elements = range(self.size)
        return self._int_elements

    def encode(self, element: tuple[int, ...]) -> int:
        value = 0
        place = 1
        for coefficient in element:
            value += (coefficient % self.prime) * place
            place *= self.prime
        return value

    def decode(self, value: int) -> tuple[int, ...]:
        coefficients = []
        for _ in range(self.degree):
            coefficients.append(value % self.prime)
            value //= self.prime
        return tuple(coefficients)

    def int_zero(self) -> int:
        return 0

    def int_one(self) -> int:
        return 1

    def int_constant(self, value: int) -> int:
        return value % self.prime

    def is_zero(self, element: tuple[int, ...]) -> bool:
        return all(coefficient == 0 for coefficient in element)

    def is_one(self, element: tuple[int, ...]) -> bool:
        return element == self.one()

    def is_square(self, element: tuple[int, ...]) -> bool:
        return element in self.squares()

    def int_is_square(self, element: int) -> bool:
        return self.int_quadratic_character(element) >= 0

    def squares(self) -> set[tuple[int, ...]]:
        if self._squares is None:
            self._squares = {self.multiply(element, element) for element in self.elements()}
        return self._squares

    def point_contributions(self) -> dict[tuple[int, ...], int]:
        if self._point_contributions is None:
            zero = self.zero()
            squares = self.squares()
            self._point_contributions = {
                element: 1 if element == zero else 2 if element in squares else 0
                for element in self.elements()
            }
        return self._point_contributions

    def add(self, lhs: tuple[int, ...], rhs: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((a + b) % self.prime for a, b in zip(lhs, rhs))

    def int_add(self, lhs: int, rhs: int) -> int:
        result = 0
        place = 1
        for _ in range(self.degree):
            result += ((lhs % self.prime) + (rhs % self.prime)) % self.prime * place
            lhs //= self.prime
            rhs //= self.prime
            place *= self.prime
        return result

    def multiply(self, lhs: tuple[int, ...], rhs: tuple[int, ...]) -> tuple[int, ...]:
        product_coefficients = [0] * (2 * self.degree - 1)
        for i, lhs_coefficient in enumerate(lhs):
            for j, rhs_coefficient in enumerate(rhs):
                product_coefficients[i + j] = (product_coefficients[i + j] + lhs_coefficient * rhs_coefficient) % self.prime

        for d in range(len(product_coefficients) - 1, self.degree - 1, -1):
            coefficient = product_coefficients[d]
            if coefficient == 0:
                continue
            for j in range(self.degree):
                index = d - self.degree + j
                product_coefficients[index] = (product_coefficients[index] - coefficient * self.modulus[j]) % self.prime

        return tuple(product_coefficients[: self.degree])

    def int_multiply(self, lhs: int, rhs: int) -> int:
        product_coefficients = [0] * (2 * self.degree - 1)
        lhs_working = lhs
        for i in range(self.degree):
            lhs_coefficient = lhs_working % self.prime
            lhs_working //= self.prime
            if lhs_coefficient == 0:
                continue
            rhs_working = rhs
            for j in range(self.degree):
                rhs_coefficient = rhs_working % self.prime
                rhs_working //= self.prime
                if rhs_coefficient:
                    product_coefficients[i + j] = (
                        product_coefficients[i + j] + lhs_coefficient * rhs_coefficient
                    ) % self.prime

        for d in range(len(product_coefficients) - 1, self.degree - 1, -1):
            coefficient = product_coefficients[d]
            if coefficient == 0:
                continue
            for j in range(self.degree):
                index = d - self.degree + j
                product_coefficients[index] = (product_coefficients[index] - coefficient * self.modulus[j]) % self.prime

        result = 0
        place = 1
        for coefficient in product_coefficients[: self.degree]:
            result += coefficient * place
            place *= self.prime
        return result

    def pow(self, base: tuple[int, ...], exponent: int) -> tuple[int, ...]:
        result = self.one()
        while exponent > 0:
            if exponent % 2 == 1:
                result = self.multiply(result, base)
            base = self.multiply(base, base)
            exponent //= 2
        return result

    def int_pow(self, base: int, exponent: int) -> int:
        result = self.int_one()
        while exponent > 0:
            if exponent % 2 == 1:
                result = self.int_multiply(result, base)
            base = self.int_multiply(base, base)
            exponent //= 2
        return result

    def int_evaluate_polynomial(self, coefficients: tuple[int, ...], x: int) -> int:
        result = self.int_zero()
        for coefficient in reversed(coefficients):
            result = self.int_add(self.int_multiply(result, x), self.int_constant(coefficient))
        return result

    def int_quadratic_character_table(self) -> tuple[int, ...]:
        if self._int_quadratic_characters is None:
            characters = [-1] * self.size
            characters[self.int_zero()] = 0
            for element in self.int_elements():
                if element != self.int_zero():
                    characters[self.int_multiply(element, element)] = 1
            self._int_quadratic_characters = tuple(characters)
        return self._int_quadratic_characters

    def int_quadratic_character(self, element: int) -> int:
        return self.int_quadratic_character_table()[element]

    def frobenius_matrix(self) -> tuple[tuple[int, ...], ...]:
        if self._frobenius_matrix is None:
            basis_images = []
            for index in range(self.degree):
                basis_element = tuple(1 if index == j else 0 for j in range(self.degree))
                basis_images.append(self.pow(basis_element, self.prime))
            self._frobenius_matrix = tuple(basis_images)
        return self._frobenius_matrix

    def apply_frobenius(self, element: tuple[int, ...]) -> tuple[int, ...]:
        result = [0] * self.degree
        for coefficient, image in zip(element, self.frobenius_matrix()):
            if coefficient == 0:
                continue
            for index, image_coefficient in enumerate(image):
                result[index] = (result[index] + coefficient * image_coefficient) % self.prime
        return tuple(result)

    @staticmethod
    def _find_irreducible_polynomial(prime: int, degree: int) -> tuple[int, ...]:
        if degree == 1:
            return (0, 1)

        for low_coefficients in product(range(prime), repeat=degree):
            coefficients = (*low_coefficients, 1)
            if _is_irreducible_by_rabin_test(coefficients, prime):
                return coefficients

        raise RuntimeError("failed to find irreducible polynomial for finite extension")


def _prime_divisors(value: int) -> tuple[int, ...]:
    divisors = []
    remaining = value
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            divisors.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1 if divisor == 2 else 2
    if remaining > 1:
        divisors.append(remaining)
    return tuple(divisors)


def _is_irreducible_by_rabin_test(polynomial: tuple[int, ...], prime: int) -> bool:
    polynomial = _poly_monic_mod(polynomial, prime)
    degree = len(polynomial) - 1
    if degree <= 0:
        return False
    if degree == 1:
        return True

    x = (0, 1)
    for divisor in _prime_divisors(degree):
        power = _poly_pow_remainder_mod(x, prime ** (degree // divisor), polynomial, prime)
        if _poly_gcd_mod(polynomial, _poly_sub_mod(power, x, prime), prime) != (1,):
            return False

    return _poly_sub_mod(_poly_pow_remainder_mod(x, prime**degree, polynomial, prime), x, prime) == (0,)


def _is_irreducible(polynomial: Polynomial) -> bool:
    if polynomial.degree <= 0:
        return False
    field = polynomial.field
    for divisor_degree in range(1, polynomial.degree // 2 + 1):
        for low_coefficients in product(range(field.prime), repeat=divisor_degree):
            divisor = Polynomial(field, (*low_coefficients, 1))
            if polynomial.remainder(divisor).is_zero():
                return False
    return True


def _rank_mod_prime(rows: list[list[int]], prime: int) -> int:
    if not rows:
        return 0

    working = [row[:] for row in rows]
    row_count = len(working)
    column_count = len(working[0])
    rank = 0
    for column in range(column_count):
        pivot = None
        for row in range(rank, row_count):
            if working[row][column] % prime != 0:
                pivot = row
                break
        if pivot is None:
            continue

        working[rank], working[pivot] = working[pivot], working[rank]
        inverse = pow(working[rank][column] % prime, -1, prime)
        working[rank] = [(entry * inverse) % prime for entry in working[rank]]
        for row in range(row_count):
            if row == rank:
                continue
            scale = working[row][column] % prime
            if scale == 0:
                continue
            working[row] = [(entry - scale * pivot_entry) % prime for entry, pivot_entry in zip(working[row], working[rank])]
        rank += 1
        if rank == row_count:
            break
    return rank


def _normal_basis_generator(extension: FiniteExtension) -> tuple[int, ...]:
    if extension.degree == 1:
        return extension.one()

    for element in extension.elements():
        if extension.is_zero(element):
            continue
        conjugates = []
        current = element
        for _ in range(extension.degree):
            conjugates.append(list(current))
            current = extension.apply_frobenius(current)
        if _rank_mod_prime(conjugates, extension.prime) == extension.degree:
            return element

    raise RuntimeError("failed to find a normal basis generator")


def _frobenius_orbit(extension: FiniteExtension, element: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    orbit = []
    current = element
    for _ in range(extension.degree):
        orbit.append(current)
        current = extension.apply_frobenius(current)
    return tuple(orbit)


def _aperiodic_necklace_representatives(alphabet_size: int, length: int) -> Iterable[tuple[int, ...]]:
    if alphabet_size < 1:
        raise ValueError("alphabet size must be positive")
    if length < 1:
        raise ValueError("length must be positive")

    word = [0] * (length + 1)

    def generate(position: int, period: int) -> Iterable[tuple[int, ...]]:
        if position > length:
            if period == length:
                yield tuple(word[1 : length + 1])
            return

        word[position] = word[position - period]
        yield from generate(position + 1, period)
        for value in range(word[position - period] + 1, alphabet_size):
            word[position] = value
            yield from generate(position + 1, position)

    yield from generate(1, 1)


def _element_from_normal_basis_word(
    extension: FiniteExtension,
    normal_basis_conjugates: tuple[tuple[int, ...], ...],
    word: tuple[int, ...],
) -> tuple[int, ...]:
    result = extension.zero()
    for coefficient, conjugate in zip(word, normal_basis_conjugates):
        if coefficient == 0:
            continue
        result = extension.add(result, tuple(coefficient * entry % extension.prime for entry in conjugate))
    return result


def _extension_polynomial_mul(
    extension: FiniteExtension,
    lhs: tuple[tuple[int, ...], ...],
    rhs: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    result = [extension.zero() for _ in range(len(lhs) + len(rhs) - 1)]
    zero = extension.zero()
    for i, lhs_coefficient in enumerate(lhs):
        if lhs_coefficient == zero:
            continue
        for j, rhs_coefficient in enumerate(rhs):
            if rhs_coefficient == zero:
                continue
            result[i + j] = extension.add(result[i + j], extension.multiply(lhs_coefficient, rhs_coefficient))
    return tuple(result)


def _extension_polynomial_product_tree(
    extension: FiniteExtension,
    factors: tuple[tuple[tuple[int, ...], ...], ...],
) -> tuple[tuple[int, ...], ...]:
    if not factors:
        return (extension.one(),)

    level = list(factors)
    while len(level) > 1:
        next_level = []
        for index in range(0, len(level), 2):
            if index + 1 == len(level):
                next_level.append(level[index])
            else:
                next_level.append(_extension_polynomial_mul(extension, level[index], level[index + 1]))
        level = next_level
    return level[0]


def _minimal_polynomial_from_frobenius_orbit(extension: FiniteExtension, element: tuple[int, ...]) -> tuple[int, ...]:
    factors = tuple(
        (
            tuple((-entry) % extension.prime for entry in conjugate),
            extension.one(),
        )
        for conjugate in _frobenius_orbit(extension, element)
    )
    coefficients = _extension_polynomial_product_tree(extension, factors)

    ground_coefficients = []
    for coefficient in coefficients:
        if any(entry != 0 for entry in coefficient[1:]):
            raise RuntimeError("minimal polynomial coefficient did not descend to the ground field")
        ground_coefficients.append(coefficient[0])
    return tuple(ground_coefficients)


def _generate_monic_irreducible_polynomials(field: PrimeField, degree: int) -> tuple[tuple[int, ...], ...]:
    if degree < 1:
        raise ValueError("degree must be positive")

    if degree == 1:
        return tuple((constant, 1) for constant in range(field.prime))

    extension = FiniteExtension(field.prime, degree)
    normal_generator = _normal_basis_generator(extension)
    normal_basis_conjugates = _frobenius_orbit(extension, normal_generator)
    polynomials = {
        _minimal_polynomial_from_frobenius_orbit(extension, _element_from_normal_basis_word(extension, normal_basis_conjugates, word))
        for word in _aperiodic_necklace_representatives(field.prime, degree)
    }
    expected_count = monic_irreducible_count(field.prime, degree)
    if len(polynomials) != expected_count:
        raise RuntimeError(f"necklace generation produced {len(polynomials)} polynomials, expected {expected_count}")
    return tuple(sorted(polynomials))


def _stream_monic_irreducible_polynomials(field: PrimeField, degree: int) -> Iterable[tuple[int, ...]]:
    if degree < 1:
        raise ValueError("degree must be positive")
    if degree == 1:
        for constant in range(field.prime):
            yield (constant, 1)
        return

    extension = FiniteExtension(field.prime, degree)
    normal_generator = _normal_basis_generator(extension)
    normal_basis_conjugates = _frobenius_orbit(extension, normal_generator)
    for word in _aperiodic_necklace_representatives(field.prime, degree):
        yield _minimal_polynomial_from_frobenius_orbit(
            extension,
            _element_from_normal_basis_word(extension, normal_basis_conjugates, word),
        )


class PointCountingContext:
    def __init__(self, field: PrimeField, polynomial_degree: int) -> None:
        if polynomial_degree < 0:
            raise ValueError("polynomial degree must be nonnegative")
        self.field = field
        self.polynomial_degree = polynomial_degree
        self._extension_cache: dict[int, FiniteExtension] = {}
        self._power_cache: dict[int, tuple[tuple[tuple[int, ...], tuple[tuple[int, ...], ...]], ...]] = {}
        self._ground_powers: Optional[tuple[tuple[int, tuple[int, ...]], ...]] = None
        self._ground_value_table: Optional[tuple[tuple[tuple[int, ...], ...], ...]] = None
        self._ground_quadratic_residues: Optional[set[int]] = None
        self._ground_point_contributions: Optional[tuple[int, ...]] = None

    def require_compatible_polynomial(self, polynomial: Polynomial) -> None:
        if polynomial.field.prime != self.field.prime:
            raise ValueError("point-counting context is over a different field")
        if polynomial.degree > self.polynomial_degree:
            raise ValueError("point-counting context polynomial degree is too small")

    def extension(self, extension_degree: int) -> FiniteExtension:
        if extension_degree < 1:
            raise ValueError("extension degree must be positive")
        extension = self._extension_cache.get(extension_degree)
        if extension is None:
            extension = FiniteExtension(self.field.prime, extension_degree)
            self._extension_cache[extension_degree] = extension
        return extension

    def extension_powers(self, extension_degree: int) -> tuple[tuple[tuple[int, ...], tuple[tuple[int, ...], ...]], ...]:
        cached = self._power_cache.get(extension_degree)
        if cached is not None:
            return cached

        extension = self.extension(extension_degree)
        one = extension.one()
        rows = []
        for x in extension.elements():
            powers = [one]
            for _ in range(self.polynomial_degree):
                powers.append(extension.multiply(powers[-1], x))
            rows.append((x, tuple(powers)))

        cached = tuple(rows)
        self._power_cache[extension_degree] = cached
        return cached

    def ground_powers(self) -> tuple[tuple[int, tuple[int, ...]], ...]:
        if self._ground_powers is not None:
            return self._ground_powers

        rows = []
        for x in range(self.field.prime):
            powers = [1]
            for _ in range(self.polynomial_degree):
                powers.append((powers[-1] * x) % self.field.prime)
            rows.append((x, tuple(powers)))

        self._ground_powers = tuple(rows)
        return self._ground_powers

    def ground_value_table(self) -> tuple[tuple[tuple[int, ...], ...], ...]:
        if self._ground_value_table is not None:
            return self._ground_value_table

        powers_by_x = self.ground_powers()
        table = []
        for degree in range(self.polynomial_degree + 1):
            values_by_coefficient = []
            for coefficient in range(self.field.prime):
                values_by_coefficient.append(
                    tuple(coefficient * powers[degree] % self.field.prime for _, powers in powers_by_x)
                )
            table.append(tuple(values_by_coefficient))

        self._ground_value_table = tuple(table)
        return self._ground_value_table

    def ground_quadratic_residues(self) -> set[int]:
        if self._ground_quadratic_residues is None:
            self._ground_quadratic_residues = {y * y % self.field.prime for y in range(self.field.prime)}
        return self._ground_quadratic_residues

    def ground_point_contributions(self) -> tuple[int, ...]:
        if self._ground_point_contributions is None:
            contributions = [0] * self.field.prime
            contributions[0] = 1
            for square in self.ground_quadratic_residues():
                if square != 0:
                    contributions[square] = 2
            self._ground_point_contributions = tuple(contributions)
        return self._ground_point_contributions


@dataclass(frozen=True)
class HyperellipticCurve:
    defining_polynomial: Polynomial
    point_counting_context: Optional[PointCountingContext] = field(default=None, repr=False, compare=False)
    _hasse_witt_matrix: Optional[tuple[tuple[int, ...], ...]] = field(default=None, init=False, repr=False, compare=False)
    _l_polynomial_coefficients_mod_p: Optional[tuple[int, ...]] = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        self._validate_model()
        if self.point_counting_context is None:
            object.__setattr__(self, "point_counting_context", PointCountingContext(self.field, self.degree))
        self.point_counting_context.require_compatible_polynomial(self.defining_polynomial)

    @property
    def field(self) -> PrimeField:
        return self.defining_polynomial.field

    @property
    def degree(self) -> int:
        return self.defining_polynomial.degree

    @property
    def genus(self) -> int:
        return (self.degree - 1) // 2

    @property
    def degree_model(self) -> str:
        return "odd" if self.degree % 2 == 1 else "even"

    def is_monic_model(self) -> bool:
        return self.defining_polynomial.is_monic()

    def point_count_over_extension(self, extension_degree: int) -> int:
        if extension_degree < 1:
            raise ValueError("extension degree must be positive")
        if extension_degree == 1:
            return self.point_count_over_ground_field()

        extension = self.point_counting_context.extension(extension_degree)
        count = self._points_at_infinity(extension)

        point_contributions = extension.point_contributions()
        for x, powers in self.point_counting_context.extension_powers(extension_degree):
            rhs = self._evaluate_defining_polynomial_from_powers(extension, powers)
            count += point_contributions[rhs]

        return count

    def point_count_over_ground_field(self) -> int:
        p = self.field.prime
        count = 1 if self.degree_model == "odd" else self._points_at_infinity_over_ground_field()
        point_contributions = self.point_counting_context.ground_point_contributions()
        coefficients = self.defining_polynomial.coefficients
        value_table = self.point_counting_context.ground_value_table()

        for x_index in range(p):
            value = sum(value_table[degree][coefficient][x_index] for degree, coefficient in enumerate(coefficients)) % p
            count += point_contributions[value]

        return count

    def l_polynomial_coefficients(self) -> list[int]:
        coefficients = self._compute_l_polynomial_coefficients(max_sparsity=None)
        if coefficients is None:
            raise RuntimeError("unlimited L-polynomial computation unexpectedly failed")
        return coefficients

    def l_polynomial_coefficients_with_sparsity_limit(self, max_sparsity: int) -> Optional[list[int]]:
        if max_sparsity < 0:
            raise ValueError("max sparsity must be nonnegative")
        if not self.passes_hasse_witt_sparsity_filter(max_sparsity):
            return None
        return self._compute_l_polynomial_coefficients(max_sparsity=max_sparsity)

    def hasse_witt_matrix(self) -> tuple[tuple[int, ...], ...]:
        if self._hasse_witt_matrix is not None:
            return self._hasse_witt_matrix

        p = self.field.prime
        h = _poly_pow_mod(self.defining_polynomial.coefficients, (p - 1) // 2, p)
        rows = []
        for i in range(1, self.genus + 1):
            row = []
            for j in range(1, self.genus + 1):
                coefficient_index = p * i - j
                row.append(h[coefficient_index] if coefficient_index < len(h) else 0)
            rows.append(tuple(row))

        matrix = tuple(rows)
        object.__setattr__(self, "_hasse_witt_matrix", matrix)
        return matrix

    def l_polynomial_coefficients_mod_p(self) -> list[int]:
        if self._l_polynomial_coefficients_mod_p is not None:
            return list(self._l_polynomial_coefficients_mod_p)

        p = self.field.prime
        coefficients = tuple(_l_polynomial_coefficients_mod_p_from_hasse_witt_matrix(self.hasse_witt_matrix(), p))
        object.__setattr__(self, "_l_polynomial_coefficients_mod_p", coefficients)
        return list(coefficients)

    def passes_hasse_witt_sparsity_filter(self, max_sparsity: int) -> bool:
        if max_sparsity < 0:
            raise ValueError("max sparsity must be nonnegative")
        sparsity_mod_p = sum(1 for coefficient in self.l_polynomial_coefficients_mod_p()[:-1] if coefficient != 0)
        return sparsity_mod_p <= max_sparsity

    def _compute_l_polynomial_coefficients(self, max_sparsity: Optional[int]) -> Optional[list[int]]:
        power_sums = [0] * (self.genus + 1)
        coefficients = [0] * (self.genus + 1)
        coefficients[0] = 1
        sparsity = 0
        q = 1

        for k in range(1, self.genus + 1):
            q *= self.field.prime
            power_sums[k] = q + 1 - self.point_count_over_extension(k)
            total = sum(coefficients[k - i] * power_sums[i] for i in range(1, k + 1))
            if total % k != 0:
                raise RuntimeError("Newton identity produced a nonintegral coefficient")
            coefficients[k] = -total // k

            if k < self.genus and coefficients[k] != 0:
                sparsity += 1
                if max_sparsity is not None and sparsity > max_sparsity:
                    return None

        return coefficients[1:]

    def _evaluate_defining_polynomial(self, extension: FiniteExtension, x: tuple[int, ...]) -> tuple[int, ...]:
        result = extension.zero()
        for i in range(self.degree, -1, -1):
            result = extension.add(
                extension.multiply(result, x),
                extension.constant(self.defining_polynomial.coefficient(i)),
            )
        return result

    def _evaluate_defining_polynomial_from_powers(
        self,
        extension: FiniteExtension,
        powers: tuple[tuple[int, ...], ...],
    ) -> tuple[int, ...]:
        result = [0] * extension.degree
        for coefficient, power in zip(self.defining_polynomial.coefficients, powers):
            if coefficient == 0:
                continue
            for i, power_coefficient in enumerate(power):
                result[i] = (result[i] + coefficient * power_coefficient) % extension.prime
        return tuple(result)

    def _points_at_infinity(self, extension: FiniteExtension) -> int:
        if self.degree_model == "odd":
            return 1

        leading_coefficient = extension.constant(self.defining_polynomial.coefficient(self.degree))
        return 2 if extension.is_square(leading_coefficient) else 0

    def _points_at_infinity_over_ground_field(self) -> int:
        leading_coefficient = self.defining_polynomial.coefficient(self.degree)
        return 2 if leading_coefficient in self.point_counting_context.ground_quadratic_residues() else 0

    def _validate_model(self) -> None:
        if self.defining_polynomial.is_zero():
            raise ValueError("defining polynomial must be nonzero")
        if self.degree < 3:
            raise ValueError("defining polynomial degree must be at least 3")
        if not self.defining_polynomial.is_squarefree():
            raise ValueError("defining polynomial must be squarefree")


class EnumerationContext:
    def __init__(
        self,
        prime: int,
        genus: int,
        sqlite_path: str | Path | None = None,
        skip_even_models_with_rational_branch_point: bool = False,
        irreducible_cache_path: str | Path | None = None,
        irreducible_memory_budget_mb: int = DEFAULT_IRREDUCIBLE_MEMORY_BUDGET_MB,
        cache_full_orbits: bool = True,
        canonicalize_branch_before_exact: bool = True,
        sqlite_write_batch_size: int = 1,
    ) -> None:
        if genus < 1:
            raise ValueError("genus must be positive")
        self._created_at = perf_counter()
        self._processed_polynomials = 0
        self._processing_seconds = 0.0
        self._sqlite_load_seconds = 0.0
        self._sqlite_write_seconds = 0.0
        self._hasse_witt_seconds = 0.0
        self._factorized_pgl2_seconds = 0.0
        self._expansion_seconds = 0.0
        self._exact_lpoly_seconds = 0.0
        self._ground_invariant_seconds = 0.0
        self._status_counts: dict[str, int] = {}
        self.field = PrimeField(prime)
        self.genus = genus
        self.binary_degree = 2 * genus + 2
        self.skip_even_models_with_rational_branch_point = skip_even_models_with_rational_branch_point
        self.cache_full_orbits = cache_full_orbits
        self.canonicalize_branch_before_exact = canonicalize_branch_before_exact
        self.sqlite_write_batch_size = sqlite_write_batch_size
        self.point_counting_context = PointCountingContext(self.field, self.binary_degree)
        self.pgl2 = _precompute_pgl2(prime)
        self.pgl2_action_matrices = tuple(_pgl2_action_matrix(matrix, self.binary_degree, prime) for matrix in self.pgl2)
        self.seen_keys: set[tuple[int, ...]] = set()
        self.canonical_key_cache: dict[tuple[int, ...], tuple[int, ...]] = {}
        self.canonical_records: dict[tuple[int, ...], CanonicalRecord] = {}
        self.index_by_rational_branch_count: dict[int, set[tuple[int, ...]]] = {}
        self.l_polynomial_mod_p_cache: dict[tuple[int, ...], list[int]] = {}
        self.exact_l_polynomial_cache: dict[tuple[tuple[int, ...], Optional[int]], Optional[list[int]]] = {}
        self._irreducible_factor_cache: dict[int, tuple[tuple[int, ...], ...]] = {}
        self.branch_orbit_cache: dict[tuple[int, ...], tuple[int, ...]] = {}
        self.branch_canonical_to_curve_key: dict[tuple[int, ...], tuple[int, ...]] = {}
        self.branch_factor_action_matrices: dict[int, tuple[tuple[tuple[int, ...], ...], ...]] = {}
        self.branch_factor_transform_cache: OrderedDict[tuple[int, tuple[int, ...]], Optional[tuple[int, tuple[int, ...]]]] = OrderedDict()
        self.branch_factor_ids: dict[tuple[int, ...], int] = {}
        self._next_branch_factor_id = 1
        self.branch_factor_character_cache: dict[tuple[int, int], tuple[int, ...]] = {}
        self.branch_product_character_cache: OrderedDict[tuple[int, tuple[int, ...]], tuple[int, ...]] = OrderedDict()
        self.branch_factor_hasse_witt_power_cache: dict[tuple[int, int], tuple[int, ...]] = {}
        self.branch_hasse_witt_product_cache: OrderedDict[tuple[int, ...], tuple[int, ...]] = OrderedDict()
        self.branch_result_cache: dict[tuple[tuple[int, ...], Optional[int]], dict[str, object]] = {}
        self.branch_orbit_cache_by_pattern: dict[tuple[int, ...], dict[tuple[int, ...], tuple[int, ...]]] = {}
        self.branch_canonical_to_curve_key_by_pattern: dict[tuple[int, ...], dict[tuple[int, ...], tuple[int, ...]]] = {}
        self.irreducible_memory_budget_bytes = _memory_budget_bytes(irreducible_memory_budget_mb)
        self._irreducible_memory_bytes = 0
        self.irreducible_cache_path = Path(irreducible_cache_path) if irreducible_cache_path is not None else default_irreducible_cache_path(prime)
        self.irreducible_cache: Optional[IrreduciblePolynomialCache] = None
        self._sage_irreducible_process: Optional[subprocess.Popen[str]] = None
        self.sqlite_path = Path(sqlite_path) if sqlite_path is not None else None
        self._sqlite_max_sparsity: Optional[int] = _parse_max_sparsity_from_sqlite_path(self.sqlite_path)
        self.sqlite_connection: Optional[sqlite3.Connection] = None
        self._sqlite_pending_writes = 0
        self._orbit_cache_columns: Optional[frozenset[str]] = None
        if self.sqlite_path is not None:
            self.open_sqlite(self.sqlite_path)

    def open_sqlite(self, sqlite_path: str | Path) -> None:
        self.sqlite_path = Path(sqlite_path)
        self._sqlite_max_sparsity = _parse_max_sparsity_from_sqlite_path(self.sqlite_path)
        self.sqlite_connection = sqlite3.connect(self.sqlite_path)
        self._orbit_cache_columns = None
        self.sqlite_connection.execute("PRAGMA journal_mode=WAL")
        self.sqlite_connection.execute("PRAGMA synchronous=NORMAL")
        self._initialize_sqlite_schema()
        started_at = perf_counter()
        self._load_sqlite_records()
        self._sqlite_load_seconds += perf_counter() - started_at

    def close_sqlite(self) -> None:
        if self.sqlite_connection is not None:
            self._commit_sqlite_writes(force=True)
            self.sqlite_connection.close()
            self.sqlite_connection = None

    def _commit_sqlite_writes(self, *, force: bool = False) -> None:
        if self.sqlite_connection is None:
            return
        if not force and self._sqlite_pending_writes < self.sqlite_write_batch_size:
            return
        if self._sqlite_pending_writes == 0 and not force:
            return
        self.sqlite_connection.commit()
        self._sqlite_pending_writes = 0

    def _mark_sqlite_write(self, *, force: bool = False) -> None:
        self._sqlite_pending_writes += 1
        self._commit_sqlite_writes(force=force)

    def close_irreducible_cache(self) -> None:
        if self.irreducible_cache is not None:
            self.irreducible_cache.close()
            self.irreducible_cache = None
        self.close_sage_irreducible_process()

    def close_sage_irreducible_process(self) -> None:
        process = self._sage_irreducible_process
        if process is None:
            return
        if process.stdin is not None:
            process.stdin.close()
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        self._sage_irreducible_process = None

    def elapsed_seconds(self) -> float:
        return perf_counter() - self._created_at

    def timing_summary(self) -> dict[str, object]:
        return {
            "elapsed_seconds": self.elapsed_seconds(),
            "processed_polynomials": self._processed_polynomials,
            "processing_seconds": self._processing_seconds,
            "sqlite_load_seconds": self._sqlite_load_seconds,
            "sqlite_write_seconds": self._sqlite_write_seconds,
            "hasse_witt_seconds": self._hasse_witt_seconds,
            "factorized_pgl2_seconds": self._factorized_pgl2_seconds,
            "expansion_seconds": self._expansion_seconds,
            "exact_lpoly_seconds": self._exact_lpoly_seconds,
            "ground_invariant_seconds": self._ground_invariant_seconds,
            "other_seconds": max(
                0.0,
                self.elapsed_seconds()
                - self._processing_seconds
                - self._sqlite_load_seconds
                - self._sqlite_write_seconds,
            ),
            "status_counts": dict(self._status_counts),
        }

    def write_enumeration_summary(
        self,
        *,
        max_sparsity: Optional[int],
        degree_model: str,
        enumeration_mode: str,
        limit: Optional[int],
        total_coefficient_vectors: int,
        processed: int,
        skipped: int,
        final_position: int,
        sparse_presentations: int,
    ) -> None:
        if self.sqlite_connection is None:
            return

        timing = self.timing_summary()
        started_at = perf_counter()
        self.sqlite_connection.execute(
            """
            INSERT OR REPLACE INTO enumeration_summary (
                id,
                prime,
                genus,
                max_sparsity,
                hasse_witt_prefilter,
                degree_model,
                enumeration_mode,
                leading_coefficient_policy,
                limit_count,
                total_coefficient_vectors,
                processed,
                skipped,
                final_position,
                sparse_presentations,
                sparse_isomorphism_classes,
                canonicalized_isomorphism_classes,
                elapsed_seconds,
                processing_seconds,
                sqlite_load_seconds,
                sqlite_write_seconds,
                other_seconds,
                status_counts
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                1,
                self.field.prime,
                self.genus,
                max_sparsity,
                int(max_sparsity is not None),
                degree_model,
                enumeration_mode,
                "odd:monic;even:monic-and-smallest-nonsquare",
                limit,
                total_coefficient_vectors,
                processed,
                skipped,
                final_position,
                sparse_presentations,
                sparse_isomorphism_classes(self, max_sparsity),
                len(self.canonical_records),
                timing["elapsed_seconds"],
                timing["processing_seconds"],
                timing["sqlite_load_seconds"],
                timing["sqlite_write_seconds"],
                timing["other_seconds"],
                json.dumps(timing["status_counts"], sort_keys=True, separators=(",", ":")),
            ),
        )
        self._mark_sqlite_write(force=True)
        self._sqlite_write_seconds += perf_counter() - started_at

    def write_progress_snapshot(
        self,
        *,
        position: int,
        processed: int,
        skipped: int,
        max_sparsity: Optional[int],
        previous_sparse_presentations: int,
        previous_sparse_isomorphism_classes: int,
        previous_canonicalized_isomorphism_classes: int,
    ) -> tuple[int, int, int]:
        if self.sqlite_connection is None:
            return (
                previous_sparse_presentations,
                previous_sparse_isomorphism_classes,
                previous_canonicalized_isomorphism_classes,
            )

        sparse_presentations = sparse_presentations_by_orbit_size(self, max_sparsity)
        sparse_classes = sparse_isomorphism_classes(self, max_sparsity)
        canonicalized_classes = len(self.canonical_records)
        timing = self.timing_summary()
        started_at = perf_counter()
        self.sqlite_connection.execute(
            """
            INSERT OR REPLACE INTO enumeration_progress (
                position,
                processed,
                skipped,
                sparse_presentations,
                sparse_isomorphism_classes,
                canonicalized_isomorphism_classes,
                delta_sparse_presentations,
                delta_sparse_isomorphism_classes,
                delta_canonicalized_isomorphism_classes,
                status_counts,
                elapsed_seconds
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                position,
                processed,
                skipped,
                sparse_presentations,
                sparse_classes,
                canonicalized_classes,
                sparse_presentations - previous_sparse_presentations,
                sparse_classes - previous_sparse_isomorphism_classes,
                canonicalized_classes - previous_canonicalized_isomorphism_classes,
                json.dumps(timing["status_counts"], sort_keys=True, separators=(",", ":")),
                timing["elapsed_seconds"],
            ),
        )
        self._mark_sqlite_write(force=True)
        self._sqlite_write_seconds += perf_counter() - started_at
        return sparse_presentations, sparse_classes, canonicalized_classes

    def _initialize_sqlite_schema(self) -> None:
        if self.sqlite_connection is None:
            raise ValueError("sqlite connection is not open")
        self.sqlite_connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS curve_cache (
                canonical_key BLOB NOT NULL,
                rational_branch_count INTEGER NOT NULL,
                coefficients TEXT NOT NULL,
                lpoly_mod_p BLOB,
                exact_lpoly BLOB,
                sparsity INTEGER,
                status TEXT NOT NULL,
                PRIMARY KEY (canonical_key)
            );

            CREATE TABLE IF NOT EXISTS sparse_curves (
                canonical_key BLOB NOT NULL,
                coefficients TEXT NOT NULL,
                branch_factors TEXT,
                branch_infinity_branch INTEGER,
                branch_leading_coefficient INTEGER,
                branch_factorization_pattern TEXT,
                lpoly TEXT NOT NULL,
                sparsity INTEGER NOT NULL,
                rational_branch_count INTEGER NOT NULL,
                PRIMARY KEY (canonical_key)
            );

            CREATE TABLE IF NOT EXISTS orbit_cache (
                rational_branch_count INTEGER NOT NULL,
                ground_point_count INTEGER,
                hasse_witt_lpoly_mod_p BLOB NOT NULL,
                orbit_key BLOB NOT NULL,
                canonical_key BLOB NOT NULL,
                PRIMARY KEY (rational_branch_count, ground_point_count, hasse_witt_lpoly_mod_p, orbit_key)
            );

            CREATE TABLE IF NOT EXISTS branch_orbit_cache (
                orbit_branch_key BLOB NOT NULL,
                canonical_branch_key BLOB NOT NULL,
                PRIMARY KEY (orbit_branch_key)
            );

            CREATE TABLE IF NOT EXISTS branch_curve_cache (
                canonical_branch_key BLOB NOT NULL,
                canonical_key BLOB NOT NULL,
                PRIMARY KEY (canonical_branch_key)
            );

            CREATE TABLE IF NOT EXISTS branch_factor_transform_cache (
                degree INTEGER NOT NULL,
                factor_key BLOB NOT NULL,
                matrix_index INTEGER NOT NULL,
                transformed_factor_key BLOB,
                PRIMARY KEY (degree, factor_key, matrix_index)
            );

            CREATE TABLE IF NOT EXISTS enumeration_summary (
                id INTEGER PRIMARY KEY CHECK (id = 1),
                prime INTEGER NOT NULL,
                genus INTEGER NOT NULL,
                max_sparsity INTEGER,
                hasse_witt_prefilter INTEGER NOT NULL,
                degree_model TEXT NOT NULL,
                enumeration_mode TEXT NOT NULL,
                leading_coefficient_policy TEXT NOT NULL,
                limit_count INTEGER,
                total_coefficient_vectors INTEGER NOT NULL,
                processed INTEGER NOT NULL,
                skipped INTEGER NOT NULL,
                final_position INTEGER NOT NULL,
                sparse_presentations INTEGER NOT NULL,
                sparse_isomorphism_classes INTEGER NOT NULL,
                canonicalized_isomorphism_classes INTEGER NOT NULL,
                elapsed_seconds REAL NOT NULL,
                processing_seconds REAL NOT NULL,
                sqlite_load_seconds REAL NOT NULL,
                sqlite_write_seconds REAL NOT NULL,
                other_seconds REAL NOT NULL,
                status_counts TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS enumeration_progress (
                position INTEGER PRIMARY KEY,
                processed INTEGER NOT NULL,
                skipped INTEGER NOT NULL,
                sparse_presentations INTEGER NOT NULL,
                sparse_isomorphism_classes INTEGER NOT NULL,
                canonicalized_isomorphism_classes INTEGER NOT NULL,
                delta_sparse_presentations INTEGER NOT NULL,
                delta_sparse_isomorphism_classes INTEGER NOT NULL,
                delta_canonicalized_isomorphism_classes INTEGER NOT NULL,
                status_counts TEXT NOT NULL,
                elapsed_seconds REAL NOT NULL
            );
            """
        )
        self._ensure_sqlite_column("sparse_curves", "branch_factors", "TEXT")
        self._ensure_sqlite_column("sparse_curves", "branch_infinity_branch", "INTEGER")
        self._ensure_sqlite_column("sparse_curves", "branch_leading_coefficient", "INTEGER")
        self._ensure_sqlite_column("sparse_curves", "branch_factorization_pattern", "TEXT")
        self._mark_sqlite_write(force=True)

    def _ensure_sqlite_column(self, table: str, column: str, column_type: str) -> None:
        if self.sqlite_connection is None:
            raise ValueError("sqlite connection is not open")
        columns = {row[1] for row in self.sqlite_connection.execute(f"PRAGMA table_info({table})").fetchall()}
        if column not in columns:
            self.sqlite_connection.execute(f"ALTER TABLE {table} ADD COLUMN {column} {column_type}")

    def _orbit_cache_has_ground_point_count(self) -> bool:
        return "ground_point_count" in self._get_orbit_cache_columns()

    def _orbit_cache_has_hasse_witt_lpoly_mod_p(self) -> bool:
        return "hasse_witt_lpoly_mod_p" in self._get_orbit_cache_columns()

    def _get_orbit_cache_columns(self) -> frozenset[str]:
        if self.sqlite_connection is None:
            return frozenset()
        if self._orbit_cache_columns is None:
            self._orbit_cache_columns = frozenset(
                row[1]
                for row in self.sqlite_connection.execute("PRAGMA table_info(orbit_cache)").fetchall()
            )
        return self._orbit_cache_columns

    def _load_sqlite_records(self) -> None:
        if self.sqlite_connection is None:
            raise ValueError("sqlite connection is not open")
        self._load_branch_factorized_caches()

        columns = {
            row[1]
            for row in self.sqlite_connection.execute("PRAGMA table_info(curve_cache)").fetchall()
        }
        has_hasse_witt_column = "hasse_witt_lpoly_mod_p" in columns
        has_max_sparsity_column = "max_sparsity" in columns
        orbit_sizes_by_canonical_key = self._load_orbit_sizes_by_canonical_key()
        rows = self.sqlite_connection.execute(
            f"""
            SELECT
                canonical_key,
                rational_branch_count,
                coefficients,
                lpoly_mod_p,
                exact_lpoly,
                sparsity,
                status
                {', max_sparsity' if has_max_sparsity_column else ''}
                {', hasse_witt_lpoly_mod_p' if has_hasse_witt_column else ''}
            FROM curve_cache
            """
        ).fetchall()

        for row in rows:
            canonical_key_value = row[0]
            rational_branch_count = row[1]
            coefficients_text = row[2]
            lpoly_mod_p_value = row[3]
            exact_lpoly_value = row[4]
            sparsity = row[5]
            status = row[6]
            next_index = 7
            max_sparsity = row[next_index] if has_max_sparsity_column else self._sqlite_max_sparsity
            next_index += 1 if has_max_sparsity_column else 0
            hasse_witt_lpoly_mod_p_value = row[next_index] if has_hasse_witt_column else None
            canonical_key = self._unpack_field_tuple_blob(canonical_key_value)
            coefficients = _unpack_int_tuple(coefficients_text)
            hasse_witt_lpoly_mod_p = _unpack_general_int_tuple_blob(hasse_witt_lpoly_mod_p_value)
            lpoly_mod_p = _unpack_general_int_tuple_blob(lpoly_mod_p_value)
            exact_lpoly = _unpack_general_int_tuple_blob(exact_lpoly_value)
            try:
                ground_point_count = self.ground_invariants(coefficients)[1]
            except ValueError:
                ground_point_count = None

            record = self._register_or_update_record(
                canonical_key=canonical_key,
                coefficients=coefficients,
                rational_branch_count=rational_branch_count,
                ground_point_count=ground_point_count,
                hasse_witt_lpoly_mod_p=hasse_witt_lpoly_mod_p or lpoly_mod_p,
                orbit_size=orbit_sizes_by_canonical_key.get(canonical_key, 1),
            )
            record.status_by_max_sparsity[max_sparsity] = status
            record.sparsity_by_max_sparsity[max_sparsity] = sparsity
            record.exact_lpoly_by_max_sparsity[max_sparsity] = exact_lpoly

            self.seen_keys.add(canonical_key)
            self.canonical_key_cache[coefficients] = canonical_key
            if lpoly_mod_p is not None:
                self.l_polynomial_mod_p_cache[canonical_key] = list(lpoly_mod_p)
            if exact_lpoly is not None:
                self.exact_l_polynomial_cache[(canonical_key, max_sparsity)] = list(exact_lpoly)
            elif status in {"rejected_hasse_witt", "rejected_exact"}:
                self.exact_l_polynomial_cache[(canonical_key, max_sparsity)] = None

    def _load_branch_factorized_caches(self) -> None:
        if self.sqlite_connection is None:
            return
        try:
            branch_orbit_rows = self.sqlite_connection.execute(
                "SELECT orbit_branch_key, canonical_branch_key FROM branch_orbit_cache"
            ).fetchall()
            branch_curve_rows = self.sqlite_connection.execute(
                "SELECT canonical_branch_key, canonical_key FROM branch_curve_cache"
            ).fetchall()
        except sqlite3.OperationalError:
            return

        self.branch_orbit_cache.update(
            {
                self._unpack_general_tuple_key(orbit_key): self._unpack_general_tuple_key(canonical_key)
                for orbit_key, canonical_key in branch_orbit_rows
            }
        )
        self.branch_canonical_to_curve_key.update(
            {
                self._unpack_general_tuple_key(branch_key): self._unpack_field_tuple_blob(curve_key)
                for branch_key, curve_key in branch_curve_rows
            }
        )

    def _load_orbit_sizes_by_canonical_key(self) -> dict[tuple[int, ...], int]:
        if self.sqlite_connection is None:
            return {}
        try:
            rows = self.sqlite_connection.execute(
                """
                SELECT canonical_key, COUNT(*)
                FROM orbit_cache
                GROUP BY canonical_key
                """
            ).fetchall()
        except sqlite3.OperationalError:
            return {}
        return {
            self._unpack_field_tuple_blob(canonical_key_value): int(orbit_size)
            for canonical_key_value, orbit_size in rows
        }

    def _pack_field_tuple_blob(self, values: tuple[int, ...]) -> bytes:
        if self.field.prime < 256:
            return bytes(values)
        return json.dumps(list(values), separators=(",", ":")).encode("ascii")

    def _unpack_field_tuple_blob(self, data: str | bytes) -> tuple[int, ...]:
        if isinstance(data, str):
            return _unpack_int_tuple(data)
        if self.field.prime < 256:
            return tuple(data)
        return tuple(json.loads(data.decode("ascii")))

    def _pack_general_tuple_key(self, values: tuple[int, ...]) -> bytes:
        return _pack_general_int_tuple_blob(values)

    def _unpack_general_tuple_key(self, data: str | bytes) -> tuple[int, ...]:
        unpacked = _unpack_general_int_tuple_blob(data)
        if unpacked is None:
            raise ValueError("expected packed tuple data")
        return unpacked

    def polynomial(self, coefficients: list[int] | tuple[int, ...]) -> Polynomial:
        return Polynomial(self.field, coefficients)

    def curve(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> HyperellipticCurve:
        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        curve = HyperellipticCurve(polynomial, point_counting_context=self.point_counting_context)
        if curve.genus != self.genus:
            raise ValueError("curve genus does not match enumeration context")
        return curve

    def rational_branch_count(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> int:
        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        self._require_compatible_polynomial(polynomial)
        return self.ground_invariants(polynomial)[0]

    def ground_invariants(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> tuple[int, int]:
        started_at = perf_counter()
        try:
            return self._ground_invariants(polynomial_or_coefficients)
        finally:
            self._ground_invariant_seconds += perf_counter() - started_at

    def _ground_invariants(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> tuple[int, int]:
        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        self._require_compatible_polynomial(polynomial)
        rational_branch_count = 1 if polynomial.degree == 2 * self.genus + 1 else 0
        ground_point_count = 1 if polynomial.degree == 2 * self.genus + 1 else 0
        if polynomial.degree == 2 * self.genus + 2:
            leading_coefficient = polynomial.coefficient(polynomial.degree)
            if leading_coefficient in self.point_counting_context.ground_quadratic_residues():
                ground_point_count = 2

        point_contributions = self.point_counting_context.ground_point_contributions()
        for x in range(self.field.prime):
            value = 0
            for coefficient in reversed(polynomial.coefficients):
                value = (value * x + coefficient) % self.field.prime
            if value == 0:
                rational_branch_count += 1
            ground_point_count += point_contributions[value]
        return rational_branch_count, ground_point_count

    def branch_ground_invariants(self, candidate: BranchDivisorCandidate) -> tuple[int, int]:
        started_at = perf_counter()
        try:
            return self._branch_ground_invariants(candidate)
        finally:
            self._ground_invariant_seconds += perf_counter() - started_at

    def _branch_ground_invariants(self, candidate: BranchDivisorCandidate) -> tuple[int, int]:
        p = self.field.prime
        rational_branch_count = 1 if candidate.infinity_branch else 0
        for factor in candidate.factors:
            if len(factor) == 2:
                rational_branch_count += 1

        ground_point_count = 1 if candidate.infinity_branch else 0
        if not candidate.infinity_branch:
            if candidate.leading_coefficient in self.point_counting_context.ground_quadratic_residues():
                ground_point_count = 2
            else:
                ground_point_count = 0

        point_contributions = self.point_counting_context.ground_point_contributions()
        for x, powers in self.point_counting_context.ground_powers():
            value = candidate.leading_coefficient % p
            for factor in candidate.factors:
                factor_value = 0
                for coefficient, power in zip(factor, powers):
                    factor_value = (factor_value + coefficient * power) % p
                value = (value * factor_value) % p
                if value == 0:
                    break
            ground_point_count += point_contributions[value]
        return rational_branch_count, ground_point_count

    def _branch_points_at_infinity(self, candidate: BranchDivisorCandidate, extension: FiniteExtension) -> int:
        if candidate.infinity_branch:
            return 1
        leading = extension.constant(candidate.leading_coefficient)
        return 2 if extension.is_square(leading) else 0

    def _branch_factor_id(self, factor: tuple[int, ...]) -> int:
        factor_id = self.branch_factor_ids.get(factor)
        if factor_id is None:
            factor_id = self._next_branch_factor_id
            self._next_branch_factor_id += 1
            self.branch_factor_ids[factor] = factor_id
        return factor_id

    def _evaluate_branch_factor_from_powers(
        self,
        extension: FiniteExtension,
        factor: tuple[int, ...],
        powers: tuple[tuple[int, ...], ...],
    ) -> tuple[int, ...]:
        result = [0] * extension.degree
        for coefficient, power in zip(factor, powers):
            if coefficient == 0:
                continue
            for i, power_coefficient in enumerate(power):
                result[i] = (result[i] + coefficient * power_coefficient) % extension.prime
        return tuple(result)

    def _ground_quadratic_character(self, value: int) -> int:
        value %= self.field.prime
        if value == 0:
            return 0
        return 1 if value in self.point_counting_context.ground_quadratic_residues() else -1

    def _extension_quadratic_character(self, extension: FiniteExtension, value: tuple[int, ...]) -> int:
        if value == extension.zero():
            return 0
        return 1 if value in extension.squares() else -1

    def _branch_factor_character_values(
        self,
        extension_degree: int,
        factor: tuple[int, ...],
    ) -> tuple[int, ...]:
        factor_id = self._branch_factor_id(factor)
        cache_key = (extension_degree, factor_id)
        cached = self.branch_factor_character_cache.get(cache_key)
        if cached is not None:
            return cached

        if extension_degree == 1:
            values = []
            for _, powers in self.point_counting_context.ground_powers():
                factor_value = 0
                for coefficient, power in zip(factor, powers):
                    factor_value = (factor_value + coefficient * power) % self.field.prime
                values.append(self._ground_quadratic_character(factor_value))
            cached = tuple(values)
        else:
            extension = self.point_counting_context.extension(extension_degree)
            cached = tuple(
                extension.int_quadratic_character(extension.int_evaluate_polynomial(factor, x))
                for x in extension.int_elements()
            )

        self.branch_factor_character_cache[cache_key] = cached
        return cached

    def _branch_product_character_values(
        self,
        extension_degree: int,
        factors: tuple[tuple[int, ...], ...],
    ) -> tuple[int, ...]:
        factor_ids = tuple(sorted(self._branch_factor_id(factor) for factor in factors))
        cache_key = (extension_degree, factor_ids)
        cached = self.branch_product_character_cache.get(cache_key)
        if cached is not None:
            self.branch_product_character_cache.move_to_end(cache_key)
            return cached

        if not factors:
            size = self.field.prime if extension_degree == 1 else self.point_counting_context.extension(extension_degree).size
            values = (1,) * size
        elif len(factors) == 1:
            values = self._branch_factor_character_values(extension_degree, factors[0])
        else:
            sorted_factors = tuple(sorted(factors, key=lambda factor: self._branch_factor_id(factor)))
            midpoint = len(sorted_factors) // 2
            left = self._branch_product_character_values(extension_degree, sorted_factors[:midpoint])
            right = self._branch_product_character_values(extension_degree, sorted_factors[midpoint:])
            values = tuple(a * b for a, b in zip(left, right))

        self.branch_product_character_cache[cache_key] = values
        self.branch_product_character_cache.move_to_end(cache_key)
        while len(self.branch_product_character_cache) > BRANCH_PRODUCT_CHARACTER_CACHE_MAX_ENTRIES:
            self.branch_product_character_cache.popitem(last=False)
        return values

    def _branch_factor_hasse_witt_power(self, factor: tuple[int, ...]) -> tuple[int, ...]:
        exponent = (self.field.prime - 1) // 2
        factor_id = self._branch_factor_id(factor)
        cache_key = (factor_id, exponent)
        cached = self.branch_factor_hasse_witt_power_cache.get(cache_key)
        if cached is None:
            cached = _poly_pow_mod(factor, exponent, self.field.prime)
            self.branch_factor_hasse_witt_power_cache[cache_key] = cached
        return cached

    def _branch_hasse_witt_product(self, factors: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
        factor_ids = tuple(sorted(self._branch_factor_id(factor) for factor in factors))
        cached = self.branch_hasse_witt_product_cache.get(factor_ids)
        if cached is not None:
            self.branch_hasse_witt_product_cache.move_to_end(factor_ids)
            return cached

        if not factors:
            product_polynomial = (1,)
        elif len(factors) == 1:
            product_polynomial = self._branch_factor_hasse_witt_power(factors[0])
        else:
            sorted_factors = tuple(sorted(factors, key=lambda factor: self._branch_factor_id(factor)))
            midpoint = len(sorted_factors) // 2
            left = self._branch_hasse_witt_product(sorted_factors[:midpoint])
            right = self._branch_hasse_witt_product(sorted_factors[midpoint:])
            product_polynomial = _poly_mul_mod(left, right, self.field.prime)

        self.branch_hasse_witt_product_cache[factor_ids] = product_polynomial
        self.branch_hasse_witt_product_cache.move_to_end(factor_ids)
        while len(self.branch_hasse_witt_product_cache) > BRANCH_HASSE_WITT_PRODUCT_CACHE_MAX_ENTRIES:
            self.branch_hasse_witt_product_cache.popitem(last=False)
        return product_polynomial

    def _branch_l_polynomial_coefficients_mod_p(self, candidate: BranchDivisorCandidate) -> list[int]:
        exponent = (self.field.prime - 1) // 2
        leading_factor = pow(candidate.leading_coefficient, exponent, self.field.prime)
        h = self._branch_hasse_witt_product(candidate.factors)
        if leading_factor != 1:
            h = _scale_polynomial_mod(h, leading_factor, self.field.prime)

        rows = []
        for i in range(1, self.genus + 1):
            row = []
            for j in range(1, self.genus + 1):
                coefficient_index = self.field.prime * i - j
                row.append(h[coefficient_index] if coefficient_index < len(h) else 0)
            rows.append(tuple(row))
        return _l_polynomial_coefficients_mod_p_from_hasse_witt_matrix(tuple(rows), self.field.prime)

    def _branch_leading_quadratic_character(self, candidate: BranchDivisorCandidate, extension_degree: int) -> int:
        if extension_degree == 1:
            return self._ground_quadratic_character(candidate.leading_coefficient)
        extension = self.point_counting_context.extension(extension_degree)
        return extension.int_quadratic_character(extension.int_constant(candidate.leading_coefficient))

    def _branch_point_count_over_extension(self, candidate: BranchDivisorCandidate, extension_degree: int) -> int:
        if extension_degree < 1:
            raise ValueError("extension degree must be positive")
        if extension_degree == 1:
            return self.branch_ground_invariants(candidate)[1]

        extension = self.point_counting_context.extension(extension_degree)
        product_characters = self._branch_product_character_values(extension_degree, candidate.factors)
        leading_character = self._branch_leading_quadratic_character(candidate, extension_degree)
        infinity_points = 1 if candidate.infinity_branch else 2 if leading_character == 1 else 0
        return infinity_points + extension.size + leading_character * sum(product_characters)

    def _branch_l_polynomial_coefficients(
        self,
        candidate: BranchDivisorCandidate,
        max_sparsity: Optional[int],
    ) -> Optional[list[int]]:
        power_sums = [0] * (self.genus + 1)
        coefficients = [0] * (self.genus + 1)
        coefficients[0] = 1
        sparsity = 0
        q = 1

        for k in range(1, self.genus + 1):
            q *= self.field.prime
            power_sums[k] = q + 1 - self._branch_point_count_over_extension(candidate, k)
            total = sum(coefficients[k - i] * power_sums[i] for i in range(1, k + 1))
            if total % k != 0:
                raise RuntimeError("Newton identity produced a nonintegral coefficient")
            coefficients[k] = -total // k

            if k < self.genus and coefficients[k] != 0:
                sparsity += 1
                if max_sparsity is not None and sparsity > max_sparsity:
                    return None

        return coefficients[1:]

    def factorization_pattern(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> tuple[int, ...]:
        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        self._require_compatible_polynomial(polynomial)
        remaining = _poly_monic_mod(polynomial.coefficients, self.field.prime)
        pattern = [1] if polynomial.degree == 2 * self.genus + 1 else []
        x_polynomial = (0, 1)
        frobenius_power = x_polynomial
        factor_degree = 1

        while 2 * factor_degree <= len(remaining) - 1:
            frobenius_power = _poly_pow_remainder_mod(frobenius_power, self.field.prime, remaining, self.field.prime)
            degree_factor = _poly_gcd_mod(remaining, _poly_sub_mod(frobenius_power, x_polynomial, self.field.prime), self.field.prime)
            if degree_factor != (1,):
                factor_count = (len(degree_factor) - 1) // factor_degree
                pattern.extend([factor_degree] * factor_count)
                remaining = _poly_exact_quotient_mod(remaining, degree_factor, self.field.prime)
                frobenius_power = _poly_remainder_mod(frobenius_power, remaining, self.field.prime)
            factor_degree += 1

        if len(remaining) > 1:
            pattern.append(len(remaining) - 1)

        return tuple(sorted(pattern))

    def _index_record(self, record: CanonicalRecord) -> None:
        self.index_by_rational_branch_count.setdefault(record.rational_branch_count, set()).add(record.canonical_key)

    def _register_or_update_record(
        self,
        canonical_key: tuple[int, ...],
        coefficients: tuple[int, ...],
        rational_branch_count: int,
        ground_point_count: Optional[int] = None,
        hasse_witt_lpoly_mod_p: Optional[tuple[int, ...]] = None,
        orbit_size: Optional[int] = None,
    ) -> CanonicalRecord:
        record = self.canonical_records.get(canonical_key)
        if record is None:
            record = CanonicalRecord(
                canonical_key=canonical_key,
                coefficients=coefficients,
                rational_branch_count=rational_branch_count,
                ground_point_count=ground_point_count,
                hasse_witt_lpoly_mod_p=hasse_witt_lpoly_mod_p,
                orbit_size=orbit_size or 1,
            )
            self.canonical_records[canonical_key] = record
            self._index_record(record)
            return record

        if record.hasse_witt_lpoly_mod_p is None and hasse_witt_lpoly_mod_p is not None:
            record.hasse_witt_lpoly_mod_p = hasse_witt_lpoly_mod_p
        if record.ground_point_count is None and ground_point_count is not None:
            record.ground_point_count = ground_point_count
        if orbit_size is not None and orbit_size > record.orbit_size:
            record.orbit_size = orbit_size
        return record

    def _result_from_record(self, record: CanonicalRecord, max_sparsity: Optional[int]) -> Optional[dict[str, object]]:
        status = record.status_by_max_sparsity.get(max_sparsity)
        if status is None:
            return None

        result: dict[str, object] = {
            "status": status,
            "canonical_key": record.canonical_key,
            "max_sparsity": max_sparsity,
            "coefficients": record.coefficients,
            "sparsity": record.sparsity_by_max_sparsity.get(max_sparsity),
        }
        lpoly_mod_p = self.l_polynomial_mod_p_cache.get(record.canonical_key)
        if lpoly_mod_p is not None:
            result["lpoly_mod_p"] = list(lpoly_mod_p)
        exact_lpoly = record.exact_lpoly_by_max_sparsity.get(max_sparsity)
        if exact_lpoly is not None:
            result["lpoly"] = list(exact_lpoly)
        return result

    def _branch_key(
        self,
        *,
        leading_coefficient: int,
        factors: Iterable[tuple[int, ...]],
        infinity_branch: bool,
    ) -> tuple[int, ...]:
        square_class = 1 if leading_coefficient in self.point_counting_context.ground_quadratic_residues() else _smallest_nonsquare(self.field.prime)
        encoded = [1 if infinity_branch else 0, square_class]
        for factor in sorted(factors, key=lambda item: (len(item), item)):
            encoded.append(len(factor) - 1)
            encoded.extend(factor)
        return tuple(encoded)

    def _transform_monic_factor(
        self,
        factor: tuple[int, ...],
        matrix_index: int,
    ) -> Optional[tuple[int, tuple[int, ...]]]:
        cache_key = (matrix_index, factor)
        if len(factor) - 1 <= BRANCH_FACTOR_ACTION_MATRIX_CACHE_MAX_DEGREE:
            cached = self.branch_factor_transform_cache.get(cache_key)
            if cache_key in self.branch_factor_transform_cache:
                self.branch_factor_transform_cache.move_to_end(cache_key)
                return cached

        persistent_cached = self._lookup_branch_factor_transform_cache(factor, matrix_index)
        if persistent_cached is not None or self._branch_factor_transform_cache_contains(factor, matrix_index):
            self._cache_branch_factor_transform(cache_key, factor, persistent_cached)
            return persistent_cached

        transformed = self._transform_factor_binary_form(factor, matrix_index)
        transformed = _poly_trim_mod(list(transformed), self.field.prime)
        if transformed == (0,):
            result = None
            self._cache_branch_factor_transform(cache_key, factor, result)
            self._insert_branch_factor_transform_cache(factor, matrix_index, result)
            return result
        leading = transformed[-1]
        monic = _poly_monic_mod(transformed, self.field.prime)
        degree = len(monic) - 1
        if degree == 0:
            result = None
            self._cache_branch_factor_transform(cache_key, factor, result)
            self._insert_branch_factor_transform_cache(factor, matrix_index, result)
            return result
        result = (leading, monic)
        self._cache_branch_factor_transform(cache_key, factor, result)
        self._insert_branch_factor_transform_cache(factor, matrix_index, result)
        return result

    def _branch_factor_transform_cache_contains(self, factor: tuple[int, ...], matrix_index: int) -> bool:
        if self.sqlite_connection is None:
            return False
        if len(factor) - 1 > BRANCH_FACTOR_ACTION_MATRIX_CACHE_MAX_DEGREE:
            return False
        row = self.sqlite_connection.execute(
            """
            SELECT 1
            FROM branch_factor_transform_cache
            WHERE degree = ? AND factor_key = ? AND matrix_index = ?
            """,
            (
                len(factor) - 1,
                self._pack_general_tuple_key(factor),
                matrix_index,
            ),
        ).fetchone()
        return row is not None

    def _lookup_branch_factor_transform_cache(
        self,
        factor: tuple[int, ...],
        matrix_index: int,
    ) -> Optional[tuple[int, tuple[int, ...]]]:
        if self.sqlite_connection is None:
            return None
        if len(factor) - 1 > BRANCH_FACTOR_ACTION_MATRIX_CACHE_MAX_DEGREE:
            return None
        row = self.sqlite_connection.execute(
            """
            SELECT transformed_factor_key
            FROM branch_factor_transform_cache
            WHERE degree = ? AND factor_key = ? AND matrix_index = ?
            """,
            (
                len(factor) - 1,
                self._pack_general_tuple_key(factor),
                matrix_index,
            ),
        ).fetchone()
        if row is None or row[0] is None:
            return None
        unpacked = self._unpack_general_tuple_key(row[0])
        return unpacked[0], tuple(unpacked[1:])

    def _insert_branch_factor_transform_cache(
        self,
        factor: tuple[int, ...],
        matrix_index: int,
        result: Optional[tuple[int, tuple[int, ...]]],
    ) -> None:
        if self.sqlite_connection is None:
            return
        if len(factor) - 1 > BRANCH_FACTOR_ACTION_MATRIX_CACHE_MAX_DEGREE:
            return
        started_at = perf_counter()
        transformed = None if result is None else self._pack_general_tuple_key((result[0], *result[1]))
        self.sqlite_connection.execute(
            """
            INSERT OR IGNORE INTO branch_factor_transform_cache (
                degree,
                factor_key,
                matrix_index,
                transformed_factor_key
            )
            VALUES (?, ?, ?, ?)
            """,
            (
                len(factor) - 1,
                self._pack_general_tuple_key(factor),
                matrix_index,
                transformed,
            ),
        )
        self._mark_sqlite_write()
        self._sqlite_write_seconds += perf_counter() - started_at

    def _transform_factor_binary_form(self, factor: tuple[int, ...], matrix_index: int) -> tuple[int, ...]:
        degree = len(factor) - 1
        action_matrices = self._branch_factor_action_matrices_for_degree(degree)
        if action_matrices is not None:
            return _apply_binary_form_action_matrix(factor, action_matrices[matrix_index], self.field.prime)
        return _transform_binary_form(factor, self.pgl2[matrix_index], self.field.prime)

    def _branch_factor_action_matrices_for_degree(
        self,
        degree: int,
    ) -> Optional[tuple[tuple[tuple[int, ...], ...], ...]]:
        if degree > BRANCH_FACTOR_ACTION_MATRIX_CACHE_MAX_DEGREE:
            return None
        cached = self.branch_factor_action_matrices.get(degree)
        if cached is None:
            cached = tuple(_pgl2_action_matrix(matrix, degree, self.field.prime) for matrix in self.pgl2)
            self.branch_factor_action_matrices[degree] = cached
        return cached

    def _cache_branch_factor_transform(
        self,
        cache_key: tuple[int, tuple[int, ...]],
        factor: tuple[int, ...],
        result: Optional[tuple[int, tuple[int, ...]]],
    ) -> None:
        if len(factor) - 1 > BRANCH_FACTOR_ACTION_MATRIX_CACHE_MAX_DEGREE:
            return
        self.branch_factor_transform_cache[cache_key] = result
        self.branch_factor_transform_cache.move_to_end(cache_key)
        while len(self.branch_factor_transform_cache) > BRANCH_FACTOR_TRANSFORM_CACHE_MAX_ENTRIES:
            self.branch_factor_transform_cache.popitem(last=False)

    def _transformed_branch_key(
        self,
        candidate: BranchDivisorCandidate,
        matrix_index: int,
    ) -> Optional[tuple[int, ...]]:
        leading_coefficient = candidate.leading_coefficient
        infinity_branch = False
        factors: list[tuple[int, ...]] = []
        matrix = self.pgl2[matrix_index]
        a, b, c, d = matrix

        for factor in candidate.factors:
            transformed = self._transform_monic_factor(factor, matrix_index)
            if transformed is None:
                if len(factor) == 2:
                    root = (-factor[0]) % self.field.prime
                    denominator = (c * root + d) % self.field.prime
                    if denominator == 0:
                        infinity_branch = True
                        continue
                return None
            factor_scalar, monic_factor = transformed
            leading_coefficient = leading_coefficient * factor_scalar % self.field.prime
            factors.append(monic_factor)

        if candidate.infinity_branch:
            if c == 0:
                infinity_branch = True
                leading_coefficient = leading_coefficient * d % self.field.prime
            else:
                finite_factor = (d * pow(c, -1, self.field.prime) % self.field.prime, 1)
                factors.append(finite_factor)
                leading_coefficient = leading_coefficient * c % self.field.prime

        finite_degree = sum(len(factor) - 1 for factor in factors)
        if infinity_branch:
            if finite_degree != 2 * self.genus + 1:
                return None
            if leading_coefficient not in self.point_counting_context.ground_quadratic_residues():
                return None
            leading_coefficient = 1
        else:
            if finite_degree != 2 * self.genus + 2:
                return None
            leading_coefficient = 1 if leading_coefficient in self.point_counting_context.ground_quadratic_residues() else _smallest_nonsquare(self.field.prime)

        return self._branch_key(
            leading_coefficient=leading_coefficient,
            factors=factors,
            infinity_branch=infinity_branch,
        )

    def _branch_canonical_key_and_orbit(self, candidate: BranchDivisorCandidate) -> Optional[tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]]:
        orbit_keys = []
        for matrix_index in range(len(self.pgl2)):
            key = self._transformed_branch_key(candidate, matrix_index)
            if key is None:
                return None
            orbit_keys.append(key)
        orbit = tuple(set(orbit_keys))
        return min(orbit), orbit

    def _lookup_branch_orbit_cache(
        self,
        branch_key: tuple[int, ...],
        factorization_pattern: Optional[tuple[int, ...]] = None,
    ) -> Optional[tuple[int, ...]]:
        if factorization_pattern is not None:
            pattern_cache = self.branch_orbit_cache_by_pattern.get(factorization_pattern)
            if pattern_cache is not None:
                cached = pattern_cache.get(branch_key)
                if cached is not None:
                    return cached
        return self.branch_orbit_cache.get(branch_key)

    def _branch_curve_key_for_canonical_branch(
        self,
        canonical_branch_key: tuple[int, ...],
        factorization_pattern: tuple[int, ...],
    ) -> Optional[tuple[int, ...]]:
        pattern_cache = self.branch_canonical_to_curve_key_by_pattern.get(factorization_pattern)
        if pattern_cache is not None:
            cached = pattern_cache.get(canonical_branch_key)
            if cached is not None:
                return cached
        return self.branch_canonical_to_curve_key.get(canonical_branch_key)

    def _normalized_binary_form_key(self, polynomial: Polynomial) -> tuple[int, ...]:
        binary_form = _affine_polynomial_to_binary_form(polynomial, self.binary_degree)
        return _normalize_binary_form_up_to_square_scalar(binary_form, self.field.prime)

    def _canonical_key_and_orbit(self, polynomial: Polynomial) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]:
        binary_form = _affine_polynomial_to_binary_form(polynomial, self.binary_degree)
        orbit = tuple(
            _normalize_binary_form_up_to_square_scalar(
                _apply_binary_form_action_matrix(binary_form, action_matrix, self.field.prime),
                self.field.prime,
            )
            for action_matrix in self.pgl2_action_matrices
        )
        key = min(orbit)
        self.canonical_key_cache[polynomial.coefficients] = key
        return key, orbit

    def _is_enumerated_orbit_key(self, orbit_key: tuple[int, ...]) -> bool:
        if orbit_key[-1] != 0:
            return True
        return orbit_key[-2] != 0 and orbit_key[-2] in self.point_counting_context.ground_quadratic_residues()

    def _enumerated_orbit_keys(self, orbit: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
        return tuple(
            orbit_key
            for orbit_key in set(orbit)
            if self._is_enumerated_orbit_key(orbit_key)
        )

    def canonical_key(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> tuple[int, ...]:
        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        self._require_compatible_polynomial(polynomial)

        cached = self.canonical_key_cache.get(polynomial.coefficients)
        if cached is not None:
            return cached

        key, _ = self._canonical_key_and_orbit(polynomial)
        return key

    def is_seen(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> bool:
        return self.canonical_key(polynomial_or_coefficients) in self.seen_keys

    def mark_seen(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> tuple[int, ...]:
        key = self.canonical_key(polynomial_or_coefficients)
        self.seen_keys.add(key)
        return key

    def is_new_isomorphism_class(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> bool:
        key = self.canonical_key(polynomial_or_coefficients)
        if key in self.seen_keys:
            return False
        self.seen_keys.add(key)
        return True

    def l_polynomial_coefficients_mod_p(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> list[int]:
        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        key = self.canonical_key(polynomial)
        cached = self.l_polynomial_mod_p_cache.get(key)
        if cached is None:
            cached = self.curve(polynomial).l_polynomial_coefficients_mod_p()
            self.l_polynomial_mod_p_cache[key] = cached
        return list(cached)

    def passes_hasse_witt_sparsity_filter(
        self,
        polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...],
        max_sparsity: int,
    ) -> bool:
        if max_sparsity < 0:
            raise ValueError("max sparsity must be nonnegative")
        sparsity_mod_p = sum(1 for coefficient in self.l_polynomial_coefficients_mod_p(polynomial_or_coefficients)[:-1] if coefficient != 0)
        return sparsity_mod_p <= max_sparsity

    def l_polynomial_coefficients(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> list[int]:
        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        key = self.canonical_key(polynomial)
        cache_key = (key, None)
        cached = self.exact_l_polynomial_cache.get(cache_key)
        if cached is None:
            cached = self.curve(polynomial).l_polynomial_coefficients()
            self.exact_l_polynomial_cache[cache_key] = cached
        return list(cached)

    def l_polynomial_coefficients_with_sparsity_limit(
        self,
        polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...],
        max_sparsity: int,
    ) -> Optional[list[int]]:
        if max_sparsity < 0:
            raise ValueError("max sparsity must be nonnegative")

        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        key = self.canonical_key(polynomial)
        cache_key = (key, max_sparsity)
        if cache_key not in self.exact_l_polynomial_cache:
            if not self.passes_hasse_witt_sparsity_filter(polynomial, max_sparsity):
                self.exact_l_polynomial_cache[cache_key] = None
            else:
                self.exact_l_polynomial_cache[cache_key] = self.curve(polynomial).l_polynomial_coefficients_with_sparsity_limit(max_sparsity)

        cached = self.exact_l_polynomial_cache[cache_key]
        return None if cached is None else list(cached)

    def process_polynomial_for_output(
        self,
        polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...],
        max_sparsity: Optional[int],
    ) -> dict[str, object]:
        started_at = perf_counter()
        try:
            result = self._process_polynomial_for_output(polynomial_or_coefficients, max_sparsity)
        finally:
            self._processing_seconds += perf_counter() - started_at

        self._processed_polynomials += 1
        status = str(result["status"])
        self._status_counts[status] = self._status_counts.get(status, 0) + 1
        return result

    def process_branch_divisor_for_output(
        self,
        candidate: BranchDivisorCandidate,
        max_sparsity: Optional[int],
    ) -> dict[str, object]:
        branch_key = self._branch_key(
            leading_coefficient=candidate.leading_coefficient,
            factors=candidate.factors,
            infinity_branch=candidate.infinity_branch,
        )
        branch_result_cache_key = (branch_key, max_sparsity)
        cached_branch_result = self.branch_result_cache.get(branch_result_cache_key)
        if cached_branch_result is not None:
            result = dict(cached_branch_result)
            result["previous_status"] = result["status"]
            result["status"] = "duplicate"
            result["matched_by"] = "exact_branch_key"
            self._processed_polynomials += 1
            self._status_counts["duplicate"] = self._status_counts.get("duplicate", 0) + 1
            return result

        started_at = perf_counter()
        try:
            result = self._process_branch_divisor_for_output(candidate, max_sparsity)
        finally:
            self._processing_seconds += perf_counter() - started_at

        self.branch_result_cache[branch_result_cache_key] = dict(result)
        self._processed_polynomials += 1
        status = str(result["status"])
        self._status_counts[status] = self._status_counts.get(status, 0) + 1
        return result

    def _process_branch_divisor_for_output(
        self,
        candidate: BranchDivisorCandidate,
        max_sparsity: Optional[int],
    ) -> dict[str, object]:
        precomputed_lpoly_mod_p: Optional[list[int]] = None

        if max_sparsity is not None:
            started_at = perf_counter()
            try:
                precomputed_lpoly_mod_p = self._branch_l_polynomial_coefficients_mod_p(candidate)
            finally:
                self._hasse_witt_seconds += perf_counter() - started_at
            sparsity_mod_p = sum(1 for coefficient in precomputed_lpoly_mod_p[:-1] if coefficient != 0)
            if sparsity_mod_p > max_sparsity:
                return {
                    "status": "rejected_hasse_witt_uncanonicalized",
                    "lpoly_mod_p": precomputed_lpoly_mod_p,
                    "branch_key": self._branch_key(
                        leading_coefficient=candidate.leading_coefficient,
                        factors=candidate.factors,
                        infinity_branch=candidate.infinity_branch,
                    ),
                }

        if not self.canonicalize_branch_before_exact:
            return self._process_branch_divisor_sparse_first(
                candidate,
                max_sparsity,
                precomputed_lpoly_mod_p,
            )

        return self._process_branch_divisor_with_factorized_canonicalization(
            candidate,
            max_sparsity,
            precomputed_lpoly_mod_p=precomputed_lpoly_mod_p,
        )

    def _process_branch_divisor_with_factorized_canonicalization(
        self,
        candidate: BranchDivisorCandidate,
        max_sparsity: Optional[int],
        precomputed_lpoly_mod_p: Optional[list[int]],
        precomputed_exact_lpoly: Optional[list[int]] = None,
    ) -> dict[str, object]:
        coefficients: Optional[tuple[int, ...]] = None
        precomputed_ground_invariants: Optional[tuple[int, int]] = None
        started_at = perf_counter()
        try:
            canonical_branch_data = self._branch_canonical_key_and_orbit(candidate)
        finally:
            self._factorized_pgl2_seconds += perf_counter() - started_at
        if canonical_branch_data is not None:
            candidate_key = self._branch_key(
                leading_coefficient=candidate.leading_coefficient,
                factors=candidate.factors,
                infinity_branch=candidate.infinity_branch,
            )
            canonical_branch_key = self._lookup_branch_orbit_cache(candidate_key, candidate.factorization_pattern)
            if canonical_branch_key is None:
                canonical_branch_key, branch_orbit = canonical_branch_data
                for orbit_key in branch_orbit:
                    self.branch_orbit_cache[orbit_key] = canonical_branch_key
                    self.branch_orbit_cache_by_pattern.setdefault(candidate.factorization_pattern, {})[orbit_key] = canonical_branch_key
                self._insert_branch_orbit_cache(branch_orbit, canonical_branch_key)

            cached_curve_key = self._branch_curve_key_for_canonical_branch(canonical_branch_key, candidate.factorization_pattern)
            if cached_curve_key is not None:
                record = self.canonical_records.get(cached_curve_key)
                if record is not None:
                    cached_result = self._result_from_record(record, max_sparsity)
                    if cached_result is not None:
                        cached_result = dict(cached_result)
                        cached_result["previous_status"] = cached_result["status"]
                        cached_result["status"] = "duplicate"
                        cached_result["matched_by"] = "factorized_branch_orbit"
                        return cached_result

            if coefficients is None:
                started_at = perf_counter()
                coefficients = polynomial_from_branch_factors(
                    candidate.factors,
                    candidate.leading_coefficient,
                    self.field.prime,
                )
                self._expansion_seconds += perf_counter() - started_at
            result = self._process_polynomial_for_output(
                coefficients,
                max_sparsity,
                precomputed_lpoly_mod_p=precomputed_lpoly_mod_p,
                precomputed_ground_invariants=precomputed_ground_invariants or self.branch_ground_invariants(candidate),
                precomputed_exact_lpoly=precomputed_exact_lpoly,
                assume_squarefree=True,
                branch_candidate=candidate,
            )
            if "canonical_key" in result:
                self._insert_branch_curve_cache(
                    canonical_branch_key,
                    tuple(result["canonical_key"]),  # type: ignore[arg-type]
                    candidate.factorization_pattern,
                )
            return result

        if coefficients is None:
            started_at = perf_counter()
            coefficients = polynomial_from_branch_factors(
                candidate.factors,
                candidate.leading_coefficient,
                self.field.prime,
            )
            self._expansion_seconds += perf_counter() - started_at
        result = self._process_polynomial_for_output(
            coefficients,
            max_sparsity,
            precomputed_lpoly_mod_p=precomputed_lpoly_mod_p,
            precomputed_ground_invariants=precomputed_ground_invariants or self.branch_ground_invariants(candidate),
            precomputed_exact_lpoly=precomputed_exact_lpoly,
            assume_squarefree=True,
            branch_candidate=candidate,
        )
        result = dict(result)
        result["factorized_matching"] = "fallback_expanded"
        return result

    def _process_branch_divisor_sparse_first(
        self,
        candidate: BranchDivisorCandidate,
        max_sparsity: Optional[int],
        precomputed_lpoly_mod_p: Optional[list[int]],
    ) -> dict[str, object]:
        started_at = perf_counter()
        try:
            exact_lpoly = self._branch_l_polynomial_coefficients(
                candidate,
                max_sparsity=max_sparsity,
            )
        finally:
            self._exact_lpoly_seconds += perf_counter() - started_at

        if exact_lpoly is None:
            return {
                "status": "rejected_exact_uncanonicalized",
                "lpoly_mod_p": precomputed_lpoly_mod_p,
                "branch_key": self._branch_key(
                    leading_coefficient=candidate.leading_coefficient,
                    factors=candidate.factors,
                    infinity_branch=candidate.infinity_branch,
                ),
                "sparsity": -1,
            }

        return self._process_branch_divisor_with_factorized_canonicalization(
            candidate,
            max_sparsity,
            precomputed_lpoly_mod_p=precomputed_lpoly_mod_p,
            precomputed_exact_lpoly=exact_lpoly,
        )

    def _process_polynomial_for_output(
        self,
        polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...],
        max_sparsity: Optional[int],
        precomputed_lpoly_mod_p: Optional[list[int]] = None,
        precomputed_ground_invariants: Optional[tuple[int, int]] = None,
        precomputed_exact_lpoly: Optional[list[int]] = None,
        assume_squarefree: bool = False,
        branch_candidate: Optional[BranchDivisorCandidate] = None,
    ) -> dict[str, object]:
        if max_sparsity is not None and max_sparsity < 0:
            raise ValueError("max sparsity must be nonnegative")

        polynomial = self._coerce_polynomial(polynomial_or_coefficients)
        self._require_compatible_polynomial(polynomial)
        if not assume_squarefree and not polynomial.is_squarefree():
            return {"status": "singular", "coefficients": polynomial.coefficients}

        rational_branch_count, ground_point_count = (
            precomputed_ground_invariants
            if precomputed_ground_invariants is not None
            else self.ground_invariants(polynomial)
        )
        if (
            self.skip_even_models_with_rational_branch_point
            and polynomial.degree == self.binary_degree
            and rational_branch_count > 0
        ):
            return {
                "status": "covered_by_odd_model",
                "coefficients": polynomial.coefficients,
                "rational_branch_count": rational_branch_count,
                "ground_point_count": ground_point_count,
            }

        if max_sparsity is not None:
            if precomputed_lpoly_mod_p is None:
                precomputed_lpoly_mod_p = self.curve(polynomial).l_polynomial_coefficients_mod_p()
            sparsity_mod_p = sum(1 for coefficient in precomputed_lpoly_mod_p[:-1] if coefficient != 0)
            if sparsity_mod_p > max_sparsity:
                return {
                    "status": "rejected_hasse_witt_uncanonicalized",
                    "coefficients": polynomial.coefficients,
                    "lpoly_mod_p": precomputed_lpoly_mod_p,
                }

        hasse_witt_lpoly_mod_p = tuple(precomputed_lpoly_mod_p) if precomputed_lpoly_mod_p is not None else None
        orbit_key = self._normalized_binary_form_key(polynomial)
        orbit_size: Optional[int] = None
        canonical_key = self._lookup_orbit_cache(
            rational_branch_count,
            ground_point_count,
            hasse_witt_lpoly_mod_p,
            orbit_key,
        )
        if canonical_key is None:
            canonical_key, orbit = self._canonical_key_and_orbit(polynomial)
            orbit_size = len(self._enumerated_orbit_keys(orbit))
            self._insert_orbit_cache(
                rational_branch_count,
                ground_point_count,
                hasse_witt_lpoly_mod_p,
                orbit,
                canonical_key,
            )
        else:
            self.canonical_key_cache[polynomial.coefficients] = canonical_key
        existing_record = self.canonical_records.get(canonical_key)
        record = self._register_or_update_record(
            canonical_key=canonical_key,
            coefficients=polynomial.coefficients,
            rational_branch_count=rational_branch_count,
            ground_point_count=ground_point_count,
            hasse_witt_lpoly_mod_p=hasse_witt_lpoly_mod_p,
            orbit_size=orbit_size,
        )

        cached_result = self._result_from_record(record, max_sparsity)
        if cached_result is not None:
            if existing_record is not None:
                cached_result = dict(cached_result)
                cached_result["previous_status"] = cached_result["status"]
                cached_result["status"] = "duplicate"
            return cached_result

        lpoly_mod_p = (
            precomputed_lpoly_mod_p
            if precomputed_lpoly_mod_p is not None
            else list(record.hasse_witt_lpoly_mod_p)
            if record.hasse_witt_lpoly_mod_p is not None
            else self.l_polynomial_coefficients_mod_p(polynomial)
        )
        self.l_polynomial_mod_p_cache[canonical_key] = lpoly_mod_p
        record.hasse_witt_lpoly_mod_p = tuple(lpoly_mod_p)
        self._index_record(record)
        sparsity_mod_p = sum(1 for coefficient in lpoly_mod_p[:-1] if coefficient != 0)
        if max_sparsity is not None and sparsity_mod_p > max_sparsity:
            record.status_by_max_sparsity[max_sparsity] = "rejected_hasse_witt"
            record.sparsity_by_max_sparsity[max_sparsity] = -1
            record.exact_lpoly_by_max_sparsity[max_sparsity] = None
            self._insert_curve_cache(
                canonical_key=canonical_key,
                rational_branch_count=rational_branch_count,
                coefficients=polynomial.coefficients,
                lpoly_mod_p=lpoly_mod_p,
                exact_lpoly=None,
                sparsity=-1,
                status="rejected_hasse_witt",
            )
            return {
                "status": "rejected_hasse_witt",
                "canonical_key": canonical_key,
                "coefficients": polynomial.coefficients,
                "lpoly_mod_p": lpoly_mod_p,
                "sparsity": -1,
            }

        exact_cache_key = (canonical_key, max_sparsity)
        if exact_cache_key not in self.exact_l_polynomial_cache:
            if precomputed_exact_lpoly is not None:
                self.exact_l_polynomial_cache[exact_cache_key] = list(precomputed_exact_lpoly)
            else:
                started_at = perf_counter()
                try:
                    if branch_candidate is not None:
                        self.exact_l_polynomial_cache[exact_cache_key] = self._branch_l_polynomial_coefficients(
                            branch_candidate,
                            max_sparsity=max_sparsity,
                        )
                    else:
                        self.exact_l_polynomial_cache[exact_cache_key] = self.curve(polynomial)._compute_l_polynomial_coefficients(
                            max_sparsity=max_sparsity,
                        )
                finally:
                    self._exact_lpoly_seconds += perf_counter() - started_at
        exact_lpoly = self.exact_l_polynomial_cache[exact_cache_key]
        if exact_lpoly is None:
            record.status_by_max_sparsity[max_sparsity] = "rejected_exact"
            record.sparsity_by_max_sparsity[max_sparsity] = -1
            record.exact_lpoly_by_max_sparsity[max_sparsity] = None
            self._insert_curve_cache(
                canonical_key=canonical_key,
                rational_branch_count=rational_branch_count,
                coefficients=polynomial.coefficients,
                lpoly_mod_p=lpoly_mod_p,
                exact_lpoly=None,
                sparsity=-1,
                status="rejected_exact",
            )
            return {
                "status": "rejected_exact",
                "canonical_key": canonical_key,
                "coefficients": polynomial.coefficients,
                "lpoly_mod_p": lpoly_mod_p,
                "sparsity": -1,
            }

        sparsity = sum(1 for coefficient in exact_lpoly[:-1] if coefficient != 0)
        record.status_by_max_sparsity[max_sparsity] = "sparse"
        record.sparsity_by_max_sparsity[max_sparsity] = sparsity
        record.exact_lpoly_by_max_sparsity[max_sparsity] = tuple(exact_lpoly)
        self._insert_curve_cache(
            canonical_key=canonical_key,
            rational_branch_count=rational_branch_count,
            coefficients=polynomial.coefficients,
            lpoly_mod_p=lpoly_mod_p,
            exact_lpoly=exact_lpoly,
            sparsity=sparsity,
            status="sparse",
        )
        self._insert_sparse_curve(
            canonical_key=canonical_key,
            coefficients=polynomial.coefficients,
            lpoly=exact_lpoly,
            sparsity=sparsity,
            rational_branch_count=rational_branch_count,
            branch_candidate=branch_candidate,
        )
        return {
            "status": "sparse",
            "canonical_key": canonical_key,
            "coefficients": polynomial.coefficients,
            "lpoly_mod_p": lpoly_mod_p,
            "lpoly": exact_lpoly,
            "sparsity": sparsity,
        }

    def process_polynomials_for_output(
        self,
        polynomials_or_coefficients: Iterable[Polynomial | list[int] | tuple[int, ...]],
        max_sparsity: Optional[int],
    ) -> dict[str, int]:
        stats: dict[str, int] = {"processed": 0}
        for polynomial_or_coefficients in polynomials_or_coefficients:
            result = self.process_polynomial_for_output(polynomial_or_coefficients, max_sparsity)
            status = str(result["status"])
            stats["processed"] += 1
            stats[status] = stats.get(status, 0) + 1
        return stats

    def _coerce_polynomial(self, polynomial_or_coefficients: Polynomial | list[int] | tuple[int, ...]) -> Polynomial:
        if isinstance(polynomial_or_coefficients, Polynomial):
            return polynomial_or_coefficients
        return self.polynomial(polynomial_or_coefficients)

    def _require_compatible_polynomial(self, polynomial: Polynomial) -> None:
        if polynomial.field.prime != self.field.prime:
            raise ValueError("polynomial is over a different field")
        if polynomial.degree not in {2 * self.genus + 1, 2 * self.genus + 2}:
            raise ValueError("polynomial degree does not match enumeration genus")

    def branch_candidate_from_coefficients(self, coefficients: tuple[int, ...]) -> BranchDivisorCandidate:
        polynomial = self.polynomial(coefficients)
        self._require_compatible_polynomial(polynomial)
        if not polynomial.is_squarefree():
            raise ValueError("polynomial is not squarefree")

        prime = self.field.prime
        leading_coefficient = polynomial.coefficient(polynomial.degree)
        remaining = _poly_monic_mod(polynomial.coefficients, prime)
        factors: list[tuple[int, ...]] = []
        for degree in range(1, polynomial.degree + 1):
            if remaining == (1,):
                break
            while len(remaining) - 1 >= degree:
                matched_factor: Optional[tuple[int, ...]] = None
                for factor in self.iter_irreducible_polynomials(degree):
                    if len(factor) > len(remaining):
                        break
                    if _poly_remainder_mod(remaining, factor, prime) == (0,):
                        matched_factor = factor
                        break
                if matched_factor is None:
                    break
                factors.append(matched_factor)
                remaining = _poly_exact_quotient_mod(remaining, matched_factor, prime)

        if remaining != (1,):
            if _is_irreducible_by_rabin_test(remaining, prime):
                factors.append(remaining)
            else:
                raise RuntimeError("failed to factor polynomial into irreducible branch factors")

        factors.sort(key=lambda factor: (len(factor), factor))
        pattern = [0] * (polynomial.degree + 1)
        for factor in factors:
            pattern[len(factor) - 1] += 1
        return BranchDivisorCandidate(
            leading_coefficient=leading_coefficient,
            factors=tuple(factors),
            infinity_branch=polynomial.degree == 2 * self.genus + 1,
            factorization_pattern=tuple(pattern),
        )

    def _sqlite_irreducible_cache(self) -> IrreduciblePolynomialCache:
        if self.irreducible_cache is None:
            self.irreducible_cache = IrreduciblePolynomialCache(self.irreducible_cache_path)
        return self.irreducible_cache

    def can_materialize_irreducibles(self, degree: int) -> bool:
        if degree in self._irreducible_factor_cache:
            return True
        estimated_bytes = estimated_irreducible_tuple_memory_bytes(self.field.prime, degree)
        return self._irreducible_memory_bytes + estimated_bytes <= self.irreducible_memory_budget_bytes

    def should_materialize_branch_irreducibles(self, degree: int) -> bool:
        if degree in self._irreducible_factor_cache:
            return True
        if self._sqlite_irreducible_cache().has_complete(self.field.prime, degree):
            return True
        return False

    def _irreducible_polynomials(self, degree: int) -> tuple[tuple[int, ...], ...]:
        cached = self._irreducible_factor_cache.get(degree)
        if cached is not None:
            return cached

        estimated_bytes = estimated_irreducible_tuple_memory_bytes(self.field.prime, degree)
        if self._irreducible_memory_bytes + estimated_bytes > self.irreducible_memory_budget_bytes:
            raise MemoryError(
                "loading irreducible polynomials would exceed the configured memory budget "
                f"({estimated_bytes / BYTES_PER_MEGABYTE:.1f} MB needed for p={self.field.prime}, "
                f"degree={degree}; "
                f"{(self.irreducible_memory_budget_bytes - self._irreducible_memory_bytes) / BYTES_PER_MEGABYTE:.1f} MB remaining)"
            )

        sqlite_cached = self._sqlite_irreducible_cache().get_complete(self.field.prime, degree)
        if sqlite_cached is not None:
            self._irreducible_factor_cache[degree] = sqlite_cached
            self._irreducible_memory_bytes += estimated_bytes
            return sqlite_cached

        cached = _generate_monic_irreducible_polynomials(self.field, degree)
        self._irreducible_factor_cache[degree] = cached
        self._irreducible_memory_bytes += estimated_bytes
        self._sqlite_irreducible_cache().store_complete(self.field.prime, degree, cached)
        return cached

    def iter_irreducible_polynomials(self, degree: int) -> Iterable[tuple[int, ...]]:
        cached = self._irreducible_factor_cache.get(degree)
        if cached is not None:
            yield from cached
            return

        sqlite_cache = self._sqlite_irreducible_cache()
        if sqlite_cache.has_complete(self.field.prime, degree):
            yield from sqlite_cache.iter_complete(self.field.prime, degree)
            return

        yield from _stream_monic_irreducible_polynomials(self.field, degree)

    def random_sage_irreducible_polynomials(
        self,
        degree: int,
        count: int,
        rng: random.Random,
    ) -> tuple[tuple[int, ...], ...]:
        if count < 0:
            raise ValueError("irreducible count must be nonnegative")
        if count == 0:
            return ()
        seed = rng.randrange(2**63)
        in_process = self._random_sage_irreducible_polynomials_in_process(degree, count, seed)
        if in_process is not None:
            return in_process
        return self._random_sage_irreducible_polynomials_from_process(degree, count, seed)

    def _random_sage_irreducible_polynomials_in_process(
        self,
        degree: int,
        count: int,
        seed: int,
    ) -> Optional[tuple[tuple[int, ...], ...]]:
        self._ensure_sage_cache_environment()
        try:
            from sage.all import GF, PolynomialRing, set_random_seed  # type: ignore[import-not-found]
        except ModuleNotFoundError:
            return None

        set_random_seed(seed)
        ring = PolynomialRing(GF(self.field.prime), "x")
        polynomials = []
        seen = set()
        while len(polynomials) < count:
            polynomial = ring.irreducible_element(degree, algorithm="random")
            leading = int(polynomial[degree]) % self.field.prime
            leading_inverse = pow(leading, -1, self.field.prime)
            coefficients = tuple(
                (int(polynomial[index]) * leading_inverse) % self.field.prime
                for index in range(degree + 1)
            )
            if coefficients not in seen:
                seen.add(coefficients)
                polynomials.append(coefficients)
        return tuple(polynomials)

    def _random_sage_irreducible_polynomials_from_process(
        self,
        degree: int,
        count: int,
        seed: int,
    ) -> tuple[tuple[int, ...], ...]:
        process = self._sage_irreducible_server()
        if process.stdin is None or process.stdout is None:
            raise RuntimeError("Sage irreducible generator process was not opened with pipes")
        request = {
            "prime": self.field.prime,
            "degree": degree,
            "count": count,
            "seed": seed,
        }
        process.stdin.write(json.dumps(request, separators=(",", ":")) + "\n")
        process.stdin.flush()
        line = process.stdout.readline()
        if not line:
            stderr = process.stderr.read() if process.stderr is not None else ""
            raise RuntimeError(f"Sage irreducible generator stopped unexpectedly: {stderr}")
        response = json.loads(line)
        if not response.get("ok"):
            raise RuntimeError(f"Sage irreducible generator failed: {response.get('error')}")
        return tuple(tuple(int(coefficient) for coefficient in polynomial) for polynomial in response["polynomials"])

    def _sage_irreducible_server(self) -> subprocess.Popen[str]:
        if self._sage_irreducible_process is not None and self._sage_irreducible_process.poll() is None:
            return self._sage_irreducible_process
        sage_executable = shutil.which("sage")
        if sage_executable is None:
            raise RuntimeError(
                "Sage is required for random generation of irreducibles missing from SQLite. "
                "Install Sage or run with degrees precomputed in the irreducible SQLite cache."
            )
        sage_environment = os.environ.copy()
        sage_environment["DOT_SAGE"] = str(self._ensure_sage_cache_environment())
        self._sage_irreducible_process = subprocess.Popen(
            [sage_executable, "-python", "-u", "-c", SAGE_RANDOM_IRREDUCIBLE_SERVER_CODE],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=sage_environment,
        )
        return self._sage_irreducible_process

    def _ensure_sage_cache_environment(self) -> Path:
        sage_cache = Path(os.environ.get("DOT_SAGE", "/private/tmp/hyperelliptic-search-sage"))
        sage_cache.mkdir(parents=True, exist_ok=True)
        os.environ["DOT_SAGE"] = str(sage_cache)
        return sage_cache

    def _lookup_orbit_cache(
        self,
        rational_branch_count: int,
        ground_point_count: int,
        hasse_witt_lpoly_mod_p: Optional[tuple[int, ...]],
        orbit_key: tuple[int, ...],
    ) -> Optional[tuple[int, ...]]:
        if self.sqlite_connection is None:
            return None

        if self._orbit_cache_has_hasse_witt_lpoly_mod_p():
            row = self.sqlite_connection.execute(
                """
                SELECT canonical_key
                FROM orbit_cache
                WHERE rational_branch_count = ?
                  AND ground_point_count = ?
                  AND hasse_witt_lpoly_mod_p = ?
                  AND orbit_key = ?
                """,
                (
                    rational_branch_count,
                    ground_point_count,
                    _pack_general_int_tuple_blob(hasse_witt_lpoly_mod_p) if hasse_witt_lpoly_mod_p is not None else b"",
                    self._pack_field_tuple_blob(orbit_key),
                ),
            ).fetchone()
        elif self._orbit_cache_has_ground_point_count():
            row = self.sqlite_connection.execute(
                """
                SELECT canonical_key
                FROM orbit_cache
                WHERE rational_branch_count = ?
                  AND ground_point_count = ?
                  AND orbit_key = ?
                """,
                (rational_branch_count, ground_point_count, self._pack_field_tuple_blob(orbit_key)),
            ).fetchone()
        else:
            row = self.sqlite_connection.execute(
                """
                SELECT canonical_key
                FROM orbit_cache
                WHERE rational_branch_count = ?
                  AND orbit_key = ?
                """,
                (rational_branch_count, self._pack_field_tuple_blob(orbit_key)),
            ).fetchone()
        return None if row is None else self._unpack_field_tuple_blob(row[0])

    def _insert_orbit_cache(
        self,
        rational_branch_count: int,
        ground_point_count: int,
        hasse_witt_lpoly_mod_p: Optional[tuple[int, ...]],
        orbit: tuple[tuple[int, ...], ...],
        canonical_key: tuple[int, ...],
    ) -> None:
        if self.sqlite_connection is None:
            return

        started_at = perf_counter()
        packed_canonical_key = self._pack_field_tuple_blob(canonical_key)
        enumerated_orbit = self._enumerated_orbit_keys(orbit)
        if not self.cache_full_orbits:
            enumerated_orbit = (canonical_key,)
        if self._orbit_cache_has_hasse_witt_lpoly_mod_p():
            packed_hasse_witt_lpoly = (
                _pack_general_int_tuple_blob(hasse_witt_lpoly_mod_p)
                if hasse_witt_lpoly_mod_p is not None
                else b""
            )
            rows = [
                (
                    rational_branch_count,
                    ground_point_count,
                    packed_hasse_witt_lpoly,
                    self._pack_field_tuple_blob(orbit_key),
                    packed_canonical_key,
                )
                for orbit_key in enumerated_orbit
            ]
            self.sqlite_connection.executemany(
                """
                INSERT OR IGNORE INTO orbit_cache (
                    rational_branch_count,
                    ground_point_count,
                    hasse_witt_lpoly_mod_p,
                    orbit_key,
                    canonical_key
                )
                VALUES (?, ?, ?, ?, ?)
                """,
                rows,
            )
        elif self._orbit_cache_has_ground_point_count():
            rows = [
                (rational_branch_count, ground_point_count, self._pack_field_tuple_blob(orbit_key), packed_canonical_key)
                for orbit_key in enumerated_orbit
            ]
            self.sqlite_connection.executemany(
                """
                INSERT OR IGNORE INTO orbit_cache (
                    rational_branch_count,
                    ground_point_count,
                    orbit_key,
                    canonical_key
                )
                VALUES (?, ?, ?, ?)
                """,
                rows,
            )
        else:
            rows = [
                (rational_branch_count, self._pack_field_tuple_blob(orbit_key), packed_canonical_key)
                for orbit_key in enumerated_orbit
            ]
            self.sqlite_connection.executemany(
                """
                INSERT OR IGNORE INTO orbit_cache (
                    rational_branch_count,
                    orbit_key,
                    canonical_key
                )
                VALUES (?, ?, ?)
                """,
                rows,
            )
        self._mark_sqlite_write()
        self._sqlite_write_seconds += perf_counter() - started_at

    def _insert_branch_orbit_cache(
        self,
        orbit_keys: Iterable[tuple[int, ...]],
        canonical_branch_key: tuple[int, ...],
    ) -> None:
        if self.sqlite_connection is None:
            return
        started_at = perf_counter()
        packed_canonical = self._pack_general_tuple_key(canonical_branch_key)
        stored_orbit_keys = tuple(orbit_keys) if self.cache_full_orbits else (canonical_branch_key,)
        self.sqlite_connection.executemany(
            """
            INSERT OR IGNORE INTO branch_orbit_cache (orbit_branch_key, canonical_branch_key)
            VALUES (?, ?)
            """,
            ((self._pack_general_tuple_key(orbit_key), packed_canonical) for orbit_key in stored_orbit_keys),
        )
        self._mark_sqlite_write()
        self._sqlite_write_seconds += perf_counter() - started_at

    def _insert_branch_curve_cache(
        self,
        canonical_branch_key: tuple[int, ...],
        canonical_key: tuple[int, ...],
        factorization_pattern: Optional[tuple[int, ...]] = None,
    ) -> None:
        self.branch_canonical_to_curve_key[canonical_branch_key] = canonical_key
        if factorization_pattern is not None:
            self.branch_canonical_to_curve_key_by_pattern.setdefault(factorization_pattern, {})[canonical_branch_key] = canonical_key
        if self.sqlite_connection is None:
            return
        started_at = perf_counter()
        self.sqlite_connection.execute(
            """
            INSERT OR REPLACE INTO branch_curve_cache (canonical_branch_key, canonical_key)
            VALUES (?, ?)
            """,
            (
                self._pack_general_tuple_key(canonical_branch_key),
                self._pack_field_tuple_blob(canonical_key),
            ),
        )
        self._mark_sqlite_write()
        self._sqlite_write_seconds += perf_counter() - started_at

    def _insert_curve_cache(
        self,
        canonical_key: tuple[int, ...],
        rational_branch_count: int,
        coefficients: tuple[int, ...],
        lpoly_mod_p: Optional[list[int]],
        exact_lpoly: Optional[list[int]],
        sparsity: Optional[int],
        status: str,
    ) -> None:
        if self.sqlite_connection is None:
            return

        started_at = perf_counter()
        self.sqlite_connection.execute(
            """
            INSERT OR REPLACE INTO curve_cache (
                canonical_key,
                rational_branch_count,
                coefficients,
                lpoly_mod_p,
                exact_lpoly,
                sparsity,
                status
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                self._pack_field_tuple_blob(canonical_key),
                rational_branch_count,
                _pack_int_tuple(coefficients),
                _pack_general_int_tuple_blob(lpoly_mod_p) if lpoly_mod_p is not None else None,
                _pack_general_int_tuple_blob(exact_lpoly) if exact_lpoly is not None else None,
                sparsity,
                status,
            ),
        )
        self._mark_sqlite_write()
        self._sqlite_write_seconds += perf_counter() - started_at

    def _insert_sparse_curve(
        self,
        canonical_key: tuple[int, ...],
        coefficients: tuple[int, ...],
        lpoly: list[int],
        sparsity: int,
        rational_branch_count: int,
        branch_candidate: Optional[BranchDivisorCandidate] = None,
    ) -> None:
        if self.sqlite_connection is None:
            return

        if branch_candidate is None:
            branch_candidate = self.branch_candidate_from_coefficients(coefficients)

        started_at = perf_counter()
        self.sqlite_connection.execute(
            """
            INSERT OR REPLACE INTO sparse_curves (
                canonical_key,
                coefficients,
                branch_factors,
                branch_infinity_branch,
                branch_leading_coefficient,
                branch_factorization_pattern,
                lpoly,
                sparsity,
                rational_branch_count
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                self._pack_field_tuple_blob(canonical_key),
                _pack_int_tuple(coefficients),
                _pack_branch_factors_text(branch_candidate.factors),
                int(branch_candidate.infinity_branch),
                branch_candidate.leading_coefficient,
                _pack_branch_pattern_text(branch_candidate.factorization_pattern),
                _pack_int_tuple(lpoly),
                sparsity,
                rational_branch_count,
            ),
        )
        self._mark_sqlite_write()
        self._sqlite_write_seconds += perf_counter() - started_at


def coefficient_vectors(
    prime: int,
    genus: int,
    degree_model: str,
    limit: Optional[int],
) -> Iterable[tuple[int, ...]]:
    degrees = []
    if degree_model in {"odd", "both"}:
        degrees.append(2 * genus + 1)
    if degree_model in {"even", "both"}:
        degrees.append(2 * genus + 2)

    produced = 0
    for degree in degrees:
        leading_coefficients = _leading_coefficients_for_degree(prime, genus, degree)
        for coefficients in product(range(prime), repeat=degree):
            for leading_coefficient in leading_coefficients:
                yield (*coefficients, leading_coefficient)
                produced += 1
                if limit is not None and produced >= limit:
                    return


def coefficient_vectors_by_support(
    prime: int,
    genus: int,
    degree_model: str,
    limit: Optional[int],
) -> Iterable[tuple[int, ...]]:
    degrees = []
    if degree_model in {"odd", "both"}:
        degrees.append(2 * genus + 1)
    if degree_model in {"even", "both"}:
        degrees.append(2 * genus + 2)

    produced = 0
    max_degree = max(degrees, default=0)
    nonzero_values = tuple(range(1, prime))
    for support_size in range(1, max_degree + 2):
        lower_support_size = support_size - 1
        for degree in degrees:
            if lower_support_size > degree:
                continue
            leading_coefficients = _leading_coefficients_for_degree(prime, genus, degree)
            for support_positions in combinations(range(degree), lower_support_size):
                for support_values in product(nonzero_values, repeat=lower_support_size):
                    coefficients = [0] * degree
                    for position, value in zip(support_positions, support_values):
                        coefficients[position] = value
                    for leading_coefficient in leading_coefficients:
                        yield (*coefficients, leading_coefficient)
                        produced += 1
                        if limit is not None and produced >= limit:
                            return


def _leading_coefficients_for_degree(prime: int, genus: int, degree: int) -> tuple[int, ...]:
    return (1,) if degree == 2 * genus + 1 else (1, _smallest_nonsquare(prime))


def coefficient_vector_at_index(prime: int, genus: int, degree_model: str, index: int) -> tuple[int, ...]:
    if index < 0:
        raise ValueError("coefficient vector index must be nonnegative")

    degrees = []
    if degree_model in {"odd", "both"}:
        degrees.append(2 * genus + 1)
    if degree_model in {"even", "both"}:
        degrees.append(2 * genus + 2)

    remaining = index
    for degree in degrees:
        leading_coefficients = _leading_coefficients_for_degree(prime, genus, degree)
        block_size = prime ** degree * len(leading_coefficients)
        if remaining >= block_size:
            remaining -= block_size
            continue

        low_index, leading_index = divmod(remaining, len(leading_coefficients))
        coefficients = []
        for exponent in range(degree - 1, -1, -1):
            place_value = prime ** exponent
            digit, low_index = divmod(low_index, place_value)
            coefficients.append(digit)
        return (*coefficients, leading_coefficients[leading_index])

    raise IndexError("coefficient vector index is outside the enumeration range")


def default_sqlite_path(prime: int, genus: int, max_sparsity: Optional[int]) -> Path:
    sparsity_label = f"s_{max_sparsity}" if max_sparsity is not None else "all"
    return Path("data_gen") / "results" / f"p{prime}_g{genus}_{sparsity_label}.sqlite"


def total_coefficient_vectors(prime: int, genus: int, degree_model: str, limit: Optional[int]) -> int:
    total = 0
    if degree_model in {"odd", "both"}:
        total += prime ** (2 * genus + 1)
    if degree_model in {"even", "both"}:
        total += 2 * prime ** (2 * genus + 2)
    return min(total, limit) if limit is not None else total


def factorization_patterns(total_degree: int, *, skip_linear: bool = False) -> Iterable[tuple[int, ...]]:
    counts = [0] * (total_degree + 1)
    minimum_degree = 2 if skip_linear else 1

    def rec(remaining: int, smallest_degree: int) -> Iterable[tuple[int, ...]]:
        if remaining == 0:
            yield tuple(counts)
            return
        for degree in range(max(smallest_degree, minimum_degree), remaining + 1):
            counts[degree] += 1
            yield from rec(remaining - degree, degree)
            counts[degree] -= 1

    yield from rec(total_degree, minimum_degree)


def pattern_to_degree_tuple(pattern: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(degree for degree, count in enumerate(pattern) for _ in range(count))


def branch_factorization_patterns(prime: int, genus: int, degree_model: str) -> Iterable[tuple[str, int, tuple[int, ...], tuple[int, ...]]]:
    if degree_model in {"odd", "both"}:
        degree = 2 * genus + 1
        for pattern in factorization_patterns(degree):
            yield "odd", degree, (1,), pattern
    if degree_model in {"even", "both"}:
        degree = 2 * genus + 2
        leading_coefficients = (1, _smallest_nonsquare(prime))
        for pattern in factorization_patterns(degree, skip_linear=True):
            yield "even", degree, leading_coefficients, pattern


def branch_pattern_presentation_count(prime: int, pattern: tuple[int, ...], leading_coefficient_count: int) -> int:
    total = leading_coefficient_count
    for degree, multiplicity in enumerate(pattern):
        if multiplicity == 0:
            continue
        total *= comb(monic_irreducible_count(prime, degree), multiplicity)
    return total


def total_branch_divisor_presentations(prime: int, genus: int, degree_model: str, limit: Optional[int]) -> int:
    total = 0
    for _, _, leading_coefficients, pattern in branch_factorization_patterns(prime, genus, degree_model):
        total += branch_pattern_presentation_count(prime, pattern, len(leading_coefficients))
        if limit is not None and total >= limit:
            return limit
    return total


def _scale_polynomial_mod(polynomial: tuple[int, ...], scalar: int, prime: int) -> tuple[int, ...]:
    return tuple((scalar * coefficient) % prime for coefficient in polynomial)


def polynomial_from_branch_factors(factors: Iterable[tuple[int, ...]], leading_coefficient: int, prime: int) -> tuple[int, ...]:
    product_polynomial = _poly_product_tree_mod(factors, prime)
    if leading_coefficient != 1:
        product_polynomial = _scale_polynomial_mod(product_polynomial, leading_coefficient, prime)
    return product_polynomial


def branch_factor_choices_for_pattern(
    context: EnumerationContext,
    pattern: tuple[int, ...],
) -> Iterable[tuple[tuple[int, ...], ...]]:
    active_degrees = [(degree, multiplicity) for degree, multiplicity in enumerate(pattern) if multiplicity > 0]

    def choices_for_degree(degree: int, multiplicity: int) -> Iterable[tuple[tuple[int, ...], ...]]:
        if multiplicity == 1:
            if degree not in context._irreducible_factor_cache and context.should_materialize_branch_irreducibles(degree):
                context._irreducible_polynomials(degree)
            for polynomial in context.iter_irreducible_polynomials(degree):
                yield (polynomial,)
            return

        if degree not in context._irreducible_factor_cache:
            if not context.should_materialize_branch_irreducibles(degree):
                raise MemoryError(
                    "branch-divisor patterns using multiple factors of an irreducible degree missing from SQLite "
                    f"are not streamable yet (degree={degree}, multiplicity={multiplicity})"
                )
            context._irreducible_polynomials(degree)
        yield from combinations(context._irreducible_factor_cache[degree], multiplicity)

    def rec(
        index: int,
        selected: tuple[tuple[int, ...], ...],
    ) -> Iterable[tuple[tuple[int, ...], ...]]:
        if index == len(active_degrees):
            yield selected
            return
        degree, multiplicity = active_degrees[index]
        for degree_choice in choices_for_degree(degree, multiplicity):
            yield from rec(index + 1, selected + degree_choice)

    yield from rec(0, ())


def branch_divisor_polynomials(
    context: EnumerationContext,
    degree_model: str,
    limit: Optional[int],
) -> Iterable[BranchDivisorCandidate]:
    produced = 0
    for model, _, leading_coefficients, pattern in branch_factorization_patterns(context.field.prime, context.genus, degree_model):
        infinity_branch = model == "odd"
        for factors in branch_factor_choices_for_pattern(context, pattern):
            for leading_coefficient in leading_coefficients:
                yield BranchDivisorCandidate(
                    leading_coefficient=leading_coefficient,
                    factors=factors,
                    infinity_branch=infinity_branch,
                    factorization_pattern=pattern,
                )
                produced += 1
                if limit is not None and produced >= limit:
                    return


def random_composition(total: int, parts: int, minimum_part: int, rng: random.Random) -> tuple[int, ...]:
    if parts < 1:
        raise ValueError("number of parts must be positive")
    if total < parts * minimum_part:
        raise ValueError("total is too small for the requested number of parts")

    remaining = total - parts * minimum_part
    if parts == 1:
        return (total,)
    cuts = sorted(rng.sample(range(remaining + parts - 1), parts - 1))
    previous = -1
    extras = []
    for cut in [*cuts, remaining + parts - 1]:
        extras.append(cut - previous - 1)
        previous = cut
    return tuple(extra + minimum_part for extra in extras)


def random_branch_factor_count(degree: int, minimum_part: int, max_factors: int, rng: random.Random) -> int:
    largest_factor_count = min(max_factors, degree // minimum_part)
    if largest_factor_count < 1:
        raise ValueError("no valid branch factor count")
    choices = list(range(1, largest_factor_count + 1))
    weights = [1 / (choice * choice) for choice in choices]
    return rng.choices(choices, weights=weights, k=1)[0]


def random_irreducible_choices_for_degree(
    context: EnumerationContext,
    degree: int,
    multiplicity: int,
    rng: random.Random,
) -> tuple[tuple[int, ...], ...]:
    if multiplicity < 1:
        return ()

    cached = context._irreducible_factor_cache.get(degree)
    if cached is not None:
        if multiplicity > len(cached):
            raise ValueError("not enough cached irreducible polynomials for requested multiplicity")
        return tuple(rng.sample(cached, multiplicity))

    sqlite_cache = context._sqlite_irreducible_cache()
    cached_count = sqlite_cache.complete_count(context.field.prime, degree)
    if cached_count is not None:
        if multiplicity > cached_count:
            raise ValueError("not enough SQLite-cached irreducible polynomials for requested multiplicity")
        positions = rng.sample(range(cached_count), multiplicity)
        return tuple(sqlite_cache.get_at_position(context.field.prime, degree, position) for position in positions)

    if multiplicity > monic_irreducible_count(context.field.prime, degree):
        raise ValueError("not enough irreducible polynomials for requested multiplicity")
    return context.random_sage_irreducible_polynomials(degree, multiplicity, rng)


def random_branch_divisor_candidate(
    context: EnumerationContext,
    degree_model: str,
    rng: random.Random,
    max_factors: int,
) -> BranchDivisorCandidate:
    model_choices = []
    if degree_model in {"odd", "both"}:
        model_choices.append("odd")
    if degree_model in {"even", "both"}:
        model_choices.append("even")
    if not model_choices:
        raise ValueError("invalid degree model")

    model = rng.choice(model_choices)
    if model == "odd":
        degree = 2 * context.genus + 1
        minimum_part = 1
        leading_coefficient = 1
        infinity_branch = True
    else:
        degree = 2 * context.genus + 2
        minimum_part = 2
        leading_coefficient = rng.choice((1, _smallest_nonsquare(context.field.prime)))
        infinity_branch = False

    factor_count = random_branch_factor_count(degree, minimum_part, max_factors, rng)
    degrees = random_composition(degree, factor_count, minimum_part, rng)
    factors_by_degree: dict[int, int] = {}
    for factor_degree in degrees:
        factors_by_degree[factor_degree] = factors_by_degree.get(factor_degree, 0) + 1

    factors: list[tuple[int, ...]] = []
    for factor_degree, multiplicity in sorted(factors_by_degree.items()):
        factors.extend(random_irreducible_choices_for_degree(context, factor_degree, multiplicity, rng))
    rng.shuffle(factors)

    pattern = [0] * (degree + 1)
    for factor_degree in degrees:
        pattern[factor_degree] += 1
    return BranchDivisorCandidate(
        leading_coefficient=leading_coefficient,
        factors=tuple(factors),
        infinity_branch=infinity_branch,
        factorization_pattern=tuple(pattern),
    )


def random_branch_divisor_polynomials(
    context: EnumerationContext,
    degree_model: str,
    steps: Optional[int],
    rng: random.Random,
    max_factors: int,
) -> Iterable[BranchDivisorCandidate]:
    if steps is not None and steps < 0:
        raise ValueError("random branch steps must be nonnegative")
    produced = 0
    while steps is None or produced < steps:
        yield random_branch_divisor_candidate(context, degree_model, rng, max_factors)
        produced += 1


def sparse_isomorphism_classes(context: EnumerationContext, max_sparsity: Optional[int]) -> int:
    return sum(
        1
        for record in context.canonical_records.values()
        if record.status_by_max_sparsity.get(max_sparsity) == "sparse"
    )


def sparse_presentations_by_orbit_size(context: EnumerationContext, max_sparsity: Optional[int]) -> int:
    return sum(
        record.orbit_size
        for record in context.canonical_records.values()
        if record.status_by_max_sparsity.get(max_sparsity) == "sparse"
    )


def should_print_progress(processed: int, total: int, interval: int) -> bool:
    if processed == total:
        return True
    if interval <= 0:
        return False
    return processed % interval == 0


def next_progress_threshold(position: int, total: int, interval: int) -> int:
    if interval <= 0:
        return total
    if total < 0:
        return ((position // interval) + 1) * interval
    return min(((position // interval) + 1) * interval, total)


def default_lexicoskip_drought(prime: int, genus: int) -> int:
    return max(2000, min(20000, 4 * prime ** min(genus, 5)))


def default_lexicoskip_initial_skip(prime: int, genus: int) -> int:
    return default_lexicoskip_drought(prime, genus)


def default_lexicoskip_max_skip(prime: int, genus: int) -> int:
    return 20 * default_lexicoskip_initial_skip(prime, genus)


def default_lexicoskip_probe_window(prime: int, genus: int) -> int:
    return max(500, min(5000, 4 * prime ** min(genus, 4)))


def is_new_isomorphism_class_result(result: dict[str, object]) -> bool:
    return result.get("status") in {"sparse", "rejected_hasse_witt", "rejected_exact"}


def progress_line(
    processed: int,
    total: int,
    skipped: int,
    context: EnumerationContext,
    max_sparsity: Optional[int],
) -> str:
    total_label = "?" if total < 0 else str(total)
    fields = [
        f"prime: {context.field.prime}",
        f"genus: {context.genus}",
        f"progress: {processed}/{total_label}",
        f"sparse_presentations: {sparse_presentations_by_orbit_size(context, max_sparsity)}",
        f"sparse_isomorphism_classes: {sparse_isomorphism_classes(context, max_sparsity)}",
        f"canonicalized_isomorphism_classes: {len(context.canonical_records)}",
    ]
    fields.append("-")
    return "\n".join(fields)


def load_irreducibles_for_branch_mode(context: EnumerationContext, max_degree: int) -> None:
    print(f"irreducible_load: 0/{max_degree}", flush=True)
    for degree in range(1, max_degree + 1):
        estimated_mb = estimated_irreducible_tuple_memory_bytes(context.field.prime, degree) / BYTES_PER_MEGABYTE
        if not context.should_materialize_branch_irreducibles(degree):
            print(
                f"irreducible_load: {degree}/{max_degree} "
                f"degree={degree} mode=streaming reason=not-in-sqlite estimated_memory_mb={estimated_mb:.1f}",
                flush=True,
            )
            continue

        started_at = perf_counter()
        from_cache = context._sqlite_irreducible_cache().has_complete(context.field.prime, degree)
        polynomials = context._irreducible_polynomials(degree)
        source = "sqlite" if from_cache else "necklace"
        print(
            f"irreducible_load: {degree}/{max_degree} "
            f"degree={degree} mode=memory source={source} count={len(polynomials)} "
            f"estimated_memory_mb={estimated_mb:.1f} seconds={perf_counter() - started_at:.6f}",
            flush=True,
        )


def precompute_irreducible_polynomials(
    prime: int,
    max_degree: int,
    cache_path: str | Path | None = None,
    memory_budget_mb: int = DEFAULT_IRREDUCIBLE_MEMORY_BUDGET_MB,
) -> list[dict[str, object]]:
    if max_degree < 1:
        raise ValueError("max degree must be positive")

    field = PrimeField(prime)
    cache = IrreduciblePolynomialCache(cache_path if cache_path is not None else default_irreducible_cache_path(prime))
    memory_budget = _memory_budget_bytes(memory_budget_mb)
    summaries: list[dict[str, object]] = []
    try:
        for degree in range(1, max_degree + 1):
            started_at = perf_counter()
            estimated_bytes = estimated_irreducible_tuple_memory_bytes(prime, degree)
            if estimated_bytes > memory_budget:
                raise MemoryError(
                    "precomputing this irreducible degree would exceed the configured memory budget "
                    f"({estimated_bytes / BYTES_PER_MEGABYTE:.1f} MB needed for p={prime}, "
                    f"degree={degree}; budget is {memory_budget_mb} MB)"
                )
            polynomials = cache.get_complete(prime, degree)
            from_cache = polynomials is not None
            if polynomials is None:
                polynomials = _generate_monic_irreducible_polynomials(field, degree)
                cache.store_complete(prime, degree, polynomials)
            summaries.append(
                {
                    "prime": prime,
                    "degree": degree,
                    "count": len(polynomials),
                    "estimated_memory_mb": estimated_bytes / BYTES_PER_MEGABYTE,
                    "from_cache": from_cache,
                    "seconds": perf_counter() - started_at,
                }
            )
    finally:
        cache.close()
    return summaries


def parse_enumeration_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Enumerate hyperelliptic curves and write sparse L-polynomial data to SQLite.")
    parser.add_argument("--p", type=int, required=True, help="Odd prime field characteristic.")
    parser.add_argument("--genus", "-g", type=int, help="Curve genus.")
    parser.add_argument("--genus-start", type=int, help="First genus for a batch run over increasing genera.")
    parser.add_argument("--genus-end", type=int, help="Last genus for a batch run over increasing genera, inclusive.")
    parser.add_argument("--genus-step", type=int, default=1, help="Genus increment for batch runs. Defaults to 1.")
    parser.add_argument(
        "--max-sparsity",
        type=int,
        help="Optional maximum allowed sparsity among a_1, ..., a_{g-1}. Omit to compute without a sparsity restriction.",
    )
    parser.add_argument("--out", type=Path, help="SQLite output path. Defaults to data_gen/results/p{p}_g{g}_s_{max}.sqlite.")
    parser.add_argument("--out-dir", type=Path, default=Path("data_gen") / "results", help="Output directory for batch runs. Defaults to data_gen/results.")
    parser.add_argument(
        "--irreducible-cache",
        type=Path,
        help="SQLite cache path for monic irreducible polynomials. Defaults to data_gen/irreducibles/irreducibles_p{p}.sqlite.",
    )
    parser.add_argument(
        "--precompute-irreducibles-up-to-degree",
        type=int,
        help="Precompute monic irreducible polynomials over F_p for degrees 1..N, then exit.",
    )
    parser.add_argument(
        "--irreducible-memory-budget-mb",
        type=int,
        default=DEFAULT_IRREDUCIBLE_MEMORY_BUDGET_MB,
        help="Maximum memory, in MB, to use when materializing irreducible polynomial tables. Defaults to 1024.",
    )
    parser.add_argument("--limit", type=int, help="Optional maximum number of branch-divisor presentations to process.")
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=1000,
        help="Print progress every N branch-divisor presentations. Use 0 to print only final output.",
    )
    parser.add_argument(
        "--enumeration-mode",
        choices=("enumerate", "random"),
        default="enumerate",
        help="Use deterministic branch-divisor enumeration or random branch-divisor sparse search.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        help="Seed for random mode. Omit for non-deterministic sampling.",
    )
    parser.add_argument(
        "--random-steps",
        type=int,
        help="Number of random branch divisors to sample. In random mode, defaults to --limit when provided; if both are omitted, runs until interrupted.",
    )
    parser.add_argument(
        "--random-max-factors",
        type=int,
        default=5,
        help="Maximum number of irreducible branch factors sampled in random mode. Defaults to 5.",
    )
    parser.add_argument(
        "--branch-random-max-factors",
        dest="random_max_factors",
        type=int,
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def _is_batch_run(args: argparse.Namespace) -> bool:
    return getattr(args, "genus_start", None) is not None or getattr(args, "genus_end", None) is not None


def _batch_genus_values(args: argparse.Namespace) -> range:
    genus_start = getattr(args, "genus_start", None)
    genus_end = getattr(args, "genus_end", None)
    genus_step = getattr(args, "genus_step", 1)
    if genus_start is None or genus_end is None:
        raise ValueError("--genus-start and --genus-end must be supplied together")
    if genus_step <= 0:
        raise ValueError("--genus-step must be positive")
    if genus_start > genus_end:
        raise ValueError("--genus-start must be at most --genus-end")
    return range(genus_start, genus_end + 1, genus_step)


def run_batch_enumeration_from_args(args: argparse.Namespace) -> None:
    if getattr(args, "genus", None) is not None:
        raise ValueError("--genus cannot be combined with --genus-start/--genus-end")
    if getattr(args, "out", None) is not None:
        raise ValueError("--out is for single-genus runs; use --out-dir for batch runs")
    if getattr(args, "precompute_irreducibles_up_to_degree", None) is not None:
        raise ValueError("--precompute-irreducibles-up-to-degree cannot be combined with genus batch runs")

    out_dir = getattr(args, "out_dir", Path("data_gen") / "results")
    out_dir.mkdir(parents=True, exist_ok=True)
    genus_values = list(_batch_genus_values(args))
    for index, genus in enumerate(genus_values, start=1):
        batch_args = argparse.Namespace(**vars(args))
        batch_args.genus = genus
        batch_args.genus_start = None
        batch_args.genus_end = None
        batch_args.out = out_dir / default_sqlite_path(args.p, genus, args.max_sparsity).name
        print(f"batch: {index}/{len(genus_values)} p={args.p} genus={genus} out={batch_args.out}", flush=True)
        run_enumeration_from_args(batch_args)


def run_enumeration_from_args(args: argparse.Namespace) -> None:
    if _is_batch_run(args):
        run_batch_enumeration_from_args(args)
        return

    irreducible_cache_path = getattr(args, "irreducible_cache", None)
    if irreducible_cache_path is None:
        irreducible_cache_path = default_irreducible_cache_path(args.p)
    irreducible_memory_budget_mb = getattr(
        args,
        "irreducible_memory_budget_mb",
        DEFAULT_IRREDUCIBLE_MEMORY_BUDGET_MB,
    )
    precompute_irreducibles_up_to_degree = getattr(args, "precompute_irreducibles_up_to_degree", None)
    if precompute_irreducibles_up_to_degree is not None:
        summaries = precompute_irreducible_polynomials(
            prime=args.p,
            max_degree=precompute_irreducibles_up_to_degree,
            cache_path=irreducible_cache_path,
            memory_budget_mb=irreducible_memory_budget_mb,
        )
        for summary in summaries:
            source = "cache" if summary["from_cache"] else "generated"
            print(
                "irreducibles: "
                f"p={summary['prime']} "
                f"degree={summary['degree']} "
                f"count={summary['count']} "
                f"estimated_memory_mb={summary['estimated_memory_mb']:.1f} "
                f"source={source} "
                f"seconds={summary['seconds']:.6f}",
                flush=True,
            )
        print(f"irreducible_cache: {irreducible_cache_path}")
        return

    if args.genus is None:
        raise ValueError("--genus is required unless --precompute-irreducibles-up-to-degree is used")

    enumeration_mode = getattr(args, "enumeration_mode", "enumerate")
    if enumeration_mode not in {"enumerate", "random"}:
        raise ValueError("enumeration mode must be 'enumerate' or 'random'")

    degree_model = "both"
    random_max_factors = getattr(args, "random_max_factors", 5)
    if random_max_factors < 1:
        raise ValueError("random max factors must be positive")

    random_steps = getattr(args, "random_steps", None)
    if enumeration_mode == "random":
        if random_steps is None:
            random_steps = args.limit
        if random_steps is not None and random_steps < 0:
            raise ValueError("random steps must be nonnegative")

    sqlite_path = args.out if args.out is not None else default_sqlite_path(args.p, args.genus, args.max_sparsity)
    sqlite_path.parent.mkdir(parents=True, exist_ok=True)

    if enumeration_mode == "enumerate":
        total = total_branch_divisor_presentations(
            prime=args.p,
            genus=args.genus,
            degree_model=degree_model,
            limit=args.limit,
        )
    elif enumeration_mode == "random":
        total = random_steps if random_steps is not None else -1
    stats: dict[str, int] = {"processed": 0, "skipped": 0, "position": 0}
    timing: Optional[dict[str, object]] = None
    last_snapshot = (0, 0, 0)
    last_snapshot_position: Optional[int] = None
    context = EnumerationContext(
        prime=args.p,
        genus=args.genus,
        sqlite_path=sqlite_path,
        skip_even_models_with_rational_branch_point=True,
        irreducible_cache_path=irreducible_cache_path,
        irreducible_memory_budget_mb=irreducible_memory_budget_mb,
        cache_full_orbits=enumeration_mode == "enumerate",
        canonicalize_branch_before_exact=enumeration_mode == "enumerate",
        sqlite_write_batch_size=SQLITE_WRITE_BATCH_SIZE,
    )
    try:
        load_irreducibles_for_branch_mode(context, 2 * args.genus + 2)

        next_progress = next_progress_threshold(0, total, args.progress_interval)
        print(progress_line(0, total, 0, context, args.max_sparsity), flush=True)
        last_snapshot = context.write_progress_snapshot(
            position=0,
            processed=0,
            skipped=0,
            max_sparsity=args.max_sparsity,
            previous_sparse_presentations=last_snapshot[0],
            previous_sparse_isomorphism_classes=last_snapshot[1],
            previous_canonicalized_isomorphism_classes=last_snapshot[2],
        )
        last_snapshot_position = 0

        if enumeration_mode == "enumerate":
            vectors = branch_divisor_polynomials(
                context,
                degree_model=degree_model,
                limit=args.limit,
            )
        else:
            vectors = random_branch_divisor_polynomials(
                context,
                degree_model=degree_model,
                steps=random_steps,
                rng=random.Random(getattr(args, "random_seed", None)),
                max_factors=random_max_factors,
            )
        for vector in vectors:
            result = context.process_branch_divisor_for_output(vector, max_sparsity=args.max_sparsity)
            status = str(result["status"])
            stats["processed"] += 1
            stats["position"] = stats["processed"]
            stats[status] = stats.get(status, 0) + 1

            if stats["position"] >= next_progress:
                print(
                    progress_line(
                        stats["position"],
                        total,
                        stats["skipped"],
                        context,
                        args.max_sparsity,
                    ),
                    flush=True,
                )
                last_snapshot = context.write_progress_snapshot(
                    position=stats["position"],
                    processed=stats["processed"],
                    skipped=stats["skipped"],
                    max_sparsity=args.max_sparsity,
                    previous_sparse_presentations=last_snapshot[0],
                    previous_sparse_isomorphism_classes=last_snapshot[1],
                    previous_canonicalized_isomorphism_classes=last_snapshot[2],
                )
                last_snapshot_position = stats["position"]
                next_progress = next_progress_threshold(stats["position"], total, args.progress_interval)

    finally:
        if last_snapshot_position != stats["position"]:
            context.write_progress_snapshot(
                position=stats["position"],
                processed=stats["processed"],
                skipped=stats["skipped"],
                max_sparsity=args.max_sparsity,
                previous_sparse_presentations=last_snapshot[0],
                previous_sparse_isomorphism_classes=last_snapshot[1],
                previous_canonicalized_isomorphism_classes=last_snapshot[2],
            )
        context.write_enumeration_summary(
            max_sparsity=args.max_sparsity,
            degree_model=degree_model,
            enumeration_mode=enumeration_mode,
            limit=random_steps if enumeration_mode == "random" else args.limit,
            total_coefficient_vectors=total,
            processed=stats["processed"],
            skipped=stats["skipped"],
            final_position=stats["position"],
            sparse_presentations=sparse_presentations_by_orbit_size(context, args.max_sparsity),
        )
        timing = context.timing_summary()
        context.close_sqlite()
        context.close_irreducible_cache()

    print(f"output: {sqlite_path}")
    print(f"stats: {stats}")
    print(f"timing: {timing}")


def main() -> None:
    run_enumeration_from_args(parse_enumeration_args())


if __name__ == "__main__":
    main()
