import argparse
import contextlib
import io
import json
import random
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from data_gen.hyperelliptic import (
    DEFAULT_IRREDUCIBLE_MEMORY_BUDGET_MB,
    EnumerationContext,
    FiniteExtension,
    HyperellipticCurve,
    PointCountingContext,
    Polynomial,
    PrimeField,
    _affine_polynomial_to_binary_form,
    _apply_binary_form_action_matrix,
    _extension_polynomial_mul,
    _extension_polynomial_product_tree,
    _generate_monic_irreducible_polynomials,
    _transform_binary_form,
    branch_divisor_polynomials,
    branch_factorization_patterns,
    coefficient_vector_at_index,
    coefficient_vectors,
    coefficient_vectors_by_support,
    default_lexicoskip_drought,
    default_lexicoskip_initial_skip,
    default_lexicoskip_max_skip,
    default_lexicoskip_probe_window,
    estimated_irreducible_tuple_memory_bytes,
    l_polynomial_coefficients_mod_p_from_branch_factors,
    monic_irreducible_count,
    polynomial_from_branch_factors,
    precompute_irreducible_polynomials,
    random_branch_divisor_polynomials,
    random_irreducible_choices_for_degree,
    run_enumeration_from_args,
    sparse_presentations_by_orbit_size,
    total_branch_divisor_presentations,
    total_coefficient_vectors,
)


class HyperellipticTests(unittest.TestCase):
    def test_coefficient_vectors_use_normalized_leading_coefficients(self):
        odd_vectors = list(coefficient_vectors(prime=5, genus=1, degree_model="odd", limit=2))
        self.assertEqual([vector[-1] for vector in odd_vectors], [1, 1])

        even_vectors = list(coefficient_vectors(prime=5, genus=1, degree_model="even", limit=4))
        self.assertEqual([vector[-1] for vector in even_vectors], [1, 2, 1, 2])

        self.assertEqual(total_coefficient_vectors(prime=5, genus=1, degree_model="odd", limit=None), 125)
        self.assertEqual(total_coefficient_vectors(prime=5, genus=1, degree_model="even", limit=None), 1250)
        self.assertEqual(total_coefficient_vectors(prime=5, genus=1, degree_model="both", limit=None), 1375)
        indexed_vectors = [
            coefficient_vector_at_index(prime=5, genus=1, degree_model="even", index=index)
            for index in range(4)
        ]
        self.assertEqual(indexed_vectors, even_vectors)

    def test_coefficient_vectors_by_support_increase_support_size(self):
        vectors = list(coefficient_vectors_by_support(prime=5, genus=1, degree_model="both", limit=8))

        self.assertEqual(vectors[0], (0, 0, 0, 1))
        self.assertEqual(vectors[1], (0, 0, 0, 0, 1))
        self.assertEqual(vectors[2], (0, 0, 0, 0, 2))
        support_sizes = [sum(1 for coefficient in vector if coefficient != 0) for vector in vectors]
        self.assertEqual(support_sizes, sorted(support_sizes))

    def test_branch_factorization_patterns_skip_even_linear_factors(self):
        patterns = list(branch_factorization_patterns(prime=3, genus=1, degree_model="both"))
        even_patterns = [pattern for model, _, _, pattern in patterns if model == "even"]

        self.assertTrue(even_patterns)
        self.assertTrue(all(pattern[1] == 0 for pattern in even_patterns))
        self.assertEqual(total_branch_divisor_presentations(prime=3, genus=1, degree_model="both", limit=None), 60)

    def test_branch_divisor_polynomials_generate_squarefree_normalized_models(self):
        context = EnumerationContext(prime=3, genus=1)
        candidates = list(branch_divisor_polynomials(context, degree_model="both", limit=10))
        vectors = [
            polynomial_from_branch_factors(
                candidate.factors,
                candidate.leading_coefficient,
                context.field.prime,
            )
            for candidate in candidates
        ]

        self.assertEqual(len(vectors), 10)
        self.assertTrue(all(vector[-1] in {1, 2} for vector in vectors))
        self.assertTrue(all(context.polynomial(vector).is_squarefree() for vector in vectors))
        context.close_irreducible_cache()

    def test_factorized_hasse_witt_matches_expanded_polynomial(self):
        context = EnumerationContext(prime=5, genus=2)
        candidates = list(branch_divisor_polynomials(context, degree_model="both", limit=12))

        for candidate in candidates:
            coefficients = polynomial_from_branch_factors(
                candidate.factors,
                candidate.leading_coefficient,
                context.field.prime,
            )
            self.assertEqual(
                l_polynomial_coefficients_mod_p_from_branch_factors(
                    candidate.factors,
                    candidate.leading_coefficient,
                    context.genus,
                    context.field.prime,
                ),
                context.curve(context.polynomial(coefficients)).l_polynomial_coefficients_mod_p(),
            )
            self.assertEqual(
                context._branch_l_polynomial_coefficients_mod_p(candidate),
                context.curve(context.polynomial(coefficients)).l_polynomial_coefficients_mod_p(),
            )
        self.assertGreater(len(context.branch_factor_hasse_witt_power_cache), 0)
        context.close_irreducible_cache()

    def test_factorized_exact_lpoly_matches_expanded_polynomial(self):
        context = EnumerationContext(prime=5, genus=2)
        candidates = list(branch_divisor_polynomials(context, degree_model="both", limit=12))

        for candidate in candidates:
            coefficients = polynomial_from_branch_factors(
                candidate.factors,
                candidate.leading_coefficient,
                context.field.prime,
            )
            self.assertEqual(
                context._branch_l_polynomial_coefficients(candidate, max_sparsity=None),
                context.curve(context.polynomial(coefficients)).l_polynomial_coefficients(),
            )
        context.close_irreducible_cache()

    def test_factorized_character_point_counts_match_expanded_polynomial(self):
        context = EnumerationContext(prime=3, genus=3)
        candidates = list(branch_divisor_polynomials(context, degree_model="both", limit=8))

        for candidate in candidates:
            coefficients = polynomial_from_branch_factors(
                candidate.factors,
                candidate.leading_coefficient,
                context.field.prime,
            )
            curve = context.curve(context.polynomial(coefficients))
            for extension_degree in range(1, context.genus + 1):
                self.assertEqual(
                    context._branch_point_count_over_extension(candidate, extension_degree),
                    curve.point_count_over_extension(extension_degree),
                )
        self.assertGreater(len(context.branch_factor_character_cache), 0)
        context.close_irreducible_cache()

    def test_exact_branch_key_duplicate_short_circuits_hasse_witt(self):
        context = EnumerationContext(prime=3, genus=2)
        candidate = next(iter(branch_divisor_polynomials(context, degree_model="both", limit=1)))
        first = context.process_branch_divisor_for_output(candidate, max_sparsity=1)
        self.assertIn("status", first)

        with mock.patch("data_gen.hyperelliptic.l_polynomial_coefficients_mod_p_from_branch_factors", side_effect=AssertionError):
            second = context.process_branch_divisor_for_output(candidate, max_sparsity=1)

        self.assertEqual(second["status"], "duplicate")
        self.assertEqual(second["matched_by"], "exact_branch_key")
        context.close_irreducible_cache()

    def test_random_sparse_first_rejects_exact_failures_before_canonicalization(self):
        context = EnumerationContext(prime=3, genus=2, canonicalize_branch_before_exact=False)
        candidate = next(iter(branch_divisor_polynomials(context, degree_model="both", limit=1)))

        with (
            mock.patch.object(context, "_branch_l_polynomial_coefficients_mod_p", return_value=[0, 0]),
            mock.patch.object(context, "_branch_l_polynomial_coefficients", return_value=None) as exact_lpoly,
            mock.patch.object(context, "_canonical_key_and_orbit", side_effect=AssertionError),
            mock.patch("data_gen.hyperelliptic.polynomial_from_branch_factors", side_effect=AssertionError),
        ):
            result = context.process_branch_divisor_for_output(candidate, max_sparsity=0)

        self.assertEqual(result["status"], "rejected_exact_uncanonicalized")
        self.assertEqual(result["sparsity"], -1)
        self.assertEqual(len(context.canonical_records), 0)
        exact_lpoly.assert_called_once_with(candidate, max_sparsity=0)
        context.close_irreducible_cache()

    def test_random_sparse_first_reuses_exact_lpoly_for_sparse_survivor(self):
        context = EnumerationContext(prime=3, genus=2, canonicalize_branch_before_exact=False)
        candidate = next(iter(branch_divisor_polynomials(context, degree_model="both", limit=1)))
        original_exact_lpoly = context._branch_l_polynomial_coefficients

        with mock.patch.object(context, "_branch_l_polynomial_coefficients", wraps=original_exact_lpoly) as exact_lpoly:
            result = context.process_branch_divisor_for_output(candidate, max_sparsity=None)

        self.assertEqual(result["status"], "sparse")
        self.assertEqual(exact_lpoly.call_count, 1)
        self.assertEqual(len(context.canonical_records), 1)
        context.close_irreducible_cache()

    def test_branch_ground_invariants_match_expanded_polynomial(self):
        context = EnumerationContext(prime=5, genus=2)
        candidates = list(branch_divisor_polynomials(context, degree_model="both", limit=12))

        for candidate in candidates:
            coefficients = polynomial_from_branch_factors(
                candidate.factors,
                candidate.leading_coefficient,
                context.field.prime,
            )
            self.assertEqual(context.branch_ground_invariants(candidate), context.ground_invariants(coefficients))
        context.close_irreducible_cache()

    def test_random_branch_divisor_polynomials_are_reproducible_squarefree_models(self):
        first_context = EnumerationContext(prime=3, genus=2)
        second_context = EnumerationContext(prime=3, genus=2)
        first_candidates = list(
            random_branch_divisor_polynomials(
                first_context,
                degree_model="both",
                steps=5,
                rng=random.Random(17),
                max_factors=3,
            )
        )
        second_candidates = list(
            random_branch_divisor_polynomials(
                second_context,
                degree_model="both",
                steps=5,
                rng=random.Random(17),
                max_factors=3,
            )
        )

        self.assertEqual(first_candidates, second_candidates)
        self.assertEqual(len(first_candidates), 5)
        for candidate in first_candidates:
            vector = polynomial_from_branch_factors(
                candidate.factors,
                candidate.leading_coefficient,
                first_context.field.prime,
            )
            self.assertIn(sum(candidate.factorization_pattern[index] * index for index in range(len(candidate.factorization_pattern))), {5, 6})
            self.assertTrue(first_context.polynomial(vector).is_squarefree())
        first_context.close_irreducible_cache()
        second_context.close_irreducible_cache()

    def test_random_missing_irreducible_uses_sage_generator(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            context = EnumerationContext(
                prime=3,
                genus=2,
                irreducible_cache_path=Path(tmpdir) / "irreducibles.sqlite",
            )
            with mock.patch.object(
                context,
                "random_sage_irreducible_polynomials",
                return_value=((1, 0, 1, 1),),
            ) as sage_generator:
                choices = random_irreducible_choices_for_degree(context, degree=3, multiplicity=1, rng=random.Random(1))

            self.assertEqual(choices, ((1, 0, 1, 1),))
            sage_generator.assert_called_once()
            context.close_irreducible_cache()

    def test_branch_divisor_duplicate_matching_avoids_expansion(self):
        context = EnumerationContext(prime=3, genus=1)
        candidates = list(branch_divisor_polynomials(context, degree_model="both", limit=20))

        first = context.process_branch_divisor_for_output(candidates[18], max_sparsity=None)
        self.assertIn(first["status"], {"sparse", "rejected_exact", "rejected_hasse_witt"})

        with mock.patch("data_gen.hyperelliptic.polynomial_from_branch_factors", side_effect=AssertionError):
            duplicate = context.process_branch_divisor_for_output(candidates[19], max_sparsity=None)

        self.assertEqual(duplicate["status"], "duplicate")
        self.assertEqual(duplicate["matched_by"], "factorized_branch_orbit")

    def test_default_lexicoskip_parameters_depend_on_prime_and_genus(self):
        self.assertEqual(default_lexicoskip_drought(5, 4), 2500)
        self.assertEqual(default_lexicoskip_initial_skip(5, 4), 2500)
        self.assertEqual(default_lexicoskip_max_skip(5, 4), 50000)
        self.assertEqual(default_lexicoskip_probe_window(5, 4), 2500)

    def test_default_irreducible_memory_budget_is_one_gb(self):
        self.assertEqual(DEFAULT_IRREDUCIBLE_MEMORY_BUDGET_MB, 1024)
        self.assertEqual(monic_irreducible_count(5, 3), 40)
        self.assertGreater(estimated_irreducible_tuple_memory_bytes(5, 3), 0)

    def test_irreducible_polynomial_sqlite_cache_reloads_complete_blocks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_path = Path(tmpdir) / "irreducibles.sqlite"
            context = EnumerationContext(prime=5, genus=2, irreducible_cache_path=cache_path)
            degree_two = context._irreducible_polynomials(2)
            context.close_irreducible_cache()

            self.assertGreater(len(degree_two), 0)
            connection = sqlite3.connect(cache_path)
            self.assertEqual(
                connection.execute(
                    "SELECT count FROM irreducible_cache_metadata WHERE prime = 5 AND degree = 2"
                ).fetchone()[0],
                len(degree_two),
            )
            connection.close()

            reloaded_context = EnumerationContext(prime=5, genus=2, irreducible_cache_path=cache_path)
            with mock.patch("data_gen.hyperelliptic._is_irreducible", side_effect=AssertionError):
                self.assertEqual(reloaded_context._irreducible_polynomials(2), degree_two)
            reloaded_context.close_irreducible_cache()

    def test_irreducible_polynomial_generator_uses_necklaces_not_trial_filtering(self):
        field = PrimeField(5)
        with mock.patch("data_gen.hyperelliptic._is_irreducible", side_effect=AssertionError):
            degree_three = _generate_monic_irreducible_polynomials(field, 3)

        self.assertEqual(len(degree_three), 40)
        self.assertEqual(len(set(degree_three)), 40)
        self.assertTrue(all(polynomial[-1] == 1 and len(polynomial) == 4 for polynomial in degree_three))

    def test_precompute_irreducible_polynomials_populates_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_path = Path(tmpdir) / "irreducibles.sqlite"
            summaries = precompute_irreducible_polynomials(prime=5, max_degree=2, cache_path=cache_path)
            self.assertEqual([summary["degree"] for summary in summaries], [1, 2])
            self.assertEqual([summary["from_cache"] for summary in summaries], [False, False])
            self.assertIn("estimated_memory_mb", summaries[0])

            cached_summaries = precompute_irreducible_polynomials(prime=5, max_degree=2, cache_path=cache_path)
            self.assertEqual([summary["from_cache"] for summary in cached_summaries], [True, True])

    def test_irreducible_materialization_respects_memory_budget(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_path = Path(tmpdir) / "irreducibles.sqlite"
            context = EnumerationContext(
                prime=5,
                genus=2,
                irreducible_cache_path=cache_path,
                irreducible_memory_budget_mb=1,
            )
            with self.assertRaises(MemoryError):
                context._irreducible_polynomials(8)
            context.close_irreducible_cache()

    def test_polynomial_coefficients_must_be_canonical(self):
        field = PrimeField(5)
        with self.assertRaises(ValueError):
            Polynomial(field, [1, 5])
        with self.assertRaises(ValueError):
            Polynomial(field, [-1, 1])

    def test_rejects_repeated_root_model(self):
        field = PrimeField(5)
        with self.assertRaises(ValueError):
            HyperellipticCurve(Polynomial(field, [1, 0, 0, 0, 0, 1]))

    def test_point_contribution_tables(self):
        field = PrimeField(5)
        context = PointCountingContext(field, polynomial_degree=3)
        self.assertEqual(context.ground_point_contributions(), (1, 2, 0, 0, 2))

        extension = FiniteExtension(5, 2)
        contributions = extension.point_contributions()
        self.assertEqual(contributions[extension.zero()], 1)
        self.assertEqual(sum(contributions.values()), extension.size)

    def test_frobenius_matrix_matches_p_power_map(self):
        extension = FiniteExtension(5, 3)
        matrix = extension.frobenius_matrix()

        for element in extension.elements()[:20]:
            self.assertEqual(extension.apply_frobenius(element), extension.pow(element, extension.prime))
        self.assertIs(extension.frobenius_matrix(), matrix)

    def test_integer_encoded_extension_arithmetic_matches_tuple_arithmetic(self):
        extension = FiniteExtension(5, 3)
        polynomial = (2, 0, 4, 1)

        for lhs in extension.elements()[:10]:
            lhs_int = extension.encode(lhs)
            expected_value = extension.zero()
            for coefficient in reversed(polynomial):
                expected_value = extension.add(extension.multiply(expected_value, lhs), extension.constant(coefficient))
            self.assertEqual(extension.decode(lhs_int), lhs)
            self.assertEqual(extension.int_quadratic_character(lhs_int), 0 if lhs == extension.zero() else 1 if extension.is_square(lhs) else -1)
            self.assertEqual(
                extension.decode(extension.int_evaluate_polynomial(polynomial, lhs_int)),
                expected_value,
            )
            for rhs in extension.elements()[5:15]:
                rhs_int = extension.encode(rhs)
                self.assertEqual(
                    extension.decode(extension.int_add(lhs_int, rhs_int)),
                    extension.add(lhs, rhs),
                )
                self.assertEqual(
                    extension.decode(extension.int_multiply(lhs_int, rhs_int)),
                    extension.multiply(lhs, rhs),
                )
            self.assertEqual(
                extension.decode(extension.int_pow(lhs_int, 7)),
                extension.pow(lhs, 7),
            )

    def test_extension_polynomial_product_tree_matches_sequential_multiplication(self):
        extension = FiniteExtension(5, 2)
        factors = tuple(
            (
                tuple((-entry) % extension.prime for entry in element),
                extension.one(),
            )
            for element in extension.elements()[1:6]
        )
        sequential = factors[0]
        for factor in factors[1:]:
            sequential = _extension_polynomial_mul(extension, sequential, factor)

        self.assertEqual(_extension_polynomial_product_tree(extension, factors), sequential)

    def test_ground_value_table(self):
        field = PrimeField(5)
        context = PointCountingContext(field, polynomial_degree=3)
        value_table = context.ground_value_table()

        self.assertEqual(value_table[0][3], (3, 3, 3, 3, 3))
        self.assertEqual(value_table[2][2], (0, 2, 3, 3, 2))

        curve = HyperellipticCurve(Polynomial(field, [1, 1, 0, 1]), point_counting_context=context)
        self.assertEqual(curve.point_count_over_ground_field(), 9)
        self.assertIs(context.ground_value_table(), value_table)

    def test_enumeration_context_combines_ground_invariants(self):
        context = EnumerationContext(prime=5, genus=2)
        coefficients = [1, 1, 0, 0, 0, 1]
        curve = context.curve(coefficients)
        self.assertEqual(context.ground_invariants(coefficients), (context.rational_branch_count(coefficients), curve.point_count_over_ground_field()))

    def test_genus_one_l_polynomial_coefficient(self):
        field = PrimeField(5)
        curve = HyperellipticCurve(Polynomial(field, [1, 1, 0, 1]))

        self.assertEqual(curve.point_count_over_ground_field(), 9)
        self.assertEqual(curve.point_count_over_extension(1), 9)
        self.assertEqual(curve.l_polynomial_coefficients(), [3])
        self.assertEqual(curve.l_polynomial_coefficients_with_sparsity_limit(0), [3])

    def test_genus_two_l_polynomial_coefficients(self):
        field = PrimeField(5)
        curve = HyperellipticCurve(Polynomial(field, [1, 1, 0, 0, 0, 1]))

        self.assertEqual(curve.point_count_over_extension(1), 6)
        self.assertEqual(curve.point_count_over_extension(2), 46)
        self.assertEqual(curve.l_polynomial_coefficients(), [0, 10])
        self.assertEqual(curve.l_polynomial_coefficients_with_sparsity_limit(0), [0, 10])

    def test_hasse_witt_matrix_and_l_polynomial_mod_p(self):
        field = PrimeField(5)
        elliptic_curve = HyperellipticCurve(Polynomial(field, [1, 1, 0, 1]))
        self.assertEqual(elliptic_curve.hasse_witt_matrix(), ((2,),))
        self.assertEqual(elliptic_curve.l_polynomial_coefficients_mod_p(), [3])

        genus_two_curve = HyperellipticCurve(Polynomial(field, [1, 1, 0, 0, 0, 1]))
        self.assertEqual(genus_two_curve.hasse_witt_matrix(), ((0, 0), (0, 0)))
        self.assertEqual(genus_two_curve.l_polynomial_coefficients_mod_p(), [0, 0])

    def test_hasse_witt_sparsity_filter(self):
        field = PrimeField(5)
        context = PointCountingContext(field, polynomial_degree=5)
        curve = HyperellipticCurve(Polynomial(field, [1, 0, 1, 0, 0, 1]), point_counting_context=context)

        self.assertEqual(curve.l_polynomial_coefficients(), [-1, 0])
        self.assertEqual(curve.l_polynomial_coefficients_mod_p(), [4, 0])
        self.assertFalse(curve.passes_hasse_witt_sparsity_filter(0))
        self.assertTrue(curve.passes_hasse_witt_sparsity_filter(1))

        context = PointCountingContext(field, polynomial_degree=5)
        curve = HyperellipticCurve(Polynomial(field, [1, 0, 1, 0, 0, 1]), point_counting_context=context)
        self.assertIsNone(curve.l_polynomial_coefficients_with_sparsity_limit(0))
        self.assertEqual(context._extension_cache, {})

    def test_sparsity_limit_stops_when_exceeded(self):
        field = PrimeField(5)
        curve = HyperellipticCurve(Polynomial(field, [1, 4, 0, 0, 0, 1]))

        self.assertIsNone(curve.l_polynomial_coefficients_with_sparsity_limit(0))

    def test_point_counting_context_can_be_shared(self):
        field = PrimeField(5)
        context = PointCountingContext(field, polynomial_degree=5)
        sparse_curve = HyperellipticCurve(Polynomial(field, [1, 1, 0, 0, 0, 1]), point_counting_context=context)
        nonsparse_curve = HyperellipticCurve(Polynomial(field, [1, 4, 0, 0, 0, 1]), point_counting_context=context)

        self.assertEqual(sparse_curve.l_polynomial_coefficients(), [0, 10])
        self.assertIsNone(nonsparse_curve.l_polynomial_coefficients_with_sparsity_limit(0))
        self.assertIn(2, context._extension_cache)
        self.assertIn(2, context._power_cache)

    def test_point_counting_context_rejects_incompatible_curves(self):
        context = PointCountingContext(PrimeField(5), polynomial_degree=3)

        with self.assertRaises(ValueError):
            HyperellipticCurve(Polynomial(PrimeField(7), [1, 1, 0, 1]), point_counting_context=context)
        with self.assertRaises(ValueError):
            HyperellipticCurve(Polynomial(PrimeField(5), [1, 1, 0, 0, 0, 1]), point_counting_context=context)

    def test_enumeration_context_canonicalizes_isomorphic_polynomials(self):
        context = EnumerationContext(prime=5, genus=2)
        f = context.polynomial([1, 1, 0, 0, 0, 1])
        translated_f = context.polynomial([3, 1, 0, 0, 0, 1])
        square_scaled_f = context.polynomial([4, 4, 0, 0, 0, 4])

        self.assertEqual(context.canonical_key(f), context.canonical_key(translated_f))
        self.assertEqual(context.canonical_key(f), context.canonical_key(square_scaled_f))
        self.assertTrue(context.is_new_isomorphism_class(f))
        self.assertFalse(context.is_new_isomorphism_class(translated_f))
        self.assertEqual(len(context.pgl2_action_matrices), len(context.pgl2))

    def test_cached_pgl2_action_matrices_match_direct_transforms(self):
        context = EnumerationContext(prime=5, genus=2)
        binary_form = _affine_polynomial_to_binary_form(context.polynomial([1, 1, 0, 0, 0, 1]), context.binary_degree)

        for matrix, action_matrix in zip(context.pgl2, context.pgl2_action_matrices):
            self.assertEqual(
                _apply_binary_form_action_matrix(binary_form, action_matrix, context.field.prime),
                _transform_binary_form(binary_form, matrix, context.field.prime),
            )

    def test_branch_factor_transform_cache_matches_direct_transform(self):
        context = EnumerationContext(prime=5, genus=2)
        factor = (1, 0, 1)
        matrix_index = 3

        cached_transform = context._transform_monic_factor(factor, matrix_index)
        direct_transform = _transform_binary_form(factor, context.pgl2[matrix_index], context.field.prime)
        direct_transform = context.polynomial(direct_transform).coefficients
        leading = direct_transform[-1]

        self.assertIsNotNone(cached_transform)
        self.assertEqual(cached_transform, (leading, Polynomial(context.field, direct_transform).monic().coefficients))
        self.assertEqual(context._transform_monic_factor(factor, matrix_index), cached_transform)
        self.assertIn((matrix_index, factor), context.branch_factor_transform_cache)

    def test_enumeration_context_reuses_curve_level_caches(self):
        context = EnumerationContext(prime=5, genus=2)
        f = context.polynomial([1, 1, 0, 0, 0, 1])
        translated_f = context.polynomial([3, 1, 0, 0, 0, 1])

        self.assertEqual(context.l_polynomial_coefficients_mod_p(f), [0, 0])
        self.assertEqual(context.l_polynomial_coefficients_mod_p(translated_f), [0, 0])
        self.assertEqual(len(context.l_polynomial_mod_p_cache), 1)

        self.assertEqual(context.l_polynomial_coefficients(f), [0, 10])
        self.assertEqual(context.l_polynomial_coefficients(translated_f), [0, 10])
        self.assertEqual(len(context.exact_l_polynomial_cache), 1)

    def test_enumeration_context_rejects_incompatible_polynomials(self):
        context = EnumerationContext(prime=5, genus=2)
        with self.assertRaises(ValueError):
            context.canonical_key(Polynomial(PrimeField(7), [1, 1, 0, 0, 0, 1]))
        with self.assertRaises(ValueError):
            context.canonical_key([1, 1, 0, 1])

    def test_exhaustive_search_uses_sqlite_orbit_lookup(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)

            result = context.process_polynomial_for_output([1, 1, 0, 0, 0, 1], max_sparsity=0)
            self.assertEqual(result["status"], "sparse")
            self.assertEqual(result["lpoly"], [0, 10])

            translated_result = context.process_polynomial_for_output([3, 1, 0, 0, 0, 1], max_sparsity=0)
            self.assertEqual(translated_result["status"], "duplicate")
            self.assertEqual(translated_result["previous_status"], "sparse")
            self.assertEqual(translated_result["lpoly"], [0, 10])
            self.assertEqual(translated_result["canonical_key"], result["canonical_key"])

            connection = sqlite3.connect(db_path)
            self.assertGreater(connection.execute("SELECT COUNT(*) FROM orbit_cache").fetchone()[0], 0)
            orbit_cache_columns = {
                row[1]
                for row in connection.execute("PRAGMA table_info(orbit_cache)").fetchall()
            }
            self.assertIn("ground_point_count", orbit_cache_columns)
            self.assertIn("hasse_witt_lpoly_mod_p", orbit_cache_columns)
            self.assertEqual(
                connection.execute("SELECT COUNT(DISTINCT ground_point_count) FROM orbit_cache").fetchone()[0],
                1,
            )
            self.assertEqual(
                connection.execute("SELECT COUNT(DISTINCT hasse_witt_lpoly_mod_p) FROM orbit_cache").fetchone()[0],
                1,
            )
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM curve_cache").fetchone()[0], 1)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM sparse_curves").fetchone()[0], 1)
            stored_orbit_keys = [
                context._unpack_field_tuple_blob(row[0])
                for row in connection.execute("SELECT orbit_key FROM orbit_cache").fetchall()
            ]
            self.assertTrue(all(context._is_enumerated_orbit_key(key) for key in stored_orbit_keys))
            _, full_orbit = context._canonical_key_and_orbit(context.polynomial([1, 1, 0, 0, 0, 1]))
            self.assertLess(len(stored_orbit_keys), len(set(full_orbit)))
            curve_cache_columns = {
                row[1]
                for row in connection.execute("PRAGMA table_info(curve_cache)").fetchall()
            }
            self.assertNotIn("ground_point_count", curve_cache_columns)
            self.assertNotIn("factorization_pattern", curve_cache_columns)
            self.assertNotIn("max_sparsity", curve_cache_columns)
            self.assertEqual(connection.execute("SELECT typeof(coefficients) FROM curve_cache").fetchone()[0], "text")
            self.assertEqual(connection.execute("SELECT typeof(canonical_key) FROM curve_cache").fetchone()[0], "blob")
            self.assertEqual(connection.execute("SELECT typeof(lpoly_mod_p) FROM curve_cache").fetchone()[0], "blob")
            self.assertEqual(connection.execute("SELECT typeof(exact_lpoly) FROM curve_cache").fetchone()[0], "blob")
            self.assertEqual(connection.execute("SELECT typeof(coefficients) FROM sparse_curves").fetchone()[0], "text")
            self.assertEqual(connection.execute("SELECT typeof(branch_factors) FROM sparse_curves").fetchone()[0], "text")
            self.assertEqual(connection.execute("SELECT typeof(branch_factorization_pattern) FROM sparse_curves").fetchone()[0], "text")
            self.assertEqual(connection.execute("SELECT typeof(lpoly) FROM sparse_curves").fetchone()[0], "text")
            self.assertEqual(connection.execute("SELECT typeof(canonical_key) FROM sparse_curves").fetchone()[0], "blob")
            branch_factors, branch_pattern = connection.execute(
                "SELECT branch_factors, branch_factorization_pattern FROM sparse_curves"
            ).fetchone()
            self.assertIsInstance(json.loads(branch_factors), list)
            self.assertIsInstance(json.loads(branch_pattern), list)
            self.assertEqual(len(context.canonical_records), 1)
            record = context.canonical_records[result["canonical_key"]]
            self.assertIn(result["canonical_key"], context.index_by_rational_branch_count[record.rational_branch_count])
            self.assertEqual(record.orbit_size, len(stored_orbit_keys))
            self.assertEqual(sparse_presentations_by_orbit_size(context, max_sparsity=0), len(stored_orbit_keys))
            connection.close()
            context.close_sqlite()

    def test_sqlite_enumeration_output_rejects_hasse_witt_failures_before_canonicalization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)

            result = context.process_polynomial_for_output([1, 0, 1, 0, 0, 1], max_sparsity=0)
            self.assertEqual(result["status"], "rejected_hasse_witt_uncanonicalized")
            self.assertEqual(result["lpoly_mod_p"], [4, 0])
            self.assertEqual(len(context.canonical_records), 0)
            stats = context.process_polynomials_for_output(([1, 0, 1, 0, 0, 1],), max_sparsity=0)
            self.assertEqual(stats, {"processed": 1, "rejected_hasse_witt_uncanonicalized": 1})

            connection = sqlite3.connect(db_path)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM curve_cache").fetchone()[0], 0)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM sparse_curves").fetchone()[0], 0)
            connection.close()
            context.close_sqlite()

    def test_enumeration_can_skip_even_models_with_rational_branch_points(self):
        context = EnumerationContext(
            prime=5,
            genus=1,
            skip_even_models_with_rational_branch_point=True,
        )

        result = context.process_polynomial_for_output([0, 1, 0, 0, 1], max_sparsity=0)

        self.assertEqual(result["status"], "covered_by_odd_model")
        self.assertGreater(result["rational_branch_count"], 0)
        self.assertEqual(len(context.canonical_records), 0)

    def test_sqlite_enumeration_output_records_exact_rejection_sparsity_as_negative_one(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)

            result = context.process_polynomial_for_output([1, 4, 0, 0, 0, 1], max_sparsity=0)
            self.assertEqual(result["status"], "rejected_exact")
            self.assertEqual(result["sparsity"], -1)

            connection = sqlite3.connect(db_path)
            self.assertEqual(connection.execute("SELECT status FROM curve_cache").fetchone()[0], "rejected_exact")
            self.assertEqual(connection.execute("SELECT sparsity FROM curve_cache").fetchone()[0], -1)
            connection.close()
            context.close_sqlite()

    def test_enumeration_output_without_sparsity_bound_computes_full_lpoly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/p5_g2_all.sqlite"
            context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)

            result = context.process_polynomial_for_output([1, 0, 1, 0, 0, 1], max_sparsity=None)
            self.assertEqual(result["status"], "sparse")
            self.assertEqual(result["lpoly"], [-1, 0])

            connection = sqlite3.connect(db_path)
            self.assertEqual(connection.execute("SELECT status FROM curve_cache").fetchone()[0], "sparse")
            self.assertEqual(connection.execute("SELECT lpoly FROM sparse_curves").fetchone()[0], "[-1,0]")
            branch_factors, branch_leading, infinity_branch = connection.execute(
                "SELECT branch_factors, branch_leading_coefficient, branch_infinity_branch FROM sparse_curves"
            ).fetchone()
            self.assertIsInstance(json.loads(branch_factors), list)
            self.assertEqual(branch_leading, 1)
            self.assertEqual(infinity_branch, 1)
            connection.close()
            context.close_sqlite()

            resumed_context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)
            translated_result = resumed_context.process_polynomial_for_output([2, 0, 1, 0, 0, 1], max_sparsity=None)
            self.assertIn(translated_result["status"], {"sparse", "duplicate"})
            resumed_context.close_sqlite()

    def test_hasse_witt_filter_skips_canonicalization_for_mod_p_rejections(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)

            result = context.process_polynomial_for_output([1, 0, 1, 0, 0, 1], max_sparsity=0)
            self.assertEqual(result["status"], "rejected_hasse_witt_uncanonicalized")
            self.assertEqual(result["lpoly_mod_p"], [4, 0])
            self.assertEqual(len(context.canonical_records), 0)

            connection = sqlite3.connect(db_path)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM curve_cache").fetchone()[0], 0)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM sparse_curves").fetchone()[0], 0)
            connection.close()
            context.close_sqlite()

    def test_hasse_witt_filter_uses_sqlite_orbit_lookup_for_survivors(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)

            result = context.process_polynomial_for_output([1, 1, 0, 0, 0, 1], max_sparsity=0)
            self.assertEqual(result["status"], "sparse")

            translated_result = context.process_polynomial_for_output([3, 1, 0, 0, 0, 1], max_sparsity=0)
            self.assertEqual(translated_result["status"], "duplicate")
            self.assertEqual(translated_result["previous_status"], "sparse")
            self.assertEqual(translated_result["canonical_key"], result["canonical_key"])

            connection = sqlite3.connect(db_path)
            self.assertGreater(connection.execute("SELECT COUNT(*) FROM orbit_cache").fetchone()[0], 0)
            self.assertEqual(
                connection.execute("SELECT COUNT(DISTINCT hasse_witt_lpoly_mod_p) FROM orbit_cache").fetchone()[0],
                1,
            )
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM curve_cache").fetchone()[0], 1)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM sparse_curves").fetchone()[0], 1)
            connection.close()
            context.close_sqlite()

    def test_sqlite_resume_loads_canonical_records_into_memory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/p5_g2_s_0.sqlite"
            context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)
            result = context.process_polynomial_for_output([1, 1, 0, 0, 0, 1], max_sparsity=0)
            context.close_sqlite()

            resumed_context = EnumerationContext(prime=5, genus=2, sqlite_path=db_path)
            self.assertEqual(len(resumed_context.canonical_records), 1)
            self.assertIn(result["canonical_key"], resumed_context.canonical_records)
            self.assertIn(result["canonical_key"], resumed_context.seen_keys)

            translated_result = resumed_context.process_polynomial_for_output([3, 1, 0, 0, 0, 1], max_sparsity=0)
            self.assertEqual(translated_result["status"], "duplicate")
            self.assertEqual(translated_result["previous_status"], "sparse")
            self.assertEqual(translated_result["lpoly"], [0, 10])

            connection = sqlite3.connect(db_path)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM curve_cache").fetchone()[0], 1)
            self.assertEqual(connection.execute("SELECT COUNT(*) FROM sparse_curves").fetchone()[0], 1)
            connection.close()
            resumed_context.close_sqlite()

    def test_enumeration_context_records_timing_summary(self):
        context = EnumerationContext(prime=5, genus=2)
        context.process_polynomial_for_output([1, 1, 0, 0, 0, 1], max_sparsity=0)
        context.process_polynomial_for_output([3, 1, 0, 0, 0, 1], max_sparsity=0)

        timing = context.timing_summary()
        self.assertEqual(timing["processed_polynomials"], 2)
        self.assertEqual(timing["status_counts"], {"sparse": 1, "duplicate": 1})
        self.assertGreaterEqual(timing["elapsed_seconds"], timing["processing_seconds"])
        self.assertGreater(timing["processing_seconds"], 0.0)
        self.assertEqual(timing["sqlite_load_seconds"], 0.0)
        self.assertEqual(timing["sqlite_write_seconds"], 0.0)

    def test_enumeration_cli_prints_progress(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "data_gen.hyperelliptic",
                    "--p",
                    "5",
                    "--genus",
                    "2",
                    "--max-sparsity",
                    "2",
                    "--limit",
                    "3",
                    "--progress-interval",
                    "1",
                    "--out",
                    db_path,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            connection = sqlite3.connect(db_path)
            summary = connection.execute(
                """
                SELECT
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
                    canonicalized_isomorphism_classes,
                    status_counts
                FROM enumeration_summary
                """
            ).fetchone()
            progress_rows = connection.execute(
                """
                SELECT
                    position,
                    processed,
                    skipped,
                    sparse_isomorphism_classes,
                    canonicalized_isomorphism_classes,
                    delta_sparse_isomorphism_classes,
                    delta_canonicalized_isomorphism_classes
                FROM enumeration_progress
                ORDER BY position
                """
            ).fetchall()
            connection.close()

            self.assertEqual(
                summary[:14],
                (
                    5,
                    2,
                    2,
                    1,
                    "both",
                    "enumerate",
                    "odd:monic;even:monic-and-smallest-nonsquare",
                    3,
                    3,
                    3,
                    0,
                    3,
                    121,
                    3,
                ),
            )
            self.assertEqual(json.loads(summary[14]), {"sparse": 3})
            self.assertEqual(progress_rows, [(0, 0, 0, 0, 0, 0, 0), (1, 1, 0, 1, 1, 1, 1), (2, 2, 0, 2, 2, 1, 1), (3, 3, 0, 3, 3, 1, 1)])

            self.assertIn("prime: 5\ngenus: 2\nprogress: 0/3\nsparse_presentations: 0", completed.stdout)
            self.assertIn("progress: 3/3", completed.stdout)
            self.assertNotIn("\nskipped:", completed.stdout)
            self.assertIn("sparse_isomorphism_classes:", completed.stdout)
            self.assertIn("canonicalized_isomorphism_classes:", completed.stdout)
            self.assertNotIn("total_isomorphism_classes:", completed.stdout)
            self.assertIn("\n-\nprime:", completed.stdout)

    def test_enumeration_cli_enumerate_mode_records_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "data_gen.hyperelliptic",
                    "--p",
                    "3",
                    "--genus",
                    "2",
                    "--max-sparsity",
                    "1",
                    "--limit",
                    "5",
                    "--progress-interval",
                    "0",
                    "--enumeration-mode",
                    "enumerate",
                    "--out",
                    db_path,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            connection = sqlite3.connect(db_path)
            enumeration_mode, processed, final_position = connection.execute(
                """
                SELECT enumeration_mode, processed, final_position
                FROM enumeration_summary
                """
            ).fetchone()
            connection.close()

        self.assertEqual(enumeration_mode, "enumerate")
        self.assertEqual(processed, 5)
        self.assertEqual(final_position, 5)

    def test_enumeration_cli_batch_run_writes_one_database_per_genus(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "data_gen.hyperelliptic",
                    "--p",
                    "3",
                    "--genus-start",
                    "1",
                    "--genus-end",
                    "2",
                    "--max-sparsity",
                    "1",
                    "--limit",
                    "1",
                    "--progress-interval",
                    "0",
                    "--out-dir",
                    tmpdir,
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            for genus in (1, 2):
                db_path = Path(tmpdir) / f"p3_g{genus}_s_1.sqlite"
                self.assertTrue(db_path.exists())
                connection = sqlite3.connect(db_path)
                prime, recorded_genus, processed = connection.execute(
                    """
                    SELECT prime, genus, processed
                    FROM enumeration_summary
                    """
                ).fetchone()
                connection.close()
                self.assertEqual((prime, recorded_genus, processed), (3, genus, 1))

        self.assertIn("batch: 1/2 p=3 genus=1", completed.stdout)
        self.assertIn("batch: 2/2 p=3 genus=2", completed.stdout)

    def test_enumeration_cli_branch_random_records_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "data_gen.hyperelliptic",
                    "--p",
                    "3",
                    "--genus",
                    "2",
                    "--max-sparsity",
                    "1",
                    "--enumeration-mode",
                    "random",
                    "--random-steps",
                    "4",
                    "--random-seed",
                    "11",
                    "--random-max-factors",
                    "3",
                    "--progress-interval",
                    "2",
                    "--out",
                    db_path,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            connection = sqlite3.connect(db_path)
            enumeration_mode, processed, total = connection.execute(
                """
                SELECT enumeration_mode, processed, total_coefficient_vectors
                FROM enumeration_summary
                """
            ).fetchone()
            sparse_count, factorized_sparse_count = connection.execute(
                """
                SELECT count(*), count(branch_factors)
                FROM sparse_curves
                """
            ).fetchone()
            connection.close()

        self.assertEqual(enumeration_mode, "random")
        self.assertEqual(processed, 4)
        self.assertEqual(total, 4)
        self.assertEqual(factorized_sparse_count, sparse_count)
        self.assertIn("progress: 4/4", completed.stdout)

    def test_random_mode_does_not_store_full_orbit_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "data_gen.hyperelliptic",
                    "--p",
                    "3",
                    "--genus",
                    "2",
                    "--enumeration-mode",
                    "random",
                    "--random-steps",
                    "1",
                    "--random-seed",
                    "11",
                    "--progress-interval",
                    "0",
                    "--out",
                    db_path,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            connection = sqlite3.connect(db_path)
            orbit_rows = connection.execute("SELECT count(*) FROM orbit_cache").fetchone()[0]
            branch_orbit_rows = connection.execute("SELECT count(*) FROM branch_orbit_cache").fetchone()[0]
            connection.close()

        self.assertLessEqual(orbit_rows, 1)
        self.assertLessEqual(branch_orbit_rows, 1)

    def test_unbounded_branch_random_writes_partial_summary_when_interrupted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            args = argparse.Namespace(
                p=3,
                genus=2,
                max_sparsity=1,
                out=Path(db_path),
                limit=None,
                progress_interval=0,
                enumeration_mode="random",
                random_steps=None,
                random_seed=11,
                random_max_factors=3,
            )
            calls = {"count": 0}
            original_process = EnumerationContext.process_branch_divisor_for_output

            def interrupt_after_one(context, candidate, max_sparsity):
                calls["count"] += 1
                if calls["count"] == 1:
                    return original_process(context, candidate, max_sparsity)
                raise KeyboardInterrupt

            with mock.patch.object(EnumerationContext, "process_branch_divisor_for_output", autospec=True, side_effect=interrupt_after_one):
                with contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaises(KeyboardInterrupt):
                        run_enumeration_from_args(args)

            connection = sqlite3.connect(db_path)
            processed, total, final_position = connection.execute(
                """
                SELECT processed, total_coefficient_vectors, final_position
                FROM enumeration_summary
                """
            ).fetchone()
            connection.close()

        self.assertEqual(processed, 1)
        self.assertEqual(total, -1)
        self.assertEqual(final_position, 1)

    def test_enumeration_cli_allows_omitted_max_sparsity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "data_gen.hyperelliptic",
                    "--p",
                    "5",
                    "--genus",
                    "2",
                    "--limit",
                    "3",
                    "--progress-interval",
                    "0",
                    "--out",
                    db_path,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            connection = sqlite3.connect(db_path)
            max_sparsity = connection.execute("SELECT max_sparsity FROM enumeration_summary").fetchone()[0]
            connection.close()

        self.assertIsNone(max_sparsity)

    def test_interrupted_enumeration_writes_partial_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = f"{tmpdir}/curves.sqlite"
            args = argparse.Namespace(
                p=5,
                genus=2,
                max_sparsity=2,
                out=Path(db_path),
                limit=10,
                progress_interval=0,
                enumeration_mode="enumerate",
                random_steps=None,
                random_seed=None,
                random_max_factors=3,
            )
            calls = {"count": 0}
            original_process = EnumerationContext.process_branch_divisor_for_output

            def interrupt_after_one(context, candidate, max_sparsity):
                calls["count"] += 1
                if calls["count"] == 1:
                    return original_process(context, candidate, max_sparsity)
                raise KeyboardInterrupt

            with mock.patch.object(EnumerationContext, "process_branch_divisor_for_output", autospec=True, side_effect=interrupt_after_one):
                with contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaises(KeyboardInterrupt):
                        run_enumeration_from_args(args)

            connection = sqlite3.connect(db_path)
            processed, final_position, status_counts = connection.execute(
                """
                SELECT processed, final_position, status_counts
                FROM enumeration_summary
                """
            ).fetchone()
            connection.close()

        self.assertEqual(processed, 1)
        self.assertEqual(final_position, 1)
        self.assertEqual(sum(json.loads(status_counts).values()), 1)


if __name__ == "__main__":
    unittest.main()
