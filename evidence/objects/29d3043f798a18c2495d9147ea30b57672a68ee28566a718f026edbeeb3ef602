#!/usr/bin/env python3
"""Exact paired exclusion for nontrivial-R83 masks 1 and 17.

The registered v47 certificate is reused only for its exact bounded-matrix
enumeration kernel and candidate reflection.  This certificate reconstructs
the new canonical points from the v45 ambient packet, applies an explicit
determinant-one inversion and saturated-unit balancing, derives and freezes
all point-specific bounds, exhausts the determinant-one search including the
``c=0`` branch, and tests the normalized quadrilateral separately at the two
centers.

The determinant-minus-one reflection lies in ``GL_2(O_L)`` but not in
``GL_2^+(O_L)`` or ``SL_2(O_L)``.  It supplies only an exact enumeration
bijection; it does not identify the two Hilbert-modular points.  All decisive
arithmetic and polygon predicates are exact.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from math import floor

import sympy as sp

import d7_ambient_points_v45 as ambient
import d7_r83_masks14_30_v47 as pair_core
import residual_orbit_k17_k38_v42 as core


MASKS = (1, 17)
ROOT_PARTNERS = {1: 49, 17: 33}
INDEX = 2
LABEL = "R83:nontrivial"
EXPONENTS = (1, 0, 0, -1, -1)
ALL_POSITIVE = (1,) * core.DEGREE

EXPECTED_CENTER = (
    Q(-603, 83),
    Q(575, 83),
    Q(4517, 166),
    Q(408, 83),
    Q(-742, 83),
    Q(-204, 83),
)
EXPECTED_Y2 = (
    Q(1497, 6889),
    Q(2438, 6889),
    Q(124, 6889),
    Q(-780, 6889),
    Q(-297, 27556),
    Q(58, 6889),
)
EXPECTED_MODULUS2 = (
    Q(231808, 6889),
    Q(-1674196, 6889),
    Q(-3166462, 6889),
    Q(327307, 6889),
    Q(1653404, 6889),
    Q(417388, 6889),
)
EXPECTED_BOUNDS = (
    (Q(67, 20), Q(217, 100), Q(221, 100), Q(18, 25), Q(67, 100), Q(43, 25)),
    (Q(171, 20), Q(113, 25), Q(39, 100), Q(22, 25), Q(11, 10), Q(329, 100)),
    (Q(476, 25), Q(693, 100), Q(24, 25), Q(199, 20), Q(1197, 100), Q(1671, 100)),
    (Q(115487, 100), Q(3388, 25), Q(39, 5), Q(917, 50), Q(3659, 100), Q(20423, 100)),
    (Q(373, 50), Q(83, 25), Q(271, 50), Q(409, 50), Q(181, 25), Q(437, 50)),
    (Q(45239, 100), Q(6497, 100), Q(2201, 50), Q(1507, 100), Q(443, 20), Q(5343, 50)),
)
EXPECTED_HISTOGRAMS = {
    1: {
        (0, 0): 1870,
        (0, 1): 5372,
        (0, 2): 1044,
        (0, 3): 344,
        (0, 4): 28,
        (0, 5): 12,
        (1, 0): 50,
        (1, 1): 16,
        (1, 2): 6,
        (1, 3): 4,
        (1, 4): 2,
    },
    17: {
        (0, 0): 1872,
        (0, 1): 5418,
        (0, 2): 1008,
        (0, 3): 340,
        (0, 4): 16,
        (0, 5): 14,
        (1, 0): 48,
        (1, 1): 18,
        (1, 2): 6,
        (1, 3): 6,
        (1, 4): 2,
    },
}
EXPECTED_DIGESTS = {
    1: "3d6a732cc36a03fe90933cdfcffd1b9d977e7ce22cb1d7ae1f3926e6dff2e1a8",
    17: "82b2d957ba62af939e5bb39b0f6bbec06b230094aeee60faa3ae54caf3e2e17b",
}
EXPECTED_AGGREGATE = "9a33ab59a829397a361196b67c9f6199be8288363a380fe506225d29496984dd"


def field_power(value, exponent: int):
    """Compute an exact power in L, including negative unit exponents."""
    if exponent < 0:
        return field_power(ambient.reduce_field(1 / value), -exponent)
    answer = sp.Integer(1)
    base = value
    while exponent:
        if exponent & 1:
            answer = ambient.reduce_l(answer * base)
        base = ambient.reduce_l(base * base)
        exponent //= 2
    return answer


def selected_pair():
    """Reconstruct and exactly reduce the two canonical v45 representatives."""
    base = {
        label: matrix
        for label, matrix, _ in ambient.module_templates(INDEX)
    }[LABEL]
    assert ambient.reduce_l(
        ambient.SQUARE_CLASS_UNITS[17] + ambient.SQUARE_CLASS_UNITS[1]
    ) == 0
    oriented = {
        mask: ambient.conjugate_by_unit(base, ambient.SQUARE_CLASS_UNITS[mask])
        for mask in MASKS
    }
    for matrix in oriented.values():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)

    reflection = sp.diag(1, -1)
    assert reflection.det() == -1
    assert ambient.sign_vector(sp.Integer(-1)) == (-1,) * core.DEGREE
    ambient.assert_integral_matrix(reflection)
    assert reflection * reflection == sp.eye(2)
    assert ambient.matrix_reduce(
        oriented[17] - reflection * oriented[1] * reflection
    ) == sp.zeros(2)

    inversion = sp.Matrix([[0, -1], [1, 0]])
    inversion_inverse = sp.Matrix([[0, 1], [-1, 0]])
    assert inversion.det() == 1
    assert inversion * inversion_inverse == sp.eye(2)
    ambient.assert_integral_matrix(inversion)
    ambient.assert_integral_matrix(inversion_inverse)
    inverted = {
        mask: ambient.matrix_multiply(
            ambient.matrix_multiply(inversion, oriented[mask]), inversion_inverse
        )
        for mask in MASKS
    }
    for matrix in inverted.values():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)
    assert ambient.matrix_reduce(
        inverted[17] - reflection * inverted[1] * reflection
    ) == sp.zeros(2)

    generators = ambient.XI_UNITS[:4] + (ambient.POSITIVE_NONSQUARE,)
    unit = sp.Integer(1)
    for generator, exponent in zip(generators, EXPONENTS):
        unit = ambient.reduce_l(unit * field_power(generator, exponent))
    assert abs(ambient.absolute_norm(unit)) == 1
    unit_inverse = ambient.reduce_field(1 / unit)
    assert ambient.reduce_l(unit * unit_inverse - 1) == 0
    conjugator = sp.diag(unit, unit_inverse)
    conjugator_inverse = sp.diag(unit_inverse, unit)
    assert ambient.reduce_l(conjugator.det() - 1) == 0
    assert ambient.reduce_l(conjugator_inverse.det() - 1) == 0
    ambient.assert_integral_matrix(conjugator)
    ambient.assert_integral_matrix(conjugator_inverse)
    balanced = {
        mask: ambient.matrix_multiply(
            ambient.matrix_multiply(conjugator, inverted[mask]), conjugator_inverse
        )
        for mask in MASKS
    }
    for matrix in balanced.values():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)
    assert ambient.matrix_reduce(
        balanced[17] - reflection * balanced[1] * reflection
    ) == sp.zeros(2)

    invariants = {}
    for mask, matrix in balanced.items():
        center, y2, modulus2 = pair_core.v46.fixed_invariants(matrix)
        invariants[mask] = tuple(
            pair_core.v46.to_core(value) for value in (center, y2, modulus2)
        )
    assert invariants[1] == (EXPECTED_CENTER, EXPECTED_Y2, EXPECTED_MODULUS2)
    assert invariants[17] == (
        core.neg(EXPECTED_CENTER), EXPECTED_Y2, EXPECTED_MODULUS2
    )
    assert core.sign_vector(EXPECTED_Y2) == ALL_POSITIVE
    assert core.sign_vector(EXPECTED_MODULUS2) == ALL_POSITIVE
    assert core.norm(EXPECTED_Y2) == Q(1, 28217344)
    assert core.norm(EXPECTED_MODULUS2) == Q(2097130411225081, 6889)
    assert all(
        abs(ambient.absolute_norm(matrix[1, 0])) == 6889
        for matrix in balanced.values()
    )
    return balanced, invariants


def outward_hundredth(value: Q) -> Q:
    """Return the next hundredth strictly outside an exact upper endpoint."""
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(y2, modulus2):
    """Derive all six point-specific outward bounds with rational intervals."""
    bounds = [[] for _ in range(6)]
    for index in range(core.DEGREE):
        y_lo, y_hi = core.sqrt_interval(
            core.interval_eval(y2, core.ROOT_INTERVALS[index])
        )
        modulus_lo, modulus_hi = core.sqrt_interval(
            core.interval_eval(modulus2, core.ROOT_INTERVALS[index])
        )
        assert y_lo > 0 and modulus_lo > 0
        denominator = core.sqrt_interval(
            (y_hi / core.M_LOWER[index], y_hi / core.M_LOWER[index])
        )[1]
        c_entry = core.sqrt_interval(
            (
                1 / (core.M_LOWER[index] * y_lo),
                1 / (core.M_LOWER[index] * y_lo),
            )
        )[1]
        numerator = core.U_UPPER[index] * denominator
        a_entry = numerator / y_lo
        d_entry = denominator + c_entry * modulus_hi
        b_entry = numerator + a_entry * modulus_hi
        for target, value in zip(
            bounds, (denominator, numerator, a_entry, b_entry, c_entry, d_entry)
        ):
            target.append(outward_hundredth(value))
    answer = tuple(tuple(values) for values in bounds)
    assert answer == EXPECTED_BOUNDS
    return answer


def entry_lists(bounds):
    """Enumerate the exact integral A/C boxes using the positive trace form."""
    a_values = core.enumerate_trace_ellipsoid(
        core.ZERO,
        sum((value * value for value in bounds[2]), Q(0)),
        lambda value: core.within_integral_box(value, bounds[2]),
    )
    c_values = core.enumerate_trace_ellipsoid(
        core.ZERO,
        sum((value * value for value in bounds[4]), Q(0)),
        lambda value: core.within_integral_box(value, bounds[4]),
    )
    assert (len(a_values), len(c_values)) == (23959, 6601)
    assert {core.neg(value) for value in a_values} == set(a_values)
    assert {core.neg(value) for value in c_values} == set(c_values)
    return a_values, c_values


def exact_square(coefficient, constant, x, y2):
    """Return ``(coefficient*x+constant)^2 + coefficient^2*y2`` in L."""
    coefficient_q = core.as_field(coefficient)
    constant_q = core.as_field(constant)
    real = core.add(core.mul(coefficient_q, x), constant_q)
    return core.add(
        core.mul(real, real),
        core.mul(core.mul(coefficient_q, coefficient_q), y2),
    )


def validate_reflection_bijection(center1, center17, y2, candidates1, candidates17):
    """Certify the exact involution, bounds identities, and complete converse."""
    assert center17 == core.neg(center1)
    assert len(candidates1) == len(candidates17) == 8748
    assert len(set(candidates1)) == len(candidates1)
    assert len(set(candidates17)) == len(candidates17)
    for gamma1, gamma17 in zip(candidates1, candidates17):
        assert pair_core.reflect_raw(gamma1) == gamma17
        assert pair_core.reflect_raw(gamma17) == gamma1
        a, b, c, d = gamma1
        ar, br, cr, dr = gamma17
        aq, bq, cq, dq = map(core.as_field, gamma1)
        arq, brq, crq, drq = map(core.as_field, gamma17)
        assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
        assert core.sub(core.mul(arq, drq), core.mul(brq, crq)) == core.ONE
        assert exact_square(c, d, center1, y2) == exact_square(
            cr, dr, center17, y2
        )
        assert exact_square(a, b, center1, y2) == exact_square(
            ar, br, center17, y2
        )
        assert (c == core.ZERO) == (cr == core.ZERO)


def classify(mask, x, candidates):
    """Run the exact polygon predicates at one center and freeze the result."""
    survivors = []
    exclusions = []
    histogram = {}
    for gamma in candidates:
        inside, failure = core.in_all_coordinate_polygons(gamma, x)
        if inside:
            survivors.append(gamma)
        else:
            assert failure is not None
            histogram[failure] = histogram.get(failure, 0) + 1
            exclusions.append((gamma, failure))
    assert survivors == []
    assert len(exclusions) == len(candidates) == 8748
    assert histogram == EXPECTED_HISTOGRAMS[mask]
    assert sum(histogram.values()) == 8748
    digest = sha256()
    for gamma, failure in sorted(exclusions):
        flat = tuple(coefficient for entry in gamma for coefficient in entry)
        digest.update(
            (",".join(map(str, flat)) + f";{failure[0]},{failure[1]}\n").encode(
                "ascii"
            )
        )
    assert digest.hexdigest() == EXPECTED_DIGESTS[mask]
    return histogram, digest.hexdigest()


def main() -> None:
    print("EXACT NONTRIVIAL-R83 MASKS {1,17} ORBIT CERTIFICATE v50")
    assert ambient.PAIR_XOR[LABEL] == 48
    assert all((LABEL, mask) in ambient.POINT_RECORDS for mask in MASKS)
    assert ROOT_PARTNERS == {mask: mask ^ 48 for mask in MASKS}
    assert {1, 49}.isdisjoint({17, 33})
    prior_fibers = (
        {14, 62}, {30, 46}, {15, 63}, {31, 47}, {6, 54}, {22, 38}
    )
    assert all(
        {1, 49}.isdisjoint(fiber) and {17, 33}.isdisjoint(fiber)
        for fiber in prior_fibers
    )

    _, invariants = selected_pair()
    center1, y2, modulus2 = invariants[1]
    center17, y2_17, modulus2_17 = invariants[17]
    assert y2_17 == y2 and modulus2_17 == modulus2
    bounds = derive_bounds(y2, modulus2)
    pair_core.install_point(y2, bounds)
    a_values, c_values = entry_lists(bounds)
    print("exact shared a/c counts:", len(a_values), len(c_values), flush=True)

    candidates1 = pair_core.bounded_candidates(
        center1,
        a_values,
        c_values,
        expected_counts=(2017, 8748, 0),
        progress_label="masks1/17",
    )
    candidates17 = [pair_core.reflect(gamma) for gamma in candidates1]
    assert all(gamma[2] != core.ZERO for gamma in candidates1)
    assert all(gamma[2] != core.ZERO for gamma in candidates17)
    validate_reflection_bijection(
        center1, center17, y2, candidates1, candidates17
    )
    print("exact reflected candidate bijection: 8748/8748", flush=True)

    result1 = classify(1, center1, candidates1)
    result17 = classify(17, center17, candidates17)
    for mask, result in ((1, result1), (17, result17)):
        print(f"mask={mask}: bottom=2017; c!=0=8748; c=0=0; survivors=0")
        print(f"mask={mask}: histogram={result[0]}")
        print(f"mask={mask}: exclusion digest={result[1]}")

    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode(
            "ascii"
        )
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: exact exhaustive exclusion for masks 1 and 17")


if __name__ == "__main__":
    main()
