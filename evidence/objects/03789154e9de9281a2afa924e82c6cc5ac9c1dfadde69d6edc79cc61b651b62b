#!/usr/bin/env python3
"""Exact paired exclusion for nontrivial-R83 masks 6 and 22.

The registered v47 paired certificate is reused only for its exact bounded-
matrix enumeration kernel and candidate reflection.  This certificate
reconstructs the two new canonical points from the v45 ambient packet,
derives and freezes their point-specific bounds, exhausts the determinant-one
search once, proves the reflected-set bijection, and evaluates the normalized
quadrilateral separately at both centers.

The determinant-minus-one reflection lies in ``GL_2(O_L)`` but not in
``GL_2^+(O_L)`` or ``SL_2(O_L)``.  It is only a search bijection and does not
identify the two Hilbert-modular points.  All decisive arithmetic is exact.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from math import floor

import sympy as sp

import d7_ambient_points_v45 as ambient
import d7_r83_masks14_30_v47 as pair_core
import residual_orbit_k17_k38_v42 as core


MASKS = (6, 22)
ROOT_PARTNERS = {6: 54, 22: 38}
INDEX = 2
LABEL = "R83:nontrivial"
EXPONENTS = (0, 1, 1, 0, 0)
ALL_POSITIVE = (1,) * core.DEGREE

EXPECTED_CENTER = (
    Q(37, 83),
    Q(-4575, 166),
    Q(-186, 83),
    Q(2830, 83),
    Q(-213, 166),
    Q(-668, 83),
)
EXPECTED_Y2 = (
    Q(2775, 27556),
    Q(10083, 27556),
    Q(1866, 6889),
    Q(-7253, 13778),
    Q(-534, 6889),
    Q(3269, 27556),
)
EXPECTED_MODULUS2 = (
    Q(308063, 6889),
    Q(993584, 6889),
    Q(-1844296, 6889),
    Q(-4470598, 6889),
    Q(1260238, 6889),
    Q(1550051, 6889),
)
EXPECTED_BOUNDS = (
    (Q(341, 100), Q(3, 4), Q(167, 100), Q(37, 25), Q(191, 100), Q(111, 100)),
    (Q(869, 100), Q(31, 20), Q(3, 10), Q(179, 100), Q(79, 25), Q(53, 25)),
    (Q(1873, 100), Q(2017, 100), Q(127, 100), Q(97, 20), Q(209, 50), Q(648, 25)),
    (Q(117363, 100), Q(4651, 100), Q(589, 100), Q(1883, 50), Q(10493, 100), Q(3291, 25)),
    (Q(367, 50), Q(967, 100), Q(717, 100), Q(199, 50), Q(253, 100), Q(339, 25)),
    (Q(22987, 50), Q(223, 10), Q(133, 4), Q(1547, 50), Q(127, 2), Q(1722, 25)),
)
EXPECTED_HISTOGRAMS = {
    6: {
        (0, 0): 1832,
        (0, 1): 5150,
        (0, 2): 1046,
        (0, 3): 350,
        (0, 4): 26,
        (0, 5): 6,
        (1, 0): 78,
        (1, 1): 24,
        (1, 3): 2,
    },
    22: {
        (0, 0): 1852,
        (0, 1): 5108,
        (0, 2): 1014,
        (0, 3): 420,
        (0, 4): 18,
        (0, 5): 10,
        (1, 0): 70,
        (1, 1): 16,
        (1, 2): 4,
        (1, 4): 2,
    },
}
EXPECTED_DIGESTS = {
    6: "756c33defcade60c47adab07c30af7b5beb7b110e2051e8dfe025055690cb523",
    22: "5f5c2adaa1366bceab1d2b43a49fffc48c8377320e9d4d522853e1586950e04a",
}
EXPECTED_AGGREGATE = "c6c2b774f58858e6f22b005eb6134de2c449c0d4eacb38a1b96a5fc745023cc1"


def selected_pair():
    """Reconstruct and exactly reduce the two canonical v45 representatives."""
    base = {
        label: matrix
        for label, matrix, _ in ambient.module_templates(INDEX)
    }[LABEL]
    assert ambient.reduce_l(
        ambient.SQUARE_CLASS_UNITS[22] + ambient.SQUARE_CLASS_UNITS[6]
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
        oriented[22] - reflection * oriented[6] * reflection
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
        inverted[22] - reflection * inverted[6] * reflection
    ) == sp.zeros(2)

    generators = ambient.XI_UNITS[:4] + (ambient.POSITIVE_NONSQUARE,)
    unit = sp.Integer(1)
    for generator, exponent in zip(generators, EXPONENTS):
        unit = ambient.reduce_l(unit * generator**exponent)
    assert unit == ambient.reduce_l(ambient.XI_UNITS[1] * ambient.XI_UNITS[2])
    assert abs(ambient.absolute_norm(unit)) == 1
    unit_inverse = ambient.reduce_field(1 / unit)
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
        balanced[22] - reflection * balanced[6] * reflection
    ) == sp.zeros(2)

    invariants = {}
    for mask, matrix in balanced.items():
        center, y2, modulus2 = pair_core.v46.fixed_invariants(matrix)
        invariants[mask] = tuple(
            pair_core.v46.to_core(value) for value in (center, y2, modulus2)
        )
    assert invariants[6] == (EXPECTED_CENTER, EXPECTED_Y2, EXPECTED_MODULUS2)
    assert invariants[22] == (
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
    """Return the next hundredth strictly outside a rational upper endpoint."""
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(y2, modulus2):
    """Derive all six point-specific outward bounds with exact intervals."""
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
    """Enumerate the exact integral A/C boxes through the positive trace form."""
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
    assert (len(a_values), len(c_values)) == (23937, 6605)
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


def validate_reflection_bijection(center6, center22, y2, candidates6, candidates22):
    """Certify the exact involution, bounds identities, and complete converse."""
    assert center22 == core.neg(center6)
    assert len(candidates6) == len(candidates22) == 8514
    assert len(set(candidates6)) == len(candidates6)
    assert len(set(candidates22)) == len(candidates22)
    for gamma6, gamma22 in zip(candidates6, candidates22):
        assert pair_core.reflect_raw(gamma6) == gamma22
        assert pair_core.reflect_raw(gamma22) == gamma6
        a, b, c, d = gamma6
        ar, br, cr, dr = gamma22
        aq, bq, cq, dq = map(core.as_field, gamma6)
        arq, brq, crq, drq = map(core.as_field, gamma22)
        assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
        assert core.sub(core.mul(arq, drq), core.mul(brq, crq)) == core.ONE
        assert exact_square(c, d, center6, y2) == exact_square(
            cr, dr, center22, y2
        )
        assert exact_square(a, b, center6, y2) == exact_square(
            ar, br, center22, y2
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
    assert len(exclusions) == len(candidates) == 8514
    assert histogram == EXPECTED_HISTOGRAMS[mask]
    assert sum(histogram.values()) == 8514
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
    print("EXACT NONTRIVIAL-R83 MASKS {6,22} ORBIT CERTIFICATE v49")
    assert ambient.PAIR_XOR[LABEL] == 48
    assert all((LABEL, mask) in ambient.POINT_RECORDS for mask in MASKS)
    assert ROOT_PARTNERS == {mask: mask ^ 48 for mask in MASKS}
    assert {6, 54}.isdisjoint({22, 38})
    prior_fibers = ({14, 62}, {30, 46}, {15, 63}, {31, 47})
    assert all(
        {6, 54}.isdisjoint(fiber) and {22, 38}.isdisjoint(fiber)
        for fiber in prior_fibers
    )

    _, invariants = selected_pair()
    center6, y2, modulus2 = invariants[6]
    center22, y2_22, modulus2_22 = invariants[22]
    assert y2_22 == y2 and modulus2_22 == modulus2
    bounds = derive_bounds(y2, modulus2)
    pair_core.install_point(y2, bounds)
    a_values, c_values = entry_lists(bounds)
    print("exact shared a/c counts:", len(a_values), len(c_values), flush=True)

    candidates6 = pair_core.bounded_candidates(
        center6,
        a_values,
        c_values,
        expected_counts=(2019, 8514, 0),
        progress_label="masks6/22",
    )
    candidates22 = [pair_core.reflect(gamma) for gamma in candidates6]
    assert all(gamma[2] != core.ZERO for gamma in candidates6)
    assert all(gamma[2] != core.ZERO for gamma in candidates22)
    validate_reflection_bijection(
        center6, center22, y2, candidates6, candidates22
    )
    print("exact reflected candidate bijection: 8514/8514", flush=True)

    result6 = classify(6, center6, candidates6)
    result22 = classify(22, center22, candidates22)
    for mask, result in ((6, result6), (22, result22)):
        print(f"mask={mask}: bottom=2019; c!=0=8514; c=0=0; survivors=0")
        print(f"mask={mask}: histogram={result[0]}")
        print(f"mask={mask}: exclusion digest={result[1]}")

    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode(
            "ascii"
        )
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: exact exhaustive exclusion for masks 6 and 22")


if __name__ == "__main__":
    main()
