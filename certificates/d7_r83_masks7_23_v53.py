#!/usr/bin/env python3
"""Exact paired exclusion for nontrivial-R83 masks 7 and 23.

The registered v47 certificate is reused only for its exact bounded-matrix
enumeration kernel and candidate reflection.  This certificate reconstructs
the two canonical points from the v45 ambient packet, certifies their free
root fibers, applies an explicit determinant-one inversion and saturated-unit
balancing, derives and freezes all point-specific bounds, exhausts the
determinant-one search including the separate ``c=0`` branch, and tests the
normalized quadrilateral separately at the two centers.

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


MASKS = (7, 23)
ROOT_PARTNERS = {7: 55, 23: 39}
INDEX = 2
LABEL = "R83:nontrivial"
EXPONENTS = (0, 0, 0, -1, -2)
ALL_POSITIVE = (1,) * core.DEGREE

EXPECTED_BALANCING_UNIT = (
    -2 * ambient.T - 4 * ambient.T**2 - ambient.T**3
    + 2 * ambient.T**4 + ambient.T**5
)
EXPECTED_BALANCING_UNIT_INVERSE = (
    1 - 25 * ambient.T + 18 * ambient.T**3 - 3 * ambient.T**5
)
EXPECTED_CENTER = (
    Q(269, 83),
    Q(677, 83),
    Q(-2106, 83),
    Q(307, 166),
    Q(501, 83),
    Q(-195, 166),
)
EXPECTED_Y2 = (
    Q(-1037, 13778),
    Q(-11997, 27556),
    Q(19989, 27556),
    Q(487, 27556),
    Q(-5233, 27556),
    Q(397, 13778),
)
EXPECTED_MODULUS2 = (
    Q(290921, 6889),
    Q(2520855, 6889),
    Q(3953621, 6889),
    Q(-3802545, 6889),
    Q(-1079522, 6889),
    Q(850036, 6889),
)
EXPECTED_BOUNDS = (
    (Q(32, 25), Q(46, 25), Q(24, 25), Q(21, 10), Q(103, 50), Q(34, 25)),
    (Q(327, 100), Q(383, 100), Q(17, 100), Q(64, 25), Q(341, 100), Q(259, 100)),
    (Q(499, 10), Q(817, 100), Q(221, 100), Q(17, 5), Q(387, 100), Q(1061, 50)),
    (Q(44057, 100), Q(5741, 50), Q(339, 100), Q(5367, 100), Q(11321, 100), Q(16081, 100)),
    (Q(391, 20), Q(98, 25), Q(312, 25), Q(14, 5), Q(117, 50), Q(111, 10)),
    (Q(8629, 50), Q(1101, 20), Q(191, 10), Q(441, 10), Q(6851, 100), Q(4207, 50)),
)
EXPECTED_HISTOGRAMS = {
    7: {
        (0, 0): 1714,
        (0, 1): 5302,
        (0, 2): 942,
        (0, 3): 396,
        (0, 4): 26,
        (0, 5): 20,
        (1, 0): 78,
        (1, 1): 18,
        (1, 2): 8,
        (1, 3): 4,
        (1, 4): 6,
    },
    23: {
        (0, 0): 1680,
        (0, 1): 5402,
        (0, 2): 912,
        (0, 3): 374,
        (0, 4): 20,
        (0, 5): 22,
        (1, 0): 80,
        (1, 1): 20,
        (1, 3): 2,
        (1, 4): 2,
    },
}
EXPECTED_DIGESTS = {
    7: "b653b69bf5888c06df2ee9a20e5f6fd8d1d06277de0b862933ba5a4bf8811ae5",
    23: "392b61bb9ec4dcfd6720de02184032cfe211f9094b3e5ed1dac5224be84f8d52",
}
EXPECTED_AGGREGATE = "aba4002201230fea6eb9e127aa6fa38a04e63194d3970adc32c95a1a87386dab"


def field_power(value, exponent: int):
    """Compute an exact power in L, including negative unit exponents."""
    if exponent < 0:
        inverse = ambient.reduce_field(1 / value)
        ambient.assert_integral(inverse)
        return field_power(inverse, -exponent)
    answer = sp.Integer(1)
    base = value
    while exponent:
        if exponent & 1:
            answer = ambient.reduce_l(answer * base)
        base = ambient.reduce_l(base * base)
        exponent //= 2
    return answer


def certify_root_pair_labels(base, oriented) -> None:
    """Replay determinant-one bridges for the two complete free-root fibers."""
    assert ambient.PAIR_XOR[LABEL] == 48
    assert ROOT_PARTNERS == {mask: mask ^ 48 for mask in MASKS}
    assert {7, 55}.isdisjoint({23, 39})
    anti = ambient.I_83_ANTI
    anti_inverse = ambient.matrix_reduce(anti.inv())
    ambient.assert_integral_matrix(anti)
    assert ambient.reduce_l(anti.det() - (ambient.T**2 - 4)) == 0
    assert ambient.matrix_reduce(
        anti * base - (ambient.alpha_trace(INDEX) * sp.eye(2) - base) * anti
    ) == sp.zeros(2)

    for mask, partner in ROOT_PARTNERS.items():
        source_unit = ambient.SQUARE_CLASS_UNITS[mask]
        target_unit = ambient.SQUARE_CLASS_UNITS[partner]
        source_diagonal = sp.diag(source_unit, 1)
        target_diagonal = sp.diag(target_unit, 1)
        source_diagonal_inverse = ambient.matrix_reduce(source_diagonal.inv())
        bridge = ambient.matrix_multiply(
            ambient.matrix_multiply(target_diagonal, anti_inverse),
            source_diagonal_inverse,
        )
        ambient.assert_integral_matrix(bridge)
        assert ambient.reduce_l(bridge.det() - 1) == 0
        root_source = ambient.matrix_reduce(
            ambient.alpha_trace(INDEX) * sp.eye(2) - oriented[mask]
        )
        target = ambient.conjugate_by_unit(base, target_unit)
        assert ambient.matrix_reduce(
            bridge * root_source - target * bridge
        ) == sp.zeros(2)


def selected_pair():
    """Reconstruct and exactly reduce the two canonical v45 representatives."""
    base = {
        label: matrix
        for label, matrix, _ in ambient.module_templates(INDEX)
    }[LABEL]
    assert ambient.reduce_l(
        ambient.SQUARE_CLASS_UNITS[23] + ambient.SQUARE_CLASS_UNITS[7]
    ) == 0
    oriented = {
        mask: ambient.conjugate_by_unit(base, ambient.SQUARE_CLASS_UNITS[mask])
        for mask in MASKS
    }
    for matrix in oriented.values():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)
    certify_root_pair_labels(base, oriented)

    reflection = sp.diag(1, -1)
    assert reflection.det() == -1
    assert ambient.sign_vector(sp.Integer(-1)) == (-1,) * core.DEGREE
    ambient.assert_integral_matrix(reflection)
    assert reflection * reflection == sp.eye(2)
    assert ambient.matrix_reduce(
        oriented[23] - reflection * oriented[7] * reflection
    ) == sp.zeros(2)

    inversion = sp.Matrix([[0, -1], [1, 0]])
    inversion_inverse = sp.Matrix([[0, 1], [-1, 0]])
    assert inversion.det() == inversion_inverse.det() == 1
    assert inversion * inversion_inverse == sp.eye(2)
    ambient.assert_integral_matrix(inversion)
    ambient.assert_integral_matrix(inversion_inverse)
    raw_invariants = {
        mask: pair_core.v46.fixed_invariants(matrix)
        for mask, matrix in oriented.items()
    }
    raw_modulus_norm = core.norm(
        tuple(pair_core.v46.to_core(value) for value in (raw_invariants[7][2],))[0]
    )
    assert 0 < raw_modulus_norm < 1
    inverted = {
        mask: ambient.matrix_multiply(
            ambient.matrix_multiply(inversion, oriented[mask]), inversion_inverse
        )
        for mask in MASKS
    }
    for mask, matrix in inverted.items():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)
        raw_center, raw_y2, raw_modulus2 = raw_invariants[mask]
        center, y2, modulus2 = pair_core.v46.fixed_invariants(matrix)
        assert ambient.reduce_field(center + raw_center / raw_modulus2) == 0
        assert ambient.reduce_field(y2 - raw_y2 / raw_modulus2**2) == 0
        assert ambient.reduce_field(modulus2 - 1 / raw_modulus2) == 0
    assert ambient.matrix_reduce(
        inverted[23] - reflection * inverted[7] * reflection
    ) == sp.zeros(2)

    generators = ambient.XI_UNITS[:4] + (ambient.POSITIVE_NONSQUARE,)
    unit = sp.Integer(1)
    for generator, exponent in zip(generators, EXPONENTS):
        unit = ambient.reduce_l(unit * field_power(generator, exponent))
    unit_inverse = ambient.reduce_field(1 / unit)
    assert ambient.reduce_l(unit - EXPECTED_BALANCING_UNIT) == 0
    assert ambient.reduce_l(unit_inverse - EXPECTED_BALANCING_UNIT_INVERSE) == 0
    assert abs(ambient.absolute_norm(unit)) == 1
    assert ambient.reduce_l(unit * unit_inverse - 1) == 0
    ambient.assert_integral(unit)
    ambient.assert_integral(unit_inverse)
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
    square = ambient.reduce_l(unit * unit)
    fourth = ambient.reduce_l(square * square)
    invariants = {}
    for mask, matrix in balanced.items():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)
        center, y2, modulus2 = pair_core.v46.fixed_invariants(matrix)
        inverted_center, inverted_y2, inverted_modulus2 = (
            pair_core.v46.fixed_invariants(inverted[mask])
        )
        assert ambient.reduce_l(center - square * inverted_center) == 0
        assert ambient.reduce_l(y2 - fourth * inverted_y2) == 0
        assert ambient.reduce_l(modulus2 - fourth * inverted_modulus2) == 0
        invariants[mask] = tuple(
            pair_core.v46.to_core(value) for value in (center, y2, modulus2)
        )
        assert abs(ambient.absolute_norm(matrix[1, 0])) == 6889
    assert ambient.matrix_reduce(
        balanced[23] - reflection * balanced[7] * reflection
    ) == sp.zeros(2)
    assert invariants[7] == (EXPECTED_CENTER, EXPECTED_Y2, EXPECTED_MODULUS2)
    assert invariants[23] == (
        core.neg(EXPECTED_CENTER), EXPECTED_Y2, EXPECTED_MODULUS2
    )
    assert core.sign_vector(EXPECTED_Y2) == ALL_POSITIVE
    assert core.sign_vector(EXPECTED_MODULUS2) == ALL_POSITIVE
    assert core.norm(EXPECTED_Y2) == Q(1, 28217344)
    assert core.norm(EXPECTED_MODULUS2) == Q(2097130411225081, 6889)
    return balanced, invariants


def outward_hundredth(value: Q) -> Q:
    """Return the next hundredth strictly outside an exact upper endpoint."""
    answer = Q(floor(100 * value) + 1, 100)
    assert answer > value
    return answer


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
    assert (len(a_values), len(c_values)) == (23889, 6609)
    assert {core.neg(value) for value in a_values} == set(a_values)
    assert {core.neg(value) for value in c_values} == set(c_values)
    return a_values, c_values


def certify_c_zero_sector(x, bounds) -> None:
    """Identify the sole c=0 bottom row and prove its determinant obstruction."""
    den_squares = core.square_bounds(bounds[0])
    den_radius = sum(den_squares, Q(0))
    zero_rows = [
        (core.ZERO, d)
        for d in core.enumerate_trace_ellipsoid(
            core.ZERO,
            den_radius,
            lambda d: core.all_le_bounds(
                core.denominator_square(core.ZERO, core.as_field(d), x),
                den_squares,
            ),
        )
    ]
    assert zero_rows == [(core.ZERO, core.ZERO)]
    c, d = zero_rows[0]
    assert c == d == core.ZERO
    # For arbitrary a,b, det((a,b),(0,0))=a*0-b*0=0, never 1.
    assert core.ZERO != core.ONE


def exact_square(coefficient, constant, x, y2):
    """Return ``(coefficient*x+constant)^2 + coefficient^2*y2`` in L."""
    coefficient_q = core.as_field(coefficient)
    constant_q = core.as_field(constant)
    real = core.add(core.mul(coefficient_q, x), constant_q)
    return core.add(
        core.mul(real, real),
        core.mul(core.mul(coefficient_q, coefficient_q), y2),
    )


def validate_reflection_bijection(center7, center23, y2, candidates7, candidates23):
    """Certify the exact involution, bounds identities, and complete converse."""
    assert center23 == core.neg(center7)
    assert len(candidates7) == len(candidates23) == 8514
    assert len(set(candidates7)) == len(candidates7)
    assert len(set(candidates23)) == len(candidates23)
    for gamma7, gamma23 in zip(candidates7, candidates23):
        assert pair_core.reflect_raw(gamma7) == gamma23
        assert pair_core.reflect_raw(gamma23) == gamma7
        a, b, c, d = gamma7
        ar, br, cr, dr = gamma23
        aq, bq, cq, dq = map(core.as_field, gamma7)
        arq, brq, crq, drq = map(core.as_field, gamma23)
        assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
        assert core.sub(core.mul(arq, drq), core.mul(brq, crq)) == core.ONE
        assert exact_square(c, d, center7, y2) == exact_square(
            cr, dr, center23, y2
        )
        assert exact_square(a, b, center7, y2) == exact_square(
            ar, br, center23, y2
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
    print("EXACT NONTRIVIAL-R83 MASKS {7,23} ORBIT CERTIFICATE v53")
    assert ambient.PAIR_XOR[LABEL] == 48
    assert all((LABEL, mask) in ambient.POINT_RECORDS for mask in MASKS)
    assert ROOT_PARTNERS == {mask: mask ^ 48 for mask in MASKS}
    assert {7, 55}.isdisjoint({23, 39})
    prior_fibers = (
        {14, 62}, {30, 46}, {15, 63}, {31, 47}, {6, 54}, {22, 38},
        {1, 49}, {17, 33}, {0, 48}, {16, 32}, {2, 50}, {18, 34},
        {3, 51}, {19, 35}, {4, 52}, {20, 36}, {5, 53}, {21, 37},
        {8, 56}, {24, 40}, {9, 57}, {25, 41}, {10, 58}, {26, 42},
        {11, 59}, {27, 43}, {12, 60}, {28, 44}, {13, 61}, {29, 45},
    )
    assert all(
        {7, 55}.isdisjoint(fiber) and {23, 39}.isdisjoint(fiber)
        for fiber in prior_fibers
    )

    _, invariants = selected_pair()
    center7, y2, modulus2 = invariants[7]
    center23, y2_23, modulus2_23 = invariants[23]
    assert y2_23 == y2 and modulus2_23 == modulus2
    bounds = derive_bounds(y2, modulus2)
    pair_core.install_point(y2, bounds)
    a_values, c_values = entry_lists(bounds)
    print("exact shared a/c counts:", len(a_values), len(c_values), flush=True)

    certify_c_zero_sector(center7, bounds)
    candidates7 = pair_core.bounded_candidates(
        center7,
        a_values,
        c_values,
        expected_counts=(1955, 8514, 0),
        progress_label="masks7/23",
    )
    assert all(gamma[2] != core.ZERO for gamma in candidates7)
    # The complete 1,955-row search is exactly 1,954 nonzero-c rows plus the
    # separately certified row (c,d)=(0,0), which cannot have determinant one.
    assert 1955 == 1954 + 1
    candidates23 = [pair_core.reflect(gamma) for gamma in candidates7]
    assert all(gamma[2] != core.ZERO for gamma in candidates23)
    validate_reflection_bijection(
        center7, center23, y2, candidates7, candidates23
    )
    print("exact c=0 sector: one bottom row (0,0), zero candidates", flush=True)
    print("exact reflected candidate bijection: 8514/8514", flush=True)

    result7 = classify(7, center7, candidates7)
    result23 = classify(23, center23, candidates23)
    for mask, result in ((7, result7), (23, result23)):
        print(
            f"mask={mask}: bottom=1955 (1954 c!=0 + 1 c=0); "
            "c!=0 candidates=8514; c=0 candidates=0; survivors=0"
        )
        print(f"mask={mask}: histogram={result[0]}")
        print(f"mask={mask}: exclusion digest={result[1]}")

    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode(
            "ascii"
        )
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: exact exhaustive exclusion for masks 7 and 23")


if __name__ == "__main__":
    main()
