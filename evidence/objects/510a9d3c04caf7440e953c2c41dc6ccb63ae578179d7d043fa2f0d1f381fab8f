#!/usr/bin/env python3
"""Exact paired exclusion for nontrivial-R83 masks 13 and 29.

The registered v47 certificate is reused only for its exact bounded-matrix
enumeration kernel and candidate reflection.  This certificate reconstructs
the two canonical points from the v45 ambient packet, certifies their free
root fibers, applies an explicit determinant-one inversion and saturated-unit
balancing, derives and freezes all point-specific bounds, exhausts the
determinant-one search including the ``c=0`` branch, and tests the normalized
quadrilateral separately at the two centers.

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


MASKS = (13, 29)
ROOT_PARTNERS = {13: 61, 29: 45}
INDEX = 2
LABEL = "R83:nontrivial"
EXPONENTS = (1, 0, 1, 0, 0)
ALL_POSITIVE = (1,) * core.DEGREE

EXPECTED_BALANCING_UNIT = (
    1 + 3 * ambient.T - 2 * ambient.T**2 - 6 * ambient.T**3
    + ambient.T**4 + 2 * ambient.T**5
)
EXPECTED_BALANCING_UNIT_INVERSE = (
    10 + 54 * ambient.T - 7 * ambient.T**2 - 37 * ambient.T**3
    + ambient.T**4 + 6 * ambient.T**5
)
EXPECTED_CENTER = (
    Q(-1846, 83),
    Q(-1437, 166),
    Q(2106, 83),
    Q(-307, 166),
    Q(-501, 83),
    Q(195, 166),
)
EXPECTED_Y2 = (
    Q(10064, 6889),
    Q(53407, 27556),
    Q(-7100, 6889),
    Q(-38025, 27556),
    Q(4561, 27556),
    Q(1568, 6889),
)
EXPECTED_MODULUS2 = (
    Q(3644951, 6889),
    Q(4890422, 6889),
    Q(-2578147, 6889),
    Q(-3551387, 6889),
    Q(467266, 6889),
    Q(593566, 6889),
)
EXPECTED_BOUNDS = (
    (Q(29, 20), Q(31, 25), Q(179, 50), Q(103, 100), Q(31, 25), Q(81, 50)),
    (Q(369, 100), Q(129, 50), Q(16, 25), Q(63, 50), Q(41, 20), Q(77, 25)),
    (Q(4409, 100), Q(607, 50), Q(59, 100), Q(347, 50), Q(129, 20), Q(1781, 100)),
    (Q(49857, 100), Q(1932, 25), Q(1267, 100), Q(658, 25), Q(6791, 100), Q(958, 5)),
    (Q(432, 25), Q(291, 50), Q(167, 50), Q(57, 10), Q(391, 100), Q(233, 25)),
    (Q(1953, 10), Q(741, 20), Q(1431, 20), Q(2163, 100), Q(411, 10), Q(401, 4)),
)
EXPECTED_HISTOGRAMS = {
    13: {
        (0, 0): 1772,
        (0, 1): 5472,
        (0, 2): 1106,
        (0, 3): 328,
        (0, 4): 18,
        (0, 5): 6,
        (1, 0): 44,
        (1, 1): 10,
        (1, 2): 2,
        (1, 3): 2,
        (1, 4): 8,
    },
    29: {
        (0, 0): 1712,
        (0, 1): 5444,
        (0, 2): 1092,
        (0, 3): 394,
        (0, 4): 18,
        (0, 5): 14,
        (1, 0): 64,
        (1, 1): 24,
        (1, 4): 6,
    },
}
EXPECTED_DIGESTS = {
    13: "d66a04fa445505dd16cf2dad98ee6c9c57fe16d90aa4bbeb0b1bf2a1c2cc1a15",
    29: "b3ce2ca80967f6acf28bf10d6bb1707f330c768313680c369dd9d6d2cb3349de",
}
EXPECTED_AGGREGATE = "20ca8ce844583d40514ca6f65db26a302957b057fb0fd23c057cb7e52dab2e68"


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
    assert {13, 61}.isdisjoint({29, 45})
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
        ambient.SQUARE_CLASS_UNITS[29] + ambient.SQUARE_CLASS_UNITS[13]
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
        oriented[29] - reflection * oriented[13] * reflection
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
        tuple(pair_core.v46.to_core(value) for value in (raw_invariants[13][2],))[0]
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
        inverted[29] - reflection * inverted[13] * reflection
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
        balanced[29] - reflection * balanced[13] * reflection
    ) == sp.zeros(2)
    assert invariants[13] == (EXPECTED_CENTER, EXPECTED_Y2, EXPECTED_MODULUS2)
    assert invariants[29] == (
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
    assert (len(a_values), len(c_values)) == (23917, 6625)
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


def validate_reflection_bijection(center13, center29, y2, candidates13, candidates29):
    """Certify the exact involution, bounds identities, and complete converse."""
    assert center29 == core.neg(center13)
    assert len(candidates13) == len(candidates29) == 8768
    assert len(set(candidates13)) == len(candidates13)
    assert len(set(candidates29)) == len(candidates29)
    for gamma13, gamma29 in zip(candidates13, candidates29):
        assert pair_core.reflect_raw(gamma13) == gamma29
        assert pair_core.reflect_raw(gamma29) == gamma13
        a, b, c, d = gamma13
        ar, br, cr, dr = gamma29
        aq, bq, cq, dq = map(core.as_field, gamma13)
        arq, brq, crq, drq = map(core.as_field, gamma29)
        assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
        assert core.sub(core.mul(arq, drq), core.mul(brq, crq)) == core.ONE
        assert exact_square(c, d, center13, y2) == exact_square(
            cr, dr, center29, y2
        )
        assert exact_square(a, b, center13, y2) == exact_square(
            ar, br, center29, y2
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
    assert len(exclusions) == len(candidates) == 8768
    assert histogram == EXPECTED_HISTOGRAMS[mask]
    assert sum(histogram.values()) == 8768
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
    print("EXACT NONTRIVIAL-R83 MASKS {13,29} ORBIT CERTIFICATE v52")
    assert ambient.PAIR_XOR[LABEL] == 48
    assert all((LABEL, mask) in ambient.POINT_RECORDS for mask in MASKS)
    assert ROOT_PARTNERS == {mask: mask ^ 48 for mask in MASKS}
    assert {13, 61}.isdisjoint({29, 45})
    prior_fibers = (
        {14, 62}, {30, 46}, {15, 63}, {31, 47}, {6, 54}, {22, 38},
        {1, 49}, {17, 33},
    )
    assert all(
        {13, 61}.isdisjoint(fiber) and {29, 45}.isdisjoint(fiber)
        for fiber in prior_fibers
    )

    _, invariants = selected_pair()
    center13, y2, modulus2 = invariants[13]
    center29, y2_29, modulus2_29 = invariants[29]
    assert y2_29 == y2 and modulus2_29 == modulus2
    bounds = derive_bounds(y2, modulus2)
    pair_core.install_point(y2, bounds)
    a_values, c_values = entry_lists(bounds)
    print("exact shared a/c counts:", len(a_values), len(c_values), flush=True)

    candidates13 = pair_core.bounded_candidates(
        center13,
        a_values,
        c_values,
        expected_counts=(2043, 8768, 0),
        progress_label="masks13/29",
    )
    candidates29 = [pair_core.reflect(gamma) for gamma in candidates13]
    assert all(gamma[2] != core.ZERO for gamma in candidates13)
    assert all(gamma[2] != core.ZERO for gamma in candidates29)
    validate_reflection_bijection(
        center13, center29, y2, candidates13, candidates29
    )
    print("exact reflected candidate bijection: 8768/8768", flush=True)

    result13 = classify(13, center13, candidates13)
    result29 = classify(29, center29, candidates29)
    for mask, result in ((13, result13), (29, result29)):
        print(f"mask={mask}: bottom=2043; c!=0=8768; c=0=0; survivors=0")
        print(f"mask={mask}: histogram={result[0]}")
        print(f"mask={mask}: exclusion digest={result[1]}")

    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode(
            "ascii"
        )
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: exact exhaustive exclusion for masks 13 and 29")


if __name__ == "__main__":
    main()
