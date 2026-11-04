#!/usr/bin/env python3
"""Exact exclusion of the nontrivial-R83 ambient point with mask 14.

The human-readable proof is in ``docs/D7_R83_MASK14_EXCLUSION.md``.  This
primary certificate reconstructs the canonical v45 point, applies explicit
determinant-one inversion and unit conjugators, derives all point-specific
outward bounds, and exhausts the determinant-one polygon search by default.

All decisive arithmetic is exact.  Rational root intervals certify every
real comparison; floating point is not used.
"""

from fractions import Fraction as Q
from hashlib import sha256
from math import floor

import sympy as sp

import d7_ambient_points_v45 as ambient
import residual_orbit_k17_k38_v42 as core


LABEL = "R83:nontrivial"
MASK = 14
INDEX = 2
EXPONENTS = (0, 0, 0, 0, -2)
ALL_POSITIVE = (1,) * core.DEGREE

EXPECTED_CENTER = (
    Q(-659, 166), Q(-4873, 166), Q(1009, 166), Q(7581, 166),
    Q(-435, 166), Q(-1169, 83),
)
EXPECTED_HEIGHT_SQUARED = (
    Q(-4185, 27556), Q(-31293, 27556), Q(-8719, 13778),
    Q(7089, 6889), Q(8025, 27556), Q(-1428, 6889),
)
EXPECTED_MODULUS_SQUARED = (
    Q(-5131685, 6889), Q(-40296014, 6889), Q(-36397872, 6889),
    Q(25562409, 6889), Q(16802613, 6889), Q(-2317685, 6889),
)
EXPECTED_BOUNDS = (
    (Q(531, 100), Q(79, 100), Q(51, 100), Q(123, 100), Q(117, 100), Q(221, 50)),
    (Q(271, 20), Q(163, 100), Q(9, 100), Q(149, 100), Q(193, 100), Q(169, 20)),
    (Q(601, 50), Q(959, 50), Q(421, 100), Q(146, 25), Q(137, 20), Q(649, 100)),
    (Q(45744, 25), Q(1223, 25), Q(89, 50), Q(1563, 50), Q(6399, 100), Q(26281, 50)),
    (Q(471, 100), Q(46, 5), Q(594, 25), Q(24, 5), Q(207, 50), Q(17, 5)),
    (Q(2867, 4), Q(469, 20), Q(251, 25), Q(642, 25), Q(3873, 100), Q(27503, 100)),
)
EXPECTED_HISTOGRAM = {
    (0, 0): 1802,
    (0, 1): 5676,
    (0, 2): 974,
    (0, 3): 468,
    (0, 4): 30,
    (0, 5): 6,
    (1, 0): 66,
    (1, 1): 12,
    (1, 2): 2,
}
EXPECTED_DIGEST = "19f1bba6f56200bac45a52810a56669c51824f22a9b3fd410db6a1fc2d69a5a1"


def to_core(value):
    polynomial = sp.Poly(ambient.reduce_field(value), ambient.T, domain=sp.QQ)
    return tuple(
        Q(int(polynomial.nth(index).p), int(polynomial.nth(index).q))
        for index in range(6)
    )


def fixed_invariants(matrix):
    trace = ambient.alpha_trace(INDEX)
    lower_left = matrix[1, 0]
    center = ambient.reduce_field(
        (matrix[0, 0] - matrix[1, 1]) / (2 * lower_left)
    )
    height_squared = ambient.reduce_field(
        (4 * ambient.PI - trace * trace) / (4 * lower_left * lower_left)
    )
    modulus_squared = ambient.reduce_l(center * center + height_squared)
    assert ambient.sign_vector(height_squared) == ALL_POSITIVE
    return center, height_squared, modulus_squared


def selected_point():
    base = {
        label: matrix
        for label, matrix, _ in ambient.module_templates(INDEX)
    }[LABEL]
    oriented = ambient.conjugate_by_unit(base, ambient.SQUARE_CLASS_UNITS[MASK])
    ambient.assert_integral_matrix(oriented)
    source_center, source_y2, source_modulus2 = fixed_invariants(oriented)

    inversion = sp.Matrix([[0, -1], [1, 0]])
    inversion_inverse = sp.Matrix([[0, 1], [-1, 0]])
    assert inversion.det() == 1
    inverted = ambient.matrix_multiply(
        ambient.matrix_multiply(inversion, oriented), inversion_inverse
    )
    ambient.assert_integral_matrix(inverted)
    inverted_center, inverted_y2, inverted_modulus2 = fixed_invariants(inverted)
    assert ambient.reduce_field(
        inverted_center + source_center / source_modulus2
    ) == 0
    assert ambient.reduce_field(
        inverted_y2 - source_y2 / (source_modulus2**2)
    ) == 0
    assert ambient.reduce_field(inverted_modulus2 - 1 / source_modulus2) == 0

    u = ambient.POSITIVE_NONSQUARE
    unit = ambient.reduce_field(u ** EXPONENTS[-1])
    unit_inverse = ambient.reduce_field(1 / unit)
    conjugator = sp.diag(unit, unit_inverse)
    conjugator_inverse = sp.diag(unit_inverse, unit)
    assert ambient.reduce_l(conjugator.det() - 1) == 0
    ambient.assert_integral_matrix(conjugator)
    ambient.assert_integral_matrix(conjugator_inverse)
    balanced = ambient.matrix_multiply(
        ambient.matrix_multiply(conjugator, inverted), conjugator_inverse
    )
    ambient.assert_integral_matrix(balanced)
    ambient.assert_characteristic(balanced, INDEX)

    center, y2, modulus2 = fixed_invariants(balanced)
    square = ambient.reduce_l(unit * unit)
    fourth = ambient.reduce_l(square * square)
    assert ambient.reduce_l(center - square * inverted_center) == 0
    assert ambient.reduce_l(y2 - fourth * inverted_y2) == 0
    assert ambient.reduce_l(modulus2 - fourth * inverted_modulus2) == 0

    center_core = to_core(center)
    y2_core = to_core(y2)
    modulus2_core = to_core(modulus2)
    assert center_core == EXPECTED_CENTER
    assert y2_core == EXPECTED_HEIGHT_SQUARED
    assert modulus2_core == EXPECTED_MODULUS_SQUARED
    assert core.norm(y2_core) == Q(1, 28217344)
    assert core.norm(modulus2_core) == Q(2097130411225081, 6889)
    assert abs(ambient.absolute_norm(balanced[1, 0])) == 6889
    return center_core, y2_core, modulus2_core


def outward_cent(value: Q) -> Q:
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(y2, modulus2):
    bounds = [[] for _ in range(6)]
    for index in range(core.DEGREE):
        y_lo, y_hi = core.sqrt_interval(
            core.interval_eval(y2, core.ROOT_INTERVALS[index])
        )
        modulus_lo, modulus_hi = core.sqrt_interval(
            core.interval_eval(modulus2, core.ROOT_INTERVALS[index])
        )
        assert y_lo > 0 and modulus_lo > 0
        den = core.sqrt_interval(
            (y_hi / core.M_LOWER[index], y_hi / core.M_LOWER[index])
        )[1]
        c_entry = core.sqrt_interval(
            (
                1 / (core.M_LOWER[index] * y_lo),
                1 / (core.M_LOWER[index] * y_lo),
            )
        )[1]
        num = core.U_UPPER[index] * den
        a_entry = num / y_lo
        d_entry = den + c_entry * modulus_hi
        b_entry = num + a_entry * modulus_hi
        for target, value in zip(
            bounds, (den, num, a_entry, b_entry, c_entry, d_entry)
        ):
            target.append(outward_cent(value))
    answer = tuple(tuple(values) for values in bounds)
    assert answer == EXPECTED_BOUNDS
    return answer


def main() -> None:
    print("EXACT NONTRIVIAL-R83 MASK-14 ORBIT CERTIFICATE v46")
    print("canonical point: label=R83:nontrivial, mask=14, partner mask=62")
    center, y2, modulus2 = selected_point()
    bounds = derive_bounds(y2, modulus2)
    labels = ("denominator", "numerator", "a", "b", "c", "d")
    for label, bound in zip(labels, bounds):
        print(f"{label} outward bounds:", bound)

    den, num, a_bound, b_bound, c_bound, d_bound = bounds
    core.Y2 = y2
    core.DEN_BOUND = den
    core.NUM_BOUND = num
    core.A_BOUND = a_bound
    core.B_BOUND = b_bound
    core.C_BOUND = c_bound
    core.D_BOUND = d_bound

    a_values = core.enumerate_trace_ellipsoid(
        core.ZERO,
        sum((value * value for value in a_bound), Q(0)),
        lambda value: core.within_integral_box(value, a_bound),
    )
    c_values = core.enumerate_trace_ellipsoid(
        core.ZERO,
        sum((value * value for value in c_bound), Q(0)),
        lambda value: core.within_integral_box(value, c_bound),
    )
    assert (len(a_values), len(c_values)) == (23923, 6605)
    print("exact separate entry counts a,c:", len(a_values), len(c_values))

    result = core.enumerate_case(INDEX, center, a_values, c_values)
    assert result == (
        2093,
        9036,
        0,
        EXPECTED_HISTOGRAM,
        EXPECTED_DIGEST,
    )
    print("first-failure histogram:", sorted(result[3].items()))
    print("PASS: exact exhaustive exclusion of nontrivial-R83 mask 14")


if __name__ == "__main__":
    main()
