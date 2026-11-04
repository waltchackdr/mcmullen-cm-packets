#!/usr/bin/env python3
"""Exact paired exclusion for nontrivial-R83 masks 15 and 31.

The v47 reflection proof is reused only as an exact bounded-enumeration core.
This certificate independently reconstructs the new point pair, derives its
point-specific outward bounds, freezes all counts, and evaluates the
quadrilateral separately at both centers.  The determinant-minus-one
reflection is a search bijection, not an SL_2 or GL_2^+ point identification.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from math import floor

import sympy as sp

import d7_ambient_points_v45 as ambient
import d7_r83_masks14_30_v47 as pair_core
import residual_orbit_k17_k38_v42 as core


MASKS = (15, 31)
ROOT_PARTNERS = {15: 63, 31: 47}
INDEX = 2
LABEL = "R83:nontrivial"
EXPONENTS = (1, 1, 1, 0, 0)
EXPECTED_CENTER = (
    Q(3219, 166), Q(6039, 166), Q(-993, 166),
    Q(-1976, 83), Q(-343, 166), Q(399, 166),
)
EXPECTED_Y2 = (
    Q(18653, 27556), Q(12115, 13778), Q(-13885, 27556),
    Q(-17105, 27556), Q(1289, 13778), Q(1475, 13778),
)
EXPECTED_MODULUS2 = (
    Q(1236236, 6889), Q(-1626587, 6889), Q(-6792283, 6889),
    Q(-1654398, 6889), Q(2923871, 6889), Q(1172479, 6889),
)
EXPECTED_BOUNDS = (
    (Q(187, 50), Q(29, 20), Q(74, 25), Q(4, 5), Q(9, 10), Q(23, 20)),
    (Q(191, 20), Q(151, 50), Q(53, 100), Q(49, 50), Q(37, 25), Q(219, 100)),
    (Q(341, 20), Q(519, 50), Q(18, 25), Q(891, 100), Q(223, 25), Q(626, 25)),
    (Q(64481, 50), Q(9041, 100), Q(1047, 100), Q(512, 25), Q(1228, 25), Q(3406, 25)),
    (Q(167, 25), Q(249, 50), Q(101, 25), Q(183, 25), Q(27, 5), Q(1311, 100)),
    (Q(50517, 100), Q(2167, 50), Q(5909, 100), Q(1683, 100), Q(2973, 100), Q(7129, 100)),
)
EXPECTED_HISTOGRAMS = {
    15: {
        (0, 0): 1802, (0, 1): 5062, (0, 2): 1006, (0, 3): 320,
        (0, 4): 20, (0, 5): 12, (1, 0): 32, (1, 1): 18,
        (1, 2): 4, (1, 4): 2,
    },
    31: {
        (0, 0): 1814, (0, 1): 5020, (0, 2): 968, (0, 3): 366,
        (0, 4): 34, (0, 5): 2, (1, 0): 52, (1, 1): 16,
        (1, 2): 2, (1, 3): 4,
    },
}
EXPECTED_DIGESTS = {
    15: "66684a7a1aa493a7e37b6b4fc6e97b3663c70e79700724eaca94e170a4d946b4",
    31: "2fa340baeab5f877f4ea61d9c2e9cc850ee54221d257396b1320dc2c048ef023",
}
EXPECTED_AGGREGATE = "974f7951d31b613dd906ba42af6ae3e790b5980625976bf1dc9eb7710ba1ed86"


def selected_pair():
    base = {
        label: matrix
        for label, matrix, _ in ambient.module_templates(INDEX)
    }[LABEL]
    assert ambient.reduce_l(
        ambient.SQUARE_CLASS_UNITS[31] + ambient.SQUARE_CLASS_UNITS[15]
    ) == 0
    oriented = {
        mask: ambient.conjugate_by_unit(base, ambient.SQUARE_CLASS_UNITS[mask])
        for mask in MASKS
    }
    reflection = sp.diag(1, -1)
    assert reflection.det() == -1
    assert ambient.sign_vector(sp.Integer(-1)) == (-1,) * 6
    assert ambient.matrix_reduce(
        oriented[31] - reflection * oriented[15] * reflection
    ) == sp.zeros(2)

    inversion = sp.Matrix([[0, -1], [1, 0]])
    inversion_inverse = sp.Matrix([[0, 1], [-1, 0]])
    assert inversion.det() == 1
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
        inverted[31] - reflection * inverted[15] * reflection
    ) == sp.zeros(2)

    generators = ambient.XI_UNITS[:4] + (ambient.POSITIVE_NONSQUARE,)
    unit = sp.Integer(1)
    for generator, exponent in zip(generators, EXPONENTS):
        unit = ambient.reduce_l(unit * generator**exponent)
    assert abs(ambient.absolute_norm(unit)) == 1
    unit_inverse = ambient.reduce_field(1 / unit)
    conjugator = sp.diag(unit, unit_inverse)
    conjugator_inverse = sp.diag(unit_inverse, unit)
    assert ambient.reduce_l(conjugator.det() - 1) == 0
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
        balanced[31] - reflection * balanced[15] * reflection
    ) == sp.zeros(2)

    invariants = {}
    for mask, matrix in balanced.items():
        center, y2, modulus2 = pair_core.v46.fixed_invariants(matrix)
        invariants[mask] = tuple(
            pair_core.v46.to_core(value) for value in (center, y2, modulus2)
        )
    assert invariants[15] == (EXPECTED_CENTER, EXPECTED_Y2, EXPECTED_MODULUS2)
    assert invariants[31] == (
        core.neg(EXPECTED_CENTER), EXPECTED_Y2, EXPECTED_MODULUS2
    )
    assert core.norm(EXPECTED_Y2) == Q(1, 28217344)
    assert core.norm(EXPECTED_MODULUS2) == Q(2097130411225081, 6889)
    assert all(abs(ambient.absolute_norm(matrix[1, 0])) == 6889 for matrix in balanced.values())
    return balanced, invariants


def outward_hundredth(value: Q) -> Q:
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
    assert (len(a_values), len(c_values)) == (24103, 6611)
    return a_values, c_values


def classify(mask, x, candidates):
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
    assert not survivors
    assert len(exclusions) == len(candidates) == 8278
    assert histogram == EXPECTED_HISTOGRAMS[mask]
    assert sum(histogram.values()) == 8278
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
    print("EXACT NONTRIVIAL-R83 MASKS {15,31} ORBIT CERTIFICATE v48")
    assert ambient.PAIR_XOR[LABEL] == 48
    assert all((LABEL, mask) in ambient.POINT_RECORDS for mask in MASKS)
    assert ROOT_PARTNERS == {mask: mask ^ 48 for mask in MASKS}
    assert {15, 63}.isdisjoint({31, 47})

    _, invariants = selected_pair()
    center15, y2, modulus2 = invariants[15]
    center31, _, _ = invariants[31]
    bounds = derive_bounds(y2, modulus2)
    pair_core.install_point(y2, bounds)
    a_values, c_values = entry_lists(bounds)
    print("exact shared a/c counts:", len(a_values), len(c_values), flush=True)
    candidates15 = pair_core.bounded_candidates(
        center15,
        a_values,
        c_values,
        expected_counts=(2045, 8278, 0),
        progress_label="masks15/31",
    )
    candidates31 = [pair_core.reflect(gamma) for gamma in candidates15]
    assert len(set(candidates15)) == len(set(candidates31)) == 8278
    result15 = classify(15, center15, candidates15)
    result31 = classify(31, center31, candidates31)
    for mask, result in ((15, result15), (31, result31)):
        print(f"mask={mask}: bottom=2045; c!=0=8278; c=0=0; survivors=0")
        print(f"mask={mask}: histogram={result[0]}")
        print(f"mask={mask}: exclusion digest={result[1]}")
    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode(
            "ascii"
        )
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: exact exhaustive exclusion for masks 15 and 31")


if __name__ == "__main__":
    main()
