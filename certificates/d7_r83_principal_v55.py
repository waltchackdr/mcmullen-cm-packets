#!/usr/bin/env python3
"""Exact primary exclusion certificate for principal-Picard proper-R83.

The certificate reconstructs every canonical regular-order representative,
uses an explicit integral determinant-one unit balancing, derives outward
rational polygon bounds, exhausts both determinant branches, and applies the
exact normalized-quadrilateral predicates.  The canonical masks are
0,...,15,32,...,47; xor 16 pairs the two CM roots of one fixed point and is
not used as a GL2+/SL2 orbit shortcut.  Every canonical center is tested.

"Principal Picard class" refers to the identity ideal class of the proper
order R83.  It is unrelated to principal polarization.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from math import floor

import sympy as sp

import d7_ambient_points_v45 as ambient
import d7_r83_masks14_30_v47 as enumeration
import d7_r83_principal_expected_v55 as frozen
import residual_orbit_k17_k38_v42 as core


INDEX = 2
LABEL = "R83:principal"
CANONICAL_MASKS = frozen.CANONICAL_MASKS
ALL_POSITIVE = (1,) * core.DEGREE
BALANCING_EXPONENTS = {
    0: (0, 0, 0, 0, 0), 1: (0, 0, 0, 0, 0),
    2: (0, -1, 0, 0, 0), 3: (0, -1, 0, 0, 0),
    4: (0, 0, 0, 0, 0), 5: (-1, 0, 0, 0, 0),
    6: (0, -1, 0, 0, 0), 7: (-1, -1, 0, 0, 0),
    8: (0, 0, 0, 0, 0), 9: (-1, 0, 0, 0, 0),
    10: (0, -1, 0, 0, 0), 11: (-1, 0, 0, 0, 0),
    12: (0, 0, 0, -1, 0), 13: (-1, 0, 0, -1, 0),
    14: (0, 0, -1, 0, 0), 15: (0, -1, -1, 0, 0),
    32: (0, 0, 0, 0, 0), 33: (0, 0, 0, 0, 0),
    34: (0, 0, 0, 0, 0), 35: (-1, 0, 1, 0, 0),
    36: (0, 0, 0, 0, 0), 37: (0, 0, 0, 0, 0),
    38: (0, 0, 0, 0, 0), 39: (0, -1, 0, 0, 0),
    40: (0, 0, 0, 0, 0), 41: (0, 0, 0, 0, 0),
    42: (0, 0, 0, 0, 0), 43: (0, -1, 0, 0, 0),
    44: (0, 0, 0, 0, 0), 45: (-1, 0, 0, 0, 0),
    46: (0, 0, 0, -1, 1), 47: (-1, 0, 0, 0, 0),
}


def to_core(value):
    polynomial = sp.Poly(ambient.reduce_field(value), ambient.T, domain=sp.QQ)
    return tuple(
        Q(int(polynomial.nth(i).p), int(polynomial.nth(i).q))
        for i in range(core.DEGREE)
    )


def field_power(value, exponent: int):
    if exponent < 0:
        inverse = ambient.reduce_field(1 / value)
        ambient.assert_integral(inverse)
        return field_power(inverse, -exponent)
    answer = sp.Integer(1)
    base = value
    n = exponent
    while n:
        if n & 1:
            answer = ambient.reduce_l(answer * base)
        base = ambient.reduce_l(base * base)
        n >>= 1
    return answer


def balancing_unit(mask: int):
    answer = sp.Integer(1)
    for generator, exponent in zip(
        ambient.XI_UNITS, BALANCING_EXPONENTS[mask]
    ):
        answer = ambient.reduce_l(answer * field_power(generator, exponent))
    ambient.assert_integral(answer)
    ambient.assert_integral(ambient.reduce_field(1 / answer))
    assert abs(ambient.absolute_norm(answer)) == 1
    return answer


def fixed_invariants(matrix):
    trace = ambient.alpha_trace(INDEX)
    c = matrix[1, 0]
    center = ambient.reduce_field((matrix[0, 0] - matrix[1, 1]) / (2 * c))
    y2 = ambient.reduce_field((4 * ambient.PI - trace * trace) / (4 * c * c))
    modulus2 = ambient.reduce_l(center * center + y2)
    assert ambient.sign_vector(y2) == ALL_POSITIVE
    assert ambient.sign_vector(modulus2) == ALL_POSITIVE
    return center, y2, modulus2


def regular_base():
    base = {
        label: matrix for label, matrix, _ in ambient.module_templates(INDEX)
    }[LABEL]
    ambient.assert_integral_matrix(base)
    ambient.assert_characteristic(base, INDEX)
    # The regular companion lattice has proper multiplier R83: the maximal
    # generator fails its integrality test in the v45 reconstruction.
    assert ambient.PAIR_XOR[LABEL] == 16
    assert ambient.ICM_COUNTS[83] == 3 and ambient.PICARD_COUNTS[83] == 2
    return base


def selected_point(mask: int, base):
    assert mask in CANONICAL_MASKS
    raw = ambient.conjugate_by_unit(base, ambient.SQUARE_CLASS_UNITS[mask])
    ambient.assert_integral_matrix(raw)
    ambient.assert_characteristic(raw, INDEX)
    raw_x, raw_y2, raw_modulus2 = fixed_invariants(raw)

    unit = balancing_unit(mask)
    inverse = ambient.reduce_field(1 / unit)
    conjugator = sp.diag(unit, inverse)
    conjugator_inverse = sp.diag(inverse, unit)
    ambient.assert_integral_matrix(conjugator)
    ambient.assert_integral_matrix(conjugator_inverse)
    assert ambient.reduce_l(conjugator.det() - 1) == 0
    balanced = ambient.matrix_multiply(
        ambient.matrix_multiply(conjugator, raw), conjugator_inverse
    )
    ambient.assert_integral_matrix(balanced)
    ambient.assert_characteristic(balanced, INDEX)
    center, y2, modulus2 = fixed_invariants(balanced)
    assert ambient.reduce_l(center - unit**2 * raw_x) == 0
    assert ambient.reduce_l(y2 - unit**4 * raw_y2) == 0
    assert ambient.reduce_l(modulus2 - unit**4 * raw_modulus2) == 0
    return to_core(center), to_core(y2), to_core(modulus2)


def outward_hundredth(value: Q) -> Q:
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(y2, modulus2):
    columns = [[] for _ in range(6)]
    for i in range(core.DEGREE):
        y_lo, y_hi = core.sqrt_interval(
            core.interval_eval(y2, core.ROOT_INTERVALS[i])
        )
        _, modulus_hi = core.sqrt_interval(
            core.interval_eval(modulus2, core.ROOT_INTERVALS[i])
        )
        assert y_lo > 0 and modulus_hi > 0
        den_hat = core.sqrt_interval(
            (y_hi / core.M_LOWER[i], y_hi / core.M_LOWER[i])
        )[1]
        c_hat = core.sqrt_interval(
            (1 / (core.M_LOWER[i] * y_lo), 1 / (core.M_LOWER[i] * y_lo))
        )[1]
        num_hat = core.U_UPPER[i] * den_hat
        a_hat = num_hat / y_lo
        d_hat = den_hat + c_hat * modulus_hi
        b_hat = num_hat + a_hat * modulus_hi
        for target, value in zip(
            columns, (den_hat, num_hat, a_hat, b_hat, c_hat, d_hat)
        ):
            target.append(outward_hundredth(value))
    return tuple(tuple(column) for column in columns)


def within_box(value, bounds):
    try:
        return core.within_integral_box(value, bounds)
    except ArithmeticError:
        field_value = core.as_field(value)
        for index, bound in enumerate(bounds):
            bound_value = core.scale(bound, core.ONE)
            if core.sign_at(core.add(bound_value, field_value), index) < 0:
                return False
            if core.sign_at(core.sub(bound_value, field_value), index) < 0:
                return False
        return True


def entry_lists(y2, bounds):
    enumeration.install_point(y2, bounds)
    a_values = core.enumerate_trace_ellipsoid(
        core.ZERO,
        sum((value * value for value in bounds[2]), Q(0)),
        lambda value: within_box(value, bounds[2]),
    )
    c_values = core.enumerate_trace_ellipsoid(
        core.ZERO,
        sum((value * value for value in bounds[4]), Q(0)),
        lambda value: within_box(value, bounds[4]),
    )
    assert {core.neg(value) for value in a_values} == set(a_values)
    assert {core.neg(value) for value in c_values} == set(c_values)
    return a_values, c_values


def polygon_result(center, candidates):
    histogram = {}
    survivors = []
    exclusions = []
    for gamma in candidates:
        inside, failure = core.in_all_coordinate_polygons(gamma, center)
        if inside:
            survivors.append(gamma)
        else:
            assert failure is not None
            histogram[failure] = histogram.get(failure, 0) + 1
            exclusions.append((gamma, failure))
    digest = sha256()
    for gamma, failure in sorted(exclusions):
        flat = tuple(coefficient for entry in gamma for coefficient in entry)
        digest.update(
            (",".join(map(str, flat)) + f";{failure[0]},{failure[1]}\n").encode(
                "ascii"
            )
        )
    assert not survivors
    assert sum(histogram.values()) == len(candidates)
    return tuple(sorted(histogram.items())), digest.hexdigest()


def parse_masks(value: str):
    if value == "all":
        return CANONICAL_MASKS
    answer = tuple(sorted({int(item) for item in value.split(",")}))
    if not answer or any(mask not in CANONICAL_MASKS for mask in answer):
        raise argparse.ArgumentTypeError(
            "use all or comma-separated masks in 0,...,15,32,...,47"
        )
    return answer


def certify_root_coverage():
    fibers = tuple(frozenset((mask, mask ^ 16)) for mask in CANONICAL_MASKS)
    assert len(fibers) == len(set(fibers)) == 32
    assert set().union(*fibers) == set(range(64))
    assert all((LABEL, mask) in ambient.POINT_RECORDS for mask in CANONICAL_MASKS)
    assert all(
        (LABEL, mask ^ 16) not in ambient.POINT_RECORDS
        for mask in CANONICAL_MASKS
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--masks", type=parse_masks, default=CANONICAL_MASKS)
    args = parser.parse_args()
    certify_root_coverage()
    base = regular_base()
    results = {}
    print("EXACT PRINCIPAL-PICARD proper-R83 ORBIT CERTIFICATE v55")
    for mask in args.masks:
        center, y2, modulus2 = selected_point(mask, base)
        bounds = derive_bounds(y2, modulus2)
        a_values, c_values = entry_lists(y2, bounds)
        assert (len(a_values), len(c_values)) == frozen.EXPECTED_ENTRY_COUNTS[mask]
        expected = frozen.EXPECTED[mask]
        candidates = enumeration.bounded_candidates(
            center,
            a_values,
            c_values,
            expected_counts=expected[:3],
            progress_label=f"principal-mask={mask:02d}",
        )
        histogram, digest = polygon_result(center, candidates)
        result = expected[:3] + (histogram, digest)
        assert result == expected[:5]
        results[mask] = result
        print(
            f"mask={mask:02d}: bottom={result[0]}; c!=0={result[1]}; "
            f"c=0={result[2]}; digest={digest}; survivors=0",
            flush=True,
        )
    assert tuple(sorted(results)) == args.masks
    if args.masks == CANONICAL_MASKS:
        assert frozen.batch_digest(results) == frozen.EXPECTED_BATCH_DIGEST
    print("PASS: exact principal-Picard proper-R83 exclusions")


if __name__ == "__main__":
    main()
