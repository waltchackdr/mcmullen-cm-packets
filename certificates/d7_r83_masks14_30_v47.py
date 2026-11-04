#!/usr/bin/env python3
"""Exact paired exclusion for the nontrivial-R83 masks 14 and 30.

This certificate subsumes the v46 mask-14 search.  It reconstructs the two
distinct canonical ambient points, enumerates the complete bounded
determinant-one set once, transfers that set by the exact reflection

    (a,b,c,d) -> (a,-b,-c,d),

and performs separate exact quadrilateral tests at the two centers.  The
reflection is only a search bijection; it is not asserted to lie in
GL_2^+(O_L) or SL_2(O_L), and polygon invariance is not assumed.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256

import sympy as sp

import d7_ambient_points_v45 as ambient
import d7_r83_mask14_v46 as v46
import residual_orbit_k17_k38_v42 as core


MASKS = (14, 30)
ROOT_PARTNERS = {14: 62, 30: 46}
EXPECTED_HISTOGRAMS = {
    14: {
        (0, 0): 1802,
        (0, 1): 5676,
        (0, 2): 974,
        (0, 3): 468,
        (0, 4): 30,
        (0, 5): 6,
        (1, 0): 66,
        (1, 1): 12,
        (1, 2): 2,
    },
    30: {
        (0, 0): 1798,
        (0, 1): 5640,
        (0, 2): 982,
        (0, 3): 500,
        (0, 4): 14,
        (0, 5): 6,
        (1, 0): 64,
        (1, 1): 20,
        (1, 2): 6,
        (1, 3): 2,
        (1, 4): 4,
    },
}
EXPECTED_DIGESTS = {
    14: "19f1bba6f56200bac45a52810a56669c51824f22a9b3fd410db6a1fc2d69a5a1",
    30: "2334a40fed0d008d44c2e892871bec20fc04b226a99df899a8927df81595faf6",
}
EXPECTED_AGGREGATE_DIGEST = (
    "8855a9cae94b13b3b18e14110fb6bf573355ce1b3c3cec7e082f9e0ba0892a8f"
)


@dataclass(frozen=True)
class Result:
    bottom_rows: int
    nonzero_candidates: int
    zero_candidates: int
    survivors: int
    histogram: tuple[tuple[tuple[int, int], int], ...]
    digest: str


def selected_pair():
    base = {
        label: matrix
        for label, matrix, _ in ambient.module_templates(v46.INDEX)
    }[v46.LABEL]
    assert ambient.reduce_l(
        ambient.SQUARE_CLASS_UNITS[30] + ambient.SQUARE_CLASS_UNITS[14]
    ) == 0
    oriented14 = ambient.conjugate_by_unit(
        base, ambient.SQUARE_CLASS_UNITS[14]
    )
    oriented = ambient.conjugate_by_unit(base, ambient.SQUARE_CLASS_UNITS[30])
    ambient.assert_integral_matrix(oriented14)
    ambient.assert_integral_matrix(oriented)

    reflection = sp.diag(1, -1)
    assert reflection.det() == -1
    assert ambient.sign_vector(sp.Integer(-1)) == (-1,) * 6
    assert ambient.matrix_reduce(
        oriented - reflection * oriented14 * reflection
    ) == sp.zeros(2)

    inversion = sp.Matrix([[0, -1], [1, 0]])
    inversion_inverse = sp.Matrix([[0, 1], [-1, 0]])
    assert inversion.det() == 1
    inverted14 = ambient.matrix_multiply(
        ambient.matrix_multiply(inversion, oriented14), inversion_inverse
    )
    inverted = ambient.matrix_multiply(
        ambient.matrix_multiply(inversion, oriented), inversion_inverse
    )
    ambient.assert_integral_matrix(inverted14)
    ambient.assert_integral_matrix(inverted)
    assert ambient.matrix_reduce(
        inverted - reflection * inverted14 * reflection
    ) == sp.zeros(2)

    unit = ambient.reduce_field(ambient.POSITIVE_NONSQUARE ** -2)
    unit_inverse = ambient.reduce_field(1 / unit)
    conjugator = sp.diag(unit, unit_inverse)
    conjugator_inverse = sp.diag(unit_inverse, unit)
    assert ambient.reduce_l(conjugator.det() - 1) == 0
    ambient.assert_integral_matrix(conjugator)
    ambient.assert_integral_matrix(conjugator_inverse)
    balanced14 = ambient.matrix_multiply(
        ambient.matrix_multiply(conjugator, inverted14), conjugator_inverse
    )
    balanced = ambient.matrix_multiply(
        ambient.matrix_multiply(conjugator, inverted), conjugator_inverse
    )
    ambient.assert_integral_matrix(balanced14)
    ambient.assert_integral_matrix(balanced)
    ambient.assert_characteristic(balanced14, v46.INDEX)
    ambient.assert_characteristic(balanced, v46.INDEX)
    assert ambient.matrix_reduce(
        balanced - reflection * balanced14 * reflection
    ) == sp.zeros(2)
    center14, y2_14, modulus2_14 = v46.fixed_invariants(balanced14)
    center, y2, modulus2 = v46.fixed_invariants(balanced)
    return (
        balanced14,
        balanced,
        v46.to_core(center14),
        v46.to_core(center),
        v46.to_core(y2_14),
        v46.to_core(y2),
        v46.to_core(modulus2_14),
        v46.to_core(modulus2),
    )


def install_point(y2, bounds) -> None:
    core.Y2 = y2
    (
        core.DEN_BOUND,
        core.NUM_BOUND,
        core.A_BOUND,
        core.B_BOUND,
        core.C_BOUND,
        core.D_BOUND,
    ) = bounds


def entry_lists(bounds):
    a_bound = bounds[2]
    c_bound = bounds[4]
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
    return a_values, c_values


def bounded_candidates(
    x,
    a_values,
    c_values,
    expected_counts=(2093, 9036, 0),
    progress_label="masks14/30",
):
    den_squares = core.square_bounds(core.DEN_BOUND)
    num_squares = core.square_bounds(core.NUM_BOUND)
    den_radius = sum(den_squares, Q(0))
    num_radius = sum(num_squares, Q(0))
    bottom_rows = []
    x_boxes = tuple(core.scaled_field_box(x, i) for i in range(core.DEGREE))
    y2_boxes = tuple(core.scaled_field_box(core.Y2, i) for i in range(core.DEGREE))
    a_boxes_cache = {
        a: tuple(core.scaled_eval_integral(a, i) for i in range(core.DEGREE))
        for a in a_values
    }

    def numerator_bounded(a, b):
        unresolved = False
        a_boxes = a_boxes_cache[a]
        for i in range(core.DEGREE):
            b_box = core.scaled_eval_integral(b, i)
            real_box = core.scaled_add(
                core.scaled_mul(a_boxes[i], x_boxes[i]), b_box
            )
            square_box = core.scaled_add(
                core.scaled_mul(real_box, real_box),
                core.scaled_mul(
                    core.scaled_mul(a_boxes[i], a_boxes[i]), y2_boxes[i]
                ),
            )
            decision = core.scaled_le(square_box, num_squares[i])
            if decision is False:
                return False
            unresolved |= decision is None
        if not unresolved:
            return True
        return core.all_le_bounds(
            core.numerator_square(core.as_field(a), core.as_field(b), x),
            num_squares,
        )

    for c_index, c in enumerate(c_values):
        cq = core.as_field(c)
        shifted_center = core.neg(core.mul(cq, x))
        c_boxes = tuple(
            core.scaled_eval_integral(c, i) for i in range(core.DEGREE)
        )

        def denominator_predicate(d):
            unresolved = False
            for i in range(core.DEGREE):
                d_box = core.scaled_eval_integral(d, i)
                real_box = core.scaled_add(
                    core.scaled_mul(c_boxes[i], x_boxes[i]), d_box
                )
                square_box = core.scaled_add(
                    core.scaled_mul(real_box, real_box),
                    core.scaled_mul(
                        core.scaled_mul(c_boxes[i], c_boxes[i]), y2_boxes[i]
                    ),
                )
                decision = core.scaled_le(square_box, den_squares[i])
                if decision is False:
                    return False
                unresolved |= decision is None
            if not unresolved:
                return True
            return core.all_le_bounds(
                core.denominator_square(cq, core.as_field(d), x), den_squares
            )

        bottom_rows.extend(
            (c, d)
            for d in core.enumerate_trace_ellipsoid(
                shifted_center, den_radius, denominator_predicate
            )
        )
        if c_index and c_index % 500 == 0:
            print(
                f"  {progress_label}: c entries {c_index}/{len(c_values)}; "
                f"bottom rows={len(bottom_rows)}",
                flush=True,
            )

    nonzero = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != core.ZERO]
    for row_index, (c, d) in enumerate(nonzero_rows):
        cq, dq = core.as_field(c), core.as_field(d)
        affine_data = core.affine_b_data(c, d)
        for a in a_values:
            b = core.affine_integral_b_fast(a, affine_data)
            if b is None or not numerator_bounded(a, b):
                continue
            aq, bq = core.as_field(a), core.as_field(b)
            assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
            nonzero.append((a, b, c, d))
        if row_index and row_index % 250 == 0:
            print(
                f"  {progress_label}: determinant rows {row_index}/{len(nonzero_rows)}; "
                f"candidates={len(nonzero)}",
                flush=True,
            )

    zero = []
    for c, d in bottom_rows:
        if c != core.ZERO or d == core.ZERO:
            continue
        dq = core.as_field(d)
        aq = core.inverse(dq)
        if not core.integral(aq):
            continue
        a = core.integral_tuple(aq)
        if not core.within_integral_box(a, core.A_BOUND):
            continue
        shifted_center = core.neg(core.mul(aq, x))

        def numerator_predicate(b):
            return numerator_bounded(a, b)

        for b in core.enumerate_trace_ellipsoid(
            shifted_center, num_radius, numerator_predicate
        ):
            assert core.mul(aq, dq) == core.ONE
            zero.append((a, b, c, d))

    assert (len(bottom_rows), len(nonzero), len(zero)) == expected_counts
    return nonzero + zero


def reflect(gamma):
    a, b, c, d = gamma
    reflected = a, tuple(-q for q in b), tuple(-q for q in c), d
    aq, bq, cq, dq = map(core.as_field, reflected)
    assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
    assert reflect_raw(reflected) == gamma
    return reflected


def reflect_raw(gamma):
    a, b, c, d = gamma
    return a, tuple(-q for q in b), tuple(-q for q in c), d


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
    assert len(exclusions) + len(survivors) == len(candidates)
    digest = sha256()
    for gamma, failure in sorted(exclusions):
        flat = tuple(coefficient for entry in gamma for coefficient in entry)
        digest.update(
            (",".join(map(str, flat)) + f";{failure[0]},{failure[1]}\n").encode(
                "ascii"
            )
        )
    result = Result(
        2093,
        9036,
        0,
        len(survivors),
        tuple(sorted(histogram.items())),
        digest.hexdigest(),
    )
    assert result.survivors == 0
    assert dict(result.histogram) == EXPECTED_HISTOGRAMS[mask]
    assert sum(count for _, count in result.histogram) == 9036
    assert result.digest == EXPECTED_DIGESTS[mask]
    return result


def main() -> None:
    print("EXACT NONTRIVIAL-R83 MASKS {14,30} ORBIT CERTIFICATE v47")
    assert ambient.PAIR_XOR[v46.LABEL] == 48
    assert all((v46.LABEL, mask) in ambient.POINT_RECORDS for mask in MASKS)
    assert ROOT_PARTNERS == {mask: mask ^ 48 for mask in MASKS}

    center14, y2, modulus2 = v46.selected_point()
    (
        _,
        _,
        checked_center14,
        center30,
        checked_y2,
        y2_30,
        checked_modulus2,
        modulus2_30,
    ) = selected_pair()
    assert checked_center14 == center14
    assert checked_y2 == y2
    assert checked_modulus2 == modulus2
    assert center30 == core.neg(center14)
    assert y2_30 == y2
    assert modulus2_30 == modulus2
    bounds = v46.derive_bounds(y2, modulus2)
    install_point(y2, bounds)
    a_values, c_values = entry_lists(bounds)
    print("exact shared a/c counts:", len(a_values), len(c_values), flush=True)

    candidates14 = bounded_candidates(center14, a_values, c_values)
    candidates30 = [reflect(gamma) for gamma in candidates14]
    assert len(set(candidates30)) == len(candidates30) == 9036
    result14 = classify(14, center14, candidates14)
    result30 = classify(30, center30, candidates30)
    for mask, result in ((14, result14), (30, result30)):
        print(
            f"mask={mask}: bottom={result.bottom_rows}; c!=0={result.nonzero_candidates}; "
            f"c=0={result.zero_candidates}; survivors={result.survivors}"
        )
        print(f"mask={mask}: histogram={dict(result.histogram)}")
        print(f"mask={mask}: exclusion digest={result.digest}")

    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode("ascii")
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE_DIGEST
    print("aggregate labelled digest:", aggregate)
    print("PASS: exact exhaustive exclusion for two distinct nontrivial-R83 points")


if __name__ == "__main__":
    main()
