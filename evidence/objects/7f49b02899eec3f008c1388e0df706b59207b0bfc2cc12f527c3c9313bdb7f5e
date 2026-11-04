#!/usr/bin/env python3
"""Independent exact exclusion of the nontrivial-R83 mask-14 point.

This standard-library verifier imports only the independent quotient-field,
Sturm, ambient-packet, and residual-orbit foundations.  It never imports or
executes a primary certificate or experiment.  It reconstructs the point from
the conductor-83 ideal lattice, orientation mask 14, determinant-one inversion,
and saturated unit exponent ``(0,0,0,0,-2)``.  It then derives certified
outward entry bounds and exhausts the determinant-one polygon search.

All decisive arithmetic is exact.  Rational Sturm intervals certify every
real sign and bound; binary floating point is not used.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import floor

import d7_ambient_points_sturm_v45 as ambient
import residual_orbit_k17_k38_sturm_v42 as orbit


Q = Fraction
Vector = tuple[int, ...]
Gamma = tuple[Vector, Vector, Vector, Vector]
ZERO_VECTOR: Vector = (0,) * orbit.DEGREE
MASK = 14
INDEX = 2
EXPONENTS = (0, 0, 0, 0, -2)

EXPECTED_CENTER = ambient.LElement(
    (Q(-659, 166), Q(-4873, 166), Q(1009, 166), Q(7581, 166), Q(-435, 166), Q(-1169, 83))
)
EXPECTED_HEIGHT_SQUARED = ambient.LElement(
    (Q(-4185, 27556), Q(-31293, 27556), Q(-8719, 13778), Q(7089, 6889), Q(8025, 27556), Q(-1428, 6889))
)
EXPECTED_MODULUS_SQUARED = ambient.LElement(
    (Q(-5131685, 6889), Q(-40296014, 6889), Q(-36397872, 6889), Q(25562409, 6889), Q(16802613, 6889), Q(-2317685, 6889))
)

EXPECTED_BOUNDS = {
    "DEN": (Q(531, 100), Q(79, 100), Q(51, 100), Q(123, 100), Q(117, 100), Q(221, 50)),
    "NUM": (Q(271, 20), Q(163, 100), Q(9, 100), Q(149, 100), Q(193, 100), Q(169, 20)),
    "A": (Q(601, 50), Q(959, 50), Q(421, 100), Q(146, 25), Q(137, 20), Q(649, 100)),
    "B": (Q(45744, 25), Q(1223, 25), Q(89, 50), Q(1563, 50), Q(6399, 100), Q(26281, 50)),
    "C": (Q(471, 100), Q(46, 5), Q(594, 25), Q(24, 5), Q(207, 50), Q(17, 5)),
    "D": (Q(2867, 4), Q(469, 20), Q(251, 25), Q(642, 25), Q(3873, 100), Q(27503, 100)),
}

EXPECTED_HISTOGRAM = (
    ((0, 0), 1802),
    ((0, 1), 5676),
    ((0, 2), 974),
    ((0, 3), 468),
    ((0, 4), 30),
    ((0, 5), 6),
    ((1, 0), 66),
    ((1, 1), 12),
    ((1, 2), 2),
)
EXPECTED_DIGEST = "19f1bba6f56200bac45a52810a56669c51824f22a9b3fd410db6a1fc2d69a5a1"


def product(values: tuple[ambient.LElement, ...]) -> ambient.LElement:
    answer = ambient.ONE
    for value in values:
        answer *= value
    return answer


def selected_point() -> tuple[ambient.Matrix, ambient.LElement, ambient.LElement, ambient.LElement]:
    """Reconstruct the balanced point independently from the v45 inputs."""
    alpha = ambient.alpha_pair(INDEX)
    maximal = ambient.multiplication_matrix(alpha)
    assert ambient.pair_norm(alpha) == ambient.PI

    gamma83 = ambient.elt(-2, -2, 2, 0, -1)
    ideal_basis = ambient.matrix(16, gamma83, ambient.elt(58, 3), 3 * gamma83)
    ideal = ambient.conjugate_by(ambient.minverse(ideal_basis), maximal)
    assert ambient.all_entries_integral(ideal)

    xi5 = ambient.T**4 - 3 * ambient.T**2 + 1
    xi11 = ambient.T**3 + ambient.T**2 - 2 * ambient.T - 1
    xi13 = ambient.T**5 - 4 * ambient.T**3 + 3 * ambient.T + 1
    xi17 = ambient.T**3 - 2 * ambient.T
    unit_basis = (xi5, xi11, xi13, xi17, -ambient.ONE, 4 - ambient.T**2)
    orientation = product(tuple(unit for bit, unit in enumerate(unit_basis) if (MASK >> bit) & 1))
    assert orientation == xi11 * xi13 * xi17
    oriented = ambient.twist(ideal, orientation)
    assert ambient.all_entries_integral(oriented)

    inversion = ambient.matrix(0, -1, 1, 0)
    assert ambient.mdet(inversion) == ambient.ONE
    inverted = ambient.conjugate_by(inversion, oriented)
    assert ambient.all_entries_integral(inverted)

    u = 4 - ambient.T**2
    balancing_unit = u ** EXPONENTS[-1]
    conjugator = ambient.matrix(balancing_unit, 0, 0, balancing_unit.inverse())
    assert ambient.mdet(conjugator) == ambient.ONE
    assert ambient.all_entries_integral(conjugator)
    assert ambient.all_entries_integral(ambient.minverse(conjugator))
    balanced = ambient.conjugate_by(conjugator, inverted)
    assert ambient.all_entries_integral(balanced)
    assert ambient.mtrace(balanced) == ambient.pair_trace(alpha)
    assert ambient.mdet(balanced) == ambient.PI

    center, height_squared = ambient.exact_fixed_invariants(balanced)
    modulus_squared = center**2 + height_squared
    assert center == EXPECTED_CENTER
    assert height_squared == EXPECTED_HEIGHT_SQUARED
    assert modulus_squared == EXPECTED_MODULUS_SQUARED
    assert ambient.sign_vector(height_squared, ambient.embedded_roots()) == (1,) * orbit.DEGREE
    assert height_squared.norm() == Q(1, 28217344)
    assert modulus_squared.norm() == Q(2097130411225081, 6889)
    return balanced, center, height_squared, modulus_squared


def outward_hundredth(value: Q) -> Q:
    """Smallest hundredth strictly larger than a rational upper endpoint."""
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(
    height_squared: ambient.LElement,
    modulus_squared: ambient.LElement,
) -> dict[str, tuple[Q, ...]]:
    """Derive the six outward bounds from exact Sturm intervals."""
    den: list[Q] = []
    num: list[Q] = []
    a_bound: list[Q] = []
    b_bound: list[Q] = []
    c_bound: list[Q] = []
    d_bound: list[Q] = []
    for embedding, root in enumerate(orbit.ROOTS):
        height = orbit.sqrt_interval(orbit.interval_eval(height_squared, root), 160)
        modulus = orbit.sqrt_interval(orbit.interval_eval(modulus_squared, root), 160)
        assert height.lo > 0 and modulus.lo > 0

        den_upper = orbit.sqrt_interval(
            orbit.Interval.point(height.hi / orbit.M_LOWER[embedding]), 160
        ).hi
        c_upper = orbit.sqrt_interval(
            orbit.Interval.point(1 / (orbit.M_LOWER[embedding] * height.lo)), 160
        ).hi
        num_upper = orbit.U_UPPER[embedding] * den_upper
        a_upper = num_upper / height.lo
        d_upper = den_upper + c_upper * modulus.hi
        b_upper = num_upper + a_upper * modulus.hi

        den.append(outward_hundredth(den_upper))
        num.append(outward_hundredth(num_upper))
        a_bound.append(outward_hundredth(a_upper))
        b_bound.append(outward_hundredth(b_upper))
        c_bound.append(outward_hundredth(c_upper))
        d_bound.append(outward_hundredth(d_upper))

    answer = {
        "DEN": tuple(den),
        "NUM": tuple(num),
        "A": tuple(a_bound),
        "B": tuple(b_bound),
        "C": tuple(c_bound),
        "D": tuple(d_bound),
    }
    assert answer == EXPECTED_BOUNDS
    return answer


def install_point_data(height_squared: ambient.LElement, bounds: dict[str, tuple[Q, ...]]) -> None:
    """Install point-specific data in the reusable independent orbit core."""
    orbit.Y2 = height_squared
    orbit.DEN_BOUND = bounds["DEN"]
    orbit.NUM_BOUND = bounds["NUM"]
    orbit.A_BOUND = bounds["A"]
    orbit.B_BOUND = bounds["B"]
    orbit.C_BOUND = bounds["C"]
    orbit.D_BOUND = bounds["D"]
    orbit.DEN_RADIUS = orbit.sum_of_squares(orbit.DEN_BOUND)
    orbit.NUM_RADIUS = orbit.sum_of_squares(orbit.NUM_BOUND)
    orbit.Y2_SCALED = tuple(
        orbit.scaled_interval(orbit.interval_eval(height_squared, root))
        for root in orbit.ROOTS
    )


@dataclass(frozen=True)
class Result:
    bottom_rows: int
    nonzero_candidates: int
    zero_candidates: int
    histogram: tuple[tuple[tuple[int, int], int], ...]
    digest: str


def enumerate_case(x: ambient.LElement, a_values: list[Vector], c_values: list[Vector]) -> Result:
    """Exhaust the coupled determinant-one and coordinate-polygon search."""
    bottom_rows: list[tuple[Vector, Vector]] = []
    x_boxes = tuple(orbit.scaled_interval(orbit.interval_eval(x, root)) for root in orbit.ROOTS)
    a_box_cache = {
        a: tuple(orbit.scaled_vector_box(a, i) for i in range(orbit.DEGREE))
        for a in a_values
    }

    for c_index, c in enumerate(c_values):
        c_element = orbit.vector_element(c)
        center = tuple(-(c_element * x).coefficients[i] for i in range(orbit.DEGREE))
        c_boxes = tuple(orbit.scaled_vector_box(c, i) for i in range(orbit.DEGREE))
        bottom_rows.extend(
            (c, d)
            for d in orbit.enumerate_trace_ball(
                center,
                orbit.DEN_RADIUS,
                lambda d, c=c, c_boxes=c_boxes: orbit.bottom_bounded(c, d, x, c_boxes, x_boxes),
            )
        )
        if c_index and c_index % 500 == 0:
            print(
                f"  mask14: c entries {c_index}/{len(c_values)}; bottom rows={len(bottom_rows)}",
                flush=True,
            )
    print(f"  mask14: coupled bottom rows={len(bottom_rows)}", flush=True)

    nonzero: list[Gamma] = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != ZERO_VECTOR]
    for row_index, (c, d) in enumerate(nonzero_rows):
        affine = orbit.determinant_one_affine_map(c, d)
        for a in a_values:
            b = affine.evaluate(a)
            if b is not None and orbit.top_bounded(a, b, x, a_box_cache[a], x_boxes):
                gamma = (a, b, c, d)
                assert orbit.vector_element(a) * orbit.vector_element(d) - orbit.vector_element(b) * orbit.vector_element(c) == orbit.ONE
                nonzero.append(gamma)
        if row_index and row_index % 250 == 0:
            print(
                f"  mask14: determinant rows {row_index}/{len(nonzero_rows)}; candidates={len(nonzero)}",
                flush=True,
            )

    zero: list[Gamma] = []
    for c, d in bottom_rows:
        if c != ZERO_VECTOR or d == ZERO_VECTOR:
            continue
        a_element = orbit.vector_element(d).inverse()
        if not a_element.integral_power_basis():
            continue
        a = orbit.integral_vector(a_element)
        if not orbit.field_box_membership(a, orbit.A_BOUND):
            continue
        center = tuple(-(a_element * x).coefficients[i] for i in range(orbit.DEGREE))
        for b in orbit.enumerate_trace_ball(
            center,
            orbit.NUM_RADIUS,
            lambda b, a=a: orbit.top_bounded(a, b, x, a_box_cache[a], x_boxes),
        ):
            assert a_element * orbit.vector_element(d) == orbit.ONE
            zero.append((a, b, c, d))

    records: list[tuple[Gamma, tuple[int, int]]] = []
    for gamma in nonzero + zero:
        failure = orbit.first_polygon_exclusion(gamma, x)
        assert failure is not None
        records.append((gamma, failure))

    histogram_data: dict[tuple[int, int], int] = {}
    for _, failure in records:
        histogram_data[failure] = histogram_data.get(failure, 0) + 1
    result = Result(
        len(bottom_rows),
        len(nonzero),
        len(zero),
        tuple(sorted(histogram_data.items())),
        orbit.canonical_digest(records),
    )
    assert result.bottom_rows == 2093
    assert result.nonzero_candidates == 9036
    assert result.zero_candidates == 0
    assert result.histogram == EXPECTED_HISTOGRAM
    assert sum(count for _, count in result.histogram) == 9036
    assert result.digest == EXPECTED_DIGEST
    return result


def main() -> None:
    assert orbit.field_discriminant() == 453789
    print("INDEPENDENT d=7 NONTRIVIAL-R83 MASK-14 STURM CHECK v46")
    print("stdlib field/Sturm/LDL core; no primary certificate or experiment imports")
    balanced, center, height_squared, modulus_squared = selected_point()
    del balanced
    bounds = derive_bounds(height_squared, modulus_squared)
    install_point_data(height_squared, bounds)
    print("point reconstruction: mask=14; inversion; saturated exponents=(0,0,0,0,-2)")
    print("center:", center.coefficients)
    print("height squared:", height_squared.coefficients)
    print("derived outward bounds:", bounds)

    a_values = orbit.enumerate_trace_ball(
        (Q(0),) * orbit.DEGREE,
        orbit.sum_of_squares(orbit.A_BOUND),
        lambda vector: orbit.field_box_membership(vector, orbit.A_BOUND),
    )
    a_values.sort()
    print("exact a-entry count:", len(a_values), flush=True)
    assert len(a_values) == 23923
    c_values = orbit.enumerate_trace_ball(
        (Q(0),) * orbit.DEGREE,
        orbit.sum_of_squares(orbit.C_BOUND),
        lambda vector: orbit.field_box_membership(vector, orbit.C_BOUND),
    )
    c_values.sort()
    print("exact c-entry count:", len(c_values), flush=True)
    assert len(c_values) == 6605

    result = enumerate_case(center, a_values, c_values)
    print(
        f"mask14: bottom rows={result.bottom_rows}; c!=0 candidates={result.nonzero_candidates}; "
        f"c=0 candidates={result.zero_candidates}; polygon survivors=0"
    )
    print("mask14: first-exclusion histogram=", dict(result.histogram))
    print("mask14: canonical exclusion SHA-256=", result.digest)
    print("PASS: independent exact exhaustive exclusion reproduced for nontrivial-R83 mask 14")


if __name__ == "__main__":
    main()
