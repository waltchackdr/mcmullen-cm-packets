#!/usr/bin/env python3
"""Independent exact nontrivial-R83 exclusions for masks 7 and 23 (v53).

This verifier uses only the independent standard-library/Sturm ambient-packet
and orbit foundations.  It reconstructs the two conductor-83 points, derives
their compact bounds by rational interval arithmetic, exhausts the positive
trace-form lattice searches, solves ``ad-bc=1`` exactly (including ``c=0``),
and applies the six polygon tests separately at both centers.  It imports no
experiment or primary certificate and uses no floating-point comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import floor

import d7_ambient_points_sturm_v45 as ambient
import residual_orbit_k17_k38_sturm_v42 as orbit


Q = Fraction
Vector = tuple[int, ...]
Gamma = tuple[Vector, Vector, Vector, Vector]
ZERO_VECTOR: Vector = (0,) * orbit.DEGREE
INDEX = 2
MASKS = (7, 23)
EXPONENTS = (0, 0, 0, -1, -2)

EXPECTED_BALANCING_UNIT = ambient.LElement((0, -2, -4, -1, 2, 1))
EXPECTED_BALANCING_UNIT_INVERSE = ambient.LElement((1, -25, 0, 18, 0, -3))
EXPECTED_CENTER_7 = ambient.LElement(
    (Q(269, 83), Q(677, 83), Q(-2106, 83), Q(307, 166), Q(501, 83), Q(-195, 166))
)
EXPECTED_HEIGHT_SQUARED = ambient.LElement(
    (
        Q(-1037, 13778),
        Q(-11997, 27556),
        Q(19989, 27556),
        Q(487, 27556),
        Q(-5233, 27556),
        Q(397, 13778),
    )
)
EXPECTED_MODULUS_SQUARED = ambient.LElement(
    (
        Q(290921, 6889),
        Q(2520855, 6889),
        Q(3953621, 6889),
        Q(-3802545, 6889),
        Q(-1079522, 6889),
        Q(850036, 6889),
    )
)

EXPECTED_BOUNDS = {
    "DEN": (Q(32, 25), Q(46, 25), Q(24, 25), Q(21, 10), Q(103, 50), Q(34, 25)),
    "NUM": (Q(327, 100), Q(383, 100), Q(17, 100), Q(64, 25), Q(341, 100), Q(259, 100)),
    "A": (Q(499, 10), Q(817, 100), Q(221, 100), Q(17, 5), Q(387, 100), Q(1061, 50)),
    "B": (Q(44057, 100), Q(5741, 50), Q(339, 100), Q(5367, 100), Q(11321, 100), Q(16081, 100)),
    "C": (Q(391, 20), Q(98, 25), Q(312, 25), Q(14, 5), Q(117, 50), Q(111, 10)),
    "D": (Q(8629, 50), Q(1101, 20), Q(191, 10), Q(441, 10), Q(6851, 100), Q(4207, 50)),
}

# a entries, c entries, coupled bottom rows, c!=0 candidates, c=0 candidates.
EXPECTED_COUNTS = (23889, 6609, 1955, 8514, 0)
EXPECTED_NONZERO_BOTTOM_ROWS = 1954
EXPECTED_ZERO_BOTTOM_ROWS = 1
EXPECTED_HISTOGRAMS = {
    7: (
        ((0, 0), 1714),
        ((0, 1), 5302),
        ((0, 2), 942),
        ((0, 3), 396),
        ((0, 4), 26),
        ((0, 5), 20),
        ((1, 0), 78),
        ((1, 1), 18),
        ((1, 2), 8),
        ((1, 3), 4),
        ((1, 4), 6),
    ),
    23: (
        ((0, 0), 1680),
        ((0, 1), 5402),
        ((0, 2), 912),
        ((0, 3), 374),
        ((0, 4), 20),
        ((0, 5), 22),
        ((1, 0), 80),
        ((1, 1), 20),
        ((1, 3), 2),
        ((1, 4), 2),
    ),
}
EXPECTED_DIGESTS = {
    7: "b653b69bf5888c06df2ee9a20e5f6fd8d1d06277de0b862933ba5a4bf8811ae5",
    23: "392b61bb9ec4dcfd6720de02184032cfe211f9094b3e5ed1dac5224be84f8d52",
}
EXPECTED_AGGREGATE = "aba4002201230fea6eb9e127aa6fa38a04e63194d3970adc32c95a1a87386dab"


def product(values: tuple[ambient.LElement, ...]) -> ambient.LElement:
    answer = ambient.ONE
    for value in values:
        answer *= value
    return answer


def orientation_units() -> tuple[ambient.LElement, ...]:
    """Return the independently frozen v45 orientation square-class basis."""
    return (
        ambient.T**4 - 3 * ambient.T**2 + 1,
        ambient.T**3 + ambient.T**2 - 2 * ambient.T - 1,
        ambient.T**5 - 4 * ambient.T**3 + 3 * ambient.T + 1,
        ambient.T**3 - 2 * ambient.T,
        -ambient.ONE,
        4 - ambient.T**2,
    )


def orientation_unit(mask: int) -> ambient.LElement:
    return product(
        tuple(unit for bit, unit in enumerate(orientation_units()) if (mask >> bit) & 1)
    )


def balancing_unit() -> ambient.LElement:
    xi5, xi11, xi13, xi17, _, u = orientation_units()
    unit = product(
        tuple(
            generator**exponent
            for generator, exponent in zip((xi5, xi11, xi13, xi17, u), EXPONENTS)
        )
    )
    assert unit == xi17.inverse() * u.inverse() ** 2
    assert unit == EXPECTED_BALANCING_UNIT
    assert unit.integral_power_basis()
    inverse = unit.inverse()
    assert inverse == EXPECTED_BALANCING_UNIT_INVERSE
    assert inverse.integral_power_basis()
    return unit


def conductor83_ideal() -> ambient.Matrix:
    alpha = ambient.alpha_pair(INDEX)
    maximal = ambient.multiplication_matrix(alpha)
    gamma83 = ambient.elt(-2, -2, 2, 0, -1)
    basis = ambient.matrix(16, gamma83, ambient.elt(58, 3), 3 * gamma83)
    ideal = ambient.conjugate_by(ambient.minverse(basis), maximal)
    assert ambient.all_entries_integral(ideal)
    return ideal


def certify_root_pair_labels(ideal: ambient.Matrix) -> None:
    """Replay the exact determinant-one free-root bridges for both labels."""
    u = 4 - ambient.T**2
    p = ambient.elt(55, 166, -18, -109, 0, 18)
    q = ambient.elt(-4, -9, -3, 5, 1, -1)
    r = ambient.elt(-603, -6206, 905, 4245, -184, -683)
    anti = ambient.matrix(p, q, r, -p)
    assert ambient.all_entries_integral(anti)
    assert ambient.mdet(anti) == -u
    assert ambient.mmul(anti, ideal) == ambient.mmul(ambient.root_conjugate(ideal), anti)

    for mask, partner in ((7, 55), (23, 39)):
        assert partner == mask ^ 48
        source = ambient.twist(ideal, orientation_unit(mask))
        target = ambient.twist(ideal, orientation_unit(partner))
        bridge = ambient.mmul(
            ambient.mmul(
                ambient.diagonal(orientation_unit(partner)), ambient.minverse(anti)
            ),
            ambient.minverse(ambient.diagonal(orientation_unit(mask))),
        )
        assert ambient.all_entries_integral(bridge)
        assert ambient.mdet(bridge) == ambient.ONE
        assert ambient.conjugate_by(bridge, ambient.root_conjugate(source)) == target


@dataclass(frozen=True)
class Point:
    matrix: ambient.Matrix
    center: ambient.LElement
    height_squared: ambient.LElement
    modulus_squared: ambient.LElement


def selected_point(mask: int, ideal: ambient.Matrix) -> Point:
    """Reconstruct one balanced canonical point from independent v45 data."""
    assert mask in MASKS
    alpha = ambient.alpha_pair(INDEX)
    assert ambient.pair_norm(alpha) == ambient.PI

    orientation = orientation_unit(mask)
    assert orientation == (orientation_unit(7) if mask == 7 else -orientation_unit(7))
    raw = ambient.twist(ideal, orientation)
    assert ambient.all_entries_integral(raw)
    raw_x, raw_y2 = ambient.exact_fixed_invariants(raw)
    raw_modulus2 = raw_x**2 + raw_y2
    assert Q(0) < raw_modulus2.norm() < Q(1)

    inversion = ambient.matrix(0, -1, 1, 0)
    assert ambient.mdet(inversion) == ambient.ONE
    inverted = ambient.conjugate_by(inversion, raw)
    assert ambient.all_entries_integral(inverted)
    inv_x, inv_y2 = ambient.exact_fixed_invariants(inverted)
    inv_modulus2 = inv_x**2 + inv_y2
    assert inv_x == -raw_x / raw_modulus2
    assert inv_y2 == raw_y2 / raw_modulus2**2
    assert inv_modulus2 == ambient.ONE / raw_modulus2

    unit = balancing_unit()
    conjugator = ambient.matrix(unit, 0, 0, unit.inverse())
    assert ambient.mdet(conjugator) == ambient.ONE
    assert ambient.all_entries_integral(conjugator)
    assert ambient.all_entries_integral(ambient.minverse(conjugator))
    balanced = ambient.conjugate_by(conjugator, inverted)
    assert ambient.all_entries_integral(balanced)
    assert ambient.mtrace(balanced) == ambient.pair_trace(alpha)
    assert ambient.mdet(balanced) == ambient.PI

    center, height_squared = ambient.exact_fixed_invariants(balanced)
    modulus_squared = center**2 + height_squared
    square = unit**2
    fourth = square**2
    assert center == square * inv_x
    assert height_squared == fourth * inv_y2
    assert modulus_squared == fourth * inv_modulus2
    assert center == (EXPECTED_CENTER_7 if mask == 7 else -EXPECTED_CENTER_7)
    assert height_squared == EXPECTED_HEIGHT_SQUARED
    assert modulus_squared == EXPECTED_MODULUS_SQUARED
    assert ambient.sign_vector(height_squared, ambient.embedded_roots()) == (1,) * orbit.DEGREE
    assert ambient.sign_vector(modulus_squared, ambient.embedded_roots()) == (1,) * orbit.DEGREE
    assert height_squared.norm() == Q(1, 28217344)
    assert modulus_squared.norm() == Q(2097130411225081, 6889)
    return Point(balanced, center, height_squared, modulus_squared)


def reflect_matrix(value: ambient.Matrix) -> ambient.Matrix:
    return ambient.matrix(value[0][0], -value[0][1], -value[1][0], value[1][1])


def outward_hundredth(value: Q) -> Q:
    """Return the smallest hundredth strictly above an exact upper endpoint."""
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(point: Point) -> dict[str, tuple[Q, ...]]:
    """Derive all six compact entry bounds with rational Sturm intervals."""
    columns: list[list[Q]] = [[] for _ in range(6)]
    for embedding, root in enumerate(orbit.ROOTS):
        height = orbit.sqrt_interval(orbit.interval_eval(point.height_squared, root), 160)
        modulus = orbit.sqrt_interval(orbit.interval_eval(point.modulus_squared, root), 160)
        assert height.lo > 0 and modulus.lo > 0
        den = orbit.sqrt_interval(
            orbit.Interval.point(height.hi / orbit.M_LOWER[embedding]), 160
        ).hi
        c_entry = orbit.sqrt_interval(
            orbit.Interval.point(1 / (orbit.M_LOWER[embedding] * height.lo)), 160
        ).hi
        num = orbit.U_UPPER[embedding] * den
        a_entry = num / height.lo
        d_entry = den + c_entry * modulus.hi
        b_entry = num + a_entry * modulus.hi
        for target, value in zip(columns, (den, num, a_entry, b_entry, c_entry, d_entry)):
            target.append(outward_hundredth(value))
    result = dict(zip(("DEN", "NUM", "A", "B", "C", "D"), map(tuple, columns)))
    assert result == EXPECTED_BOUNDS
    return result


def install_point_data(point: Point, bounds: dict[str, tuple[Q, ...]]) -> None:
    orbit.Y2 = point.height_squared
    orbit.DEN_BOUND = bounds["DEN"]
    orbit.NUM_BOUND = bounds["NUM"]
    orbit.A_BOUND = bounds["A"]
    orbit.B_BOUND = bounds["B"]
    orbit.C_BOUND = bounds["C"]
    orbit.D_BOUND = bounds["D"]
    orbit.DEN_RADIUS = orbit.sum_of_squares(orbit.DEN_BOUND)
    orbit.NUM_RADIUS = orbit.sum_of_squares(orbit.NUM_BOUND)
    orbit.Y2_SCALED = tuple(
        orbit.scaled_interval(orbit.interval_eval(point.height_squared, root))
        for root in orbit.ROOTS
    )


@dataclass(frozen=True)
class CandidateSet:
    bottom_rows: int
    nonzero_bottom_rows: int
    zero_bottom_rows: int
    nonzero: tuple[Gamma, ...]
    zero: tuple[Gamma, ...]


def enumerate_candidates(
    point: Point, a_values: list[Vector], c_values: list[Vector]
) -> CandidateSet:
    """Exhaust every bounded determinant-one matrix at mask 7."""
    x = point.center
    x_boxes = tuple(
        orbit.scaled_interval(orbit.interval_eval(x, root)) for root in orbit.ROOTS
    )
    a_box_cache = {
        a: tuple(orbit.scaled_vector_box(a, i) for i in range(orbit.DEGREE))
        for a in a_values
    }

    bottom_rows: list[tuple[Vector, Vector]] = []
    for c_index, c in enumerate(c_values):
        center = tuple(
            -(orbit.vector_element(c) * x).coefficients[i] for i in range(orbit.DEGREE)
        )
        c_boxes = tuple(orbit.scaled_vector_box(c, i) for i in range(orbit.DEGREE))
        bottom_rows.extend(
            (c, d)
            for d in orbit.enumerate_trace_ball(
                center,
                orbit.DEN_RADIUS,
                lambda d, c=c, c_boxes=c_boxes: orbit.bottom_bounded(
                    c, d, x, c_boxes, x_boxes
                ),
            )
        )
        if c_index and c_index % 500 == 0:
            print(
                f"  masks7/23: c entries {c_index}/{len(c_values)}; "
                f"bottom rows={len(bottom_rows)}",
                flush=True,
            )

    nonzero_rows = [(c, d) for c, d in bottom_rows if c != ZERO_VECTOR]
    zero_rows = [(c, d) for c, d in bottom_rows if c == ZERO_VECTOR]
    assert len(bottom_rows) == EXPECTED_COUNTS[2]
    assert len(nonzero_rows) == EXPECTED_NONZERO_BOTTOM_ROWS
    assert len(zero_rows) == EXPECTED_ZERO_BOTTOM_ROWS
    assert zero_rows == [(ZERO_VECTOR, ZERO_VECTOR)]

    nonzero: list[Gamma] = []
    for row_index, (c, d) in enumerate(nonzero_rows):
        affine = orbit.determinant_one_affine_map(c, d)
        for a in a_values:
            b = affine.evaluate(a)
            if b is not None and orbit.top_bounded(a, b, x, a_box_cache[a], x_boxes):
                gamma = (a, b, c, d)
                ae, be, ce, de = map(orbit.vector_element, gamma)
                assert ae * de - be * ce == orbit.ONE
                nonzero.append(gamma)
        if row_index and row_index % 250 == 0:
            print(
                f"  masks7/23: determinant rows {row_index}/{len(nonzero_rows)}; "
                f"candidates={len(nonzero)}",
                flush=True,
            )

    zero: list[Gamma] = []
    for c, d in zero_rows:
        assert c == ZERO_VECTOR
        if d == ZERO_VECTOR:
            # Then ad-bc=0 for every a,b, so determinant one is impossible.
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
            gamma = (a, b, c, d)
            ae, be, ce, de = map(orbit.vector_element, gamma)
            assert ce == orbit.ZERO and ae * de - be * ce == orbit.ONE
            zero.append(gamma)

    result = CandidateSet(
        len(bottom_rows), len(nonzero_rows), len(zero_rows), tuple(nonzero), tuple(zero)
    )
    assert (len(result.nonzero), len(result.zero)) == EXPECTED_COUNTS[3:]
    return result


def negate_vector(value: Vector) -> Vector:
    return tuple(-coefficient for coefficient in value)


def reflected_gamma(gamma: Gamma) -> Gamma:
    a, b, c, d = gamma
    return a, negate_vector(b), negate_vector(c), d


def validate_reflection(
    base: Point, reflected: Point, candidates: tuple[Gamma, ...]
) -> tuple[Gamma, ...]:
    """Prove the exact bounded-candidate bijection, not a quotient identity."""
    assert reflected.center == -base.center
    assert reflected.height_squared == base.height_squared
    assert reflected.modulus_squared == base.modulus_squared
    targets = tuple(reflected_gamma(gamma) for gamma in candidates)
    assert len(set(candidates)) == len(candidates)
    assert len(set(targets)) == len(targets)
    for source, target in zip(candidates, targets):
        assert reflected_gamma(target) == source
        a, b, c, d = map(orbit.vector_element, source)
        ar, br, cr, dr = map(orbit.vector_element, target)
        assert a * d - b * c == orbit.ONE
        assert ar * dr - br * cr == orbit.ONE
        assert (c * base.center + d) ** 2 + c**2 * base.height_squared == (
            (cr * reflected.center + dr) ** 2 + cr**2 * reflected.height_squared
        )
        assert (a * base.center + b) ** 2 + a**2 * base.height_squared == (
            (ar * reflected.center + br) ** 2 + ar**2 * reflected.height_squared
        )
    return targets


@dataclass(frozen=True)
class PolygonResult:
    survivors: tuple[Gamma, ...]
    histogram: tuple[tuple[tuple[int, int], int], ...]
    digest: str


def polygon_result(mask: int, point: Point, candidates: tuple[Gamma, ...]) -> PolygonResult:
    """Run the exact polygon predicate separately at one point center."""
    survivors: list[Gamma] = []
    records: list[tuple[Gamma, tuple[int, int]]] = []
    histogram: dict[tuple[int, int], int] = {}
    for gamma in candidates:
        failure = orbit.first_polygon_exclusion(gamma, point.center)
        if failure is None:
            survivors.append(gamma)
        else:
            records.append((gamma, failure))
            histogram[failure] = histogram.get(failure, 0) + 1
    result = PolygonResult(
        tuple(survivors), tuple(sorted(histogram.items())), orbit.canonical_digest(records)
    )
    assert result.survivors == ()
    assert result.histogram == EXPECTED_HISTOGRAMS[mask]
    assert sum(count for _, count in result.histogram) == sum(EXPECTED_COUNTS[3:])
    assert result.digest == EXPECTED_DIGESTS[mask]
    return result


def certify_disjoint_labels() -> None:
    new_fibers = (frozenset((7, 55)), frozenset((23, 39)))
    assert new_fibers[0].isdisjoint(new_fibers[1])
    prior_masks = (
        14,
        30,
        15,
        31,
        6,
        22,
        1,
        17,
        0,
        16,
        2,
        18,
        3,
        19,
        4,
        20,
        5,
        21,
        8,
        24,
        9,
        25,
        10,
        26,
        11,
        27,
        12,
        28,
        13,
        29,
    )
    prior_fibers = tuple(frozenset((mask, mask ^ 48)) for mask in prior_masks)
    assert all(new.isdisjoint(old) for new in new_fibers for old in prior_fibers)


def main() -> None:
    assert orbit.field_discriminant() == 453789
    certify_disjoint_labels()
    ideal = conductor83_ideal()
    certify_root_pair_labels(ideal)
    print("INDEPENDENT d=7 NONTRIVIAL-R83 MASKS 7/23 STURM CHECK v53")
    print("stdlib field/Sturm/trace-LDL core; no primary certificate or experiment imports")
    print("canonical root pairs: {7,55} and {23,39}")
    print("saturated balancing exponents:", EXPONENTS)

    point7 = selected_point(7, ideal)
    point23 = selected_point(23, ideal)
    reflection = ambient.matrix(-1, 0, 0, 1)
    assert ambient.all_entries_integral(reflection)
    assert ambient.minverse(reflection) == reflection
    assert ambient.mdet(reflection) == -ambient.ONE
    assert ambient.sign_vector(
        ambient.mdet(reflection), ambient.embedded_roots()
    ) == (-1,) * orbit.DEGREE
    assert point23.matrix == reflect_matrix(point7.matrix)
    assert point23.matrix == ambient.conjugate_by(reflection, point7.matrix)
    print("diag(-1,1) lies in GL2(O_L), not GL2+(O_L) or SL2(O_L)")

    bounds7 = derive_bounds(point7)
    assert derive_bounds(point23) == bounds7 == EXPECTED_BOUNDS
    install_point_data(point7, bounds7)
    print("shared exact outward bounds:", bounds7)

    a_values = orbit.enumerate_trace_ball(
        (Q(0),) * orbit.DEGREE,
        orbit.sum_of_squares(orbit.A_BOUND),
        lambda vector: orbit.field_box_membership(vector, orbit.A_BOUND),
    )
    c_values = orbit.enumerate_trace_ball(
        (Q(0),) * orbit.DEGREE,
        orbit.sum_of_squares(orbit.C_BOUND),
        lambda vector: orbit.field_box_membership(vector, orbit.C_BOUND),
    )
    a_values.sort()
    c_values.sort()
    assert (len(a_values), len(c_values)) == EXPECTED_COUNTS[:2]
    assert {negate_vector(value) for value in a_values} == set(a_values)
    assert {negate_vector(value) for value in c_values} == set(c_values)
    print("single exact a/c enumeration:", len(a_values), len(c_values), flush=True)

    bounded7 = enumerate_candidates(point7, a_values, c_values)
    candidates7 = bounded7.nonzero + bounded7.zero
    candidates23 = validate_reflection(point7, point23, candidates7)
    print(f"exact reflected candidate bijection: {len(candidates7)}/{len(candidates23)}", flush=True)

    result7 = polygon_result(7, point7, candidates7)
    result23 = polygon_result(23, point23, candidates23)
    for mask, result in ((7, result7), (23, result23)):
        print(
            f"mask{mask}: bottom rows={bounded7.bottom_rows} "
            f"({bounded7.nonzero_bottom_rows} c!=0 + {bounded7.zero_bottom_rows} c=0); "
            f"c!=0 candidates={len(bounded7.nonzero)}; "
            f"c=0 candidates={len(bounded7.zero)}; polygon survivors=0"
        )
        print(f"mask{mask}: first-exclusion histogram={dict(result.histogram)}")
        print(f"mask{mask}: canonical exclusion SHA-256={result.digest}")

    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode("ascii")
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: independent exact exhaustive exclusions for distinct masks 7 and 23")


if __name__ == "__main__":
    main()
