#!/usr/bin/env python3
"""Independent exact nontrivial-R83 exclusions for masks 13 and 29 (v52).

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
MASKS = (13, 29)
EXPONENTS = (1, 0, 1, 0, 0)

EXPECTED_BALANCING_UNIT = ambient.LElement((1, 3, -2, -6, 1, 2))
EXPECTED_BALANCING_UNIT_INVERSE = ambient.LElement((10, 54, -7, -37, 1, 6))
EXPECTED_CENTER_13 = ambient.LElement(
    (Q(-1846, 83), Q(-1437, 166), Q(2106, 83), Q(-307, 166), Q(-501, 83), Q(195, 166))
)
EXPECTED_HEIGHT_SQUARED = ambient.LElement(
    (
        Q(10064, 6889),
        Q(53407, 27556),
        Q(-7100, 6889),
        Q(-38025, 27556),
        Q(4561, 27556),
        Q(1568, 6889),
    )
)
EXPECTED_MODULUS_SQUARED = ambient.LElement(
    (
        Q(3644951, 6889),
        Q(4890422, 6889),
        Q(-2578147, 6889),
        Q(-3551387, 6889),
        Q(467266, 6889),
        Q(593566, 6889),
    )
)

EXPECTED_BOUNDS = {
    "DEN": (Q(29, 20), Q(31, 25), Q(179, 50), Q(103, 100), Q(31, 25), Q(81, 50)),
    "NUM": (Q(369, 100), Q(129, 50), Q(16, 25), Q(63, 50), Q(41, 20), Q(77, 25)),
    "A": (Q(4409, 100), Q(607, 50), Q(59, 100), Q(347, 50), Q(129, 20), Q(1781, 100)),
    "B": (Q(49857, 100), Q(1932, 25), Q(1267, 100), Q(658, 25), Q(6791, 100), Q(958, 5)),
    "C": (Q(432, 25), Q(291, 50), Q(167, 50), Q(57, 10), Q(391, 100), Q(233, 25)),
    "D": (Q(1953, 10), Q(741, 20), Q(1431, 20), Q(2163, 100), Q(411, 10), Q(401, 4)),
}

EXPECTED_COUNTS = (23917, 6625, 2043, 8768, 0)
EXPECTED_HISTOGRAMS = {
    13: (
        ((0, 0), 1772),
        ((0, 1), 5472),
        ((0, 2), 1106),
        ((0, 3), 328),
        ((0, 4), 18),
        ((0, 5), 6),
        ((1, 0), 44),
        ((1, 1), 10),
        ((1, 2), 2),
        ((1, 3), 2),
        ((1, 4), 8),
    ),
    29: (
        ((0, 0), 1712),
        ((0, 1), 5444),
        ((0, 2), 1092),
        ((0, 3), 394),
        ((0, 4), 18),
        ((0, 5), 14),
        ((1, 0), 64),
        ((1, 1), 24),
        ((1, 4), 6),
    ),
}
EXPECTED_DIGESTS = {
    13: "d66a04fa445505dd16cf2dad98ee6c9c57fe16d90aa4bbeb0b1bf2a1c2cc1a15",
    29: "b3ce2ca80967f6acf28bf10d6bb1707f330c768313680c369dd9d6d2cb3349de",
}
EXPECTED_AGGREGATE = "20ca8ce844583d40514ca6f65db26a302957b057fb0fd23c057cb7e52dab2e68"


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
    assert unit == xi5 * xi13 == EXPECTED_BALANCING_UNIT
    assert unit.inverse() == EXPECTED_BALANCING_UNIT_INVERSE
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

    for mask, partner in ((13, 61), (29, 45)):
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
    assert orientation == (orientation_unit(13) if mask == 13 else -orientation_unit(13))
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
    assert center == (EXPECTED_CENTER_13 if mask == 13 else -EXPECTED_CENTER_13)
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
    nonzero: tuple[Gamma, ...]
    zero: tuple[Gamma, ...]


def enumerate_candidates(
    point: Point, a_values: list[Vector], c_values: list[Vector]
) -> CandidateSet:
    """Exhaust every bounded determinant-one matrix at mask 13."""
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
                f"  masks13/29: c entries {c_index}/{len(c_values)}; "
                f"bottom rows={len(bottom_rows)}",
                flush=True,
            )

    nonzero: list[Gamma] = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != ZERO_VECTOR]
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
                f"  masks13/29: determinant rows {row_index}/{len(nonzero_rows)}; "
                f"candidates={len(nonzero)}",
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
            gamma = (a, b, c, d)
            ae, be, ce, de = map(orbit.vector_element, gamma)
            assert ce == orbit.ZERO and ae * de - be * ce == orbit.ONE
            zero.append(gamma)

    result = CandidateSet(len(bottom_rows), tuple(nonzero), tuple(zero))
    assert (result.bottom_rows, len(result.nonzero), len(result.zero)) == EXPECTED_COUNTS[2:]
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
    assert sum(count for _, count in result.histogram) == EXPECTED_COUNTS[3] + EXPECTED_COUNTS[4]
    assert result.digest == EXPECTED_DIGESTS[mask]
    return result


def main() -> None:
    assert orbit.field_discriminant() == 453789
    assert 13 ^ 48 == 61 and 29 ^ 48 == 45
    assert {13, 61}.isdisjoint({29, 45})

    ideal = conductor83_ideal()
    certify_root_pair_labels(ideal)
    print("INDEPENDENT d=7 NONTRIVIAL-R83 MASKS 13/29 STURM CHECK v52")
    print("stdlib field/Sturm/trace-LDL core; no primary certificate or experiment imports")
    print("canonical root pairs: {13,61} and {29,45}")
    print("saturated balancing exponents:", EXPONENTS)

    point13 = selected_point(13, ideal)
    point29 = selected_point(29, ideal)
    reflection = ambient.matrix(-1, 0, 0, 1)
    assert ambient.all_entries_integral(reflection)
    assert ambient.minverse(reflection) == reflection
    assert ambient.mdet(reflection) == -ambient.ONE
    assert ambient.sign_vector(
        ambient.mdet(reflection), ambient.embedded_roots()
    ) == (-1,) * orbit.DEGREE
    assert point29.matrix == reflect_matrix(point13.matrix)
    assert point29.matrix == ambient.conjugate_by(reflection, point13.matrix)
    print("diag(-1,1) lies in GL2(O_L), not GL2+(O_L) or SL2(O_L)")

    bounds13 = derive_bounds(point13)
    assert derive_bounds(point29) == bounds13 == EXPECTED_BOUNDS
    install_point_data(point13, bounds13)
    print("shared exact outward bounds:", bounds13)

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

    bounded13 = enumerate_candidates(point13, a_values, c_values)
    candidates13 = bounded13.nonzero + bounded13.zero
    candidates29 = validate_reflection(point13, point29, candidates13)
    print(f"exact reflected candidate bijection: {len(candidates13)}/{len(candidates29)}", flush=True)

    result13 = polygon_result(13, point13, candidates13)
    result29 = polygon_result(29, point29, candidates29)
    for mask, result in ((13, result13), (29, result29)):
        print(
            f"mask{mask}: bottom rows={bounded13.bottom_rows}; "
            f"c!=0 candidates={len(bounded13.nonzero)}; "
            f"c=0 candidates={len(bounded13.zero)}; polygon survivors=0"
        )
        print(f"mask{mask}: first-exclusion histogram={dict(result.histogram)}")
        print(f"mask{mask}: canonical exclusion SHA-256={result.digest}")

    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode("ascii")
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: independent exact exhaustive exclusions for distinct masks 13 and 29")


if __name__ == "__main__":
    main()
