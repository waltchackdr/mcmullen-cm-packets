#!/usr/bin/env python3
"""Independent exact exclusions of the nontrivial-R83 masks 1 and 17.

This verifier imports only the independent standard-library/Sturm ambient and
orbit foundations.  It reconstructs both canonical points from the
conductor-83 ideal lattice.  Their balanced matrices are conjugate by the
determinant-minus-one reflection ``diag(-1,1)``: the centers are negatives,
while the height, modulus, and all outward entry bounds agree.

The reflection normalizes ``SL_2(O_L)`` and gives an exact involution on the
finite bounded candidate set.  It is not in ``GL_2^+(O_L)`` or ``SL_2(O_L)``
and does not identify the two points.  Candidates are enumerated once at mask
1, reflected exactly, and tested separately against the six coordinate
quadrilaterals at the two distinct centers.  All decisive arithmetic is exact;
rational Sturm intervals certify every sign and bound.
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
MASKS = (1, 17)
EXPONENTS = (1, 0, 0, -1, -1)

EXPECTED_CENTER_1 = ambient.LElement(
    (Q(-603, 83), Q(575, 83), Q(4517, 166), Q(408, 83), Q(-742, 83), Q(-204, 83))
)
EXPECTED_HEIGHT_SQUARED = ambient.LElement(
    (Q(1497, 6889), Q(2438, 6889), Q(124, 6889), Q(-780, 6889), Q(-297, 27556), Q(58, 6889))
)
EXPECTED_MODULUS_SQUARED = ambient.LElement(
    (Q(231808, 6889), Q(-1674196, 6889), Q(-3166462, 6889), Q(327307, 6889), Q(1653404, 6889), Q(417388, 6889))
)

EXPECTED_BOUNDS = {
    "DEN": (Q(67, 20), Q(217, 100), Q(221, 100), Q(18, 25), Q(67, 100), Q(43, 25)),
    "NUM": (Q(171, 20), Q(113, 25), Q(39, 100), Q(22, 25), Q(11, 10), Q(329, 100)),
    "A": (Q(476, 25), Q(693, 100), Q(24, 25), Q(199, 20), Q(1197, 100), Q(1671, 100)),
    "B": (Q(115487, 100), Q(3388, 25), Q(39, 5), Q(917, 50), Q(3659, 100), Q(20423, 100)),
    "C": (Q(373, 50), Q(83, 25), Q(271, 50), Q(409, 50), Q(181, 25), Q(437, 50)),
    "D": (Q(45239, 100), Q(6497, 100), Q(2201, 50), Q(1507, 100), Q(443, 20), Q(5343, 50)),
}

EXPECTED_HISTOGRAMS = {
    1: (
        ((0, 0), 1870),
        ((0, 1), 5372),
        ((0, 2), 1044),
        ((0, 3), 344),
        ((0, 4), 28),
        ((0, 5), 12),
        ((1, 0), 50),
        ((1, 1), 16),
        ((1, 2), 6),
        ((1, 3), 4),
        ((1, 4), 2),
    ),
    17: (
        ((0, 0), 1872),
        ((0, 1), 5418),
        ((0, 2), 1008),
        ((0, 3), 340),
        ((0, 4), 16),
        ((0, 5), 14),
        ((1, 0), 48),
        ((1, 1), 18),
        ((1, 2), 6),
        ((1, 3), 6),
        ((1, 4), 2),
    ),
}
EXPECTED_DIGESTS = {
    1: "3d6a732cc36a03fe90933cdfcffd1b9d977e7ce22cb1d7ae1f3926e6dff2e1a8",
    17: "82b2d957ba62af939e5bb39b0f6bbec06b230094aeee60faa3ae54caf3e2e17b",
}
EXPECTED_AGGREGATE = "9a33ab59a829397a361196b67c9f6199be8288363a380fe506225d29496984dd"
EXPECTED_BALANCING_UNIT = ambient.LElement((1, 1, -3, -3, 1, 1))


def product(values: tuple[ambient.LElement, ...]) -> ambient.LElement:
    answer = ambient.ONE
    for value in values:
        answer *= value
    return answer


def orientation_units() -> tuple[ambient.LElement, ...]:
    """Return the independent v45 six-generator square-class basis."""
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
        tuple(
            unit
            for bit, unit in enumerate(orientation_units())
            if (mask >> bit) & 1
        )
    )


def balancing_unit() -> ambient.LElement:
    xi5, xi11, xi13, xi17, _, u = orientation_units()
    saturated_generators = (xi5, xi11, xi13, xi17, u)
    return product(
        tuple(
            generator**exponent
            for generator, exponent in zip(saturated_generators, EXPONENTS)
        )
    )


def certify_root_pair_labels() -> None:
    """Recheck the exact I83 xor-48 root intertwiners for both labels.

    The independently certified v45 packet theorem proves that these root
    pairings are free and complete.  Here the two relevant determinant-one
    intertwiners are replayed so the label bridge is locally inspectable.
    """
    alpha = ambient.alpha_pair(INDEX)
    maximal = ambient.multiplication_matrix(alpha)
    gamma83 = ambient.elt(-2, -2, 2, 0, -1)
    ideal_basis = ambient.matrix(16, gamma83, ambient.elt(58, 3), 3 * gamma83)
    ideal = ambient.conjugate_by(ambient.minverse(ideal_basis), maximal)
    u = 4 - ambient.T**2
    p = ambient.elt(55, 166, -18, -109, 0, 18)
    q = ambient.elt(-4, -9, -3, 5, 1, -1)
    r = ambient.elt(-603, -6206, 905, 4245, -184, -683)
    anti = ambient.matrix(p, q, r, -p)
    assert ambient.all_entries_integral(anti)
    assert ambient.mdet(anti) == -u
    assert ambient.mmul(anti, ideal) == ambient.mmul(
        ambient.root_conjugate(ideal), anti
    )

    for mask, partner_mask in ((1, 49), (17, 33)):
        source = ambient.twist(ideal, orientation_unit(mask))
        partner = ambient.twist(ideal, orientation_unit(partner_mask))
        intertwiner = ambient.mmul(
            ambient.mmul(
                ambient.diagonal(orientation_unit(partner_mask)),
                ambient.minverse(anti),
            ),
            ambient.minverse(ambient.diagonal(orientation_unit(mask))),
        )
        assert ambient.all_entries_integral(intertwiner)
        assert ambient.mdet(intertwiner) == ambient.ONE
        assert ambient.conjugate_by(
            intertwiner, ambient.root_conjugate(source)
        ) == partner


def reflect_matrix(value: ambient.Matrix) -> ambient.Matrix:
    """Conjugate by diag(-1,1), whose determinant is -1."""
    return ambient.matrix(value[0][0], -value[0][1], -value[1][0], value[1][1])


def selected_point(
    mask: int,
) -> tuple[ambient.Matrix, ambient.LElement, ambient.LElement, ambient.LElement]:
    """Reconstruct one balanced canonical point from independent v45 data."""
    assert mask in MASKS
    alpha = ambient.alpha_pair(INDEX)
    maximal = ambient.multiplication_matrix(alpha)
    assert ambient.pair_norm(alpha) == ambient.PI

    gamma83 = ambient.elt(-2, -2, 2, 0, -1)
    ideal_basis = ambient.matrix(16, gamma83, ambient.elt(58, 3), 3 * gamma83)
    ideal = ambient.conjugate_by(ambient.minverse(ideal_basis), maximal)
    assert ambient.all_entries_integral(ideal)

    xi5, _, _, _, _, _ = orientation_units()
    expected_orientation = xi5
    orientation = orientation_unit(mask)
    assert orientation == (expected_orientation if mask == 1 else -expected_orientation)
    oriented = ambient.twist(ideal, orientation)
    assert ambient.all_entries_integral(oriented)

    inversion = ambient.matrix(0, -1, 1, 0)
    assert ambient.mdet(inversion) == ambient.ONE
    inverted = ambient.conjugate_by(inversion, oriented)
    assert ambient.all_entries_integral(inverted)

    unit = balancing_unit()
    xi5, _, _, xi17, _, u = orientation_units()
    assert unit == xi5 * xi17.inverse() * u.inverse()
    assert unit == EXPECTED_BALANCING_UNIT
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
    expected_center = EXPECTED_CENTER_1 if mask == 1 else -EXPECTED_CENTER_1
    assert center == expected_center
    assert height_squared == EXPECTED_HEIGHT_SQUARED
    assert modulus_squared == EXPECTED_MODULUS_SQUARED
    assert ambient.sign_vector(height_squared, ambient.embedded_roots()) == (1,) * orbit.DEGREE
    assert ambient.sign_vector(modulus_squared, ambient.embedded_roots()) == (1,) * orbit.DEGREE
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
    """Derive all six outward entry bounds from exact Sturm intervals."""
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
        for target, value in zip(
            (den, num, a_bound, b_bound, c_bound, d_bound),
            (den_upper, num_upper, a_upper, b_upper, c_upper, d_upper),
        ):
            target.append(outward_hundredth(value))

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


def install_point_data(
    height_squared: ambient.LElement,
    bounds: dict[str, tuple[Q, ...]],
) -> None:
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
class CandidateSet:
    bottom_rows: int
    nonzero: tuple[Gamma, ...]
    zero: tuple[Gamma, ...]


@dataclass(frozen=True)
class PolygonResult:
    survivors: tuple[Gamma, ...]
    histogram: tuple[tuple[tuple[int, int], int], ...]
    digest: str


def enumerate_bounded_candidates(
    x: ambient.LElement,
    a_values: list[Vector],
    c_values: list[Vector],
) -> CandidateSet:
    """Enumerate every mask-1 bounded determinant-one matrix exactly once."""
    bottom_rows: list[tuple[Vector, Vector]] = []
    x_boxes = tuple(orbit.scaled_interval(orbit.interval_eval(x, root)) for root in orbit.ROOTS)
    a_box_cache = {
        a: tuple(orbit.scaled_vector_box(a, i) for i in range(orbit.DEGREE))
        for a in a_values
    }
    for c_index, c in enumerate(c_values):
        center = tuple(
            -(orbit.vector_element(c) * x).coefficients[i]
            for i in range(orbit.DEGREE)
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
                f"  masks1/17: c entries {c_index}/{len(c_values)}; "
                f"bottom rows={len(bottom_rows)}",
                flush=True,
            )
    print(f"  masks1/17: coupled bottom rows={len(bottom_rows)}", flush=True)

    nonzero: list[Gamma] = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != ZERO_VECTOR]
    for row_index, (c, d) in enumerate(nonzero_rows):
        affine = orbit.determinant_one_affine_map(c, d)
        for a in a_values:
            b = affine.evaluate(a)
            if b is not None and orbit.top_bounded(a, b, x, a_box_cache[a], x_boxes):
                gamma = (a, b, c, d)
                assert (
                    orbit.vector_element(a) * orbit.vector_element(d)
                    - orbit.vector_element(b) * orbit.vector_element(c)
                    == orbit.ONE
                )
                nonzero.append(gamma)
        if row_index and row_index % 250 == 0:
            print(
                f"  masks1/17: determinant rows {row_index}/{len(nonzero_rows)}; "
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
        center = tuple(
            -(a_element * x).coefficients[i] for i in range(orbit.DEGREE)
        )
        for b in orbit.enumerate_trace_ball(
            center,
            orbit.NUM_RADIUS,
            lambda b, a=a: orbit.top_bounded(a, b, x, a_box_cache[a], x_boxes),
        ):
            assert a_element * orbit.vector_element(d) == orbit.ONE
            zero.append((a, b, c, d))

    result = CandidateSet(len(bottom_rows), tuple(nonzero), tuple(zero))
    assert result.bottom_rows == 2017
    assert len(result.nonzero) == 8748
    assert len(result.zero) == 0
    return result


def negate_vector(value: Vector) -> Vector:
    return tuple(-coefficient for coefficient in value)


def reflected_gamma(gamma: Gamma) -> Gamma:
    """Exact candidate involution induced by diag(-1,1)."""
    a, b, c, d = gamma
    return a, negate_vector(b), negate_vector(c), d


def validate_reflection_bijection(
    x1: ambient.LElement,
    x17: ambient.LElement,
    height_squared: ambient.LElement,
    candidates1: tuple[Gamma, ...],
    candidates17: tuple[Gamma, ...],
) -> None:
    """Check the exact candidate involution and its complete converse.

    If ``gamma17`` is any bounded integral determinant-one matrix at ``x17``,
    its reflection is again integral and determinant one.  The two norm
    identities below, applied with the roles reversed, put that reflection
    inside the same symmetric A/B/C/D bounds at ``x1``.  Exhaustiveness of
    the mask-1 enumeration therefore puts it in ``candidates1``.  Since the
    reflection is an involution, every mask-17 candidate is in the displayed
    image.  The argument includes ``c=0``, which reflection preserves.
    """
    assert x17 == -x1
    assert len(candidates1) == len(candidates17)
    assert len(set(candidates1)) == len(candidates1)
    assert len(set(candidates17)) == len(candidates17)
    for gamma1, gamma17 in zip(candidates1, candidates17):
        assert reflected_gamma(gamma1) == gamma17
        assert reflected_gamma(gamma17) == gamma1
        a, b, c, d = map(orbit.vector_element, gamma1)
        ar, br, cr, dr = map(orbit.vector_element, gamma17)
        assert a * d - b * c == orbit.ONE
        assert ar * dr - br * cr == orbit.ONE
        assert (c * x1 + d) ** 2 + c**2 * height_squared == (
            cr * x17 + dr
        ) ** 2 + cr**2 * height_squared
        assert (a * x1 + b) ** 2 + a**2 * height_squared == (
            ar * x17 + br
        ) ** 2 + ar**2 * height_squared


def polygon_result(mask: int, x: ambient.LElement, candidates: tuple[Gamma, ...]) -> PolygonResult:
    """Apply the exact polygon predicates separately at one center."""
    survivors: list[Gamma] = []
    records: list[tuple[Gamma, tuple[int, int]]] = []
    histogram_data: dict[tuple[int, int], int] = {}
    for gamma in candidates:
        failure = orbit.first_polygon_exclusion(gamma, x)
        if failure is None:
            survivors.append(gamma)
            continue
        records.append((gamma, failure))
        histogram_data[failure] = histogram_data.get(failure, 0) + 1
    result = PolygonResult(
        tuple(survivors),
        tuple(sorted(histogram_data.items())),
        orbit.canonical_digest(records),
    )
    assert result.survivors == ()
    assert result.histogram == EXPECTED_HISTOGRAMS[mask]
    assert sum(count for _, count in result.histogram) == 8748
    assert result.digest == EXPECTED_DIGESTS[mask]
    return result


def main() -> None:
    assert orbit.field_discriminant() == 453789
    assert 1 ^ 48 == 49 and 17 ^ 48 == 33
    assert {1, 49}.isdisjoint({17, 33})
    prior_fibers = (
        {14, 62},
        {30, 46},
        {15, 63},
        {31, 47},
        {6, 54},
        {22, 38},
    )
    assert all(
        {1, 49}.isdisjoint(fiber) and {17, 33}.isdisjoint(fiber)
        for fiber in prior_fibers
    )
    certify_root_pair_labels()
    print("INDEPENDENT d=7 NONTRIVIAL-R83 MASKS 1/17 STURM CHECK v50")
    print("stdlib field/Sturm/LDL core; no primary certificate or experiment imports")
    print("canonical root pairs: mask 1 xor 48 = 49; mask 17 xor 48 = 33")
    print("saturated balancing exponents:", EXPONENTS)

    points = {mask: selected_point(mask) for mask in MASKS}
    matrix1, center1, height1, modulus1 = points[1]
    matrix17, center17, height17, modulus17 = points[17]
    reflection = ambient.matrix(-1, 0, 0, 1)
    assert ambient.mdet(reflection) == -ambient.ONE
    assert ambient.sign_vector(
        ambient.mdet(reflection), ambient.embedded_roots()
    ) == (-1,) * orbit.DEGREE
    assert ambient.all_entries_integral(reflection)
    assert ambient.minverse(reflection) == reflection
    assert matrix17 == reflect_matrix(matrix1)
    assert matrix17 == ambient.conjugate_by(reflection, matrix1)
    assert center17 == -center1
    assert height17 == height1
    assert modulus17 == modulus1
    print(
        "matrix reflection: diag(-1,1) lies in GL2(O_L) but not "
        "GL2+(O_L) or SL2(O_L)"
    )

    bounds1 = derive_bounds(height1, modulus1)
    bounds17 = derive_bounds(height17, modulus17)
    assert bounds17 == bounds1 == EXPECTED_BOUNDS
    install_point_data(height1, bounds1)
    print("shared exact outward bounds:", bounds1)

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
    assert (len(a_values), len(c_values)) == (23959, 6601)
    assert {negate_vector(value) for value in a_values} == set(a_values)
    assert {negate_vector(value) for value in c_values} == set(c_values)
    print("single exact a/c enumeration:", len(a_values), len(c_values), flush=True)

    bounded1 = enumerate_bounded_candidates(center1, a_values, c_values)
    candidates1 = bounded1.nonzero + bounded1.zero
    candidates17 = tuple(reflected_gamma(gamma) for gamma in candidates1)
    validate_reflection_bijection(
        center1, center17, height1, candidates1, candidates17
    )
    print("exact reflected candidate bijection: 8748/8748", flush=True)

    result1 = polygon_result(1, center1, candidates1)
    result17 = polygon_result(17, center17, candidates17)
    for mask, result in ((1, result1), (17, result17)):
        print(
            f"mask{mask}: bottom rows={bounded1.bottom_rows}; "
            f"c!=0 candidates={len(bounded1.nonzero)}; "
            f"c=0 candidates={len(bounded1.zero)}; polygon survivors=0"
        )
        print(f"mask{mask}: first-exclusion histogram={dict(result.histogram)}")
        print(f"mask{mask}: canonical exclusion SHA-256={result.digest}")
    aggregate = sha256(
        "".join(f"{mask}:{EXPECTED_DIGESTS[mask]}\n" for mask in MASKS).encode(
            "ascii"
        )
    ).hexdigest()
    assert aggregate == EXPECTED_AGGREGATE
    print("aggregate labelled digest:", aggregate)
    print("PASS: independent exact exhaustive exclusions for distinct masks 1 and 17")


if __name__ == "__main__":
    main()
