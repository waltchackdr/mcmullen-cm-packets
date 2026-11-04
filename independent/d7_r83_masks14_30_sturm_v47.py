#!/usr/bin/env python3
"""Independent exact exclusions of the nontrivial-R83 masks 14 and 30.

This verifier imports only the independent standard-library/Sturm ambient and
orbit foundations.  It reconstructs both canonical points from the
conductor-83 ideal lattice.  Their balanced matrices are conjugate by the
determinant-minus-one reflection ``diag(-1,1)``: the centers are negatives,
while the height, modulus, and all outward entry bounds agree.

The reflection normalizes ``SL_2(O_L)`` and gives an exact involution on the
finite bounded candidate set.  It is not asserted to lie in ``SL_2(O_L)`` or
to identify the two ambient points.  Candidates are enumerated once at mask
14, reflected exactly, and tested separately against the six coordinate
polygons at the two distinct centers.  All decisive arithmetic is exact;
rational Sturm intervals certify every sign and bound.
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
INDEX = 2
MASKS = (14, 30)
EXPONENTS = (0, 0, 0, 0, -2)

EXPECTED_CENTER_14 = ambient.LElement(
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

EXPECTED_HISTOGRAMS = {
    14: (
        ((0, 0), 1802),
        ((0, 1), 5676),
        ((0, 2), 974),
        ((0, 3), 468),
        ((0, 4), 30),
        ((0, 5), 6),
        ((1, 0), 66),
        ((1, 1), 12),
        ((1, 2), 2),
    ),
    30: (
        ((0, 0), 1798),
        ((0, 1), 5640),
        ((0, 2), 982),
        ((0, 3), 500),
        ((0, 4), 14),
        ((0, 5), 6),
        ((1, 0), 64),
        ((1, 1), 20),
        ((1, 2), 6),
        ((1, 3), 2),
        ((1, 4), 4),
    ),
}
EXPECTED_DIGESTS = {
    14: "19f1bba6f56200bac45a52810a56669c51824f22a9b3fd410db6a1fc2d69a5a1",
    30: "2334a40fed0d008d44c2e892871bec20fc04b226a99df899a8927df81595faf6",
}


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


def certify_root_pair_labels() -> None:
    """Recheck the exact I83 xor-48 root intertwiners for both labels.

    The independently certified v45 packet theorem proves that these root
    pairings are free and complete.  Here we replay the two relevant
    determinant-one intertwiners so the point-label bridge is inspectable in
    this certificate as well.
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

    for mask, partner_mask in ((14, 62), (30, 46)):
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

    xi5, xi11, xi13, xi17, _, u = orientation_units()
    del xi5
    orientation = orientation_unit(mask)
    expected_orientation = xi11 * xi13 * xi17
    assert orientation == (expected_orientation if mask == 14 else -expected_orientation)
    oriented = ambient.twist(ideal, orientation)
    assert ambient.all_entries_integral(oriented)

    inversion = ambient.matrix(0, -1, 1, 0)
    assert ambient.mdet(inversion) == ambient.ONE
    inverted = ambient.conjugate_by(inversion, oriented)
    assert ambient.all_entries_integral(inverted)

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
    expected_center = EXPECTED_CENTER_14 if mask == 14 else -EXPECTED_CENTER_14
    assert center == expected_center
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
    """Enumerate every mask-14 bounded determinant-one matrix exactly once."""
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
                f"  masks14/30: c entries {c_index}/{len(c_values)}; "
                f"bottom rows={len(bottom_rows)}",
                flush=True,
            )
    print(f"  masks14/30: coupled bottom rows={len(bottom_rows)}", flush=True)

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
                f"  masks14/30: determinant rows {row_index}/{len(nonzero_rows)}; "
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
    assert result.bottom_rows == 2093
    assert len(result.nonzero) == 9036
    assert len(result.zero) == 0
    return result


def negate_vector(value: Vector) -> Vector:
    return tuple(-coefficient for coefficient in value)


def reflected_gamma(gamma: Gamma) -> Gamma:
    """Exact candidate involution induced by diag(-1,1)."""
    a, b, c, d = gamma
    return a, negate_vector(b), negate_vector(c), d


def validate_reflection_bijection(
    x14: ambient.LElement,
    x30: ambient.LElement,
    height_squared: ambient.LElement,
    candidates14: tuple[Gamma, ...],
    candidates30: tuple[Gamma, ...],
) -> None:
    """Check determinant and bound identities for every reflected candidate.

    Completeness in the reverse direction is structural.  If ``gamma30`` is
    any bounded integral determinant-one matrix at ``x30``, its reflection is
    again integral and determinant one.  The two asserted norm identities,
    applied with the roles reversed, put that reflection inside the same
    symmetric A/B/C/D bounds at ``x14``.  The exhaustive mask-14 enumeration
    must therefore contain it.  Since reflection is an involution, every
    mask-30 candidate is in the displayed image.  The argument also covers
    ``c=0``, because reflection preserves the equation ``c=0``.
    """
    assert x30 == -x14
    assert len(candidates14) == len(candidates30)
    assert len(set(candidates14)) == len(candidates14)
    assert len(set(candidates30)) == len(candidates30)
    for gamma14, gamma30 in zip(candidates14, candidates30):
        assert reflected_gamma(gamma14) == gamma30
        assert reflected_gamma(gamma30) == gamma14
        a, b, c, d = map(orbit.vector_element, gamma14)
        ar, br, cr, dr = map(orbit.vector_element, gamma30)
        assert a * d - b * c == orbit.ONE
        assert ar * dr - br * cr == orbit.ONE
        assert (c * x14 + d) ** 2 + c**2 * height_squared == (
            cr * x30 + dr
        ) ** 2 + cr**2 * height_squared
        assert (a * x14 + b) ** 2 + a**2 * height_squared == (
            ar * x30 + br
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
    assert sum(count for _, count in result.histogram) == 9036
    assert result.digest == EXPECTED_DIGESTS[mask]
    return result


def main() -> None:
    assert orbit.field_discriminant() == 453789
    assert 14 ^ 48 == 62 and 30 ^ 48 == 46
    assert {14, 62}.isdisjoint({30, 46})
    certify_root_pair_labels()
    print("INDEPENDENT d=7 NONTRIVIAL-R83 MASKS 14/30 STURM CHECK v47")
    print("stdlib field/Sturm/LDL core; no primary certificate or experiment imports")
    print("canonical root pairs: mask 14 xor 48 = 62; mask 30 xor 48 = 46")

    points = {mask: selected_point(mask) for mask in MASKS}
    matrix14, center14, height14, modulus14 = points[14]
    matrix30, center30, height30, modulus30 = points[30]
    reflection = ambient.matrix(-1, 0, 0, 1)
    assert ambient.mdet(reflection) == -ambient.ONE
    assert ambient.sign_vector(
        ambient.mdet(reflection), ambient.embedded_roots()
    ) == (-1,) * orbit.DEGREE
    assert ambient.all_entries_integral(reflection)
    assert ambient.minverse(reflection) == reflection
    assert matrix30 == reflect_matrix(matrix14)
    assert matrix30 == ambient.conjugate_by(reflection, matrix14)
    assert center30 == -center14
    assert height30 == height14
    assert modulus30 == modulus14
    print(
        "matrix reflection: diag(-1,1) lies in GL2(O_L) but not "
        "GL2+(O_L) or SL2(O_L)"
    )

    bounds14 = derive_bounds(height14, modulus14)
    bounds30 = derive_bounds(height30, modulus30)
    assert bounds30 == bounds14 == EXPECTED_BOUNDS
    install_point_data(height14, bounds14)
    print("shared exact outward bounds:", bounds14)

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
    assert (len(a_values), len(c_values)) == (23923, 6605)
    assert {negate_vector(value) for value in a_values} == set(a_values)
    assert {negate_vector(value) for value in c_values} == set(c_values)
    print("single exact a/c enumeration:", len(a_values), len(c_values), flush=True)

    bounded14 = enumerate_bounded_candidates(center14, a_values, c_values)
    candidates14 = bounded14.nonzero + bounded14.zero
    candidates30 = tuple(reflected_gamma(gamma) for gamma in candidates14)
    validate_reflection_bijection(
        center14, center30, height14, candidates14, candidates30
    )
    print("exact reflected candidate bijection: 9036/9036", flush=True)

    result14 = polygon_result(14, center14, candidates14)
    result30 = polygon_result(30, center30, candidates30)
    for mask, result in ((14, result14), (30, result30)):
        print(
            f"mask{mask}: bottom rows={bounded14.bottom_rows}; "
            f"c!=0 candidates={len(bounded14.nonzero)}; "
            f"c=0 candidates={len(bounded14.zero)}; polygon survivors=0"
        )
        print(f"mask{mask}: first-exclusion histogram={dict(result.histogram)}")
        print(f"mask{mask}: canonical exclusion SHA-256={result.digest}")
    print("PASS: independent exact exhaustive exclusions for distinct masks 14 and 30")


if __name__ == "__main__":
    main()
