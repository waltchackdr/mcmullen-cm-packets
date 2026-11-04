#!/usr/bin/env python3
"""Independent exact ten-pair nontrivial-R83 exclusion batch (v51).

This standard-library/Sturm verifier imports only the independent v45 ambient
packet reconstruction and the independent v42 orbit/polygon foundation.  It
reconstructs twenty canonical conductor-83 points, proves the exact reflected
candidate bijection inside each pair, and applies the polygon predicates at
the twenty distinct centers.  No experiment or primary certificate is
imported, and no floating-point arithmetic enters a decisive comparison.
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
BASE_MASKS = (0, 2, 3, 4, 5, 8, 9, 10, 11, 12)
MASK_PAIRS = tuple((mask, mask + 16) for mask in BASE_MASKS)

# Frozen saturated-basis coordinates in (xi5,xi11,xi13,xi17,4-t^2).
EXPECTED_EXPONENTS = {
    0: (0, 0, -1, -1, -2),
    2: (1, 1, 1, 0, 1),
    3: (1, 1, 0, 0, 0),
    4: (0, 0, 0, 0, -1),
    5: (1, -1, 0, -1, -2),
    8: (0, 0, 0, 0, -1),
    9: (1, 0, 0, 0, -1),
    10: (0, 1, 0, 0, -1),
    11: (1, 0, 0, 0, -1),
    12: (0, 0, 0, -1, -2),
}

# a entries, c entries, coupled bottom rows, c!=0 candidates, c=0 candidates.
EXPECTED_COUNTS = {
    0: (24035, 6633, 2023, 8136, 0),
    2: (24137, 6607, 1967, 8334, 0),
    3: (23993, 6597, 1937, 7794, 0),
    4: (24059, 6599, 2015, 8702, 10),
    5: (23959, 6619, 2083, 8878, 0),
    8: (23967, 6613, 2021, 8750, 6),
    9: (24001, 6607, 2085, 9146, 0),
    10: (24039, 6611, 2017, 9436, 6),
    11: (23973, 6617, 2013, 8360, 0),
    12: (24013, 6635, 2015, 8508, 4),
}

EXPECTED_HISTOGRAMS = {
    0: (((0, 0), 1810), ((0, 1), 4964), ((0, 2), 876), ((0, 3), 346), ((0, 4), 20), ((0, 5), 14), ((1, 0), 78), ((1, 1), 12), ((1, 2), 12), ((1, 4), 4)),
    16: (((0, 0), 1736), ((0, 1), 5054), ((0, 2), 886), ((0, 3), 346), ((0, 4), 24), ((0, 5), 10), ((1, 0), 50), ((1, 1), 18), ((1, 2), 12)),
    2: (((0, 0), 1658), ((0, 1), 5194), ((0, 2), 996), ((0, 3), 360), ((0, 4), 18), ((0, 5), 12), ((1, 0), 66), ((1, 1), 14), ((1, 2), 6), ((1, 3), 8), ((1, 4), 2)),
    18: (((0, 0), 1694), ((0, 1), 5220), ((0, 2), 1006), ((0, 3), 296), ((0, 4), 28), ((0, 5), 10), ((1, 0), 50), ((1, 1), 22), ((1, 2), 4), ((1, 3), 2), ((1, 4), 2)),
    3: (((0, 0), 1644), ((0, 1), 4848), ((0, 2), 900), ((0, 3), 282), ((0, 4), 18), ((0, 5), 14), ((1, 0), 62), ((1, 1), 18), ((1, 2), 6), ((1, 4), 2)),
    19: (((0, 0), 1630), ((0, 1), 4874), ((0, 2), 896), ((0, 3), 288), ((0, 4), 16), ((0, 5), 8), ((1, 0), 60), ((1, 1), 14), ((1, 2), 2), ((1, 3), 4), ((1, 5), 2)),
    4: (((0, 0), 1806), ((0, 1), 5434), ((0, 2), 928), ((0, 3), 394), ((0, 4), 38), ((0, 5), 10), ((1, 0), 68), ((1, 1), 26), ((1, 2), 2), ((1, 3), 2), ((1, 4), 4)),
    20: (((0, 0), 1808), ((0, 1), 5352), ((0, 2), 980), ((0, 3), 420), ((0, 4), 26), ((0, 5), 12), ((1, 0), 72), ((1, 1), 24), ((1, 2), 8), ((1, 3), 6), ((1, 4), 4)),
    5: (((0, 0), 1774), ((0, 1), 5618), ((0, 2), 972), ((0, 3), 392), ((0, 4), 26), ((0, 5), 20), ((1, 0), 52), ((1, 1), 16), ((1, 2), 6), ((1, 4), 2)),
    21: (((0, 0), 1870), ((0, 1), 5484), ((0, 2), 1022), ((0, 3), 392), ((0, 4), 24), ((0, 5), 16), ((1, 0), 46), ((1, 1), 6), ((1, 2), 8), ((1, 3), 4), ((1, 4), 6)),
    8: (((0, 0), 1838), ((0, 1), 5528), ((0, 2), 886), ((0, 3), 386), ((0, 4), 28), ((0, 5), 8), ((1, 0), 54), ((1, 1), 18), ((1, 2), 2), ((1, 4), 6), ((1, 5), 2)),
    24: (((0, 0), 1798), ((0, 1), 5550), ((0, 2), 900), ((0, 3), 402), ((0, 4), 24), ((0, 5), 10), ((1, 0), 50), ((1, 1), 20), ((1, 4), 2)),
    9: (((0, 0), 1752), ((0, 1), 5858), ((0, 2), 1042), ((0, 3), 348), ((0, 4), 26), ((0, 5), 14), ((1, 0), 82), ((1, 1), 14), ((1, 2), 2), ((1, 3), 8)),
    25: (((0, 0), 1806), ((0, 1), 5822), ((0, 2), 1052), ((0, 3), 328), ((0, 4), 14), ((0, 5), 16), ((1, 0), 76), ((1, 1), 22), ((1, 2), 4), ((1, 3), 4), ((1, 4), 2)),
    10: (((0, 0), 1964), ((0, 1), 5880), ((0, 2), 1066), ((0, 3), 404), ((0, 4), 36), ((0, 5), 4), ((1, 0), 66), ((1, 1), 16), ((1, 2), 4), ((1, 3), 2)),
    26: (((0, 0), 1970), ((0, 1), 5884), ((0, 2), 1046), ((0, 3), 392), ((0, 4), 36), ((0, 5), 16), ((1, 0), 76), ((1, 1), 12), ((1, 2), 2), ((1, 3), 4), ((1, 4), 4)),
    11: (((0, 0), 1658), ((0, 1), 5182), ((0, 2), 986), ((0, 3), 408), ((0, 4), 24), ((0, 5), 12), ((1, 0), 66), ((1, 1), 18), ((1, 2), 2), ((1, 3), 2), ((1, 4), 2)),
    27: (((0, 0), 1792), ((0, 1), 5110), ((0, 2), 944), ((0, 3), 402), ((0, 4), 26), ((0, 5), 8), ((1, 0), 62), ((1, 1), 10), ((1, 2), 2), ((1, 3), 4)),
    12: (((0, 0), 1760), ((0, 1), 5214), ((0, 2), 1012), ((0, 3), 364), ((0, 4), 34), ((0, 5), 20), ((1, 0), 82), ((1, 1), 18), ((1, 2), 2), ((1, 4), 6)),
    28: (((0, 0), 1852), ((0, 1), 5196), ((0, 2), 982), ((0, 3), 364), ((0, 4), 34), ((0, 5), 10), ((1, 0), 56), ((1, 1), 12), ((1, 2), 2), ((1, 4), 4)),
}

EXPECTED_DIGESTS = {
    0: "45804d2344a890e78e1346763d937d6ad1432cc50ca4490b1e502406685f065d",
    16: "bf73238c70cdf49152d44e9fee300041080e96be78c598d45cdf9b38c55cfe92",
    2: "b24eee6e97452507710752855ba850802909c6bb07fce8590ef91e39f354fe32",
    18: "af37bb608248a19682e7a348b2c8cc20dd861f014138877f21daff0b5d814628",
    3: "de62d03cd5d4c5f403d2495661f8229545070a104bb8ceeb46acd7f77f1b8fc2",
    19: "b80ae682698dd57828230fee4c6b82cc75615fc5e8eda375864ef2840b7b519a",
    4: "46710d59df6c922304a528bba2f38e9affc857e8ad436294765defff0635a4e0",
    20: "2101d71d1d4f392c6bf360651f092e2e90b2a201633506b22816b05cf049c822",
    5: "491b99f68f4b8cb95ea9650a355e737dfb67934f505388c62a7b36e3cf44938e",
    21: "484642f37fbccdbd1ca08e9f090d455e0a4c7417e1fdd8b0d4c7ca57bf0a966e",
    8: "dcc56539ddcc7ca63d1b52401bac336528397525b66981519a887cab063c5d22",
    24: "8ce35e562295148abab3aae0025ff13f118608542097f877adcf7961c5bcc042",
    9: "fddabcc85df1115417de6efe4908b691c695eb06809a1af4d6526d24b7bc05f5",
    25: "99fcae815bea2acc753ed956237497a9e8f3a83e308bfacd732b83ef78a7146d",
    10: "5465692d11f69ff0bf3848eba868528787ef10a74c48aac38fd0b616bd0503d6",
    26: "a00e7cc9b32857d6731c1cd177926ed6f46ad109212a67d7c7cdf99642d87b04",
    11: "81454a94c45e6a015f79431b32d2947cc08f2f1d68107a3af47831f2f131c55d",
    27: "06feaeb73141073f4483109e0f3e385fb96853cc27fce9cd650a2c3158d6a0c1",
    12: "e87d1de9b7f8a963be479adedb8962472fc969b68bc4c051c4f36457f9d57a0b",
    28: "70c63bd13e29a08d4b5d56ee159955cfc8fe118ea64d76d5ccc79cafcf9f8ea0",
}

EXPECTED_AGGREGATES = {
    0: "4ccbb78db598145e486f6ce920ce2b242414c20fd1db41d0a48b24a46756f7e6",
    2: "b1260e07bd41712965921d1e2270784e0e68e3a032ca2eca6ddb475911d393e1",
    3: "3205f62662e42a3c4b267f242533242535e463789ce27f27ae11ee12ac868cb4",
    4: "eb8317988da0d674f56fa88997d3adde2eeaa3749ea2839fa94e6b2cb187d656",
    5: "393f2e41d0fdc003649301b07ad0673c86e7dcf31223facf82d7ad1f93a5781b",
    8: "850196c198f67497ebec722b6296026ae43f926fd793d0a9694ad1c693ba137d",
    9: "cdd94a0e3e1871145f93793cbc0bbf1f6facf3aa96bcd3c02b0fe3df1279a7b3",
    10: "d45e2712f974189abdaf4ca4179c969072e9c16a52fe27ccbea0a02a15ebac34",
    11: "b4af0ce6777f6fe9a3cc2ebba3226fa1a8ec8f8f2ba23322ccf3014ccf8e451d",
    12: "3d655c9b51c04c70655bf9d52b1b23be5780bebe5607e2e5250094c7ec25dd47",
}


def product(values: tuple[ambient.LElement, ...]) -> ambient.LElement:
    answer = ambient.ONE
    for value in values:
        answer *= value
    return answer


def orientation_units() -> tuple[ambient.LElement, ...]:
    return (
        ambient.T**4 - 3 * ambient.T**2 + 1,
        ambient.T**3 + ambient.T**2 - 2 * ambient.T - 1,
        ambient.T**5 - 4 * ambient.T**3 + 3 * ambient.T + 1,
        ambient.T**3 - 2 * ambient.T,
        -ambient.ONE,
        4 - ambient.T**2,
    )


def orientation_unit(mask: int) -> ambient.LElement:
    return product(tuple(unit for bit, unit in enumerate(orientation_units()) if (mask >> bit) & 1))


def balancing_unit(base_mask: int) -> ambient.LElement:
    xi5, xi11, xi13, xi17, _, u = orientation_units()
    return product(tuple(generator**exponent for generator, exponent in zip((xi5, xi11, xi13, xi17, u), EXPECTED_EXPONENTS[base_mask])))


def conductor83_ideal() -> ambient.Matrix:
    maximal = ambient.multiplication_matrix(ambient.alpha_pair(INDEX))
    gamma = ambient.elt(-2, -2, 2, 0, -1)
    basis = ambient.matrix(16, gamma, ambient.elt(58, 3), 3 * gamma)
    ideal = ambient.conjugate_by(ambient.minverse(basis), maximal)
    assert ambient.all_entries_integral(ideal)
    return ideal


def certify_root_labels(ideal: ambient.Matrix) -> None:
    u = 4 - ambient.T**2
    p = ambient.elt(55, 166, -18, -109, 0, 18)
    q = ambient.elt(-4, -9, -3, 5, 1, -1)
    r = ambient.elt(-603, -6206, 905, 4245, -184, -683)
    anti = ambient.matrix(p, q, r, -p)
    assert ambient.all_entries_integral(anti)
    assert ambient.mdet(anti) == -u
    assert ambient.mmul(anti, ideal) == ambient.mmul(ambient.root_conjugate(ideal), anti)
    for mask in tuple(value for pair in MASK_PAIRS for value in pair):
        partner = mask ^ 48
        source = ambient.twist(ideal, orientation_unit(mask))
        target = ambient.twist(ideal, orientation_unit(partner))
        bridge = ambient.mmul(
            ambient.mmul(ambient.diagonal(orientation_unit(partner)), ambient.minverse(anti)),
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
    base_mask = mask % 16
    assert base_mask in BASE_MASKS and mask in (base_mask, base_mask + 16)
    alpha = ambient.alpha_pair(INDEX)
    assert ambient.pair_norm(alpha) == ambient.PI
    raw = ambient.twist(ideal, orientation_unit(mask))
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

    unit = balancing_unit(base_mask)
    assert unit.integral_power_basis() and unit.inverse().integral_power_basis()
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
    assert ambient.sign_vector(height_squared, ambient.embedded_roots()) == (1,) * orbit.DEGREE
    assert ambient.sign_vector(modulus_squared, ambient.embedded_roots()) == (1,) * orbit.DEGREE
    assert height_squared.norm() > 0 and modulus_squared.norm() > 0
    return Point(balanced, center, height_squared, modulus_squared)


def reflect_matrix(value: ambient.Matrix) -> ambient.Matrix:
    return ambient.matrix(value[0][0], -value[0][1], -value[1][0], value[1][1])


def outward_hundredth(value: Q) -> Q:
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(point: Point) -> dict[str, tuple[Q, ...]]:
    columns: list[list[Q]] = [[] for _ in range(6)]
    for embedding, root in enumerate(orbit.ROOTS):
        height = orbit.sqrt_interval(orbit.interval_eval(point.height_squared, root), 160)
        modulus = orbit.sqrt_interval(orbit.interval_eval(point.modulus_squared, root), 160)
        assert height.lo > 0 and modulus.lo > 0
        den = orbit.sqrt_interval(orbit.Interval.point(height.hi / orbit.M_LOWER[embedding]), 160).hi
        c_entry = orbit.sqrt_interval(orbit.Interval.point(1 / (orbit.M_LOWER[embedding] * height.lo)), 160).hi
        num = orbit.U_UPPER[embedding] * den
        a_entry = num / height.lo
        d_entry = den + c_entry * modulus.hi
        b_entry = num + a_entry * modulus.hi
        for target, value in zip(columns, (den, num, a_entry, b_entry, c_entry, d_entry)):
            target.append(outward_hundredth(value))
    return dict(zip(("DEN", "NUM", "A", "B", "C", "D"), map(tuple, columns)))


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
    orbit.Y2_SCALED = tuple(orbit.scaled_interval(orbit.interval_eval(point.height_squared, root)) for root in orbit.ROOTS)


@dataclass(frozen=True)
class CandidateSet:
    bottom_rows: int
    nonzero: tuple[Gamma, ...]
    zero: tuple[Gamma, ...]


def enumerate_candidates(base_mask: int, point: Point, a_values: list[Vector], c_values: list[Vector]) -> CandidateSet:
    x = point.center
    x_boxes = tuple(orbit.scaled_interval(orbit.interval_eval(x, root)) for root in orbit.ROOTS)
    a_box_cache = {a: tuple(orbit.scaled_vector_box(a, i) for i in range(orbit.DEGREE)) for a in a_values}
    bottom_rows: list[tuple[Vector, Vector]] = []
    for c_index, c in enumerate(c_values):
        center = tuple(-(orbit.vector_element(c) * x).coefficients[i] for i in range(orbit.DEGREE))
        c_boxes = tuple(orbit.scaled_vector_box(c, i) for i in range(orbit.DEGREE))
        bottom_rows.extend(
            (c, d)
            for d in orbit.enumerate_trace_ball(
                center,
                orbit.DEN_RADIUS,
                lambda d, c=c, c_boxes=c_boxes: orbit.bottom_bounded(c, d, x, c_boxes, x_boxes),
            )
        )
        if c_index and c_index % 1000 == 0:
            print(f"  pair {base_mask}/{base_mask + 16}: c {c_index}/{len(c_values)}, rows {len(bottom_rows)}", flush=True)

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
        if row_index and row_index % 500 == 0:
            print(f"  pair {base_mask}/{base_mask + 16}: determinant rows {row_index}/{len(nonzero_rows)}, candidates {len(nonzero)}", flush=True)

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
    expected = EXPECTED_COUNTS[base_mask]
    assert (result.bottom_rows, len(result.nonzero), len(result.zero)) == expected[2:]
    return result


def negate_vector(value: Vector) -> Vector:
    return tuple(-coefficient for coefficient in value)


def reflected_gamma(gamma: Gamma) -> Gamma:
    a, b, c, d = gamma
    return a, negate_vector(b), negate_vector(c), d


def validate_reflection(base: Point, reflected: Point, candidates: tuple[Gamma, ...]) -> tuple[Gamma, ...]:
    assert reflected.center == -base.center
    assert reflected.height_squared == base.height_squared
    assert reflected.modulus_squared == base.modulus_squared
    reflected_candidates = tuple(reflected_gamma(gamma) for gamma in candidates)
    assert len(set(candidates)) == len(candidates)
    assert len(set(reflected_candidates)) == len(reflected_candidates)
    for source, target in zip(candidates, reflected_candidates):
        assert reflected_gamma(target) == source
        a, b, c, d = map(orbit.vector_element, source)
        ar, br, cr, dr = map(orbit.vector_element, target)
        assert a * d - b * c == orbit.ONE
        assert ar * dr - br * cr == orbit.ONE
        assert (c * base.center + d) ** 2 + c**2 * base.height_squared == (cr * reflected.center + dr) ** 2 + cr**2 * reflected.height_squared
        assert (a * base.center + b) ** 2 + a**2 * base.height_squared == (ar * reflected.center + br) ** 2 + ar**2 * reflected.height_squared
    return reflected_candidates


def polygon_result(mask: int, point: Point, candidates: tuple[Gamma, ...]) -> str:
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
    frozen = tuple(sorted(histogram.items()))
    digest = orbit.canonical_digest(records)
    assert survivors == []
    assert frozen == EXPECTED_HISTOGRAMS[mask]
    assert sum(histogram.values()) == len(candidates)
    assert digest == EXPECTED_DIGESTS[mask]
    return digest


def certify_disjoint_labels() -> None:
    new_fibers = tuple(frozenset((mask, mask ^ 48)) for pair in MASK_PAIRS for mask in pair)
    assert len(new_fibers) == 20 and len(set(new_fibers)) == 20
    assert all(len(fiber) == 2 for fiber in new_fibers)
    prior = tuple(map(frozenset, ({14, 62}, {30, 46}, {15, 63}, {31, 47}, {6, 54}, {22, 38}, {1, 49}, {17, 33})))
    assert all(left.isdisjoint(right) for left in new_fibers for right in prior)
    assert all(new_fibers[i].isdisjoint(new_fibers[j]) for i in range(len(new_fibers)) for j in range(i))


def main() -> None:
    assert orbit.field_discriminant() == 453789
    certify_disjoint_labels()
    ideal = conductor83_ideal()
    certify_root_labels(ideal)
    reflection = ambient.matrix(-1, 0, 0, 1)
    assert ambient.all_entries_integral(reflection)
    assert ambient.minverse(reflection) == reflection
    assert ambient.mdet(reflection) == -ambient.ONE
    assert ambient.sign_vector(ambient.mdet(reflection), ambient.embedded_roots()) == (-1,) * orbit.DEGREE
    print("INDEPENDENT d=7 TEN-PAIR NONTRIVIAL-R83 STURM CHECK v51")
    print("stdlib field/Sturm/trace-LDL core; no primary certificate or experiment imports")
    print("diag(-1,1) is in GL2(O_L), not GL2+(O_L) or SL2(O_L)")

    for pair_index, (base_mask, reflected_mask) in enumerate(MASK_PAIRS, 1):
        print(f"PAIR {pair_index}/10: masks {base_mask}/{reflected_mask}; root fibers {{{base_mask},{base_mask ^ 48}}}, {{{reflected_mask},{reflected_mask ^ 48}}}", flush=True)
        print("  balancing exponents:", EXPECTED_EXPONENTS[base_mask], flush=True)
        base = selected_point(base_mask, ideal)
        reflected = selected_point(reflected_mask, ideal)
        assert reflected.matrix == reflect_matrix(base.matrix)
        assert reflected.matrix == ambient.conjugate_by(reflection, base.matrix)
        bounds = derive_bounds(base)
        assert derive_bounds(reflected) == bounds
        install_point_data(base, bounds)
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
        expected = EXPECTED_COUNTS[base_mask]
        assert (len(a_values), len(c_values)) == expected[:2]
        assert {negate_vector(value) for value in a_values} == set(a_values)
        assert {negate_vector(value) for value in c_values} == set(c_values)
        if base_mask == 3:
            endpoint = (-7, 0, 0, 0, 0, 0)
            assert bounds["C"][1] == Q(7)
            endpoint_value = orbit.vector_element(endpoint)
            assert orbit.sign_at(ambient.LElement(7) + endpoint_value, 1) == 0
            assert orbit.sign_at(ambient.LElement(7) - endpoint_value, 1) == 1
            # The equality at embedding 1 is accepted exactly; a later
            # embedding rejects this vector from the full six-place box.
            assert not orbit.field_box_membership(endpoint, bounds["C"])
            print("  exact endpoint equality resolved: (-7,0,0,0,0,0) at C[1]=7", flush=True)
        print(f"  exact a/c lists: {len(a_values)}/{len(c_values)}", flush=True)

        bounded = enumerate_candidates(base_mask, base, a_values, c_values)
        candidates = bounded.nonzero + bounded.zero
        reflected_candidates = validate_reflection(base, reflected, candidates)
        digest_base = polygon_result(base_mask, base, candidates)
        digest_reflected = polygon_result(reflected_mask, reflected, reflected_candidates)
        aggregate = sha256(
            f"{base_mask}:{digest_base}\n{reflected_mask}:{digest_reflected}\n".encode("ascii")
        ).hexdigest()
        assert aggregate == EXPECTED_AGGREGATES[base_mask]
        print(
            f"  rows={bounded.bottom_rows}; c!=0={len(bounded.nonzero)}; c=0={len(bounded.zero)}; survivors=0/0",
            flush=True,
        )
        print(f"  digests={digest_base}, {digest_reflected}; aggregate={aggregate}", flush=True)

    assert {base: EXPECTED_COUNTS[base][4] for base in (4, 8, 10, 12)} == {4: 10, 8: 6, 10: 6, 12: 4}
    assert all(EXPECTED_COUNTS[base][4] == 0 for base in (0, 2, 3, 5, 9, 11))
    print("PASS: independent exact exhaustive exclusions for all twenty distinct points")


if __name__ == "__main__":
    main()
