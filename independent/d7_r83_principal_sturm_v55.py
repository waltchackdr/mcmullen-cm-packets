#!/usr/bin/env python3
"""Independent exact principal-Picard proper-R83 orbit search (v55).

Exact independent verifier; theorem status is governed by the claim ledger
and its clean-replay gate. This verifier imports only the retained
standard-library / Sturm independent arithmetic foundations. It does not
import SymPy, a primary certificate, or an experiment.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from math import floor
from pathlib import Path
import sys

import d7_ambient_points_sturm_v45 as ambient
import residual_orbit_k17_k38_sturm_v42 as orbit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "certificates"))
import d7_r83_principal_expected_v55 as frozen  # noqa: E402


Vector = tuple[int, ...]
Gamma = tuple[Vector, Vector, Vector, Vector]
ZERO_VECTOR: Vector = (0,) * orbit.DEGREE
INDEX = 2
CANONICAL_MASKS = tuple(range(16)) + tuple(range(32, 48))
assert CANONICAL_MASKS == frozen.CANONICAL_MASKS

# Exact determinant-one balancing witnesses in the independent cyclotomic-unit
# basis (xi_5, xi_11, xi_13, xi_17, xi_19). They are witnesses, not numerical
# approximations; every unit, inverse, transformation law, and bound is checked.
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


def product(values: tuple[ambient.LElement, ...]) -> ambient.LElement:
    answer = ambient.ONE
    for value in values:
        answer *= value
    return answer


def square_class_basis() -> tuple[ambient.LElement, ...]:
    return (
        ambient.T**4 - 3 * ambient.T**2 + 1,
        ambient.T**3 + ambient.T**2 - 2 * ambient.T - 1,
        ambient.T**5 - 4 * ambient.T**3 + 3 * ambient.T + 1,
        ambient.T**3 - 2 * ambient.T,
        -ambient.ONE,
        4 - ambient.T**2,
    )


def cyclotomic_units() -> tuple[ambient.LElement, ...]:
    xi5, xi11, xi13, xi17, _, positive_nonsquare = square_class_basis()
    # The exact cyclotomic relation used in the primary power basis gives
    # xi_19=t. Recheck it independently rather than importing that assertion.
    xi19 = positive_nonsquare**-3 / (xi5 * xi11 * xi13 * xi17)
    assert xi19 == ambient.T
    return xi5, xi11, xi13, xi17, xi19


def square_class_unit(mask: int) -> ambient.LElement:
    return product(tuple(
        unit for bit, unit in enumerate(square_class_basis())
        if (mask >> bit) & 1
    ))


def balancing_unit(mask: int) -> ambient.LElement:
    return product(tuple(
        generator**exponent
        for generator, exponent in zip(
            cyclotomic_units(), BALANCING_EXPONENTS[mask]
        )
    ))


def regular_action() -> ambient.Matrix:
    alpha = ambient.alpha_pair(INDEX)
    trace = ambient.pair_trace(alpha)
    assert ambient.pair_norm(alpha) == ambient.PI
    value = ambient.companion(trace, ambient.PI)
    assert ambient.all_entries_integral(value)
    assert (ambient.mtrace(value), ambient.mdet(value)) == (trace, ambient.PI)

    # Certify the proper regular-order label independently: the maximal
    # generator z is not integral on this lattice, and the discriminant has
    # conductor norm 83^2.
    maximal = ambient.multiplication_matrix(alpha)
    z_matrix = ambient.matrix(0, -1, 1, ambient.T)
    bridge = ambient.column_intertwiner(maximal)
    assert ambient.mmul(bridge, value) == ambient.mmul(maximal, bridge)
    assert not ambient.all_entries_integral(
        ambient.conjugate_by(ambient.minverse(bridge), z_matrix)
    )
    discriminant = trace**2 - 4 * ambient.PI
    assert abs(discriminant.norm()) == 83**2
    return value


def certify_root_labels(base: ambient.Matrix) -> None:
    trace = ambient.mtrace(base)
    anti = ambient.matrix(1, trace, 0, -1)
    assert ambient.all_entries_integral(anti)
    assert ambient.mdet(anti) == -ambient.ONE
    assert ambient.mmul(anti, base) == ambient.mmul(
        ambient.root_conjugate(base), anti
    )
    raw_masks = set()
    for mask in CANONICAL_MASKS:
        partner = mask ^ 16
        assert partner not in CANONICAL_MASKS
        source = ambient.twist(base, square_class_unit(mask))
        target = ambient.twist(base, square_class_unit(partner))
        bridge = ambient.mmul(
            ambient.mmul(
                ambient.diagonal(square_class_unit(partner)),
                ambient.minverse(anti),
            ),
            ambient.minverse(ambient.diagonal(square_class_unit(mask))),
        )
        assert ambient.all_entries_integral(bridge)
        assert ambient.mdet(bridge) == ambient.ONE
        assert ambient.conjugate_by(
            bridge, ambient.root_conjugate(source)
        ) == target
        raw_masks.update((mask, partner))
    assert raw_masks == set(range(64))


@dataclass(frozen=True)
class Point:
    matrix: ambient.Matrix
    center: ambient.LElement
    height_squared: ambient.LElement
    modulus_squared: ambient.LElement


def selected_point(mask: int, base: ambient.Matrix) -> Point:
    assert mask in CANONICAL_MASKS
    raw = ambient.twist(base, square_class_unit(mask))
    assert ambient.all_entries_integral(raw)
    raw_x, raw_y2 = ambient.exact_fixed_invariants(raw)
    raw_modulus2 = raw_x**2 + raw_y2

    unit = balancing_unit(mask)
    assert unit.integral_power_basis()
    assert unit.inverse().integral_power_basis()
    conjugator = ambient.matrix(unit, 0, 0, unit.inverse())
    assert ambient.all_entries_integral(conjugator)
    assert ambient.all_entries_integral(ambient.minverse(conjugator))
    assert ambient.mdet(conjugator) == ambient.ONE
    balanced = ambient.conjugate_by(conjugator, raw)
    assert ambient.all_entries_integral(balanced)
    assert (ambient.mtrace(balanced), ambient.mdet(balanced)) == (
        ambient.mtrace(base), ambient.PI
    )

    center, height_squared = ambient.exact_fixed_invariants(balanced)
    modulus_squared = center**2 + height_squared
    assert center == unit**2 * raw_x
    assert height_squared == unit**4 * raw_y2
    assert modulus_squared == unit**4 * raw_modulus2
    roots = ambient.embedded_roots()
    assert ambient.sign_vector(height_squared, roots) == (1,) * orbit.DEGREE
    assert ambient.sign_vector(modulus_squared, roots) == (1,) * orbit.DEGREE
    return Point(balanced, center, height_squared, modulus_squared)


def outward_hundredth(value: Q) -> Q:
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(point: Point) -> dict[str, tuple[Q, ...]]:
    columns: list[list[Q]] = [[] for _ in range(6)]
    for embedding, root in enumerate(orbit.ROOTS):
        height = orbit.sqrt_interval(
            orbit.interval_eval(point.height_squared, root), 160
        )
        modulus = orbit.sqrt_interval(
            orbit.interval_eval(point.modulus_squared, root), 160
        )
        assert height.lo > 0 and modulus.lo > 0
        # Each stored entry rounds outward from its own pre-rounding hatted
        # rational endpoint. D and B do not reuse separately rounded entries.
        den_hat = orbit.sqrt_interval(
            orbit.Interval.point(height.hi / orbit.M_LOWER[embedding]), 160
        ).hi
        c_hat = orbit.sqrt_interval(
            orbit.Interval.point(
                1 / (orbit.M_LOWER[embedding] * height.lo)
            ),
            160,
        ).hi
        num_hat = orbit.U_UPPER[embedding] * den_hat
        a_hat = num_hat / height.lo
        d_hat = den_hat + c_hat * modulus.hi
        b_hat = num_hat + a_hat * modulus.hi
        for target, value in zip(
            columns, (den_hat, num_hat, a_hat, b_hat, c_hat, d_hat)
        ):
            target.append(outward_hundredth(value))
    return dict(
        zip(("DEN", "NUM", "A", "B", "C", "D"), map(tuple, columns))
    )


def install_point_data(
    point: Point, bounds: dict[str, tuple[Q, ...]]
) -> None:
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


def negate_vector(value: Vector) -> Vector:
    return tuple(-coefficient for coefficient in value)


def entry_lists(point: Point, bounds: dict[str, tuple[Q, ...]]):
    install_point_data(point, bounds)
    a_values = orbit.enumerate_trace_ball(
        (Q(0),) * orbit.DEGREE,
        orbit.sum_of_squares(orbit.A_BOUND),
        lambda value: orbit.field_box_membership(value, orbit.A_BOUND),
    )
    c_values = orbit.enumerate_trace_ball(
        (Q(0),) * orbit.DEGREE,
        orbit.sum_of_squares(orbit.C_BOUND),
        lambda value: orbit.field_box_membership(value, orbit.C_BOUND),
    )
    a_values.sort()
    c_values.sort()
    assert {negate_vector(value) for value in a_values} == set(a_values)
    assert {negate_vector(value) for value in c_values} == set(c_values)
    return a_values, c_values


@dataclass(frozen=True)
class CandidateSet:
    bottom_rows: int
    nonzero: tuple[Gamma, ...]
    zero: tuple[Gamma, ...]
    zero_bottom_rows: int
    zero_d_rows: int
    integral_inverse_rows: int
    a_admissible_rows: int


def enumerate_candidates(
    mask: int, point: Point, a_values: list[Vector], c_values: list[Vector]
) -> CandidateSet:
    x = point.center
    x_boxes = tuple(
        orbit.scaled_interval(orbit.interval_eval(x, root))
        for root in orbit.ROOTS
    )
    a_boxes = {
        a: tuple(orbit.scaled_vector_box(a, i) for i in range(orbit.DEGREE))
        for a in a_values
    }
    bottom_rows: list[tuple[Vector, Vector]] = []
    for c_index, c in enumerate(c_values):
        center = tuple(
            -(orbit.vector_element(c) * x).coefficients[i]
            for i in range(orbit.DEGREE)
        )
        c_boxes = tuple(
            orbit.scaled_vector_box(c, i) for i in range(orbit.DEGREE)
        )
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
        if c_index and c_index % 25 == 0:
            print(
                f"  mask={mask:02d}: c rows {c_index}/{len(c_values)}; "
                f"bottom={len(bottom_rows)}",
                flush=True,
            )

    nonzero: list[Gamma] = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != ZERO_VECTOR]
    for row_index, (c, d) in enumerate(nonzero_rows):
        affine = orbit.determinant_one_affine_map(c, d)
        for a in a_values:
            b = affine.evaluate(a)
            if b is not None and orbit.top_bounded(
                a, b, x, a_boxes[a], x_boxes
            ):
                gamma = (a, b, c, d)
                ae, be, ce, de = map(orbit.vector_element, gamma)
                assert ae * de - be * ce == orbit.ONE
                nonzero.append(gamma)
        if row_index and row_index % 500 == 0:
            print(
                f"  mask={mask:02d}: determinant rows "
                f"{row_index}/{len(nonzero_rows)}; candidates={len(nonzero)}",
                flush=True,
            )

    zero: list[Gamma] = []
    zero_bottom_rows = 0
    zero_d_rows = 0
    integral_inverse_rows = 0
    a_admissible_rows = 0
    a_value_set = set(a_values)
    for c, d in bottom_rows:
        if c != ZERO_VECTOR:
            continue
        zero_bottom_rows += 1
        if d == ZERO_VECTOR:
            # For c=d=0 the determinant is identically zero, so this row
            # cannot contribute a determinant-one matrix.
            zero_d_rows += 1
            continue
        a_element = orbit.vector_element(d).inverse()
        if not a_element.integral_power_basis():
            continue
        integral_inverse_rows += 1
        a = orbit.integral_vector(a_element)
        if not orbit.field_box_membership(a, orbit.A_BOUND):
            continue
        assert a in a_value_set
        a_admissible_rows += 1
        center = tuple(
            -(a_element * x).coefficients[i] for i in range(orbit.DEGREE)
        )
        for b in orbit.enumerate_trace_ball(
            center,
            orbit.NUM_RADIUS,
            lambda b, a=a: orbit.top_bounded(
                a, b, x, a_boxes[a], x_boxes
            ),
        ):
            gamma = (a, b, c, d)
            ae, be, ce, de = map(orbit.vector_element, gamma)
            assert ce == orbit.ZERO and ae * de - be * ce == orbit.ONE
            zero.append(gamma)
    assert zero_d_rows == 1
    assert len(set(nonzero)) == len(nonzero)
    assert len(set(zero)) == len(zero)
    assert set(nonzero).isdisjoint(zero)
    return CandidateSet(
        len(bottom_rows),
        tuple(nonzero),
        tuple(zero),
        zero_bottom_rows,
        zero_d_rows,
        integral_inverse_rows,
        a_admissible_rows,
    )


def polygon_result(point: Point, candidates: tuple[Gamma, ...]):
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
    assert sum(histogram.values()) + len(survivors) == len(candidates)
    return tuple(sorted(histogram.items())), orbit.canonical_digest(records), tuple(survivors)


def parse_masks(value: str) -> tuple[int, ...]:
    if value == "all":
        return CANONICAL_MASKS
    answer = tuple(sorted({int(item) for item in value.split(",")}))
    if not answer or any(mask not in CANONICAL_MASKS for mask in answer):
        raise argparse.ArgumentTypeError(
            "use 'all' or comma-separated masks in 0,...,15,32,...,47"
        )
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masks", type=parse_masks, default=CANONICAL_MASKS)
    parser.add_argument(
        "--mode", choices=("bounds", "count", "enumerate"), default="bounds"
    )
    args = parser.parse_args()

    assert ambient.field_discriminant() == 453789
    base = regular_action()
    certify_root_labels(base)
    print("INDEPENDENT PRINCIPAL-PICARD proper-R83 STURM SEARCH v55")
    print("stdlib field/Sturm/trace-LDL core; no primary imports")
    print("root involution xor16 certified by integral determinant-one bridges")

    results = {}
    for mask in args.masks:
        point = selected_point(mask, base)
        bounds = derive_bounds(point)
        print(
            f"mask={mask:02d} exponents={BALANCING_EXPONENTS[mask]}",
            flush=True,
        )
        if args.mode == "bounds":
            continue
        a_values, c_values = entry_lists(point, bounds)
        assert (len(a_values), len(c_values)) == frozen.EXPECTED_ENTRY_COUNTS[mask]
        print(
            f"mask={mask:02d} exact a/c counts="
            f"{len(a_values)}/{len(c_values)}",
            flush=True,
        )
        if args.mode == "count":
            continue
        bounded = enumerate_candidates(mask, point, a_values, c_values)
        print(
            f"mask={mask:02d} c=0 audit="
            f"{bounded.zero_bottom_rows}/{bounded.zero_d_rows}/"
            f"{bounded.integral_inverse_rows}/{bounded.a_admissible_rows}/"
            f"{len(bounded.zero)}",
            flush=True,
        )
        candidates = bounded.nonzero + bounded.zero
        histogram, digest, survivors = polygon_result(point, candidates)
        result = (
            bounded.bottom_rows,
            len(bounded.nonzero),
            len(bounded.zero),
            histogram,
            digest,
            (
                bounded.zero_bottom_rows,
                bounded.zero_d_rows,
                bounded.integral_inverse_rows,
                bounded.a_admissible_rows,
                len(bounded.zero),
            ),
        )
        assert survivors == ()
        assert result == frozen.EXPECTED[mask]
        results[mask] = result
        print(
            f"mask={mask:02d} exact polygon result="
            f"{result[:-1] + (survivors,)}",
            flush=True,
        )

    if args.mode == "enumerate":
        assert tuple(sorted(results)) == args.masks
        if args.masks == CANONICAL_MASKS:
            assert frozen.batch_digest(results) == frozen.EXPECTED_BATCH_DIGEST
        print("PASS: independent exact principal-Picard proper-R83 exclusions")
    else:
        print(f"SCOPE: independent exact {args.mode} stage only")


if __name__ == "__main__":
    main()
