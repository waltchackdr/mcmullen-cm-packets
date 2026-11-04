#!/usr/bin/env python3
"""Independent exact exclusion of the residual imprimitive pair k={1,22}.

This is a standard-library-only verifier.  It imports only the quotient-field
and Sturm primitives from ``quadrilateral_k17_stdlib_v41.py``; it does not
import or execute any file in ``experiments/`` or ``certificates/``.

The finite lattice search is reconstructed here from the positive trace form.
Its rational LDL decomposition is computed at run time.  Final polygon tests
use the cross-multiplied L-valued side form

    |az+b|^2 - 2 h Re((az+b) conjugate(cz+d))
      + (h^2-r^2)|cz+d|^2,

whose sign is the sign of ``|gamma*z-h|^2-r^2`` because the omitted
denominator is totally positive.  Thus no division or floating-point
arithmetic occurs in the decisive exclusions.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import ceil, floor, isqrt, lcm
from typing import Callable, Iterable

from quadrilateral_k17_stdlib_v41 import (
    DEGREE,
    EMBEDDINGS,
    Interval,
    LElement,
    ONE,
    T,
    ZERO,
    conjugate_by,
    elt,
    field_discriminant,
    fixed_point_invariants,
    interval_eval,
    isolate_real_roots,
    matrix,
    mdet,
    mmul,
    power_trace,
    sqrt_interval,
)


Q = Fraction
Vector = tuple[int, ...]
RationalVector = tuple[Q, ...]
Gamma = tuple[Vector, Vector, Vector, Vector]
ZERO_VECTOR: Vector = (0,) * DEGREE
T_POWERS = tuple(T**j for j in range(DEGREE))


def exact_roots(bits: int = 156) -> tuple[Interval, ...]:
    """Return roots in the embedding order 1,5,11,13,17,19."""
    ascending = isolate_real_roots(bits)
    assert len(ascending) == DEGREE
    by_label = {
        label: Interval(lo, hi)
        for label, (lo, hi) in zip((19, 17, 13, 11, 5, 1), ascending)
    }
    return tuple(by_label[label] for label in EMBEDDINGS)


ROOTS = exact_roots()


def sign_at(value: LElement, embedding: int) -> int:
    if value == ZERO:
        return 0
    return interval_eval(value, ROOTS[embedding]).strict_sign()


def signs(value: LElement) -> tuple[int, ...]:
    return tuple(sign_at(value, i) for i in range(DEGREE))


# Power-basis embedding boxes are evaluated once.  Linear combinations of
# these exact rational intervals are the fast filters in the lattice loops.
BASIS_BOXES = tuple(
    tuple(interval_eval(T_POWERS[j], root) for j in range(DEGREE))
    for root in ROOTS
)

# The trace balls contain many points that are rejected by their coordinate
# boxes.  Store the same certified root intervals with one fixed dyadic
# denominator for this hot filter.  Endpoints are rounded outward, so a fast
# rejection can never discard a genuine point; an unresolved endpoint falls
# back to the exact L/Sturm comparison below.
SCALE = 1 << 132


def floor_scaled(value: Q) -> int:
    return value.numerator * SCALE // value.denominator


def ceil_scaled(value: Q) -> int:
    return -((-value.numerator * SCALE) // value.denominator)


def scaled_interval(value: Interval) -> tuple[int, int]:
    return floor_scaled(value.lo), ceil_scaled(value.hi)


BASIS_SCALED = tuple(
    tuple(scaled_interval(value) for value in embedding_values)
    for embedding_values in BASIS_BOXES
)


def scaled_vector_box(vector: Vector, embedding: int) -> tuple[int, int]:
    lo = hi = 0
    for coefficient, (basis_lo, basis_hi) in zip(vector, BASIS_SCALED[embedding]):
        if coefficient >= 0:
            lo += coefficient * basis_lo
            hi += coefficient * basis_hi
        else:
            lo += coefficient * basis_hi
            hi += coefficient * basis_lo
    return lo, hi


def scaled_add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] + right[0], left[1] + right[1]


def scaled_mul(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products) // SCALE, -((-max(products)) // SCALE)


def scaled_le_rational(value: tuple[int, int], bound: Q) -> bool | None:
    numerator = bound.numerator * SCALE
    denominator = bound.denominator
    if value[1] * denominator <= numerator:
        return True
    if value[0] * denominator > numerator:
        return False
    return None


def vector_element(vector: Vector) -> LElement:
    return LElement(vector)


def integral_vector(value: LElement) -> Vector:
    assert value.integral_power_basis()
    return tuple(coefficient.numerator for coefficient in value.coefficients)


def field_box_membership(vector: Vector, bounds: tuple[Q, ...]) -> bool:
    """Certify |sigma_i(vector)| <= bounds[i] for every i."""
    value: LElement | None = None
    for i, bound in enumerate(bounds):
        lo, hi = scaled_vector_box(vector, i)
        numerator = bound.numerator * SCALE
        denominator = bound.denominator
        if lo * denominator > numerator or hi * denominator < -numerator:
            return False
        # Usually the interval lies strictly within the requested box.  If an
        # endpoint is unresolved, decide the exact algebraic comparison.
        if hi * denominator >= numerator or lo * denominator <= -numerator:
            if value is None:
                value = vector_element(vector)
            right = LElement(bound) - value
            left = LElement(bound) + value
            if sign_at(right, i) < 0 or sign_at(left, i) < 0:
                return False
    return True


def ldl(matrix_data: tuple[tuple[Q, ...], ...]) -> tuple[tuple[tuple[Q, ...], ...], tuple[Q, ...]]:
    """Compute G=L D L^T exactly, with unit lower-triangular L."""
    n = len(matrix_data)
    lower = [[Q(0) for _ in range(n)] for _ in range(n)]
    diagonal = [Q(0) for _ in range(n)]
    for j in range(n):
        lower[j][j] = Q(1)
        diagonal[j] = matrix_data[j][j] - sum(
            lower[j][k] * lower[j][k] * diagonal[k] for k in range(j)
        )
        assert diagonal[j] > 0
        for i in range(j + 1, n):
            numerator = matrix_data[i][j] - sum(
                lower[i][k] * lower[j][k] * diagonal[k] for k in range(j)
            )
            lower[i][j] = numerator / diagonal[j]
    reconstructed = tuple(
        tuple(
            sum(lower[i][k] * diagonal[k] * lower[j][k] for k in range(n))
            for j in range(n)
        )
        for i in range(n)
    )
    assert reconstructed == matrix_data
    return tuple(tuple(row) for row in lower), tuple(diagonal)


TRACE_GRAM = tuple(
    tuple((T ** (i + j)).trace() for j in range(DEGREE))
    for i in range(DEGREE)
)
TRACE_LOWER, TRACE_DIAGONAL = ldl(TRACE_GRAM)


def ceil_sqrt_rational(value: Q) -> int:
    assert value >= 0
    candidate = isqrt(value.numerator // value.denominator)
    while Q(candidate * candidate) < value:
        candidate += 1
    return candidate


def enumerate_trace_ball(
    center: RationalVector,
    radius_squared: Q,
    predicate: Callable[[Vector], bool] | None = None,
) -> list[Vector]:
    """Enumerate every v in Z^6 with Tr((v-center)^2)<=radius_squared.

    Since ``TRACE_GRAM = L D L^T``, the quadratic form is the sum of
    ``D_i ((v_i-center_i)+sum_{j>i} L[j,i](v_j-center_j))^2``.  Descending
    recursion therefore gives exact finite integer intervals at every level.
    """
    coordinates = [0] * DEGREE
    output: list[Vector] = []

    def descend(index: int, remaining: Q) -> None:
        if index < 0:
            vector = tuple(coordinates)
            if predicate is None or predicate(vector):
                output.append(vector)
            return
        already_fixed = sum(
            TRACE_LOWER[j][index] * (Q(coordinates[j]) - center[j])
            for j in range(index + 1, DEGREE)
        )
        target = center[index] - already_fixed
        span = ceil_sqrt_rational(remaining / TRACE_DIAGONAL[index])
        for integer in range(floor(target) - span, ceil(target) + span + 1):
            cost = TRACE_DIAGONAL[index] * (Q(integer) - target) ** 2
            if cost <= remaining:
                coordinates[index] = integer
                descend(index - 1, remaining - cost)

    descend(DEGREE - 1, radius_squared)
    return output


def normalized_vertices():
    raw = (
        matrix(
            elt(0, 3, -1, -4, 0, 1),
            elt(0, -10, -3, 6, 1, -1),
            elt(-1, 0, 0, 1),
            elt(0, -6, 1, 5, 0, -1),
        ),
        matrix(0, -1, 1, T),
        matrix(
            elt(0, 2, -4, -4, 1, 1),
            elt(0, -5, -4, 1, 1),
            elt(-1, 5, 4, -1, -1),
            elt(0, -1, 4, 4, -1, -1),
        ),
        matrix(
            elt(-1, -13, -6, 7, 2, -1),
            elt(0, 6, -5, -9, 1, 2),
            elt(1, -5, 2, 5, 0, -1),
            elt(1, 14, 6, -7, -2, 1),
        ),
    )
    assert all(mdet(vertex) == ONE for vertex in raw)
    assert mmul(mmul(mmul(raw[0], raw[1]), raw[2]), raw[3]) == matrix(1, 0, 0, 1)
    normalizer = matrix(T, 0, 0, 1)
    return tuple(conjugate_by(normalizer, vertex) for vertex in raw)


VERTICES = normalized_vertices()
VERTEX_DATA = tuple(
    (center, height2, center * center + height2)
    for center, height2 in map(fixed_point_invariants, VERTICES)
)
assert all(signs(height2) == (1,) * DEGREE for _, height2, _ in VERTEX_DATA)


@dataclass(frozen=True)
class Side:
    center: LElement
    radius_squared: LElement
    interior_signs: tuple[int, ...]


def construct_sides() -> tuple[Side, ...]:
    answer: list[Side] = []
    for side_index in range(4):
        x0, y20, modulus20 = VERTEX_DATA[side_index]
        x1, _, modulus21 = VERTEX_DATA[(side_index + 1) % 4]
        center = (modulus21 - modulus20) / (2 * (x1 - x0))
        radius2 = (x0 - center) ** 2 + y20

        def equation(vertex_index: int) -> LElement:
            x, _, modulus2 = VERTEX_DATA[vertex_index]
            return modulus2 - 2 * center * x + center * center - radius2

        assert equation(side_index) == ZERO
        assert equation((side_index + 1) % 4) == ZERO
        first = signs(equation((side_index + 2) % 4))
        second = signs(equation((side_index + 3) % 4))
        assert 0 not in first and first == second
        answer.append(Side(center, radius2, first))
    return tuple(answer)


SIDES = construct_sides()


# Strictly outward rational bounds.  The side arcs are hyperbolic geodesics:
# Im(z) has its minimum at an endpoint of each compact arc, while |z|^2 is
# affine in Re(z) along a circle orthogonal to R and hence has endpoint
# extrema.  The maximum principle then reduces both polygon bounds to the
# four exact vertex checks below.
M_LOWER = tuple(Q(value, 100000) for value in (4006, 13933, 8395, 17094, 20761, 6663))
U_UPPER = tuple(Q(value, 100000) for value in (255285, 208594, 17707, 121724, 165248, 191115))

for i, root in enumerate(ROOTS):
    for _, height2, modulus2 in VERTEX_DATA:
        assert M_LOWER[i] < sqrt_interval(interval_eval(height2, root), 160).lo
        assert U_UPPER[i] > sqrt_interval(interval_eval(modulus2, root), 160).hi


PI = T**4 + T**3 - 4 * T**2 - 2 * T + 5
C1 = power_trace(1) + power_trace(2) + power_trace(5)
C22 = power_trace(22) + power_trace(23) + power_trace(26)
Y2 = PI - C1 * C1 / 4
assert C1 == elt(-2, 6, 1, -5, 0, 1)
assert C22 == -C1
assert Y2 == elt(Q(4), Q(21, 4), Q(-11, 4), Q(-13, 4), Q(1, 2), Q(1, 2))
assert PI.norm() == 43 and signs(PI) == (1,) * DEGREE
assert signs(Y2) == (1,) * DEGREE


# Separate entry bounds obtained from D=|cz+d|^2<=y/m and
# N=|az+b|^2<=U^2 D.  Here |Re(z)|<=|z|=sqrt(pi).
DEN_BOUND = tuple(Q(value, 100) for value in (530, 283, 461, 51, 103, 331))
NUM_BOUND = tuple(Q(value, 100) for value in (1353, 590, 82, 62, 170, 632))
A_BOUND = tuple(Q(value, 100) for value in (1204, 530, 46, 1403, 776, 868))
C_BOUND = tuple(Q(value, 100) for value in (472, 254, 259, 1153, 470, 454))
D_BOUND = tuple(Q(value, 100) for value in (1899, 566, 1051, 2419, 371, 675))
B_BOUND = tuple(Q(value, 100) for value in (4846, 1180, 187, 2944, 613, 1289))


def certify_derived_bounds() -> None:
    for i, root in enumerate(ROOTS):
        y = sqrt_interval(interval_eval(Y2, root), 160)
        sqrt_pi = sqrt_interval(interval_eval(PI, root), 160)
        den = sqrt_interval(y / Interval.point(M_LOWER[i]), 160)
        c_abs = sqrt_interval(
            Interval.point(1) / (Interval.point(M_LOWER[i]) * y), 160
        )
        num = Interval.point(U_UPPER[i]) * den
        a_abs = num / y
        d_abs = den + c_abs * sqrt_pi
        b_abs = num + a_abs * sqrt_pi
        assert DEN_BOUND[i] > den.hi
        assert C_BOUND[i] > c_abs.hi
        assert NUM_BOUND[i] > num.hi
        assert A_BOUND[i] > a_abs.hi
        assert D_BOUND[i] > d_abs.hi
        assert B_BOUND[i] > b_abs.hi


certify_derived_bounds()


def sum_of_squares(bounds: tuple[Q, ...]) -> Q:
    return sum((bound * bound for bound in bounds), Q(0))


DEN_RADIUS = sum_of_squares(DEN_BOUND)
NUM_RADIUS = sum_of_squares(NUM_BOUND)


Y2_SCALED = tuple(scaled_interval(interval_eval(Y2, root)) for root in ROOTS)


def scaled_complex_square(
    real: tuple[int, int], coefficient: tuple[int, int], y2: tuple[int, int]
) -> tuple[int, int]:
    return scaled_add(scaled_mul(real, real), scaled_mul(scaled_mul(coefficient, coefficient), y2))


def certified_le(value: LElement, rational_bound: Q, embedding: int) -> bool:
    difference = LElement(rational_bound) - value
    return difference == ZERO or sign_at(difference, embedding) > 0


def bottom_bounded(
    c: Vector,
    d: Vector,
    x: LElement,
    c_boxes: tuple[tuple[int, int], ...],
    x_boxes: tuple[tuple[int, int], ...],
) -> bool:
    c_element = vector_element(c)
    d_element = vector_element(d)
    for i, bound in enumerate(DEN_BOUND):
        c_box = c_boxes[i]
        d_box = scaled_vector_box(d, i)
        quick = scaled_complex_square(
            scaled_add(scaled_mul(c_box, x_boxes[i]), d_box), c_box, Y2_SCALED[i]
        )
        bound2 = bound * bound
        decision = scaled_le_rational(quick, bound2)
        if decision is False:
            return False
        if decision is None:
            exact = (c_element * x + d_element) ** 2 + c_element**2 * Y2
            if not certified_le(exact, bound2, i):
                return False
    return True


def top_bounded(
    a: Vector,
    b: Vector,
    x: LElement,
    a_boxes: tuple[tuple[int, int], ...],
    x_boxes: tuple[tuple[int, int], ...],
) -> bool:
    a_element = vector_element(a)
    b_element = vector_element(b)
    for i, bound in enumerate(NUM_BOUND):
        a_box = a_boxes[i]
        b_box = scaled_vector_box(b, i)
        quick = scaled_complex_square(
            scaled_add(scaled_mul(a_box, x_boxes[i]), b_box), a_box, Y2_SCALED[i]
        )
        bound2 = bound * bound
        decision = scaled_le_rational(quick, bound2)
        if decision is False:
            return False
        if decision is None:
            exact = (a_element * x + b_element) ** 2 + a_element**2 * Y2
            if not certified_le(exact, bound2, i):
                return False
    return True


@dataclass(frozen=True)
class IntegralAffineMap:
    denominator: int
    rows: tuple[tuple[int, ...], ...]
    offsets: tuple[int, ...]

    def evaluate(self, vector: Vector) -> Vector | None:
        output: list[int] = []
        for row, offset in zip(self.rows, self.offsets):
            numerator = offset
            for coefficient, entry in zip(vector, row):
                numerator += coefficient * entry
            if numerator % self.denominator:
                return None
            output.append(numerator // self.denominator)
        return tuple(output)


def determinant_one_affine_map(c: Vector, d: Vector) -> IntegralAffineMap:
    """Return the coefficient map a -> b=(ad-1)/c."""
    inverse_c = vector_element(c).inverse()
    multiplier = vector_element(d) * inverse_c
    offset = -inverse_c
    columns = tuple(power * multiplier for power in T_POWERS)
    denominator = 1
    for value in (offset,) + columns:
        for coefficient in value.coefficients:
            denominator = lcm(denominator, coefficient.denominator)
    rows = tuple(
        tuple(int(columns[column].coefficients[row] * denominator) for column in range(DEGREE))
        for row in range(DEGREE)
    )
    offsets = tuple(int(coefficient * denominator) for coefficient in offset.coefficients)
    return IntegralAffineMap(denominator, rows, offsets)


def side_numerator(gamma: Gamma, x: LElement, side: Side) -> LElement:
    """Cross-multiplied side equation at gamma*(x+i*sqrt(Y2))."""
    a, b, c, d = map(vector_element, gamma)
    top_real = a * x + b
    bottom_real = c * x + d
    denominator = bottom_real**2 + c**2 * Y2
    mixed_real = top_real * bottom_real + a * c * Y2
    top_norm = top_real**2 + a**2 * Y2
    return (
        top_norm
        - 2 * side.center * mixed_real
        + (side.center**2 - side.radius_squared) * denominator
    )


def first_polygon_exclusion(gamma: Gamma, x: LElement) -> tuple[int, int] | None:
    for side_index, side in enumerate(SIDES):
        value = side_numerator(gamma, x, side)
        if value == ZERO:
            continue
        for embedding, wanted in enumerate(side.interior_signs):
            if sign_at(value, embedding) != wanted:
                return side_index, embedding
    return None


def canonical_digest(records: Iterable[tuple[Gamma, tuple[int, int]]]) -> str:
    serialized: list[bytes] = []
    for gamma, failure in sorted(records):
        coefficients = ",".join(str(coefficient) for entry in gamma for coefficient in entry)
        serialized.append(f"{coefficients};{failure[0]},{failure[1]}\n".encode("ascii"))
    digest = sha256()
    for record in serialized:
        digest.update(record)
    return digest.hexdigest()


@dataclass(frozen=True)
class Result:
    bottom_rows: int
    nonzero_candidates: int
    zero_candidates: int
    histogram: tuple[tuple[tuple[int, int], int], ...]
    digest: str


EXPECTED_HISTOGRAM = (
    ((0, 0), 1932),
    ((0, 1), 5502),
    ((0, 2), 966),
    ((0, 3), 546),
    ((1, 0), 84),
    ((1, 1), 42),
)
EXPECTED_DIGESTS = {
    1: "fac4b968a1d907fd25d4877ae7373a930a1ecb23ee6685926a26e9d686cfe22d",
    22: "6a08e191a62fdc81807e1c6447f83c7a0afa78a8ee170aa4edcf92752b73951a",
}


def enumerate_case(k: int, x: LElement, a_values: list[Vector], c_values: list[Vector]) -> Result:
    bottom_rows: list[tuple[Vector, Vector]] = []
    x_boxes = tuple(scaled_interval(interval_eval(x, root)) for root in ROOTS)
    a_box_cache = {
        a: tuple(scaled_vector_box(a, i) for i in range(DEGREE))
        for a in a_values
    }
    for c_index, c in enumerate(c_values):
        center = tuple(-(vector_element(c) * x).coefficients[i] for i in range(DEGREE))
        c_boxes = tuple(scaled_vector_box(c, i) for i in range(DEGREE))
        bottom_rows.extend(
            (c, d)
            for d in enumerate_trace_ball(
                center,
                DEN_RADIUS,
                lambda d, c=c, c_boxes=c_boxes: bottom_bounded(c, d, x, c_boxes, x_boxes),
            )
        )
        if c_index and c_index % 100 == 0:
            print(f"  k={k}: c entries {c_index}/{len(c_values)}; bottom rows={len(bottom_rows)}", flush=True)

    print(f"  k={k}: coupled bottom rows={len(bottom_rows)}", flush=True)
    nonzero: list[Gamma] = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != ZERO_VECTOR]
    for row_index, (c, d) in enumerate(nonzero_rows):
        affine = determinant_one_affine_map(c, d)
        for a in a_values:
            b = affine.evaluate(a)
            if b is not None and top_bounded(a, b, x, a_box_cache[a], x_boxes):
                gamma = (a, b, c, d)
                assert vector_element(a) * vector_element(d) - vector_element(b) * vector_element(c) == ONE
                nonzero.append(gamma)
        if row_index and row_index % 250 == 0:
            print(f"  k={k}: determinant rows {row_index}/{len(nonzero_rows)}; candidates={len(nonzero)}", flush=True)

    zero: list[Gamma] = []
    for c, d in bottom_rows:
        if c != ZERO_VECTOR or d == ZERO_VECTOR:
            continue
        a_element = vector_element(d).inverse()
        if not a_element.integral_power_basis():
            continue
        a = integral_vector(a_element)
        if not field_box_membership(a, A_BOUND):
            continue
        center = tuple(-(a_element * x).coefficients[i] for i in range(DEGREE))
        for b in enumerate_trace_ball(
            center,
            NUM_RADIUS,
            lambda b, a=a: top_bounded(a, b, x, a_box_cache[a], x_boxes),
        ):
            assert a_element * vector_element(d) == ONE
            zero.append((a, b, c, d))

    records: list[tuple[Gamma, tuple[int, int]]] = []
    for gamma in nonzero + zero:
        failure = first_polygon_exclusion(gamma, x)
        assert failure is not None
        records.append((gamma, failure))

    histogram_data: dict[tuple[int, int], int] = {}
    for _, failure in records:
        histogram_data[failure] = histogram_data.get(failure, 0) + 1
    histogram = tuple(sorted(histogram_data.items()))
    digest = canonical_digest(records)
    result = Result(len(bottom_rows), len(nonzero), len(zero), histogram, digest)

    assert result.bottom_rows == 1849
    assert result.nonzero_candidates == 9072
    assert result.zero_candidates == 0
    assert result.histogram == EXPECTED_HISTOGRAM
    assert sum(count for _, count in result.histogram) == 9072
    assert result.digest == EXPECTED_DIGESTS[k]
    return result


def main() -> None:
    assert field_discriminant() == 453789
    assert TRACE_GRAM[0][0] == DEGREE
    print("INDEPENDENT RESIDUAL ORBIT k={1,22} STURM/LDL CHECK")
    print("stdlib field/Sturm primitives; dynamically reconstructed exact trace LDL")
    print("cross-multiplied L-valued side predicates; no binary floating point")

    a_values = enumerate_trace_ball(
        (Q(0),) * DEGREE,
        sum_of_squares(A_BOUND),
        lambda vector: field_box_membership(vector, A_BOUND),
    )
    c_values = enumerate_trace_ball(
        (Q(0),) * DEGREE,
        sum_of_squares(C_BOUND),
        lambda vector: field_box_membership(vector, C_BOUND),
    )
    d_values = enumerate_trace_ball(
        (Q(0),) * DEGREE,
        sum_of_squares(D_BOUND),
        lambda vector: field_box_membership(vector, D_BOUND),
    )
    a_values.sort()
    c_values.sort()
    d_values.sort()
    separate_counts = (len(a_values), len(c_values), len(d_values))
    print("separate exact entry counts a,c,d:", separate_counts)
    assert separate_counts == (2641, 715, 65009)

    result1 = enumerate_case(1, -C1 / 2, a_values, c_values)
    result22 = enumerate_case(22, -C22 / 2, a_values, c_values)
    for k, result in ((1, result1), (22, result22)):
        print(
            f"k={k}: bottom rows={result.bottom_rows}; "
            f"c!=0 candidates={result.nonzero_candidates}; "
            f"c=0 candidates={result.zero_candidates}; polygon survivors=0"
        )
        print(f"k={k}: first-exclusion histogram={dict(result.histogram)}")
        print(f"k={k}: canonical exclusion SHA-256={result.digest}")
    assert result1.histogram == result22.histogram
    print("PASS: independent exact exhaustive exclusion reproduced for k=1 and k=22")


if __name__ == "__main__":
    main()
