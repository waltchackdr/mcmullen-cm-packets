#!/usr/bin/env python3
"""Exact exclusion certificate for the residual primitive pair k={17,38}.

The human-readable completeness proof is in
``docs/RESIDUAL_17_38_EXCLUSION.md``.  No binary floating-point value
participates in an assertion.

The arithmetic core uses the integral power basis of O_L=Z[t], the exact
positive trace form Tr_{L/Q}(x^2), and a rational LDL decomposition to
enumerate every lattice point in the required archimedean boxes.  The final
quadrilateral predicates are exact elements of L: for z in H, both Re(z) and
|z|^2 lie in L for the CM points considered here, so every geodesic-circle
side test has a certified sign at all six real embeddings.
"""

from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from math import ceil, floor, isqrt, lcm
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
t = sp.symbols("t")
f = sp.Poly(t**6 + t**5 - 6*t**4 - 6*t**3 + 8*t**2 + 8*t + 1,
            t, domain=sp.QQ)
DEGREE = 6
EMBEDDINGS = (1, 5, 11, 13, 17, 19)
ZERO = (Q(0),) * DEGREE
ONE = (Q(1),) + (Q(0),) * (DEGREE - 1)
T = (Q(0), Q(1), Q(0), Q(0), Q(0), Q(0))

# t^6 = -1 - 8t - 8t^2 + 6t^3 + 6t^4 - t^5.
REDUCTION = (Q(-1), Q(-8), Q(-8), Q(6), Q(6), Q(-1))


def element(*coefficients: object) -> tuple[Q, ...]:
    out = [Q(x) for x in coefficients]
    out.extend([Q(0)] * (DEGREE - len(out)))
    if len(out) != DEGREE:
        raise ValueError("element has degree at least six")
    return tuple(out)


def add(x: tuple[Q, ...], y: tuple[Q, ...]) -> tuple[Q, ...]:
    return tuple(a + b for a, b in zip(x, y))


def neg(x: tuple[Q, ...]) -> tuple[Q, ...]:
    return tuple(-a for a in x)


def sub(x: tuple[Q, ...], y: tuple[Q, ...]) -> tuple[Q, ...]:
    return add(x, neg(y))


def scale(q: object, x: tuple[Q, ...]) -> tuple[Q, ...]:
    scalar = Q(q)
    return tuple(scalar * a for a in x)


def mul(x: tuple[Q, ...], y: tuple[Q, ...]) -> tuple[Q, ...]:
    out = [Q(0)] * 11
    for i, xi in enumerate(x):
        if xi:
            for j, yj in enumerate(y):
                if yj:
                    out[i + j] += xi * yj
    for degree in range(10, 5, -1):
        q = out[degree]
        if not q:
            continue
        out[degree] = Q(0)
        shift = degree - 6
        for j, rj in enumerate(REDUCTION):
            out[shift + j] += q * rj
    return tuple(out[:DEGREE])


def power(x: tuple[Q, ...], exponent: int) -> tuple[Q, ...]:
    answer = ONE
    base = x
    n = exponent
    while n:
        if n & 1:
            answer = mul(answer, base)
        base = mul(base, base)
        n >>= 1
    return answer


def to_sympy(x: tuple[Q, ...]) -> sp.Poly:
    return sp.Poly(sum(sp.Rational(a.numerator, a.denominator) * t**i
                       for i, a in enumerate(x)), t, domain=sp.QQ)


def from_sympy(x: sp.Poly) -> tuple[Q, ...]:
    p = x.rem(f)
    return tuple(Q(int(p.nth(i).p), int(p.nth(i).q)) for i in range(DEGREE))


@lru_cache(maxsize=None)
def inverse(x: tuple[Q, ...]) -> tuple[Q, ...]:
    if x == ZERO:
        raise ZeroDivisionError
    return from_sympy(sp.invert(to_sympy(x), f))


def div(x: tuple[Q, ...], y: tuple[Q, ...]) -> tuple[Q, ...]:
    return mul(x, inverse(y))


def integral(x: tuple[Q, ...]) -> bool:
    return all(a.denominator == 1 for a in x)


def integral_tuple(x: tuple[Q, ...]) -> tuple[int, ...]:
    if not integral(x):
        raise ArithmeticError("nonintegral field element")
    return tuple(a.numerator for a in x)


def as_field(x: tuple[int, ...]) -> tuple[Q, ...]:
    return tuple(Q(a) for a in x)


def trace(x: tuple[Q, ...]) -> Q:
    # Newton sums for f, equivalently the first row of the trace Gram matrix.
    return sum(q * a for q, a in zip((6, -1, 13, -1, 41, -1), x))


def norm(x: tuple[Q, ...]) -> Q:
    return Q(sp.resultant(f.as_expr(), to_sympy(x).as_expr(), t))


# Exact rational isolating intervals, ordered as 2*cos(r*pi/21) for
# r=1,5,11,13,17,19.  SymPy's intervals routine uses exact Sturm arithmetic.
ROOT_INTERVALS_RAW = [ab for ab, _ in sp.polys.polytools.intervals(
    f, eps=sp.Rational(1, 10)**45)]
ROOT_INTERVALS_RAW.sort(key=lambda ab: ab[0], reverse=True)
assert len(ROOT_INTERVALS_RAW) == DEGREE
ROOT_INTERVALS = tuple(
    (Q(int(a.p), int(a.q)), Q(int(b.p), int(b.q)))
    for a, b in ROOT_INTERVALS_RAW
)


def interval_add(x: tuple[Q, Q], y: tuple[Q, Q]) -> tuple[Q, Q]:
    return x[0] + y[0], x[1] + y[1]


def interval_neg(x: tuple[Q, Q]) -> tuple[Q, Q]:
    return -x[1], -x[0]


def interval_mul(x: tuple[Q, Q], y: tuple[Q, Q]) -> tuple[Q, Q]:
    products = (x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
    return min(products), max(products)


def interval_div(x: tuple[Q, Q], y: tuple[Q, Q]) -> tuple[Q, Q]:
    if y[0] <= 0 <= y[1]:
        raise ZeroDivisionError("interval denominator meets zero")
    reciprocal = (min(1 / y[0], 1 / y[1]), max(1 / y[0], 1 / y[1]))
    return interval_mul(x, reciprocal)


def interval_eval(x: tuple[Q, ...], root: tuple[Q, Q]) -> tuple[Q, Q]:
    value = (Q(0), Q(0))
    for coefficient in reversed(x):
        value = interval_add(interval_mul(value, root), (coefficient, coefficient))
    return value


def sqrt_interval(x: tuple[Q, Q], bits: int = 150) -> tuple[Q, Q]:
    if x[0] <= 0:
        raise ArithmeticError(f"square-root interval is not positive: {x}")
    denominator = 1 << bits

    def lower(q: Q) -> Q:
        n = q.numerator * denominator * denominator // q.denominator
        return Q(isqrt(n), denominator)

    def upper(q: Q) -> Q:
        n = q.numerator * denominator * denominator // q.denominator
        r = isqrt(n)
        if Q(r * r, denominator * denominator) < q:
            r += 1
        return Q(r, denominator)

    return lower(x[0]), upper(x[1])


def sign_at(x: tuple[Q, ...], embedding_index: int) -> int:
    if x == ZERO:
        return 0
    lo, hi = interval_eval(x, ROOT_INTERVALS[embedding_index])
    if lo > 0:
        return 1
    if hi < 0:
        return -1
    raise ArithmeticError(("uncertified sign", x, EMBEDDINGS[embedding_index], lo, hi))


def sign_vector(x: tuple[Q, ...]) -> tuple[int, ...]:
    return tuple(sign_at(x, i) for i in range(DEGREE))


# Fixed-scale outward boxes for fast exact filtering in the hot lattice loops.
# Every endpoint below is an integer divided by SCALE.
SCALE_BITS = 120
SCALE = 1 << SCALE_BITS


def floor_scaled(q: Q) -> int:
    return q.numerator * SCALE // q.denominator


def ceil_scaled(q: Q) -> int:
    return -((-q.numerator * SCALE) // q.denominator)


BASIS_BOXES: tuple[tuple[tuple[int, int], ...], ...] = tuple(
    tuple((floor_scaled(interval_eval(power(T, j), root)[0]),
           ceil_scaled(interval_eval(power(T, j), root)[1]))
          for j in range(DEGREE))
    for root in ROOT_INTERVALS
)


def scaled_eval_integral(x: tuple[int, ...], i: int) -> tuple[int, int]:
    lo = hi = 0
    for coefficient, (basis_lo, basis_hi) in zip(x, BASIS_BOXES[i]):
        if coefficient >= 0:
            lo += coefficient * basis_lo
            hi += coefficient * basis_hi
        else:
            lo += coefficient * basis_hi
            hi += coefficient * basis_lo
    return lo, hi


def scaled_field_box(x: tuple[Q, ...], i: int) -> tuple[int, int]:
    lo, hi = interval_eval(x, ROOT_INTERVALS[i])
    return floor_scaled(lo), ceil_scaled(hi)


def scaled_add(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return x[0] + y[0], x[1] + y[1]


def scaled_mul(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    products = (x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
    lo = min(products) // SCALE
    hi = -((-max(products)) // SCALE)
    return lo, hi


def scaled_le(box: tuple[int, int], bound: Q) -> bool | None:
    threshold_num = bound.numerator * SCALE
    threshold_den = bound.denominator
    if box[1] * threshold_den <= threshold_num:
        return True
    if box[0] * threshold_den > threshold_num:
        return False
    return None


def within_integral_box(x: tuple[int, ...], bounds: tuple[Q, ...]) -> bool:
    for i, bound in enumerate(bounds):
        lo, hi = scaled_eval_integral(x, i)
        threshold_num = bound.numerator * SCALE
        threshold_den = bound.denominator
        if lo * threshold_den > threshold_num or hi * threshold_den < -threshold_num:
            return False
        if lo * threshold_den <= -threshold_num <= hi * threshold_den:
            raise ArithmeticError(("left box boundary unresolved", x, i, bound))
        if lo * threshold_den <= threshold_num <= hi * threshold_den:
            raise ArithmeticError(("right box boundary unresolved", x, i, bound))
    return True


# Trace Gram matrix G_ij=Tr(t^(i+j)) and the exact decomposition G=L D L^T.
G = (
    (Q(6), Q(-1), Q(13), Q(-1), Q(41), Q(-1)),
    (Q(-1), Q(13), Q(-1), Q(41), Q(-1), Q(139)),
    (Q(13), Q(-1), Q(41), Q(-1), Q(139), Q(6)),
    (Q(-1), Q(41), Q(-1), Q(139), Q(6), Q(489)),
    (Q(41), Q(-1), Q(139), Q(6), Q(489), Q(62)),
    (Q(-1), Q(139), Q(6), Q(489), Q(62), Q(1763)),
)
L = (
    (Q(1), Q(0), Q(0), Q(0), Q(0), Q(0)),
    (Q(-1, 6), Q(1), Q(0), Q(0), Q(0), Q(0)),
    (Q(13, 6), Q(1, 11), Q(1), Q(0), Q(0), Q(0)),
    (Q(-1, 6), Q(35, 11), Q(-1, 5), Q(1), Q(0), Q(0)),
    (Q(41, 6), Q(5, 11), Q(39, 10), Q(1, 2), Q(1), Q(0)),
    (Q(-1, 6), Q(119, 11), Q(-7, 20), Q(11, 2), Q(0), Q(1)),
)
DTRACE = (Q(6), Q(77, 6), Q(140, 11), Q(42, 5), Q(21, 2), Q(21, 4))


def matmul_q(a: tuple[tuple[Q, ...], ...], b: tuple[tuple[Q, ...], ...]) -> tuple[tuple[Q, ...], ...]:
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(DEGREE))
                       for j in range(DEGREE)) for i in range(DEGREE))


def transpose_q(a: tuple[tuple[Q, ...], ...]) -> tuple[tuple[Q, ...], ...]:
    return tuple(tuple(a[j][i] for j in range(DEGREE)) for i in range(DEGREE))


LD = tuple(tuple(L[i][j] * DTRACE[j] for j in range(DEGREE)) for i in range(DEGREE))
assert matmul_q(LD, transpose_q(L)) == G
assert int(sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row]
                      for row in G]).det()) == 453789
assert int(sp.discriminant(f.as_expr(), t)) == 453789


def ceil_sqrt_fraction(q: Q) -> int:
    if q < 0:
        raise ValueError(q)
    n, d = q.numerator, q.denominator
    r = isqrt(n // d)
    while Q(r * r) < q:
        r += 1
    return r


def enumerate_trace_ellipsoid(center: tuple[Q, ...], radius: Q,
                               predicate=None) -> list[tuple[int, ...]]:
    """Enumerate all v in Z^6 with Tr((v-center)^2)<=radius exactly."""
    coefficients = [0] * DEGREE
    answer: list[tuple[int, ...]] = []

    def visit(i: int, remaining: Q) -> None:
        if i < 0:
            candidate = tuple(coefficients)
            if predicate is None or predicate(candidate):
                answer.append(candidate)
            return
        tail = sum(L[j][i] * (Q(coefficients[j]) - center[j])
                   for j in range(i + 1, DEGREE))
        target = center[i] - tail
        span = ceil_sqrt_fraction(remaining / DTRACE[i])
        lo = floor(target) - span - 1
        hi = ceil(target) + span + 1
        for value in range(lo, hi + 1):
            term = DTRACE[i] * (Q(value) - target)**2
            if term <= remaining:
                coefficients[i] = value
                visit(i - 1, remaining - term)

    visit(DEGREE - 1, radius)
    return answer


def power_trace(n: int) -> tuple[Q, ...]:
    if n == 0:
        return scale(2, ONE)
    previous, current = scale(2, ONE), T
    for _ in range(1, n):
        previous, current = current, sub(mul(T, current), previous)
    return current


# Exact quadrilateral matrices in McMullen's cyclic A,B,C,D order.
Matrix = tuple[tuple[tuple[Q, ...], tuple[Q, ...]],
               tuple[tuple[Q, ...], tuple[Q, ...]]]


def matrix(a: tuple[Q, ...], b: tuple[Q, ...],
           c: tuple[Q, ...], d: tuple[Q, ...]) -> Matrix:
    return ((a, b), (c, d))


def mtrace(a: Matrix) -> tuple[Q, ...]:
    return add(a[0][0], a[1][1])


def mdet(a: Matrix) -> tuple[Q, ...]:
    return sub(mul(a[0][0], a[1][1]), mul(a[0][1], a[1][0]))


def mmul(a: Matrix, b: Matrix) -> Matrix:
    return ((add(mul(a[0][0], b[0][0]), mul(a[0][1], b[1][0])),
             add(mul(a[0][0], b[0][1]), mul(a[0][1], b[1][1]))),
            (add(mul(a[1][0], b[0][0]), mul(a[1][1], b[1][0])),
             add(mul(a[1][0], b[0][1]), mul(a[1][1], b[1][1]))))


I2 = matrix(ONE, ZERO, ZERO, ONE)


def normalized_by_t(a: Matrix) -> Matrix:
    return matrix(a[0][0], mul(T, a[0][1]),
                  mul(inverse(T), a[1][0]), a[1][1])


A = matrix(
    element(0, 3, -1, -4, 0, 1),
    element(0, -10, -3, 6, 1, -1),
    element(-1, 0, 0, 1),
    element(0, -6, 1, 5, 0, -1),
)
B = matrix(ZERO, neg(ONE), ONE, T)
C = matrix(
    element(0, 2, -4, -4, 1, 1),
    element(0, -5, -4, 1, 1),
    element(-1, 5, 4, -1, -1),
    element(0, -1, 4, 4, -1, -1),
)
DVERTEX = matrix(
    element(-1, -13, -6, 7, 2, -1),
    element(0, 6, -5, -9, 1, 2),
    element(1, -5, 2, 5, 0, -1),
    element(1, 14, 6, -7, -2, 1),
)
RAW_VERTICES = (A, B, C, DVERTEX)
VERTEX_NAMES = ("A", "B", "C", "D")
VERTICES = tuple(normalized_by_t(a) for a in RAW_VERTICES)
assert all(mdet(a) == ONE for a in RAW_VERTICES + VERTICES)
assert mmul(mmul(mmul(A, B), C), DVERTEX) == I2


def fixed_invariants(a: Matrix) -> tuple[tuple[Q, ...], tuple[Q, ...], tuple[Q, ...]]:
    aa, cc, dd = a[0][0], a[1][0], a[1][1]
    x = div(sub(aa, dd), scale(2, cc))
    y2 = div(sub(scale(4, ONE), mul(mtrace(a), mtrace(a))),
             scale(4, mul(cc, cc)))
    modulus2 = add(mul(x, x), y2)
    return x, y2, modulus2


VERTEX_INVARIANTS = tuple(fixed_invariants(a) for a in VERTICES)
assert all(sign_vector(y2) == (1,) * DEGREE for _, y2, _ in VERTEX_INVARIANTS)


# Each side is represented by (real circle center, radius squared).  The two
# nonincident vertices must lie strictly on the same side in every embedding;
# this certifies the cyclic convex polygon encoded by A,B,C,D.
SIDES: list[tuple[tuple[Q, ...], tuple[Q, ...], tuple[int, ...]]] = []
for j in range(4):
    x1, y21, m21 = VERTEX_INVARIANTS[j]
    x2, _, m22 = VERTEX_INVARIANTS[(j + 1) % 4]
    circle_center = div(sub(m22, m21), scale(2, sub(x2, x1)))
    radius2 = add(mul(sub(x1, circle_center), sub(x1, circle_center)), y21)

    def side_value(point_index: int) -> tuple[Q, ...]:
        xp, _, mp = VERTEX_INVARIANTS[point_index]
        return sub(add(sub(mp, scale(2, mul(circle_center, xp))),
                       mul(circle_center, circle_center)), radius2)

    assert side_value(j) == ZERO and side_value((j + 1) % 4) == ZERO
    interior1 = sign_vector(side_value((j + 2) % 4))
    interior2 = sign_vector(side_value((j + 3) % 4))
    assert 0 not in interior1 and interior1 == interior2
    SIDES.append((circle_center, radius2, interior1))


# Rational lower/upper bounds for the six normalized coordinate polygons.
M_LOWER = tuple(Q(x, 100000) for x in (4006, 13933, 8395, 17094, 20761, 6663))
U_UPPER = tuple(Q(x, 100000) for x in (255285, 208594, 17707, 121724, 165248, 191115))

for i in range(DEGREE):
    for _, y2, modulus2 in VERTEX_INVARIANTS:
        y_box = sqrt_interval(interval_eval(y2, ROOT_INTERVALS[i]))
        modulus_box = sqrt_interval(interval_eval(modulus2, ROOT_INTERVALS[i]))
        assert M_LOWER[i] < y_box[0]
        assert U_UPPER[i] > modulus_box[1]


pi = element(5, -2, -4, 1, 1)
c17 = add(add(power_trace(17), power_trace(18)), power_trace(21))
c38 = add(add(power_trace(38), power_trace(39)), power_trace(42))
assert c17 == element(-4, 3, 4, -1, -1)
assert c38 == neg(c17)
assert norm(pi) == 43
assert sign_vector(pi) == (1,) * DEGREE
Y2 = sub(pi, scale(Q(1, 4), mul(c17, c17)))
assert Y2 == scale(Q(1, 4), add(T, scale(2, ONE)))
assert sign_vector(Y2) == (1,) * DEGREE


# Certified outward bounds.  The shorter B and N arrays are the coupled
# denominator/numerator bounds; A,C,DENTRY,BENTRY are the derived separate
# entry bounds used only where needed.
DEN_BOUND = tuple(Q(x, 100) for x in (500, 259, 285, 182, 120, 150))
NUM_BOUND = tuple(Q(x, 100) for x in (1274, 540, 51, 222, 198, 287))
A_BOUND = tuple(Q(x, 100) for x in (1278, 580, 75, 393, 669, 1919))
C_BOUND = tuple(Q(x, 100) for x in (501, 278, 419, 323, 405, 1004))
D_BOUND = tuple(Q(x, 100) for x in (1952, 568, 1240, 844, 351, 909))
B_BOUND = tuple(Q(x, 100) for x in (4983, 1185, 220, 1027, 579, 1738))

for i in range(DEGREE):
    y_lo, y_hi = sqrt_interval(interval_eval(Y2, ROOT_INTERVALS[i]))
    pi_lo, pi_hi = sqrt_interval(interval_eval(pi, ROOT_INTERVALS[i]))
    b_actual_upper = sqrt_interval((y_hi / M_LOWER[i], y_hi / M_LOWER[i]))[1]
    c_actual_upper = sqrt_interval((1 / (M_LOWER[i] * y_lo),
                                    1 / (M_LOWER[i] * y_lo)))[1]
    n_actual_upper = U_UPPER[i] * b_actual_upper
    a_actual_upper = n_actual_upper / y_lo
    d_actual_upper = b_actual_upper + c_actual_upper * pi_hi
    bentry_actual_upper = n_actual_upper + a_actual_upper * pi_hi
    assert DEN_BOUND[i] > b_actual_upper
    assert NUM_BOUND[i] > n_actual_upper
    assert A_BOUND[i] > a_actual_upper
    assert C_BOUND[i] > c_actual_upper
    assert D_BOUND[i] > d_actual_upper
    assert B_BOUND[i] > bentry_actual_upper


def comparison_le_at(x: tuple[Q, ...], bound: Q, i: int) -> bool:
    lo, hi = interval_eval(x, ROOT_INTERVALS[i])
    if hi <= bound:
        return True
    if lo > bound:
        return False
    raise ArithmeticError(("comparison unresolved", x, bound, EMBEDDINGS[i], lo, hi))


def all_le_bounds(x: tuple[Q, ...], bounds: tuple[Q, ...]) -> bool:
    return all(comparison_le_at(x, bound, i) for i, bound in enumerate(bounds))


def square_bounds(bounds: tuple[Q, ...]) -> tuple[Q, ...]:
    return tuple(x * x for x in bounds)


def denominator_square(c: tuple[Q, ...], d: tuple[Q, ...], x: tuple[Q, ...]) -> tuple[Q, ...]:
    real = add(mul(c, x), d)
    return add(mul(real, real), mul(mul(c, c), Y2))


def numerator_square(a: tuple[Q, ...], b: tuple[Q, ...], x: tuple[Q, ...]) -> tuple[Q, ...]:
    real = add(mul(a, x), b)
    return add(mul(real, real), mul(mul(a, a), Y2))


def affine_b_data(c: tuple[int, ...], d: tuple[int, ...]):
    """Integer-scaled affine map a -> (ad-1)/c for a fixed bottom row."""
    inv_c = inverse(as_field(c))
    multiplier = mul(as_field(d), inv_c)
    offset = neg(inv_c)
    columns = [mul(power(T, j), multiplier) for j in range(DEGREE)]
    denominator = 1
    for q in offset:
        denominator = lcm(denominator, q.denominator)
    for column in columns:
        for q in column:
            denominator = lcm(denominator, q.denominator)
    rows = tuple(tuple(int(columns[j][i] * denominator) for j in range(DEGREE))
                 for i in range(DEGREE))
    offsets = tuple(int(q * denominator) for q in offset)
    return denominator, rows, offsets


def affine_integral_b_fast(a: tuple[int, ...], data) -> tuple[int, ...] | None:
    denominator, rows, offsets = data
    answer = []
    for row, offset in zip(rows, offsets):
        value = offset + sum(q * coefficient for q, coefficient in zip(row, a))
        if value % denominator:
            return None
        answer.append(value // denominator)
    return tuple(answer)


def image_invariants(gamma: tuple[tuple[int, ...], tuple[int, ...],
                                  tuple[int, ...], tuple[int, ...]],
                     x: tuple[Q, ...]) -> tuple[tuple[Q, ...], tuple[Q, ...]]:
    aq, bq, cq, dq = map(as_field, gamma)
    cx_d = add(mul(cq, x), dq)
    ax_b = add(mul(aq, x), bq)
    denominator = add(mul(cx_d, cx_d), mul(mul(cq, cq), Y2))
    real_numerator = add(mul(ax_b, cx_d), mul(mul(aq, cq), Y2))
    modulus_numerator = add(mul(ax_b, ax_b), mul(mul(aq, aq), Y2))
    return div(real_numerator, denominator), div(modulus_numerator, denominator)


def in_all_coordinate_polygons(gamma, x: tuple[Q, ...]) -> tuple[bool, tuple[int, int] | None]:
    real, modulus2 = image_invariants(gamma, x)
    for side_index, (circle_center, radius2, interior_signs) in enumerate(SIDES):
        value = sub(add(sub(modulus2, scale(2, mul(circle_center, real))),
                        mul(circle_center, circle_center)), radius2)
        if value == ZERO:
            continue
        for i, wanted in enumerate(interior_signs):
            if sign_at(value, i) != wanted:
                return False, (side_index, i)
    return True, None


def enumerate_case(k: int, x: tuple[Q, ...], a_values: list[tuple[int, ...]],
                   c_values: list[tuple[int, ...]]):
    den_squares = square_bounds(DEN_BOUND)
    num_squares = square_bounds(NUM_BOUND)
    den_radius = sum(den_squares, Q(0))
    num_radius = sum(num_squares, Q(0))
    bottom_rows: list[tuple[tuple[int, ...], tuple[int, ...]]] = []

    x_boxes = tuple(scaled_field_box(x, i) for i in range(DEGREE))
    y2_boxes = tuple(scaled_field_box(Y2, i) for i in range(DEGREE))
    a_boxes_cache = {
        a: tuple(scaled_eval_integral(a, i) for i in range(DEGREE))
        for a in a_values
    }

    def numerator_bounded_scaled(a: tuple[int, ...], b: tuple[int, ...]) -> bool:
        unresolved = False
        a_boxes = a_boxes_cache.get(a)
        if a_boxes is None:
            a_boxes = tuple(scaled_eval_integral(a, i) for i in range(DEGREE))
        for i in range(DEGREE):
            b_box = scaled_eval_integral(b, i)
            real_box = scaled_add(scaled_mul(a_boxes[i], x_boxes[i]), b_box)
            square_box = scaled_add(scaled_mul(real_box, real_box),
                                    scaled_mul(scaled_mul(a_boxes[i], a_boxes[i]),
                                               y2_boxes[i]))
            decision = scaled_le(square_box, num_squares[i])
            if decision is False:
                return False
            unresolved |= decision is None
        if not unresolved:
            return True
        return all_le_bounds(numerator_square(as_field(a), as_field(b), x), num_squares)

    for c_index, c in enumerate(c_values):
        cq = as_field(c)
        shifted_center = neg(mul(cq, x))
        c_boxes = tuple(scaled_eval_integral(c, i) for i in range(DEGREE))

        def denominator_predicate(d: tuple[int, ...]) -> bool:
            unresolved = False
            for i in range(DEGREE):
                d_box = scaled_eval_integral(d, i)
                real_box = scaled_add(scaled_mul(c_boxes[i], x_boxes[i]), d_box)
                square_box = scaled_add(scaled_mul(real_box, real_box),
                                        scaled_mul(scaled_mul(c_boxes[i], c_boxes[i]),
                                                   y2_boxes[i]))
                decision = scaled_le(square_box, den_squares[i])
                if decision is False:
                    return False
                unresolved |= decision is None
            if not unresolved:
                return True
            return all_le_bounds(denominator_square(cq, as_field(d), x), den_squares)

        for d in enumerate_trace_ellipsoid(shifted_center, den_radius,
                                           denominator_predicate):
            bottom_rows.append((c, d))
        if c_index and c_index % 100 == 0:
            print(f"  k={k}: processed c entries {c_index}/{len(c_values)}; "
                  f"bottom rows={len(bottom_rows)}", flush=True)

    print(f"  k={k}: exact coupled bottom-row enumeration complete: {len(bottom_rows)}",
          flush=True)

    candidates = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != (0,) * DEGREE]
    for row_index, (c, d) in enumerate(nonzero_rows):
        if c == (0,) * DEGREE:
            continue
        cq, dq = as_field(c), as_field(d)
        affine_data = affine_b_data(c, d)
        for a in a_values:
            b = affine_integral_b_fast(a, affine_data)
            if b is None:
                continue
            if numerator_bounded_scaled(a, b):
                aq, bq = as_field(a), as_field(b)
                assert sub(mul(aq, dq), mul(bq, cq)) == ONE
                candidates.append((a, b, c, d))
        if row_index and row_index % 250 == 0:
            print(f"  k={k}: determinant scan rows {row_index}/{len(nonzero_rows)}; "
                  f"bounded candidates={len(candidates)}", flush=True)

    zero_candidates = []
    for c, d in bottom_rows:
        if c != (0,) * DEGREE:
            continue
        if d == (0,) * DEGREE:
            continue
        dq = as_field(d)
        aq = inverse(dq)
        if not integral(aq):
            continue
        a = integral_tuple(aq)
        if not within_integral_box(a, A_BOUND):
            continue
        shifted_center = neg(mul(aq, x))

        def numerator_predicate(b: tuple[int, ...]) -> bool:
            return numerator_bounded_scaled(a, b)

        for b in enumerate_trace_ellipsoid(shifted_center, num_radius,
                                           numerator_predicate):
            gamma = (a, b, c, d)
            assert mul(aq, dq) == ONE
            zero_candidates.append(gamma)

    failure_histogram: dict[tuple[int, int], int] = {}
    polygon_survivors = []
    exclusion_records = []
    for gamma in candidates + zero_candidates:
        inside, failure = in_all_coordinate_polygons(gamma, x)
        if inside:
            polygon_survivors.append(gamma)
        else:
            assert failure is not None
            failure_histogram[failure] = failure_histogram.get(failure, 0) + 1
            exclusion_records.append((gamma, failure))

    assert not polygon_survivors
    digest = sha256()
    for gamma, failure in sorted(exclusion_records):
        flat = tuple(coefficient for entry in gamma for coefficient in entry)
        digest.update((",".join(map(str, flat)) +
                       f";{failure[0]},{failure[1]}\n").encode("ascii"))
    digest_hex = digest.hexdigest()
    print(f"k={k}: bottom rows={len(bottom_rows)}; c!=0 candidates={len(candidates)}; "
          f"c=0 candidates={len(zero_candidates)}; polygon survivors=0; "
          f"exclusion digest={digest_hex}")
    return (len(bottom_rows), len(candidates), len(zero_candidates),
            failure_histogram, digest_hex)


def main() -> None:
    print("EXACT RESIDUAL ORBIT CERTIFICATE v42")
    print("root intervals: exact SymPy isolation; trace-form determinant=453789")
    print("polygon lower bounds m_i:", M_LOWER)
    print("polygon upper bounds U_i:", U_UPPER)

    a_radius = sum((x * x for x in A_BOUND), Q(0))
    c_radius = sum((x * x for x in C_BOUND), Q(0))
    d_radius = sum((x * x for x in D_BOUND), Q(0))
    a_values = enumerate_trace_ellipsoid(
        ZERO, a_radius, lambda x: within_integral_box(x, A_BOUND))
    c_values = enumerate_trace_ellipsoid(
        ZERO, c_radius, lambda x: within_integral_box(x, C_BOUND))
    d_values = enumerate_trace_ellipsoid(
        ZERO, d_radius, lambda x: within_integral_box(x, D_BOUND))
    print("exact separate entry counts a,c,d:", len(a_values), len(c_values), len(d_values))
    assert (len(a_values), len(c_values), len(d_values)) == (2673, 731, 35181)

    result17 = enumerate_case(17, scale(Q(-1, 2), c17), a_values, c_values)
    result38 = enumerate_case(38, scale(Q(-1, 2), c38), a_values, c_values)
    print("first-failure histogram k=17:", sorted(result17[3].items()))
    print("first-failure histogram k=38:", sorted(result38[3].items()))
    expected_histogram = {
        (0, 0): 2392,
        (0, 1): 5826,
        (0, 2): 1048,
        (0, 3): 504,
        (1, 0): 42,
    }
    assert result17[:3] == result38[:3] == (2307, 9756, 56)
    assert result17[3] == result38[3] == expected_histogram
    assert sum(expected_histogram.values()) == 9756 + 56
    assert result17[4] == "964a7c3d8b2a7740e1c9bd72cd8057e0449371dac87dc3c19807349835e98aad"
    assert result38[4] == "174961ec5c5f06fe176b3c02d372901048dd4e4f90ac78eb260dd67bf0189d10"
    print("PASS: exact exhaustive exclusion for both primitive residual indices")


if __name__ == "__main__":
    main()
