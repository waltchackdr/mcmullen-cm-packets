#!/usr/bin/env python3
"""Independent exact reconstruction of McMullen's quadrilateral and k=17.

Independence boundary
---------------------
This program uses only the Python standard library.  It does not import or
execute repository certificate/experiment code.  The mathematical inputs are
the displayed definitions in the manuscript and the four matrices printed in
Section 4 of McMullen, "Triangle groups and Hilbert modular varieties".

The script deliberately stops short of a full Hilbert-modular orbit test for
k=17.  Its output is an exact arithmetic/fixed-point/orientation fingerprint
which can be consumed by such a test.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt


Q = Fraction
DEGREE = 6
# f(T) = T^6 + T^5 - 6T^4 - 6T^3 + 8T^2 + 8T + 1,
# stored in ascending order.
MODULUS = (Q(1), Q(8), Q(8), Q(-6), Q(-6), Q(1), Q(1))
EMBEDDINGS = (1, 5, 11, 13, 17, 19)


def _trim(poly: list[Q]) -> list[Q]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def _reduce(coefficients: list[Q]) -> tuple[Q, ...]:
    a = list(coefficients)
    if len(a) < DEGREE:
        a.extend([Q(0)] * (DEGREE - len(a)))
    for degree in range(len(a) - 1, DEGREE - 1, -1):
        lead = a[degree]
        if lead:
            shift = degree - DEGREE
            for j in range(DEGREE):
                a[shift + j] -= lead * MODULUS[j]
        a[degree] = Q(0)
    return tuple(a[:DEGREE])


def determinant(matrix: list[list[Q]]) -> Q:
    """Fraction-preserving Gaussian determinant (dimensions here are <= 11)."""
    a = [row[:] for row in matrix]
    n = len(a)
    value = Q(1)
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return Q(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            value = -value
        p = a[col][col]
        value *= p
        for j in range(col, n):
            a[col][j] /= p
        for r in range(col + 1, n):
            factor = a[r][col]
            if factor:
                for j in range(col, n):
                    a[r][j] -= factor * a[col][j]
    return value


def solve(matrix: list[list[Q]], rhs: list[Q]) -> list[Q]:
    a = [row[:] + [b] for row, b in zip(matrix, rhs)]
    n = len(a)
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            raise ZeroDivisionError("singular linear system")
        a[col], a[pivot] = a[pivot], a[col]
        p = a[col][col]
        a[col] = [x / p for x in a[col]]
        for r in range(n):
            if r == col:
                continue
            factor = a[r][col]
            if factor:
                a[r] = [x - factor * y for x, y in zip(a[r], a[col])]
    return [a[i][-1] for i in range(n)]


@dataclass(frozen=True)
class LElement:
    """An element of Q[t]/(f), represented in the power basis."""

    coefficients: tuple[Q, ...]

    def __init__(self, coefficients: object = 0):
        if isinstance(coefficients, LElement):
            c = coefficients.coefficients
        elif isinstance(coefficients, (int, Q)):
            c = (Q(coefficients),)
        else:
            c = tuple(Q(x) for x in coefficients)  # type: ignore[arg-type]
        object.__setattr__(self, "coefficients", _reduce(list(c)))

    @staticmethod
    def coerce(value: object) -> "LElement":
        return value if isinstance(value, LElement) else LElement(value)

    def __add__(self, other: object) -> "LElement":
        b = self.coerce(other)
        return LElement([x + y for x, y in zip(self.coefficients, b.coefficients)])

    __radd__ = __add__

    def __neg__(self) -> "LElement":
        return LElement([-x for x in self.coefficients])

    def __sub__(self, other: object) -> "LElement":
        return self + (-self.coerce(other))

    def __rsub__(self, other: object) -> "LElement":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "LElement":
        b = self.coerce(other)
        product = [Q(0)] * (2 * DEGREE - 1)
        for i, x in enumerate(self.coefficients):
            for j, y in enumerate(b.coefficients):
                product[i + j] += x * y
        return LElement(product)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "LElement":
        if exponent < 0:
            return (self.inverse()) ** (-exponent)
        answer = ONE
        base = self
        n = exponent
        while n:
            if n & 1:
                answer *= base
            base *= base
            n >>= 1
        return answer

    def multiplication_matrix(self) -> list[list[Q]]:
        columns = [(self * (T**j)).coefficients for j in range(DEGREE)]
        return [[columns[j][i] for j in range(DEGREE)] for i in range(DEGREE)]

    def inverse(self) -> "LElement":
        vector = solve(self.multiplication_matrix(), [Q(1)] + [Q(0)] * 5)
        return LElement(vector)

    def __truediv__(self, other: object) -> "LElement":
        return self * self.coerce(other).inverse()

    def __rtruediv__(self, other: object) -> "LElement":
        return self.coerce(other) / self

    def norm(self) -> Q:
        return determinant(self.multiplication_matrix())

    def trace(self) -> Q:
        m = self.multiplication_matrix()
        return sum((m[i][i] for i in range(DEGREE)), Q(0))

    def integral_power_basis(self) -> bool:
        return all(x.denominator == 1 for x in self.coefficients)

    def __str__(self) -> str:
        terms: list[str] = []
        for degree, coefficient in enumerate(self.coefficients):
            if not coefficient:
                continue
            terms.append(f"({coefficient})*t^{degree}")
        return " + ".join(terms) if terms else "0"


ZERO = LElement(0)
ONE = LElement(1)
T = LElement((0, 1))


def elt(*coefficients: int) -> LElement:
    return LElement(coefficients)


Matrix = tuple[tuple[LElement, LElement], tuple[LElement, LElement]]


def matrix(a: object, b: object, c: object, d: object) -> Matrix:
    return (
        (LElement.coerce(a), LElement.coerce(b)),
        (LElement.coerce(c), LElement.coerce(d)),
    )


IDENTITY = matrix(1, 0, 0, 1)
NEG_IDENTITY = matrix(-1, 0, 0, -1)


def madd(a: Matrix, b: Matrix) -> Matrix:
    return matrix(*(a[i][j] + b[i][j] for i in range(2) for j in range(2)))


def msub(a: Matrix, b: Matrix) -> Matrix:
    return matrix(*(a[i][j] - b[i][j] for i in range(2) for j in range(2)))


def mmul(a: Matrix, b: Matrix) -> Matrix:
    return (
        (
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ),
        (
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ),
    )


def mpow(a: Matrix, exponent: int) -> Matrix:
    answer = IDENTITY
    base = a
    n = exponent
    while n:
        if n & 1:
            answer = mmul(answer, base)
        base = mmul(base, base)
        n >>= 1
    return answer


def mtrace(a: Matrix) -> LElement:
    return a[0][0] + a[1][1]


def mdet(a: Matrix) -> LElement:
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def mscale(x: LElement, a: Matrix) -> Matrix:
    return matrix(*(x * a[i][j] for i in range(2) for j in range(2)))


def minverse(a: Matrix) -> Matrix:
    determinant_element = mdet(a)
    return matrix(
        a[1][1] / determinant_element,
        -a[0][1] / determinant_element,
        -a[1][0] / determinant_element,
        a[0][0] / determinant_element,
    )


def conjugate_by(h: Matrix, a: Matrix) -> Matrix:
    return mmul(mmul(h, a), minverse(h))


def all_entries_integral(a: Matrix) -> bool:
    return all(a[i][j].integral_power_basis() for i in range(2) for j in range(2))


def poly_derivative(a: list[Q]) -> list[Q]:
    return _trim([Q(i) * a[i] for i in range(1, len(a))])


def poly_divmod(a: list[Q], b: list[Q]) -> tuple[list[Q], list[Q]]:
    dividend = _trim(a[:])
    divisor = _trim(b[:])
    if divisor == [0]:
        raise ZeroDivisionError
    quotient = [Q(0)] * max(1, len(dividend) - len(divisor) + 1)
    while len(dividend) >= len(divisor) and dividend != [0]:
        shift = len(dividend) - len(divisor)
        factor = dividend[-1] / divisor[-1]
        quotient[shift] += factor
        for j, coefficient in enumerate(divisor):
            dividend[shift + j] -= factor * coefficient
        _trim(dividend)
    return _trim(quotient), _trim(dividend)


def poly_eval(a: list[Q] | tuple[Q, ...], x: Q) -> Q:
    value = Q(0)
    for coefficient in reversed(a):
        value = value * x + coefficient
    return value


def sturm_sequence(a: list[Q]) -> list[list[Q]]:
    sequence = [_trim(a[:]), poly_derivative(a)]
    while sequence[-1] != [0]:
        _, remainder = poly_divmod(sequence[-2], sequence[-1])
        if remainder == [0]:
            break
        sequence.append([-x for x in remainder])
    return sequence


def sign(value: Q) -> int:
    return (value > 0) - (value < 0)


def variations(sequence: list[list[Q]], x: Q) -> int:
    signs = [sign(poly_eval(p, x)) for p in sequence]
    signs = [s for s in signs if s]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def isolate_real_roots(bits: int = 72) -> list[tuple[Q, Q]]:
    sequence = sturm_sequence(list(MODULUS))

    def count(a: Q, b: Q) -> int:
        return variations(sequence, a) - variations(sequence, b)

    pending = [(Q(-2), Q(2), count(Q(-2), Q(2)))]
    roots: list[tuple[Q, Q]] = []
    target = Q(1, 1 << bits)
    while pending:
        a, b, n = pending.pop()
        if n == 0:
            continue
        if n == 1 and b - a <= target:
            roots.append((a, b))
            continue
        midpoint = (a + b) / 2
        left_count = count(a, midpoint)
        pending.append((midpoint, b, n - left_count))
        pending.append((a, midpoint, left_count))
    return sorted(roots)


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def point(value: Q | int) -> "Interval":
        return Interval(Q(value), Q(value))

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: "Interval") -> "Interval":
        return self + (-other)

    def __mul__(self, other: "Interval") -> "Interval":
        products = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(products), max(products))

    def reciprocal(self) -> "Interval":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("interval contains zero")
        return Interval(min(1 / self.lo, 1 / self.hi), max(1 / self.lo, 1 / self.hi))

    def __truediv__(self, other: "Interval") -> "Interval":
        return self * other.reciprocal()

    def strict_sign(self) -> int:
        if self.lo > 0:
            return 1
        if self.hi < 0:
            return -1
        raise ArithmeticError(f"interval does not certify a sign: {self}")

    def absolute(self) -> "Interval":
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return Interval(Q(0), max(-self.lo, self.hi))

    def midpoint(self) -> Q:
        return (self.lo + self.hi) / 2


def interval_eval(element: LElement, t_interval: Interval) -> Interval:
    value = Interval.point(0)
    for coefficient in reversed(element.coefficients):
        value = value * t_interval + Interval.point(coefficient)
    return value


def sqrt_interval(value: Interval, bits: int = 100) -> Interval:
    if value.lo <= 0:
        raise ArithmeticError("square-root interval is not certified positive")
    scale = 1 << bits

    def lower(q: Q) -> Q:
        floor_scaled_square = q.numerator * scale * scale // q.denominator
        return Q(isqrt(floor_scaled_square), scale)

    def upper(q: Q) -> Q:
        floor_scaled_square = q.numerator * scale * scale // q.denominator
        return Q(isqrt(floor_scaled_square) + 1, scale)

    return Interval(lower(value.lo), upper(value.hi))


def embedded_roots() -> dict[int, Interval]:
    ascending = isolate_real_roots()
    if len(ascending) != DEGREE:
        raise AssertionError(f"expected six real roots, found {len(ascending)}")
    # 2*cos(k*pi/21) decreases for 1 <= k < 21.
    ascending_labels = (19, 17, 13, 11, 5, 1)
    return {k: Interval(a, b) for k, (a, b) in zip(ascending_labels, ascending)}


def sign_vector(element: LElement, roots: dict[int, Interval]) -> tuple[int, ...]:
    return tuple(interval_eval(element, roots[k]).strict_sign() for k in EMBEDDINGS)


def fixed_point_box(a: Matrix, root: Interval) -> tuple[Interval, Interval]:
    aa = interval_eval(a[0][0], root)
    cc = interval_eval(a[1][0], root)
    dd = interval_eval(a[1][1], root)
    tr = aa + dd
    two = Interval.point(2)
    four = Interval.point(4)
    center = (aa - dd) / (two * cc)
    radicand = four - tr * tr
    height = sqrt_interval(radicand) / (two * cc.absolute())
    if height.lo <= 0:
        raise AssertionError("upper fixed-point height is not positive")
    return center, height


def fixed_point_invariants(a: Matrix) -> tuple[LElement, LElement]:
    """Return exact (real center, squared upper height) in L."""
    aa, cc, dd = a[0][0], a[1][0], a[1][1]
    center = (aa - dd) / (2 * cc)
    height_squared = (4 - mtrace(a) ** 2) / (4 * cc**2)
    return center, height_squared


def orientation_box(a: Matrix, b: Matrix, c: Matrix, root: Interval) -> Interval:
    ar, ai = fixed_point_box(a, root)
    br, bi = fixed_point_box(b, root)
    cr, ci = fixed_point_box(c, root)
    return (br - ar) * (ci - ai) - (bi - ai) * (cr - ar)


def power_trace(n: int) -> LElement:
    """Return zeta^n + zeta^-n as a polynomial in t=zeta+zeta^-1."""
    if n == 0:
        return LElement(2)
    a, b = LElement(2), T
    for _ in range(1, n):
        a, b = b, T * b - a
    return b


def sine_quotient(n: int) -> LElement:
    """Return (zeta^n-zeta^-n)/(zeta-zeta^-1)."""
    if n == 0:
        return ZERO
    a, b = ZERO, ONE
    for _ in range(1, n):
        a, b = b, T * b - a
    return b


def alpha_at_matrix(a: Matrix, k: int) -> Matrix:
    return madd(madd(mpow(a, k), mpow(a, k + 1)), mpow(a, k + 4))


def companion(trace_element: LElement, determinant_element: LElement) -> Matrix:
    return matrix(0, -determinant_element, 1, trace_element)


def column_intertwiner(a: Matrix) -> Matrix:
    """P=(e1,a*e1), for comparison with the companion matrix."""
    return matrix(1, a[0][0], 0, a[1][0])


def field_discriminant() -> Q:
    gram: list[list[Q]] = []
    for i in range(DEGREE):
        gram.append([(T ** (i + j)).trace() for j in range(DEGREE)])
    return determinant(gram)


def fmt_float(value: Q) -> str:
    # Human-readable output only.  No assertion depends on this conversion.
    return f"{float(value): .12f}"


def main() -> None:
    roots = embedded_roots()

    # Section 4 of McMullen, with t=2*cos(pi/21).
    A = matrix(
        elt(0, 3, -1, -4, 0, 1),
        elt(0, -10, -3, 6, 1, -1),
        elt(-1, 0, 0, 1),
        elt(0, -6, 1, 5, 0, -1),
    )
    B = matrix(0, -1, 1, T)
    C = matrix(
        elt(0, 2, -4, -4, 1, 1),
        elt(0, -5, -4, 1, 1),
        elt(-1, 5, 4, -1, -1),
        elt(0, -1, 4, 4, -1, -1),
    )
    D = matrix(
        elt(-1, -13, -6, 7, 2, -1),
        elt(0, 6, -5, -9, 1, 2),
        elt(1, -5, 2, 5, 0, -1),
        elt(1, 14, 6, -7, -2, 1),
    )
    vertices = {"A": A, "B": B, "C": C, "D": D}

    assert field_discriminant() == 453789
    assert T.norm() == 1
    assert all(mdet(v) == ONE for v in vertices.values())
    assert mtrace(A) == T**3 - 3 * T
    assert all(mtrace(vertices[name]) == T for name in ("B", "C", "D"))
    assert mpow(A, 7) == NEG_IDENTITY
    assert all(mpow(vertices[name], 21) == NEG_IDENTITY for name in ("B", "C", "D"))
    assert mmul(mmul(mmul(A, B), C), D) == IDENTITY

    h = matrix(T, 0, 0, 1)
    normalized = {name: conjugate_by(h, value) for name, value in vertices.items()}
    assert all(all_entries_integral(value) for value in normalized.values())
    assert all(mdet(value) == ONE for value in normalized.values())

    raw_orientation = tuple(
        orientation_box(A, B, C, roots[k]).strict_sign() for k in EMBEDDINGS
    )
    relative_orientation = tuple(x * raw_orientation[0] for x in raw_orientation)
    assert raw_orientation == (-1, -1, 1, 1, 1, 1)
    assert relative_orientation == (1, 1, -1, -1, -1, -1)
    assert sign_vector(T, roots) == relative_orientation

    # The B--D rational intertwiner printed conceptually in the manuscript.
    p_bd = column_intertwiner(normalized["D"])
    assert mmul(p_bd, B) == mmul(normalized["D"], p_bd)
    q_bd = mmul(p_bd, minverse(h))
    assert mmul(q_bd, normalized["B"]) == mmul(normalized["D"], q_bd)
    det_p_bd = mdet(p_bd)
    det_q_bd = mdet(q_bd)
    assert abs(det_p_bd.norm()) == 13693
    assert sign_vector(det_p_bd, roots) == relative_orientation
    assert sign_vector(det_q_bd, roots) == (1, 1, 1, 1, 1, 1)

    # Corrected ell=43 packet, reconstructed directly from alpha_k.
    pi = T**4 + T**3 - 4 * T**2 - 2 * T + 5
    k = 17
    c17 = power_trace(k) + power_trace(k + 1) + power_trace(k + 4)
    delta17 = c17**2 - 4 * pi
    assert pi.norm() == 43
    assert sign_vector(pi, roots) == (1, 1, 1, 1, 1, 1)
    assert delta17.norm() == 1
    assert sign_vector(delta17, roots) == (-1, -1, -1, -1, -1, -1)

    # alpha_17=x17+y17*zeta.  The sign of y17 selects alpha versus conjugate
    # on the upper companion eigenline because sin(r*pi/21)>0 for our r's.
    y17 = sine_quotient(k) + sine_quotient(k + 1) + sine_quotient(k + 4)
    x17 = (c17 - y17 * T) / 2
    assert x17.integral_power_basis() and y17.integral_power_basis()
    assert y17.norm() == 1
    assert delta17 == y17**2 * (T**2 - 4)
    assert (T**2 - 4).norm() == 1
    selected_signs = sign_vector(y17, roots)
    phi17 = tuple(r if s > 0 else (-r) % 42 for r, s in zip(EMBEDDINGS, selected_signs))
    galois_group = tuple(a for a in range(42) if gcd(a, 42) == 1)
    phi17_stabilizer = tuple(
        a for a in galois_group if {(a * exponent) % 42 for exponent in phi17} == set(phi17)
    )
    assert phi17_stabilizer == (1,)
    k_signature = (
        sum(exponent % 3 == 1 for exponent in phi17),
        sum(exponent % 3 == 2 for exponent in phi17),
    )
    assert k_signature == (3, 3)

    # sqrt(-3)=1+2*zeta^14 = X17+Y17*alpha17.
    q13 = sine_quotient(13)
    q14 = sine_quotient(14)
    capital_y17 = (2 * q14) / y17
    capital_x17 = 1 - 2 * q13 - capital_y17 * x17
    assert capital_x17 + capital_y17 * x17 == 1 - 2 * q13
    assert capital_y17 * y17 == 2 * q14
    assert capital_y17.norm() == -1728
    hermitian_determinant = LElement(field_discriminant()) * LElement(capital_y17.norm())
    assert hermitian_determinant == LElement(-(2**6) * (3**6) * (7**5))

    companion_center = -c17 / 2
    companion_height_squared = pi - c17**2 / 4
    assert companion_height_squared == (T + 2) / 4
    normalized_fixed_invariants = {
        name: fixed_point_invariants(value) for name, value in normalized.items()
    }
    for center, height_squared in normalized_fixed_invariants.values():
        assert sign_vector(height_squared, roots) == (1, 1, 1, 1, 1, 1)
        assert (center, height_squared) != (companion_center, companion_height_squared)

    companion17 = companion(c17, pi)
    vertex_fingerprints: dict[str, tuple[tuple[int, ...], tuple[int, ...], int]] = {}
    for name in ("B", "C", "D"):
        value = normalized[name]
        t17 = alpha_at_matrix(value, k)
        assert mtrace(t17) == c17
        assert mdet(t17) == pi
        p0 = column_intertwiner(t17)
        assert mmul(p0, companion17) == mmul(t17, p0)
        direct_signs = sign_vector(mdet(p0), roots)

        conjugate_t17 = msub(mscale(c17, IDENTITY), t17)
        conjugate_p0 = column_intertwiner(conjugate_t17)
        assert mmul(conjugate_p0, companion17) == mmul(conjugate_t17, conjugate_p0)
        conjugate_signs = sign_vector(mdet(conjugate_p0), roots)
        determinant_norm = abs(mdet(p0).norm())
        vertex_fingerprints[name] = (direct_signs, conjugate_signs, int(determinant_norm))

    # Expected-output comparison is intentionally last: all preceding values
    # were reconstructed from definitions, not used to guide the computation.
    expected_fingerprints = {
        "B": ((1, 1, 1, -1, -1, 1), (-1, -1, -1, 1, 1, -1), 1),
        "C": ((1, 1, -1, 1, 1, -1), (-1, -1, 1, -1, -1, 1), 127),
        "D": ((1, 1, 1, -1, -1, 1), (-1, -1, -1, 1, 1, -1), 13693),
    }
    assert vertex_fingerprints == expected_fingerprints

    print("INDEPENDENT QUADRILATERAL + k=17 STANDARD-LIBRARY CHECK")
    print(f"disc(Z[t])={field_discriminant()} ; Norm(t)={T.norm()}")
    print("raw Omega(A,B,C) signs:", raw_orientation)
    print("relative component signs:", relative_orientation)
    print("B--D rational orientation determinant norm:", abs(det_p_bd.norm()))
    print("B--D normalized rational intertwiner signs:", sign_vector(det_q_bd, roots))
    print("c_17 =", c17)
    print("Delta_17 =", delta17)
    print("Norm(pi)=", pi.norm(), "; Norm(Delta_17)=", delta17.norm())
    print("Phi_17 =", phi17, "; stabilizer=", phi17_stabilizer, "; K-signature=", k_signature)
    print("Norm(Y_17 for sqrt(-3)=X_17+Y_17*alpha_17) =", capital_y17.norm())
    print("Hermitian determinant integer =", hermitian_determinant.coefficients[0])
    print("exact fixed-point invariants (center, height^2):")
    print("  k=17 companion:", companion_center, ";", companion_height_squared)
    for name in ("A", "B", "C", "D"):
        center, height_squared = normalized_fixed_invariants[name]
        print(f"  normalized {name}:", center, ";", height_squared)
    print("k=17 vertex orientation fingerprints:")
    for name, fingerprint in vertex_fingerprints.items():
        print(f"  {name}: direct={fingerprint[0]} conjugate={fingerprint[1]} |Norm(det P0)|={fingerprint[2]}")
    print("certified companion fixed-point boxes (displayed as decimal midpoints only):")
    for embedding in EMBEDDINGS:
        root = roots[embedding]
        c_box = interval_eval(c17, root)
        pi_box = interval_eval(pi, root)
        center = -c_box / Interval.point(2)
        height = sqrt_interval(Interval.point(4) * pi_box - c_box * c_box) / Interval.point(2)
        print(
            f"  sigma_{embedding:02d}: t={fmt_float(root.midpoint())} "
            f"w17={fmt_float(center.midpoint())}+{fmt_float(height.midpoint())}i"
        )
    print("PASS: exact identities and every asserted sign used rational arithmetic/intervals")
    print("SCOPE: this is not a full Hilbert-modular orbit-membership decision for k=17")


if __name__ == "__main__":
    main()
