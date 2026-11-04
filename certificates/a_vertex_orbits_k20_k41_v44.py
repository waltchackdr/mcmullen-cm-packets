#!/usr/bin/env python3
"""Exact A-vertex orbit certificate for maximal-order indices k={20,41}.

The certificate exhibits integral intertwiners from the two companion
matrices to holomorphically correct degree-43 endomorphisms of the standard
normalized A vertex.  Their common determinant is a totally positive unit.
An explicit cyclotomic unit gives the integral centralizer adjustment to
determinant one; no existential unit-norm step remains.

All polynomial reductions, resultants, root isolations, and signs are exact.
No binary floating-point value participates in an assertion.
"""

from __future__ import annotations

import sympy as sp


t = sp.symbols("t")
f = sp.Poly(
    t**6 + t**5 - 6*t**4 - 6*t**3 + 8*t**2 + 8*t + 1,
    t,
    domain=sp.QQ,
)
I = sp.eye(2)
BASE_EMBEDDINGS = (1, 5, 11, 13, 17, 19)


def polynomial(value: sp.Expr) -> sp.Poly:
    return sp.Poly(sp.expand(value), t, domain=sp.QQ)


def red(value: sp.Expr) -> sp.Expr:
    return polynomial(value).rem(f).as_expr()


def inverse(value: sp.Expr) -> sp.Expr:
    return sp.invert(polynomial(value), f).as_expr()


def matrix_red(value: sp.Matrix) -> sp.Matrix:
    return value.applyfunc(red)


def matrix_mul(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    return matrix_red(left * right)


def matrix_power(value: sp.Matrix, exponent: int) -> sp.Matrix:
    result = I
    base = value
    n = exponent
    while n:
        if n & 1:
            result = matrix_mul(result, base)
        base = matrix_mul(base, base)
        n >>= 1
    return result


def assert_zero_matrix(value: sp.Matrix) -> None:
    assert all(red(entry) == 0 for entry in value)


def assert_integral_matrix(value: sp.Matrix) -> None:
    for entry in value:
        sp.Poly(red(entry), t, domain=sp.ZZ)


def norm(value: sp.Expr) -> sp.Rational:
    return sp.Rational(sp.resultant(f.as_expr(), polynomial(value).as_expr(), t))


ROOT_INTERVALS = [
    interval
    for interval, _ in sp.polys.polytools.intervals(
        f, eps=sp.Rational(1, 10)**40
    )
]
ROOT_INTERVALS.sort(key=lambda interval: interval[0], reverse=True)
assert len(ROOT_INTERVALS) == 6


def interval_eval(value: sp.Expr, interval) -> tuple[sp.Rational, sp.Rational]:
    p = polynomial(value)
    lo, hi = map(sp.Rational, interval)
    lower = upper = p.all_coeffs()[0]
    for coefficient in p.all_coeffs()[1:]:
        products = (lower * lo, lower * hi, upper * lo, upper * hi)
        lower = min(products) + coefficient
        upper = max(products) + coefficient
    return lower, upper


def signs(value: sp.Expr) -> tuple[int, ...]:
    result = []
    for interval in ROOT_INTERVALS:
        lower, upper = interval_eval(red(value), interval)
        assert lower > 0 or upper < 0
        result.append(1 if lower > 0 else -1)
    return tuple(result)


def generator_cm_type(value: sp.Matrix) -> tuple[int, ...]:
    """CM type of the full M-action defined by zeta_42 mapping to value.

    This predicate is valid for the cyclotomic generator, whose trace is t.
    It must not be applied to a single Hecke matrix alpha_k(value): that
    matrix's lower-left signs do not by themselves identify the full M-action.
    """
    orientation = signs(value[1, 0])
    return tuple(
        embedding if sign > 0 else 42 - embedding
        for embedding, sign in zip(BASE_EMBEDDINGS, orientation)
    )


def alpha_k(value: sp.Matrix, k: int) -> sp.Matrix:
    alpha_zero = matrix_red(I + value + matrix_power(value, 4))
    return matrix_mul(matrix_power(value, k), alpha_zero)


def pair_mul(
    left: tuple[sp.Expr, sp.Expr], right: tuple[sp.Expr, sp.Expr]
) -> tuple[sp.Expr, sp.Expr]:
    """Multiply a+bz pairs using z^2=tz-1."""
    a, b = left
    c, d = right
    return red(a*c - b*d), red(a*d + b*c + t*b*d)


def pair_power(exponent: int) -> tuple[sp.Expr, sp.Expr]:
    result = (sp.Integer(1), sp.Integer(0))
    base = (sp.Integer(0), sp.Integer(1))
    n = exponent
    while n:
        if n & 1:
            result = pair_mul(result, base)
        base = pair_mul(base, base)
        n >>= 1
    return result


def pair_sub(
    left: tuple[sp.Expr, sp.Expr], right: tuple[sp.Expr, sp.Expr]
) -> tuple[sp.Expr, sp.Expr]:
    return red(left[0] - right[0]), red(left[1] - right[1])


def alpha_pair(k: int) -> tuple[sp.Expr, sp.Expr]:
    alpha_zero = (red(2 - t**2), red(t**3 - 2*t + 1))
    return pair_mul(pair_power(k), alpha_zero)


def pair_norm(value: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    a, b = value
    return red(a*a + t*a*b + b*b)


def main() -> None:
    print("EXACT A-VERTEX ORBIT CERTIFICATE v44 FOR k={20,41}")

    a_printed = sp.Matrix(
        [
            [t**5 - 4*t**3 - t**2 + 3*t,
             -t**5 + t**4 + 6*t**3 - 3*t**2 - 10*t],
            [t**3 - 1,
             -t**5 + 5*t**3 + t**2 - 6*t],
        ]
    )
    t_inverse = inverse(t)
    assert abs(norm(t)) == 1
    sp.Poly(t_inverse, t, domain=sp.ZZ)
    h = sp.diag(t, 1)
    a_normalized = matrix_mul(matrix_mul(h, a_printed), sp.diag(t_inverse, 1))
    assert_integral_matrix(a_normalized)

    z_a = sp.Matrix(
        [
            [-t**5 + 7*t**3 - 11*t - 1,
             -t**5 + 8*t**3 - 15*t - 2],
            [-t**4 + 6*t**2 - 7,
             t**5 - 7*t**3 + 12*t + 1],
        ]
    )
    assert_integral_matrix(z_a)
    assert_zero_matrix(matrix_power(z_a, 2) - t*z_a + I)
    assert_zero_matrix(matrix_power(z_a, 3) - a_normalized)
    assert_zero_matrix(matrix_power(z_a, 21) + I)
    assert_zero_matrix(matrix_power(z_a, 42) - I)
    assert red(z_a.trace() - t) == 0
    assert red(z_a.det() - 1) == 0
    assert abs(norm(t**2 - 4)) == 1
    assert signs(4 - t**2) == (1, 1, 1, 1, 1, 1)
    phi_a = generator_cm_type(z_a)
    assert phi_a == (1, 5, 31, 29, 17, 19)
    z_a_inverse = matrix_red(t*I - z_a)
    assert_zero_matrix(matrix_mul(z_a, z_a_inverse) - I)
    phi_a_conjugate = generator_cm_type(z_a_inverse)
    assert phi_a_conjugate == (41, 37, 11, 13, 25, 23)

    pi = red(t**4 + t**3 - 4*t**2 - 2*t + 5)
    assert abs(norm(pi)) == 43
    c20 = red(-t**3 + 2*t - 2)
    c41 = red(t**3 - 2*t + 2)
    assert c41 == red(-c20)

    raw20 = alpha_k(z_a, 20)
    t20 = matrix_red(c20*I - raw20)
    t41 = alpha_k(z_a, 41)
    assert_zero_matrix(t20 - alpha_k(z_a_inverse, 20))
    cases = ((20, c20, t20), (41, c41, t41))
    for k, trace, hecke in cases:
        assert_integral_matrix(hecke)
        assert red(hecke.trace() - trace) == 0
        assert red(hecke.det() - pi) == 0
        assert abs(norm(trace**2 - 4*pi)) == 1
        assert signs(trace**2 - 4*pi) == (-1, -1, -1, -1, -1, -1)
        assert hecke[1, 0] != 0

    common_upper_right = red(3*t**5 - 17*t**3 - 2*t**2 + 23*t + 3)
    q20 = sp.Matrix(
        [
            [0, common_upper_right],
            [1, -2*t**5 + 11*t**3 + t**2 - 15*t - 2],
        ]
    )
    q41 = sp.Matrix(
        [
            [0, common_upper_right],
            [1, -2*t**5 + 12*t**3 + t**2 - 17*t],
        ]
    )
    determinant = red(-3*t**5 + 17*t**3 + 2*t**2 - 23*t - 3)
    assert_integral_matrix(q20)
    assert_integral_matrix(q41)
    assert red(q20.det() - determinant) == 0
    assert red(q41.det() - determinant) == 0
    assert norm(determinant) == 1
    assert signs(determinant) == (1, 1, 1, 1, 1, 1)
    determinant_inverse = inverse(determinant)
    sp.Poly(determinant_inverse, t, domain=sp.ZZ)
    assert signs(determinant_inverse) == (1, 1, 1, 1, 1, 1)

    one_pair = (sp.Integer(1), sp.Integer(0))
    epsilon = pair_mul(
        pair_sub(one_pair, pair_power(22)),
        pair_sub(one_pair, pair_power(26)),
    )
    assert pair_norm(epsilon) == red(determinant_inverse)

    for k, trace, hecke, intertwiner in (
        (20, c20, t20, q20),
        (41, c41, t41, q41),
    ):
        companion = sp.Matrix([[0, -pi], [1, trace]])
        assert_integral_matrix(companion)
        assert_zero_matrix(intertwiner*companion - hecke*intertwiner)
        alpha_a, alpha_b = alpha_pair(k)
        assert abs(norm(alpha_b)) == 1
        centralizer_y = red(epsilon[1] * inverse(alpha_b))
        centralizer_x = red(epsilon[0] - centralizer_y*alpha_a)
        sp.Poly(centralizer_x, t, domain=sp.ZZ)
        sp.Poly(centralizer_y, t, domain=sp.ZZ)
        centralizer = matrix_red(centralizer_x*I + centralizer_y*companion)
        assert_integral_matrix(centralizer)
        assert_zero_matrix(centralizer*companion - companion*centralizer)
        assert red(centralizer.det() - determinant_inverse) == 0
        special = matrix_mul(intertwiner, centralizer)
        assert_integral_matrix(special)
        assert red(special.det() - 1) == 0
        assert_zero_matrix(special*companion - hecke*special)
        print(
            f"k={k}: {'conjugate' if k == 20 else 'direct'} A realization; "
            f"Q A_c = T Q; explicit G in SL2(O_L); "
            f"full O_M CM type="
            f"{phi_a_conjugate if k == 20 else phi_a}"
        )
        print(f"  epsilon = x + y alpha_{k}: x={centralizer_x}; y={centralizer_y}")

    print("A cube root generates the maximal order: |Norm(t^2-4)|=1")
    print("explicit epsilon=(1-zeta_21^11)(1-zeta_21^13) has norm det(Q)^-1")
    print("PASS: k=20 and k=41 lie in the normalized A SL2(O_L)-orbit")


if __name__ == "__main__":
    main()
