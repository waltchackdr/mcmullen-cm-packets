#!/usr/bin/env python3
"""Independent exact A-vertex realization for the residual pair {20,41}.

This verifier imports only the standard-library quotient-field and Sturm
primitives from ``quadrilateral_k17_stdlib_v41.py``.  It does not import or
execute a certificate or experiment.  All matrix identities and integrality
claims are exact; every real sign is certified by rational Sturm intervals.

The script first verifies explicit totally-positive-unit determinant
intertwiners, then reconstructs a specific cyclotomic centralizer unit whose
relative norm cancels that determinant.  The resulting matrices are checked
directly in SL_2(O_L); no general norm-surjectivity theorem is invoked.
"""

from __future__ import annotations

from quadrilateral_k17_stdlib_v41 import (
    EMBEDDINGS,
    IDENTITY,
    NEG_IDENTITY,
    ONE,
    T,
    ZERO,
    LElement,
    Matrix,
    all_entries_integral,
    companion,
    conjugate_by,
    elt,
    embedded_roots,
    field_discriminant,
    fixed_point_invariants,
    madd,
    mdet,
    minverse,
    mmul,
    mpow,
    mscale,
    msub,
    mtrace,
    matrix,
    power_trace,
    sign_vector,
)


Pair = tuple[LElement, LElement]
POSITIVE_SIGNS = (1,) * len(EMBEDDINGS)


def pair_multiply(left: Pair, right: Pair) -> Pair:
    """Multiply A+B*z pairs using z^2-t*z+1=0."""
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c + T * b * d


def z_power(exponent: int) -> Pair:
    answer = (ONE, ZERO)
    base = (ZERO, ONE)
    for _ in range(exponent):
        answer = pair_multiply(answer, base)
    return answer


def alpha_pair(index: int) -> Pair:
    """Return z^k(1+z+z^4) as an exact A+B*z pair."""
    z4 = z_power(4)
    seed = (ONE + z4[0], ONE + z4[1])
    return pair_multiply(z_power(index), seed)


def pair_conjugate(value: Pair) -> Pair:
    """Complex conjugation z -> t-z on A+B*z."""
    a, b = value
    return a + T * b, -b


def pair_trace(value: Pair) -> LElement:
    a, b = value
    return 2 * a + T * b


def pair_norm(value: Pair) -> LElement:
    a, b = value
    return a * a + T * a * b + b * b


def evaluate_pair(value: Pair, generator: Matrix) -> Matrix:
    a, b = value
    return madd(mscale(a, IDENTITY), mscale(b, generator))


def upper_generator_type(generator: Matrix, roots) -> tuple[int, ...]:
    """Recover the CM type of a zeta generator at the upper fixed point.

    The generator has trace t and determinant one.  At every printed
    embedding r in (1,5,11,13,17,19), sin(r*pi/21)>0.  Hence the sign of the
    lower-left entry says whether its upper eigenvalue is zeta^r or zeta^-r.
    The signs are exact Sturm decisions.
    """
    assert mtrace(generator) == T
    assert mdet(generator) == ONE
    lower_signs = sign_vector(generator[1][0], roots)
    return tuple(
        r if sign > 0 else (-r) % 42
        for r, sign in zip(EMBEDDINGS, lower_signs)
    )


def main() -> None:
    roots = embedded_roots()
    assert field_discriminant() == 453789

    # McMullen's printed order-14 A generator and its h=diag(t,1)
    # normalization.  The expected matrix is asserted only after reconstructing
    # it by exact conjugation.
    a_printed = matrix(
        elt(0, 3, -1, -4, 0, 1),
        elt(0, -10, -3, 6, 1, -1),
        elt(-1, 0, 0, 1),
        elt(0, -6, 1, 5, 0, -1),
    )
    assert mdet(a_printed) == ONE
    assert mtrace(a_printed) == T**3 - 3 * T
    assert mpow(a_printed, 7) == NEG_IDENTITY

    h = matrix(T, 0, 0, 1)
    a_t = conjugate_by(h, a_printed)
    expected_a_t = matrix(
        elt(0, 3, -1, -4, 0, 1),
        elt(1, 8, -2, -9, 0, 2),
        elt(8, 8, -5, -6, 1, 1),
        elt(0, -6, 1, 5, 0, -1),
    )
    assert a_t == expected_a_t
    assert all_entries_integral(a_t)
    assert mdet(a_t) == ONE
    assert mpow(a_t, 7) == NEG_IDENTITY

    # In O_L[z]/(z^2-tz+1), z^3=-t+(t^2-1)z.  Solving
    # z=x+y*z^3 gives the following exact matrix cube root of A_t.
    z3 = z_power(3)
    assert z3 == (-T, T**2 - 1)
    coefficient_y = ONE / z3[1]
    coefficient_x = -z3[0] * coefficient_y
    z_a = madd(mscale(coefficient_x, IDENTITY), mscale(coefficient_y, a_t))
    expected_z_a = matrix(
        elt(-1, -11, 0, 7, 0, -1),
        elt(-2, -15, 0, 8, 0, -1),
        elt(-7, 0, 6, 0, -1),
        elt(1, 12, 0, -7, 0, 1),
    )
    assert z_a == expected_z_a
    assert all_entries_integral(z_a)
    assert msub(madd(mpow(z_a, 2), IDENTITY), mscale(T, z_a)) == matrix(0, 0, 0, 0)
    assert mpow(z_a, 3) == a_t
    assert mtrace(z_a) == T
    assert mdet(z_a) == ONE
    assert mpow(z_a, 21) == NEG_IDENTITY
    assert fixed_point_invariants(z_a) == fixed_point_invariants(a_t)

    z_a_inverse = msub(mscale(T, IDENTITY), z_a)
    assert z_a_inverse == minverse(z_a)
    assert mmul(z_a, z_a_inverse) == IDENTITY
    assert all_entries_integral(z_a_inverse)
    assert fixed_point_invariants(z_a_inverse) == fixed_point_invariants(z_a)

    # The order discriminant is t^2-4.  Its absolute norm is one, so its
    # discriminant ideal is the unit ideal: the integral quadratic order
    # O_L[z_a] admits no proper integral overorder and is the maximal O_M.
    relative_discriminant = T**2 - 4
    assert relative_discriminant.integral_power_basis()
    assert relative_discriminant.norm() == 1
    assert sign_vector(relative_discriminant, roots) == (-1,) * len(EMBEDDINGS)

    # Reconstruct alpha_20 and alpha_41 from z^k(1+z+z^4), rather than taking
    # their coefficient pairs or traces as inputs.
    alpha20 = alpha_pair(20)
    alpha41 = alpha_pair(41)
    assert alpha20 == (-ONE, 2 - T**2)
    assert alpha41 == (ONE, T**2 - 2)
    assert alpha41 == (-alpha20[0], -alpha20[1])

    pi = T**4 + T**3 - 4 * T**2 - 2 * T + 5
    c20 = -T**3 + 2 * T - 2
    c41 = -c20
    assert pair_trace(alpha20) == c20
    assert pair_trace(alpha41) == c41
    assert c20 == power_trace(20) + power_trace(21) + power_trace(24)
    assert c41 == power_trace(41) + power_trace(42) + power_trace(45)
    assert pair_norm(alpha20) == pi
    assert pair_norm(alpha41) == pi
    assert alpha20[1].norm() == 1 and alpha41[1].norm() == 1
    delta20 = c20**2 - 4 * pi
    delta41 = c41**2 - 4 * pi
    assert delta20 == alpha20[1] ** 2 * relative_discriminant
    assert delta41 == alpha41[1] ** 2 * relative_discriminant
    assert delta20.norm() == 1 and delta41.norm() == 1
    assert pi.norm() == 43
    assert sign_vector(pi, roots) == POSITIVE_SIGNS

    # The assignment is deliberately asymmetric.  k=41 uses z->Z_A, while
    # k=20 uses z->Z_A^{-1}; looking only at one Hecke matrix's lower-left
    # signs would not certify this full M-action or its CM type.
    t41 = evaluate_pair(alpha41, z_a)
    t20 = evaluate_pair(alpha20, z_a_inverse)
    assert t20 == evaluate_pair(pair_conjugate(alpha20), z_a)
    expected_t20 = matrix(
        elt(0, 17, -1, -12, 0, 2),
        elt(3, 23, -2, -17, 0, 3),
        elt(15, 8, -11, -6, 2, 1),
        elt(-2, -15, 1, 11, 0, -2),
    )
    expected_t41 = matrix(
        elt(2, 15, -1, -11, 0, 2),
        elt(3, 23, -2, -17, 0, 3),
        elt(15, 8, -11, -6, 2, 1),
        elt(0, -17, 1, 12, 0, -2),
    )
    assert t20 == expected_t20
    assert t41 == expected_t41
    assert all_entries_integral(t20) and all_entries_integral(t41)
    assert (mtrace(t20), mdet(t20)) == (c20, pi)
    assert (mtrace(t41), mdet(t41)) == (c41, pi)

    companion20 = companion(c20, pi)
    companion41 = companion(c41, pi)

    # Derive Q=(e_2,T e_2) from each T, then compare it with a fully explicit
    # integral matrix.  Cayley-Hamilton gives Q A_c=T Q; the verifier checks
    # the identity entry by entry rather than invoking that argument blindly.
    q20 = matrix(0, t20[0][1], 1, t20[1][1])
    q41 = matrix(0, t41[0][1], 1, t41[1][1])
    expected_q20 = matrix(
        0,
        elt(3, 23, -2, -17, 0, 3),
        1,
        elt(-2, -15, 1, 11, 0, -2),
    )
    expected_q41 = matrix(
        0,
        elt(3, 23, -2, -17, 0, 3),
        1,
        elt(0, -17, 1, 12, 0, -2),
    )
    assert q20 == expected_q20
    assert q41 == expected_q41
    assert mmul(q20, companion20) == mmul(t20, q20)
    assert mmul(q41, companion41) == mmul(t41, q41)

    common_determinant = elt(-3, -23, 2, 17, 0, -3)
    assert mdet(q20) == common_determinant
    assert mdet(q41) == common_determinant
    assert common_determinant.norm() == 1
    assert common_determinant.inverse().integral_power_basis()
    assert sign_vector(common_determinant, roots) == POSITIVE_SIGNS
    assert all_entries_integral(q20) and all_entries_integral(q41)
    assert all_entries_integral(minverse(q20)) and all_entries_integral(minverse(q41))

    # Instead of appealing to a general norm-surjectivity theorem, reconstruct
    # the concrete centralizer unit
    #   epsilon=(1-zeta_21^11)(1-zeta_21^13)
    #          =(1-z^22)(1-z^26).
    z22 = z_power(22)
    z26 = z_power(26)
    epsilon = pair_multiply(
        (ONE - z22[0], -z22[1]),
        (ONE - z26[0], -z26[1]),
    )
    assert epsilon == (
        elt(0, 2, 3, -1, -1, 0),
        elt(2, 3, -3, -4, 1, 1),
    )
    epsilon_norm = pair_norm(epsilon)
    assert epsilon_norm == common_determinant.inverse()
    assert epsilon_norm == elt(3, 4, -3, -4, 1, 1)
    assert epsilon_norm.norm() == 1
    assert sign_vector(epsilon_norm, roots) == POSITIVE_SIGNS

    centralizers: dict[int, Matrix] = {}
    sl_intertwiners: dict[int, Matrix] = {}
    expected_x = elt(1, 4, 6, 2, -2, -1)
    expected_y20 = elt(1, 2, 3, 3, -1, -1)
    expected_y41 = -expected_y20
    data = (
        (20, alpha20, companion20, t20, q20, expected_y20),
        (41, alpha41, companion41, t41, q41, expected_y41),
    )
    for index, alpha, companion_matrix, target, q_matrix, expected_y in data:
        # Since the z-coefficient b_k has norm one, it is an integral unit.
        # Thus epsilon=x_k+y_k*alpha_k has integral x_k,y_k.
        y_coefficient = epsilon[1] / alpha[1]
        x_coefficient = epsilon[0] - y_coefficient * alpha[0]
        assert x_coefficient == expected_x
        assert y_coefficient == expected_y
        assert x_coefficient.integral_power_basis()
        assert y_coefficient.integral_power_basis()
        assert (
            x_coefficient + y_coefficient * alpha[0],
            y_coefficient * alpha[1],
        ) == epsilon

        c_epsilon = madd(
            mscale(x_coefficient, IDENTITY),
            mscale(y_coefficient, companion_matrix),
        )
        assert all_entries_integral(c_epsilon)
        assert mmul(c_epsilon, companion_matrix) == mmul(companion_matrix, c_epsilon)
        assert mdet(c_epsilon) == epsilon_norm
        centralizers[index] = c_epsilon

        g_matrix = mmul(q_matrix, c_epsilon)
        assert all_entries_integral(g_matrix)
        assert mdet(g_matrix) == ONE
        assert mmul(g_matrix, companion_matrix) == mmul(target, g_matrix)
        assert all_entries_integral(minverse(g_matrix))
        sl_intertwiners[index] = g_matrix

    expected_c20 = matrix(
        elt(1, 4, 6, 2, -2, -1),
        elt(-6, -15, -8, 3, 4, 1),
        elt(1, 2, 3, 3, -1, -1),
        elt(-2, -6, -5, -1, 2, 1),
    )
    expected_c41 = matrix(
        elt(1, 4, 6, 2, -2, -1),
        elt(6, 15, 8, -3, -4, -1),
        elt(-1, -2, -3, -3, 1, 1),
        elt(-2, -6, -5, -1, 2, 1),
    )
    expected_g20 = matrix(
        elt(-1, -5, -3, 1, 1, 0),
        elt(2, 10, 0, -3, 0, 0),
        elt(1, 2, 0, 0, 0, 0),
        elt(-6, -11, 7, 6, -2, -1),
    )
    expected_g41 = matrix(
        elt(1, 5, 3, -1, -1, 0),
        elt(2, 10, 0, -3, 0, 0),
        elt(-2, -4, 1, 1, 0, 0),
        elt(0, -6, 0, 2, 0, 0),
    )
    assert centralizers == {20: expected_c20, 41: expected_c41}
    assert sl_intertwiners == {20: expected_g20, 41: expected_g41}

    # Certify the *generator* embeddings.  This avoids the invalid shortcut of
    # inferring the full M CM type from signs of the single Hecke matrix T_k.
    phi_z_a = upper_generator_type(z_a, roots)
    phi_z_a_inverse = upper_generator_type(z_a_inverse, roots)
    expected_phi41 = (1, 5, 31, 29, 17, 19)
    expected_phi20 = (41, 37, 11, 13, 25, 23)
    assert phi_z_a == expected_phi41
    assert phi_z_a_inverse == expected_phi20
    assert set(phi_z_a_inverse) == {(-exponent) % 42 for exponent in phi_z_a}
    assert set(phi_z_a).isdisjoint(phi_z_a_inverse)
    assert set(phi_z_a) | set(phi_z_a_inverse) == {
        1, 5, 11, 13, 17, 19, 23, 25, 29, 31, 37, 41
    }

    print("INDEPENDENT A-VERTEX k={20,41} STURM CHECK")
    print("Z_A integral; Z_A^2-t Z_A+I=0; Z_A^3=A_t; det(Z_A)=1")
    print("Norm(t^2-4)=1: O_L[Z_A] has maximal multiplier order O_M")
    print("c20 =", c20)
    print("c41 =", c41)
    print("Norm(pi)=", pi.norm(), "; both Hecke determinants=pi")
    print("det(Q20)=det(Q41)=", common_determinant)
    print("Norm(det Q)=", common_determinant.norm())
    print("signs(det Q)=", sign_vector(common_determinant, roots))
    print("Phi20 via Z_A^{-1} =", phi_z_a_inverse)
    print("Phi41 via Z_A      =", phi_z_a)
    print("Norm_{M/L}(epsilon)=", epsilon_norm, "=det(Q)^-1")
    print("det(G20)=det(G41)=1; both G_k are integral SL_2 intertwiners")
    print("PASS: exact A-vertex SL_2(O_L) intertwiners for k=20 and k=41")


if __name__ == "__main__":
    main()
