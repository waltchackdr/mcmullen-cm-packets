#!/usr/bin/env python3
"""New exact reconstruction of the manuscript's k=11 vertex bridge.

This source is newly written, not a recovered legacy certificate. The primary
path starts in Q[z]/Phi_42, then extracts coordinates in L+Lz. Geometry of
McMullen's vertex remains an explicitly separate human-proof dependency.
Registration and claim-ledger promotion require the root agent's review.
"""

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O/-OO")
    source = Path(__file__)
    before = sha256(source.read_bytes()).hexdigest()
    z, t = sp.symbols("z t")
    cyclotomic = sp.Poly(sp.cyclotomic_poly(42, z), z, domain=sp.QQ)
    minimal = sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1, t, domain=sp.QQ)
    assert cyclotomic.degree() == 12 and cyclotomic.is_irreducible
    assert minimal.degree() == 6 and minimal.is_irreducible

    def rm(value):
        return sp.Poly(value, z, domain=sp.QQ).rem(cyclotomic).as_expr()

    def rl(value):
        return sp.Poly(value, t, domain=sp.QQ).rem(minimal).as_expr()

    def at_m(value, image):
        answer = sp.Integer(0)
        for coefficient in sp.Poly(value, z, domain=sp.QQ).all_coeffs():
            answer = rm(answer*image + coefficient)
        return answer

    z_inverse = sp.invert(sp.Poly(z, z, domain=sp.QQ), cyclotomic).as_expr()
    real_generator = rm(z + z_inverse)
    real_powers = [sp.Integer(1)]
    for _ in range(6):
        real_powers.append(rm(real_powers[-1]*real_generator))
    assert rm(sum(minimal.nth(j)*real_powers[j] for j in range(7))) == 0
    basis = real_powers[:6] + [rm(z*p) for p in real_powers[:6]]
    coordinate_matrix = sp.Matrix([[sp.Poly(p, z).nth(i) for p in basis] for i in range(12)])
    assert coordinate_matrix.det() != 0
    coordinate_inverse = coordinate_matrix.inv()

    def split(value):
        poly = sp.Poly(rm(value), z, domain=sp.QQ)
        coefficients = coordinate_inverse * sp.Matrix([poly.nth(i) for i in range(12)])
        return tuple(sum(coefficients[6*j+i]*t**i for i in range(6)) for j in range(2))

    alpha0 = rm(1+z+z**4)
    alpha11 = rm(z**11 * alpha0)
    conjugate0 = at_m(alpha0, z_inverse)
    conjugate11 = at_m(alpha11, z_inverse)
    x_from_alpha, y_from_alpha = split(alpha11)
    pi, pi_nonreal = split(rm(alpha0*conjugate0))
    trace, trace_nonreal = split(rm(alpha11+conjugate11))
    norm11, norm_nonreal = split(rm(alpha11*conjugate11))
    assert pi_nonreal == trace_nonreal == norm_nonreal == 0
    assert norm11 == pi == rl(t**4+t**3-4*t**2-2*t+5)
    assert trace == rl(2*t**5+t**4-12*t**3-6*t**2+16*t+7)
    x = rl((t**2-2)*(t**3-t**2-6*t-2))
    y = rl(t*(t**2-2)*(t**2+t-1))
    assert (x_from_alpha, y_from_alpha) == (x, y)
    assert rl(2*x+t*y) == trace
    assert rl(x*x+t*x*y+y*y) == pi

    def mr(matrix):
        return matrix.applyfunc(rl)

    def mm(left, right):
        return mr(left*right)

    def integral(value):
        return all(c.q == 1 for c in sp.Poly(rl(value), t, domain=sp.QQ).all_coeffs())

    t_inverse = sp.invert(sp.Poly(t, t, domain=sp.QQ), minimal).as_expr()
    assert integral(t_inverse) and rl(t*t_inverse) == 1
    companion = sp.Matrix([[0, -pi], [1, trace]])
    raw_b = sp.Matrix([[0, -1], [1, t]])
    h, h_inverse = sp.diag(t, 1), sp.diag(t_inverse, 1)
    normalized_b = mm(mm(h, raw_b), h_inverse)
    assert normalized_b == sp.Matrix([[0, -t], [t_inverse, t]])
    g = mr(sp.Matrix([
        [2*t**5-t**4-9*t**3+3*t**2+7*t+1, (t**2-t-1)*(t**3-4*t-1)],
        [-2*(t**5-5*t**3-t**2+5*t+3), -(t-2)*(3*t**4+7*t**3-4*t**2-16*t-7)],
    ]))
    assert rl(g.det()) == 1
    g_inverse = mr(g.adjugate())
    assert mm(g, g_inverse) == mm(g_inverse, g) == sp.eye(2)
    transformed = mm(mm(g, companion), g_inverse)
    affine = mr(x*sp.eye(2)+y*normalized_b)
    assert transformed == affine
    assert mm(transformed, normalized_b) == mm(normalized_b, transformed)
    for matrix in (companion, raw_b, h, h_inverse, normalized_b, g, g_inverse, transformed):
        assert all(integral(value) for value in matrix)
    assert integral(x) and integral(y) and y != 0
    assert rl(normalized_b.det()) == 1 and rl(sp.trace(normalized_b)) == t
    assert rl(companion.det()) == pi and rl(sp.trace(companion)) == trace

    norm_pi = sp.resultant(minimal.as_expr(), pi, t)
    norm_y = sp.resultant(minimal.as_expr(), y, t)
    norm_t = sp.resultant(minimal.as_expr(), t, t)
    assert norm_pi == 43 and norm_y != 0 and norm_t == 1
    discriminant = rl(trace**2-4*pi)
    assert discriminant == rl(y*y*(t*t-4))
    isolated = minimal.intervals(eps=sp.Rational(1, 10)**12)
    assert len(isolated) == 6 and all(multiplicity == 1 for _, multiplicity in isolated)
    assert all(-2 < interval[0] < interval[1] < 2 for interval, _ in isolated)
    # At each of these six real roots, y is nonzero (norm_y!=0), so the
    # discriminant is strictly negative. Also
    # pi=(x+t*y/2)^2+(4-t^2)*y^2/4 is strictly positive.
    assert rl(4*pi-(2*x+t*y)**2-(4-t*t)*y*y) == 0

    def coefficients(value):
        poly = sp.Poly(rl(value), t, domain=sp.QQ)
        return [str(poly.nth(i)) for i in range(6)]

    def matrix_coefficients(matrix):
        return [[coefficients(matrix[i, j]) for j in range(2)] for i in range(2)]

    mathematical_data = {
        "alpha11_pair": [coefficients(x_from_alpha), coefficients(y_from_alpha)],
        "trace11": coefficients(trace), "pi": coefficients(pi),
        "x": coefficients(x), "y": coefficients(y),
        "G": matrix_coefficients(g), "G_inverse": matrix_coefficients(g_inverse),
        "B_t": matrix_coefficients(normalized_b),
        "companion": matrix_coefficients(companion),
        "conjugated_companion": matrix_coefficients(transformed),
        "discriminant": coefficients(discriminant),
        "norm_pi": str(norm_pi), "norm_y": str(norm_y), "norm_t": str(norm_t),
        "det_G": "1", "total_ellipticity": True, "pi_totally_positive": True,
    }
    mathematical_digest = sha256(json.dumps(mathematical_data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert sha256(source.read_bytes()).hexdigest() == before
    print(json.dumps({
        "status": "PASS_NEW_EXACT_K11_BRIDGE_PRIMARY",
        "scope": "Displayed arithmetic vertex bridge; geometry remains a named human-proof dependency; no legacy source recovery",
        "mathematical_data": mathematical_data, "mathematical_digest": mathematical_digest,
        "source_sha256": {source.name: before}, "sympy_version": sp.__version__,
        "real_root_intervals": [[str(a), str(b)] for (a, b), _ in isolated],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
