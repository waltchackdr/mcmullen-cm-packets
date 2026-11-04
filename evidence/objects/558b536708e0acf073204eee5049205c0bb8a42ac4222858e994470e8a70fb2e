#!/usr/bin/env python3
"""Separate exact reconstruction of the displayed k=11 vertex bridge.

Uses only the standard library and the retained independent field/Sturm
module. Packet elements are derived by quadratic-pair recurrence, without
importing the primary certificate or reading its outputs. This is new
source, not a recovery of the absent legacy vertex certificates.
"""

from hashlib import sha256
import json
from pathlib import Path
import sys


FOUNDATION = Path(__file__).with_name("quadrilateral_k17_stdlib_v41.py")
EXPECTED_FOUNDATION = "c5bd45f334a028820a071b5c758acfb8e48afc1e8454b18a3174dd04ca58a38c"


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O/-OO")
    source = Path(__file__)
    before = sha256(source.read_bytes()).hexdigest()
    assert sha256(FOUNDATION.read_bytes()).hexdigest() == EXPECTED_FOUNDATION
    sys.dont_write_bytecode = True
    import quadrilateral_k17_stdlib_v41 as f
    assert Path(f.__file__).resolve() == FOUNDATION.resolve()
    t, one, zero = f.T, f.ONE, f.ZERO

    # A pair (a,b) denotes a+b*z, where z^2=t*z-1. Applying z is the
    # linear recurrence (a,b) -> (-b,a+t*b), distinct from the primary
    # full cyclotomic quotient and twelve-dimensional coordinate solve.
    powers = [(one, zero)]
    for _ in range(42):
        a, b = powers[-1]
        powers.append((-b, a+t*b))
    assert powers[21] == (-one, zero) and powers[42] == (one, zero)
    a0 = sum((powers[j][0] for j in (0, 1, 4)), zero)
    b0 = sum((powers[j][1] for j in (0, 1, 4)), zero)
    a11 = sum((powers[j][0] for j in (11, 12, 15)), zero)
    b11 = sum((powers[j][1] for j in (11, 12, 15)), zero)
    pi = a0*a0+t*a0*b0+b0*b0
    trace = 2*a11+t*b11
    assert a11*a11+t*a11*b11+b11*b11 == pi
    assert pi == t**4+t**3-4*t**2-2*t+5
    assert trace == 2*t**5+t**4-12*t**3-6*t**2+16*t+7

    # Derive the affine coefficients from alpha11; compare with the human
    # display only after the packet recurrence has produced them.
    x, y = a11, b11
    assert x == (t*t-2)*(t**3-t*t-6*t-2)
    assert y == t*(t*t-2)*(t*t+t-1)
    b_vertex = f.matrix(0, -t, t.inverse(), t)
    raw_b = f.matrix(0, -1, 1, t)
    h = f.matrix(t, 0, 0, 1)
    assert f.conjugate_by(h, raw_b) == b_vertex
    companion = f.matrix(0, -pi, 1, trace)

    # These are the displayed witness entries; all operations take place
    # in the separate standard-library field implementation.
    g = f.matrix(
        2*t**5-t**4-9*t**3+3*t*t+7*t+1,
        (t*t-t-1)*(t**3-4*t-1),
        -2*(t**5-5*t**3-t*t+5*t+3),
        -(t-2)*(3*t**4+7*t**3-4*t*t-16*t-7),
    )
    assert f.mdet(g) == one
    g_inverse = f.minverse(g)
    assert f.mmul(g, g_inverse) == f.mmul(g_inverse, g) == f.IDENTITY
    affine = f.madd(f.mscale(x, f.IDENTITY), f.mscale(y, b_vertex))
    # Intertwining is checked first, rather than deriving the affine side
    # from the primary conjugacy calculation.
    assert f.mmul(g, companion) == f.mmul(affine, g)
    transformed = f.conjugate_by(g, companion)
    assert transformed == affine
    assert f.mmul(transformed, b_vertex) == f.mmul(b_vertex, transformed)
    for matrix in (g, g_inverse, b_vertex, raw_b, companion, h, f.minverse(h), transformed):
        assert f.all_entries_integral(matrix)
    assert x.integral_power_basis() and y.integral_power_basis()
    assert f.mdet(companion) == pi and f.mtrace(companion) == trace
    assert f.mdet(b_vertex) == one and f.mtrace(b_vertex) == t

    norm_pi, norm_y, norm_t = pi.norm(), y.norm(), t.norm()
    assert norm_pi == 43 and norm_y != 0 and norm_t == 1
    discriminant = trace*trace-4*pi
    assert discriminant == y*y*(t*t-4)
    roots = f.embedded_roots()
    assert len(roots) == 6
    assert all(-2 < interval.lo < interval.hi < 2 for interval in roots.values())
    assert f.sign_vector(4-t*t, roots) == (1,)*6
    # Independent direct exact interval signs cross-check the structural
    # ellipticity argument. Unresolved intervals fail, never pass by tolerance.
    assert f.sign_vector(discriminant, roots) == (-1,)*6
    assert f.sign_vector(pi, roots) == (1,)*6
    assert y != zero and (2*x+t*y)**2+(4-t*t)*y*y == 4*pi

    def coefficients(value):
        return list(map(str, value.coefficients))

    def matrix_coefficients(matrix):
        return [[coefficients(value) for value in row] for row in matrix]

    mathematical_data = {
        "alpha11_pair": [coefficients(a11), coefficients(b11)],
        "trace11": coefficients(trace), "pi": coefficients(pi),
        "x": coefficients(x), "y": coefficients(y),
        "G": matrix_coefficients(g), "G_inverse": matrix_coefficients(g_inverse),
        "B_t": matrix_coefficients(b_vertex),
        "companion": matrix_coefficients(companion),
        "conjugated_companion": matrix_coefficients(transformed),
        "discriminant": coefficients(discriminant),
        "norm_pi": str(norm_pi), "norm_y": str(norm_y), "norm_t": str(norm_t),
        "det_G": "1", "total_ellipticity": True, "pi_totally_positive": True,
    }
    mathematical_digest = sha256(json.dumps(mathematical_data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert sha256(source.read_bytes()).hexdigest() == before
    assert sha256(FOUNDATION.read_bytes()).hexdigest() == EXPECTED_FOUNDATION
    print(json.dumps({
        "status": "PASS_NEW_EXACT_K11_BRIDGE_STDLIB",
        "scope": "Displayed arithmetic vertex bridge; geometry remains a named human-proof dependency; no legacy source recovery",
        "mathematical_data": mathematical_data, "mathematical_digest": mathematical_digest,
        "source_sha256": {source.name: before, FOUNDATION.name: EXPECTED_FOUNDATION},
        "real_root_intervals": [[str(roots[k].lo), str(roots[k].hi)] for k in f.EMBEDDINGS],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
