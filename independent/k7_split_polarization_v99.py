#!/usr/bin/env python3
"""Exact independent cyclotomic reconstruction of the K7 polarization.

Standard library only. Signs use the exact locations of cyclotomic angles,
not numerical roots or the primary polynomial interval implementation.
This is an arithmetic regression supporting a separately reviewed human
polarization proof, not a new McMullen curve-point certificate.
"""
from fractions import Fraction
from math import gcd
import json


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def divide(a, b):
    a = list(a)
    q = [0] * max(1, len(a)-len(b)+1)
    assert b[-1] == 1
    while len(a) >= len(b) and a != [0]:
        k = len(a)-len(b)
        q[k] = a[-1]
        for j, v in enumerate(b):
            a[j+k] -= q[k]*v
        trim(a)
    return trim(q), trim(a)


def cyclotomic(n):
    polys = {}
    for k in range(1, n+1):
        if n % k:
            continue
        p = [-1] + [0]*(k-1) + [1]
        for d, f in polys.items():
            if k % d == 0:
                p, rem = divide(p, f)
                assert rem == [0]
        polys[k] = p
    return polys[n]


PHI = cyclotomic(42)
DEG = len(PHI)-1


def reduce(a):
    _, rem = divide(trim(list(a)), PHI)
    return tuple(rem + [0]*(DEG-len(rem)))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, k):
    return tuple(k*x for x in a)


def mul(a, b):
    p = [0]*(2*DEG-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p[i+j] += x*y
    return reduce(p)


def zpower(k):
    k %= 42
    return reduce([0]*k+[1])


ZERO = reduce([0])
ONE = reduce([1])


def evaluate(coefficients, x):
    y = ZERO
    for c in reversed(coefficients):
        y = add(mul(y, x), scale(ONE, c))
    return y


def conjugate(a, r):
    y = ZERO
    for i, c in enumerate(a):
        y = add(y, scale(zpower(i*r), c))
    return y


def determinant(matrix):
    rows = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for i in range(len(rows)):
        pivot = next(j for j in range(i, len(rows)) if rows[j][i])
        if pivot != i:
            rows[i], rows[pivot] = rows[pivot], rows[i]
            result = -result
        v = rows[i][i]
        result *= v
        for j in range(i+1, len(rows)):
            q = rows[j][i]/v
            rows[j] = [a-q*b for a, b in zip(rows[j], rows[i])]
    assert result.denominator == 1
    return int(result)


def real_trace(x, representatives):
    y = ZERO
    for r in representatives:
        y = add(y, conjugate(x, r))
    assert not any(y[1:])
    return y[0]


def real_norm(x, representatives):
    y = ONE
    for r in representatives:
        y = mul(y, conjugate(x, r))
    assert not any(y[1:])
    return y[0]


def chi3(r):
    return 1 if r % 3 == 1 else -1


def chi7(r):
    return 1 if r % 7 in (1, 2, 4) else -1


def cosine_sign(r):
    r %= 42
    assert gcd(r, 42) == 1
    # cos(2*pi*r/42)>0 exactly when its distance from 42Z is <42/4.
    return 1 if 4*min(r, 42-r) < 42 else -1


def check():
    units = tuple(r for r in range(42) if gcd(r, 42) == 1)
    reps = tuple(r for r in units if r < 21)
    cm_type = {1, 5, 11, 13, 19, 25}
    s3 = add(ONE, scale(zpower(14), 2))
    s7 = ZERO
    for a in range(1, 7):
        s7 = add(s7, scale(zpower(6*a), chi7(a)))
    assert mul(s3, s3) == scale(ONE, -3)
    assert mul(s7, s7) == scale(ONE, -7)
    sqrt21 = scale(mul(s3, s7), -1)
    assert mul(sqrt21, sqrt21) == scale(ONE, 21)
    ts = {a: add(zpower(a), zpower(-a)) for a in (1, 13, 17)}
    t = ts[1]
    assert evaluate([1, 8, 8, -6, -6, 1, 1], t) == ZERO
    u = mul(mul(t, ts[13]), ts[17])
    b = mul(sqrt21, u)
    assert b == evaluate([3, 4, -13, -15, 3, 4], t)
    assert sqrt21 == evaluate([-3, 12, 8, -10, -2, 2], t)
    assert conjugate(b, 41) == b
    sign_rows = []
    for r in reps:
        assert conjugate(s3, r) == scale(s3, chi3(r))
        assert conjugate(s7, r) == scale(s7, chi7(r))
        sign = chi3(r)*chi7(r)
        for a in (1, 13, 17):
            sign *= cosine_sign(r*a)
        selection = 1 if r in cm_type else -1
        assert sign == chi7(r)*selection
        sign_rows.append({"embedding": r, "b_sign": sign,
                          "selection_sign": selection, "chi7": chi7(r)})
    assert sum(chi7(r) == 1 for r in cm_type) == 3
    norms = {str(a): real_norm(v, reps) for a, v in ts.items()}
    assert set(norms.values()) == {1}
    norm_b = real_norm(b, reps)
    assert norm_b == -21**3
    tpowers = [ONE]
    for _ in range(10):
        tpowers.append(mul(tpowers[-1], t))
    trace_matrix = [[real_trace(tpowers[i+j], reps) for j in range(6)] for i in range(6)]
    hermitian = [[real_trace(mul(b, tpowers[i+j]), reps) for j in range(6)] for i in range(6)]
    disc = determinant(trace_matrix)
    det = determinant(hermitian)
    assert disc == 3**3*7**5
    assert det == disc*norm_b == -(3**3*7**4)**2
    return {"status": "PASS_EXACT_K7_SPLIT_POLARIZATION_ARITHMETIC",
            "method": "cyclotomic integer quotient, exact Galois products and trace determinants; elementary angle signs",
            "cyclotomic_polynomial_ascending": PHI,
            "b_cyclotomic_coordinates": list(b), "signs": sign_rows,
            "unit_norms": norms, "norm_b": norm_b, "trace_discriminant": disc,
            "trace_hermitian_matrix": hermitian, "hermitian_determinant": det,
            "negative_determinant_square_root": 3**3*7**4,
            "identity_embedding_convention": "s3=i*sqrt3 and s7=i*sqrt7; positivity of s7 follows from sin(2pi/7)>sin(pi/7) and sin(4pi/7)>0",
            "geometry_claim": False, "principal_polarization_claim": False,
            "scope": "Arithmetic support for human Riemann-form proof and theorem-hypothesis review; not proof of Markman's theorem"}


if __name__ == "__main__":
    assert __debug__
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_CPU, (15, 16))
    signal.alarm(25)
    print(json.dumps(check(), indent=2, sort_keys=True))
