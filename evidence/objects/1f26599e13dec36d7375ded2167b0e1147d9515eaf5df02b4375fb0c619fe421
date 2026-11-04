#!/usr/bin/env python3
"""Exact small regression supporting the human known-CM comparison, C35.

Registered exact regression, not a reconstruction of the geometric inputs.
The three CM types below are INPUTS from existing certificates. This script
does not reconstruct vertex geometry, integral lattices, or polarizations.
It imports no repository mathematical implementation and uses no numerics.
"""

from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import json


if not __debug__:
    raise RuntimeError("Run without Python optimization: assertions are required.")

PHI_A = frozenset((1, 5, 17, 19, 29, 31))
PHI_B = frozenset((1, 5, 23, 25, 29, 31))
PHI_C = frozenset((1, 5, 11, 13, 17, 19))
S21 = frozenset((1, 2, 4, 5, 8, 10))
J7_TYPE = frozenset((1, 2, 3))


def units(n):
    return tuple(a for a in range(n) if gcd(a, n) == 1)


def multiply_type(values, multiplier, modulus):
    return frozenset(multiplier * a % modulus for a in values)


def stabilizer(values, modulus):
    return tuple(a for a in units(modulus)
                 if multiply_type(values, a, modulus) == values)


# Z[r]/(1+r+...+r^6); all coefficient tuples have length six.
def reduce_poly(coefficients):
    coefficients = list(coefficients)
    if len(coefficients) < 6:
        coefficients += [0] * (6 - len(coefficients))
    for degree in range(len(coefficients) - 1, 5, -1):
        leading = coefficients[degree]
        for shift in range(7):
            coefficients[degree - 6 + shift] -= leading
    assert all(a == 0 for a in coefficients[6:])
    return tuple(coefficients[:6])


ZERO = reduce_poly((0,))
ONE = reduce_poly((1,))
RING_R = reduce_poly((0, 1))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def negate(a):
    return tuple(-x for x in a)


def multiply(a, b):
    coefficients = [0] * 11
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coefficients[i + j] += x * y
    return reduce_poly(coefficients)


# Matrices are row-major tuples of four ring elements.
IDENTITY = (ONE, ZERO, ZERO, ONE)
ZERO_MATRIX = (ZERO,) * 4


def matrix_multiply(a, b):
    return tuple(add(multiply(a[2*i], b[j]),
                     multiply(a[2*i + 1], b[2+j]))
                 for i in range(2) for j in range(2))


def matrix_power(a, exponent):
    result = IDENTITY
    for _ in range(exponent):
        result = matrix_multiply(result, a)
    return result


def matrix_add(a, b):
    return tuple(add(x, y) for x, y in zip(a, b))


def matrix_scale(a, scalar):
    return tuple(tuple(scalar * x for x in value) for value in a)


def main():
    for values in (PHI_A, PHI_B, PHI_C):
        assert len(values) == 6
        assert values <= set(units(42))
        assert values.isdisjoint(multiply_type(values, -1, 42))
        assert values | multiply_type(values, -1, 42) == set(units(42))

    assert {a % 21 for a in PHI_B} == S21
    b_map = {j: (21 + 22*j) % 42 for j in sorted(S21)}
    assert frozenset(b_map.values()) == PHI_B
    assert multiply_type(PHI_B, 19, 42) == PHI_C
    assert 19 * 31 % 42 == 1
    assert stabilizer(S21, 21) == (1,)
    assert stabilizer(PHI_B, 42) == stabilizer(PHI_C, 42) == (1,)

    psi_a = frozenset((1, 3, 5))
    assert Counter(a % 7 for a in PHI_A) == Counter({1: 2, 3: 2, 5: 2})
    assert multiply_type(J7_TYPE, 5, 7) == psi_a
    assert frozenset(a for a in units(42) if a % 7 in psi_a) == PHI_A
    assert tuple(a for a in units(42) if a % 7 == 1) == (1, 29)
    assert stabilizer(PHI_A, 42) == (1, 29)
    assert stabilizer(J7_TYPE, 7) == (1,)

    r = (RING_R, ZERO, ZERO, RING_R)
    u = (ZERO, negate(ONE), ONE, negate(ONE))
    assert matrix_power(r, 7) == IDENTITY
    assert matrix_add(matrix_add(matrix_power(u, 2), u), IDENTITY) == ZERO_MATRIX
    assert matrix_power(u, 3) == IDENTITY
    assert matrix_multiply(r, u) == matrix_multiply(u, r)
    t = matrix_scale(matrix_multiply(matrix_power(r, 2), matrix_power(u, 2)), -1)
    r2 = multiply(RING_R, RING_R)
    assert t == (r2, negate(r2), r2, ZERO)
    for exponent, expected in (
        (6, matrix_power(r, 5)), (14, u), (18, r),
        (21, matrix_scale(IDENTITY, -1)), (42, IDENTITY)
    ):
        assert matrix_power(t, exponent) == expected

    # Phi_42(X)=X^12+X^11-X^9-X^8+X^6-X^4-X^3+X+1.
    phi42 = {12: 1, 11: 1, 9: -1, 8: -1, 6: 1, 4: -1, 3: -1, 1: 1, 0: 1}
    value = ZERO_MATRIX
    for exponent, coefficient in phi42.items():
        value = matrix_add(value, matrix_scale(matrix_power(t, exponent), coefficient))
    assert value == ZERO_MATRIX

    a_map = {f"j={j},b={b}": (21 + 12*j + 28*b) % 42
             for j in sorted(J7_TYPE) for b in (1, 2)}
    assert frozenset(a_map.values()) == PHI_A
    assert all(gcd(a, 42) == 1 for a in a_map.values())

    # Bind the existing input files without importing or executing them.
    root = Path(__file__).resolve().parents[1]
    input_paths = (
        "certificates/a_vertex_orbits_k20_k41_v44.py",
        "certificates/a_vertex_certificate_v31.py",
        "certificates/vertex_saturation_v35.py",
        "certificates/four_vertex_audit_v40.py",
    )
    hashes = {}
    for path in input_paths:
        source = root / path
        if not source.is_file():
            raise FileNotFoundError(source)
        hashes[path] = sha256(source.read_bytes()).hexdigest()

    print(json.dumps({
        "status": "PASS_EXACT_CM_COMPARISON_REGRESSION",
        "scope": "comparison of supplied geometric CM types; not geometry or polarization",
        "types_mod42": {"A": sorted(PHI_A), "B": sorted(PHI_B), "C": sorted(PHI_C)},
        "X21_type_mod21": sorted(S21),
        "X21_z42_character_map": b_map,
        "B_to_C_multiplier_mod42": 19,
        "A_restricted_type_mod7": sorted(psi_a),
        "J7_to_A_multiplier_mod7": 5,
        "A_product_character_map": a_map,
        "stabilizers_mod42": {"A": list(stabilizer(PHI_A, 42)),
                                "B": list(stabilizer(PHI_B, 42)),
                                "C": list(stabilizer(PHI_C, 42))},
        "product_matrix_identities": ["T^6=R^5", "T^14=U", "T^18=R",
                                      "T^21=-I", "T^42=I", "Phi42(T)=0"],
        "input_source_sha256": hashes,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
