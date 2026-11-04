#!/usr/bin/env python3
"""Independent quadratic-model reconstruction of the d=3 arithmetic manifest.

Only the retained standard-library field implementation is imported.  No
degree-12 checker, primary certificate, experiment, or build artifact supplies
definitions.  The manifest is opened only after the complete expected matrix
list and finite arithmetic witnesses have been reconstructed.

The class-number bound, local lattice classification, unit saturation and
SL2 correspondence remain the human proof reviewed in
docs/REFEREE_D3_COMPLETENESS_V58.md.  This is not a geometry certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

import quadrilateral_k17_stdlib_v41 as field

if not __debug__:
    raise RuntimeError("Assertions must be enabled; do not use -O/-OO")

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "certificates/d3_matrix_manifest_v60.jsonl"
EXPECTED_MANIFEST_SHA256 = "fc87260544e8d8418e7bc6fe3ca80ef7d9715a82fcdaea9b2a9716016e00004e"
T, ONE, ZERO = field.T, field.ONE, field.ZERO
PI = T**4 + T**3 - 4*T**2 - 2*T + 5


def coefficients(value):
    assert value.integral_power_basis()
    return [int(c) for c in value.coefficients]


def multiply(left, right):
    a, b = left
    c, d = right
    return a*c - b*d, a*d + b*c + T*b*d


def relative_norm(value):
    a, b = value
    return a*a + T*a*b + b*b


def evaluate_mod(coeffs, root, prime):
    result = 0
    for coefficient in reversed(coeffs):
        assert Fraction(coefficient).denominator == 1
        result = (result*root + int(coefficient)) % prime
    return result


def conductor(b):
    size = abs(b.norm())
    assert size.denominator == 1
    size = int(size)
    assert size in (7, 43, 49)
    prime, exponent = (7, 2) if size == 49 else (size, 1)
    roots = [r for r in range(prime)
             if evaluate_mod(field.MODULUS, r, prime) == 0
             and evaluate_mod(b.coefficients, r, prime) == 0]
    assert len(roots) == 1
    if exponent == 2:
        ratio = b/(T-1)**2
        assert ratio.integral_power_basis() and ratio.inverse().integral_power_basis()
        assert abs(ratio.norm()) == 1 and (T-1).norm() == 7
    return prime, roots[0], exponent


def serialize(row):
    return json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n"


def reconstruct():
    roots = field.embedded_roots()
    assert PI.norm() == 43
    assert field.sign_vector(PI, roots) == (1,)*6
    positive = 4 - T*T
    assert positive.norm() == 1
    assert field.sign_vector(positive, roots) == (1,)*6
    assert evaluate_mod(field.MODULUS, 10, 41) == 0
    residue = evaluate_mod(positive.coefficients, 10, 41)
    assert pow(residue, 20, 41) == 40
    derivative = field.LElement([i*c for i, c in enumerate(field.MODULUS)][1:])
    assert -derivative.norm() == 453789
    assert Fraction(factorial(6)**2*453789, 6**12) < Fraction(52, 5)**2
    assert Fraction(factorial(12)*453789, 12**12)*Fraction(200, 157)**6 < 105

    powers = [(ONE, ZERO)]
    for _ in range(45):
        a, b = powers[-1]
        powers.append((-b, a+T*b))
    assert powers[42] == (ONE, ZERO)
    assert all(powers[j] != (ONE, ZERO) for j in (1, 2, 3, 6, 7, 14, 21))
    alpha = [tuple(sum((powers[k+j][i] for j in (0, 1, 4)), ZERO) for i in (0, 1))
             for k in range(42)]
    traces = [2*a+T*b for a, b in alpha]
    assert len(set(traces)) == 42
    assert all(relative_norm(value) == PI for value in alpha)
    d3 = []
    for k, (a, b) in enumerate(alpha):
        assert b != ZERO
        assert traces[k]**2 - 4*PI == b*b*(T*T-4)
        phi = [r if sign > 0 else 42-r
               for r, sign in zip(field.EMBEDDINGS, field.sign_vector(b, roots))]
        if sum(r % 3 == 1 for r in phi) == 3:
            d3.append(k)
    selected = [k for k in d3 if abs(alpha[k][1].norm()) != 1]
    assert selected == [0, 4, 8, 12, 15, 18, 21, 25, 29, 33, 36, 39]
    conductors = {k: conductor(alpha[k][1]) for k in selected}

    generators = []
    for a in (5, 11, 13, 17):
        total = tuple(sum((powers[(2*j) % 42][i] for j in range(a)), ZERO)
                      for i in (0, 1))
        real, imaginary = multiply(powers[(1-a) % 42], total)
        assert imaginary == ZERO and abs(real.norm()) == 1
        assert field.sign_vector(real, roots) == tuple((-1)**((a*r)//21) for r in field.EMBEDDINGS)
        generators.append(real)
    generators.extend((-ONE, positive))
    units, inverses, signs = [], [], []
    for mask in range(64):
        unit = ONE
        for bit, generator in enumerate(generators):
            if mask & (1 << bit):
                unit *= generator
        inverse = unit.inverse()
        assert unit.integral_power_basis() and inverse.integral_power_basis()
        assert unit*inverse == ONE
        units.append(unit)
        inverses.append(inverse)
        signs.append(field.sign_vector(unit, roots))
    assert len(set(signs[:32])) == 32
    assert len(set(units)) == 64

    # Separate quadratic-model checks of the class-number generators.
    norm_witnesses = []
    for terms, expected in ((((0, 1), (6, 1), (12, 1), (24, 1)), 64),
                            (((0, 1), (1, 1), (3, -1)), 7),
                            (((0, 1), (1, 1), (4, 1)), 43)):
        value = tuple(sum((c*powers[e][i] for e, c in terms), ZERO) for i in (0, 1))
        assert relative_norm(value).norm() == expected
        norm_witnesses.append(dict(a=coefficients(value[0]), b=coefficients(value[1]), norm=expected))
    epsilon = (ONE-powers[2][0], -powers[2][1])
    epsilon_bar = (epsilon[0]+T*epsilon[1], -epsilon[1])
    assert relative_norm(epsilon) == positive
    assert multiply(powers[23], epsilon_bar) == epsilon

    matrices = []
    for k in selected:
        p, r, exponent = conductors[k]
        a, b = alpha[k]
        for j in range(exponent+1):
            if j == 0:
                template = (a, -b, b, a+T*b)
            elif j == exponent:
                template = (ZERO, -PI, ONE, traces[k])
            else:
                assert (p, r, exponent, j) == (7, 1, 2, 1)
                template = (a, -(T-1)*b, b/(T-1), a+T*b)
            for mask in range(32 if j == 0 else 64):
                aa, bb, cc, dd = template
                matrix = (aa, units[mask]*bb, inverses[mask]*cc, dd)
                aa, bb, cc, dd = matrix
                assert aa+dd == traces[k] and aa*dd-bb*cc == PI
                matrices.append(dict(trace_index=k, conductor=[p, r, exponent],
                                     multiplier="OM" if j == 0 else [p, r, j],
                                     orientation_mask=mask, unit=coefficients(units[mask]),
                                     unit_signs=list(signs[mask]),
                                     matrix_row_major_L_coefficients=[coefficients(v) for v in matrix]))
    expected_labels = {(k, j, m) for k in selected for j in range(conductors[k][2]+1)
                       for m in range(32 if j == 0 else 64)}
    labels = [(row["trace_index"], 0 if row["multiplier"] == "OM" else row["multiplier"][2],
               row["orientation_mask"]) for row in matrices]
    assert len(labels) == len(set(labels)) == 1280 and set(labels) == expected_labels
    assert sum(row["unit_signs"] == [1]*6 for row in matrices) == 40
    raw = "".join(serialize(row) for row in matrices).encode("ascii")
    assert sha256(raw).hexdigest() == EXPECTED_MANIFEST_SHA256
    summary = dict(arithmetic_d3_indices=d3, selected_nonmaximal_indices=selected,
                   conductors={str(k): list(v) for k, v in conductors.items()},
                   counts_per_trace=dict(Counter(str(row["trace_index"]) for row in matrices)),
                   matrix_count=1280, positive_sector_count=40,
                   class_number_norm_witnesses=norm_witnesses)
    return raw, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    expected, summary = reconstruct()
    observed = args.manifest.read_bytes()
    assert observed == expected, "Manifest differs from independently reconstructed canonical matrix list"
    report = dict(status="PASS_INDEPENDENT_D3_QUADRATIC_MANIFEST", **summary,
                  manifest_sha256=sha256(observed).hexdigest(),
                  implementation_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  field_sha256=sha256(Path(field.__file__).read_bytes()).hexdigest(),
                  scope="Exact finite inputs and 1280 labelled representatives; class completeness uses the reviewed human proof; no ambient-point deduplication or geometry")
    print(json.dumps(report, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
