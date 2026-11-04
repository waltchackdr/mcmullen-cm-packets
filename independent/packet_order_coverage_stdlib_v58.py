#!/usr/bin/env python3
"""Independent exact arithmetic coverage of all 42 packet traces.

Uses only the retained standard-library quotient field, not the SymPy
discovery/certificate path. All coefficients and order groups are computed
from definitions before reading the discovery JSON as a comparison target.
No module-class count, global point count, or geometric decision is checked.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path

import quadrilateral_k17_stdlib_v41 as field


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COMPARISON = ROOT / "experiments/packet_scope/trace_coverage_v58.json"
T, ONE, ZERO = field.T, field.ONE, field.ZERO
PI = T**4 + T**3 - 4*T**2 - 2*T + 5
EXPECTED_KEYS = ((1, 0, 0), (7, 1, 1), (7, 1, 2), (41, 12, 1),
                 (41, 22, 1), (43, 20, 1), (43, 31, 1), (83, 36, 1), (127, 94, 1))


def coefficients(value):
    assert value.integral_power_basis()
    return [int(coefficient) for coefficient in value.coefficients]


def z_powers():
    """z*(a+bz)=-b+(a+tb)z in the displayed quadratic extension."""
    values = [(ONE, ZERO)]
    for _ in range(45):
        a, b = values[-1]
        values.append((-b, a + T*b))
    assert values[21] == (-ONE, ZERO)
    assert values[42] == (ONE, ZERO)
    assert all(values[d] != (ONE, ZERO) for d in (1, 2, 3, 6, 7, 14, 21))
    return values


def norm_integer(value):
    norm = value.norm()
    assert norm.denominator == 1
    return int(norm)


def evaluate_mod(coeffs, residue, prime):
    answer = 0
    for coefficient in reversed(coeffs):
        assert coefficient.denominator == 1
        answer = (answer*residue + int(coefficient)) % prime
    return answer


def factor_integer(value):
    assert value > 1
    factors, divisor = [], 2
    while divisor*divisor <= value:
        exponent = 0
        while value % divisor == 0:
            value //= divisor
            exponent += 1
        if exponent:
            factors.append((divisor, exponent))
        divisor += 1
    if value > 1:
        factors.append((value, 1))
    return factors


def conductor_key(b):
    """Norm-p ideals are pinned by a quotient map, not just by their norm."""
    norm = abs(norm_integer(b))
    if norm == 1:
        assert b.inverse().integral_power_basis()
        return (1, 0, 0)
    factors = factor_integer(norm)
    assert len(factors) == 1, (norm, factors)
    prime, exponent = factors[0]
    roots = [r for r in range(prime)
             if evaluate_mod(field.MODULUS, r, prime) == 0
             and evaluate_mod(b.coefficients, r, prime) == 0]
    assert len(roots) == 1, (prime, roots)
    # For exponent one, (b) lies in the surjective evaluation-map kernel
    # (p,t-r), and both additive lattice indices are p. Exponent two is
    # certified below by a unit ratio to the square of a norm-p generator.
    assert exponent == 1 or (prime, roots[0], exponent) == (7, 1, 2)
    return prime, roots[0], exponent


def monomial(value, exponent=0):
    return {} if value == ZERO else {exponent: value}


def ladd(left, right):
    result = dict(left)
    for exponent, value in right.items():
        result[exponent] = result.get(exponent, ZERO) + value
        if result[exponent] == ZERO:
            del result[exponent]
    return result


def lscale(value, scalar, shift=0):
    return {exponent + shift: coefficient*scalar for exponent, coefficient in value.items()
            if coefficient*scalar != ZERO}


def lmul(left, right):
    result = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            result = ladd(result, monomial(left_coefficient*right_coefficient,
                                          left_exponent + right_exponent))
    return result


def universal_maximal_invariants(a, b):
    """Laurent-polynomial identities for an arbitrary nonzero variable u."""
    aa, bb = monomial(a), monomial(-b, 1)
    cc, dd = monomial(b, -1), monomial(a + T*b)
    trace = ladd(aa, dd)
    determinant = ladd(lmul(aa, dd), lscale(lmul(bb, cc), -ONE))
    assert determinant == monomial(PI)
    inverse_c = monomial(b.inverse(), 1)
    x = lscale(lmul(ladd(aa, lscale(dd, -ONE)), inverse_c), ONE/2)
    discriminant_negative = ladd(lscale(determinant, 4), lscale(lmul(trace, trace), -ONE))
    y2 = lscale(lmul(discriminant_negative, lmul(inverse_c, inverse_c)), ONE/4)
    assert x == monomial(-T/2, 1)
    assert y2 == monomial((4 - T**2)/4, 2)
    assert ladd(lmul(x, x), y2) == monomial(ONE, 2)


def reconstruct():
    assert __debug__, "Assertions must be enabled"
    assert norm_integer(PI) == 43
    assert norm_integer(T**2 - 4) == 1
    roots = field.embedded_roots()
    assert field.sign_vector(4 - T**2, roots) == (1,)*6
    assert field.sign_vector(PI, roots) == (1,)*6
    powers = z_powers()
    alpha = [tuple(sum((powers[k+j][i] for j in (0, 1, 4)), ZERO) for i in range(2))
             for k in range(42)]
    conjugates = [(a + T*b, -b) for a, b in alpha]
    constructed_fiber = alpha + conjugates
    assert len(set(constructed_fiber)) == 84
    trace_fibers = defaultdict(list)
    for a, b in constructed_fiber:
        assert a.integral_power_basis() and b.integral_power_basis()
        assert a*a + T*a*b + b*b == PI
        conjugate = (a + T*b, -b)
        assert conjugate in constructed_fiber and conjugate != (a, b)
        assert (conjugate[0] + T*conjugate[1], -conjugate[1]) == (a, b)
        trace_fibers[2*a + T*b].append((a, b))
    assert len(trace_fibers) == 42 and all(len(values) == 2 for values in trace_fibers.values())
    for k, value in enumerate(alpha):
        a, b = value
        assert set(trace_fibers[2*a + T*b]) == {value, conjugates[k]}
    rows, groups = [], defaultdict(list)
    for k, (a, b) in enumerate(alpha):
        assert a.integral_power_basis() and b.integral_power_basis() and b != ZERO
        assert alpha[(k + 1) % 42] == (-b, a + T*b)
        assert alpha[(k + 21) % 42] == (-a, -b)
        trace = 2*a + T*b
        assert a*a + T*a*b + b*b == PI
        discriminant = trace**2 - 4*PI
        assert discriminant == b*b*(T**2 - 4)
        conductor_norm = abs(norm_integer(b))
        assert abs(norm_integer(discriminant)) == conductor_norm**2
        key = conductor_key(b)
        assert conductor_norm == (1 if key[0] == 1 else key[0]**key[2])
        universal_maximal_invariants(a, b)
        groups[key].append(k)
        rows.append(dict(trace_index=k, a=coefficients(a), b=coefficients(b),
                         trace_coefficients=coefficients(trace), discriminant_coefficients=coefficients(discriminant),
                         conductor_norm=conductor_norm, conductor_key=list(key),
                         opposite_trace_index=(k + 21) % 42))
    assert len({tuple(row["trace_coefficients"]) for row in rows}) == 42
    assert tuple(sorted(groups)) == EXPECTED_KEYS

    # R=OL+OL*alpha=OL+b*OL*z. For u=A+Bz the two conductor tests are
    # B=0 mod b and A+tB=0 mod b. This integral unimodular change of
    # coordinates proves that they are equivalent to A=B=0 mod b.
    conductor_tests = field.matrix(0, 1, 1, T)
    inverse_tests = field.matrix(-T, 1, 1, 0)
    assert field.mdet(conductor_tests) == -ONE
    assert field.mmul(conductor_tests, inverse_tests) == field.IDENTITY
    assert field.all_entries_integral(conductor_tests) and field.all_entries_integral(inverse_tests)

    # A computed norm-seven generator, separate from the primary's t-1
    # generator, identifies the squared prime by exact ideal association.
    seven_generator = alpha[groups[(7, 1, 1)][0]][1]
    assert abs(norm_integer(seven_generator)) == 7
    assert evaluate_mod(seven_generator.coefficients, 1, 7) == 0
    square_witnesses = []
    for k in groups[(7, 1, 2)]:
        ratio = alpha[k][1] / seven_generator**2
        inverse = ratio.inverse()
        assert ratio.integral_power_basis() and inverse.integral_power_basis()
        assert ratio*inverse == ONE and abs(norm_integer(ratio)) == 1
        square_witnesses.append(dict(trace_index=k, ratio=coefficients(ratio), inverse=coefficients(inverse)))

    # Equal-norm conductors at 41 and 43 must remain different embedded ideals.
    for prime, residues in ((41, (12, 22)), (43, (20, 31))):
        for residue, other in (residues, tuple(reversed(residues))):
            b = alpha[groups[(prime, residue, 1)][0]][1]
            assert evaluate_mod(b.coefficients, residue, prime) == 0
            assert evaluate_mod(b.coefficients, other, prime) != 0

    output_groups = []
    for key, indices in sorted(groups.items()):
        representative = indices[0]
        a0, b0 = alpha[representative]
        original = field.matrix(a0, -b0, b0, a0 + T*b0)
        witnesses = []
        for k in indices:
            a, b = alpha[k]
            unit = b/b0
            inverse = unit.inverse()
            shift = a - unit*a0
            inverse_shift = -inverse*shift
            assert all(value.integral_power_basis() for value in (unit, inverse, shift, inverse_shift))
            assert unit*inverse == ONE and abs(norm_integer(unit)) == 1
            assert (shift + unit*a0, unit*b0) == (a, b)
            assert (inverse_shift + inverse*a, inverse*b) == (a0, b0)
            action = field.matrix(a, -b, b, a + T*b)
            transported = field.madd(field.mscale(shift, field.IDENTITY), field.mscale(unit, original))
            assert action == transported
            assert field.madd(field.mscale(inverse_shift, field.IDENTITY), field.mscale(inverse, action)) == original
            witnesses.append(dict(trace_index=k, unit_ratio=coefficients(unit),
                                  inverse_unit=coefficients(inverse), affine_shift=coefficients(shift)))
            rows[k]["equal_embedded_order_representative"] = representative
        output_groups.append(dict(conductor_key=list(key), representative=representative,
                                  trace_indices=indices, affine_unit_witnesses=witnesses))
    return dict(traces=rows, embedded_order_groups=output_groups,
                constructed_norm_pi_elements=[dict(trace_index=k, conjugate=conjugate,
                                                  a=coefficients(value[0]), b=coefficients(value[1]))
                                             for conjugate, values in ((False, alpha), (True, conjugates))
                                             for k, value in enumerate(values)],
                square_prime_generator=coefficients(seven_generator), square_prime_witnesses=square_witnesses,
                common_OM_formula=dict(center="-u*t/2", height_squared="u^2*(4-t^2)/4", modulus_squared="u^2"))


def compare_discovery(derived, path):
    """Read the target only after reconstruct() has completed independently."""
    target = json.loads(path.read_text())
    assert target["common_OM_formula"] == derived["common_OM_formula"]
    reported = target["traces"]
    assert len(reported) == 42 and {row["trace_index"] for row in reported} == set(range(42))
    reported_by_index = {row["trace_index"]: row for row in reported}
    for row in derived["traces"]:
        other = reported_by_index[row["trace_index"]]
        for key in ("trace_coefficients", "opposite_trace_index", "equal_embedded_order_representative"):
            assert row[key] == other[key], (row["trace_index"], key)
        conductor = other["conductor"]
        assert conductor["norm"] == row["conductor_norm"]
        if row["conductor_norm"] == 1:
            assert conductor["ideal"] == "O_L"
        else:
            assert [conductor["prime"], conductor["root_mod_prime"], conductor["exponent"]] == row["conductor_key"]
            assert conductor["generator_coefficients"] == row["b"]
    reported_groups = target["embedded_order_groups"]
    assert len(reported_groups) == 9
    reduced_groups = [{key: group[key] for key in ("conductor_key", "representative", "trace_indices", "affine_unit_witnesses")}
                      for group in reported_groups]
    assert sorted(reduced_groups, key=lambda group: group["conductor_key"]) == derived["embedded_order_groups"]
    # Deliberately do not read claims from known_ambient_points, provisional
    # orientation counts, module coverage, or the geometry status fields.


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O/-OO")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compare-json", type=Path, default=DEFAULT_COMPARISON)
    args = parser.parse_args()
    derived = reconstruct()
    compare_discovery(derived, args.compare_json)
    encoded = json.dumps(derived, sort_keys=True, separators=(",", ":")).encode("ascii")
    report = dict(status="INDEPENDENT_42_TRACE_ORDER_ARITHMETIC_PASS", trace_count=42,
                  constructed_norm_pi_elements=84, conjugation_orbits=42,
                  embedded_order_count=9, conductor_keys=[list(key) for key in EXPECTED_KEYS],
                  order_groups=[dict(conductor_key=g["conductor_key"], trace_indices=g["trace_indices"])
                                for g in derived["embedded_order_groups"]],
                  arithmetic_payload_sha256=sha256(encoded).hexdigest(),
                  compared_json_sha256=sha256(args.compare_json.read_bytes()).hexdigest(),
                  independent_field_sha256=sha256(Path(field.__file__).read_bytes()).hexdigest(),
                  scope="Exact 84-element/42-trace/order/affine/fixed-invariant identities; norm-fiber completeness remains the human proof, not recovery of a missing historical file; no module-class/global-point counts or geometry decisions")
    print(json.dumps(report, sort_keys=True, separators=(",", ":")))
    print("PASS: 84 distinct constructed norm-pi elements in 42 conjugation/trace pairs; nine distinct embedded orders; exact unit/inverse affine witnesses; universal maximal invariant formulas")


if __name__ == "__main__":
    main()
