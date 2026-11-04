#!/usr/bin/env python3
"""Exact finite inputs and representatives for the reviewed d=3 completeness proof.

No repository mathematics imports, floating arithmetic, class-group black box,
or geometry search. Human completeness is proved in docs/D3_COMPLETENESS_AUDIT_2026-09-18.md.
The cyclotomic degree-12 model reconstructs coefficients in the relative
basis; integer determinants and rational root brackets supply exact witnesses.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import factorial, gcd
from pathlib import Path
import argparse
import json

if not __debug__:
    raise RuntimeError("Assertions must be enabled")


def determinant(rows):
    a = [list(row) for row in rows]
    n, sign, previous = len(a), 1, 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j]*a[k][k] - a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        previous = a[k][k]
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign*a[-1][-1]


def inverse_matrix(rows):
    n = len(rows)
    a = [[Q(v) for v in row] + [Q(i == j) for j in range(n)]
         for i, row in enumerate(rows)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if a[i][k])
        a[k], a[pivot] = a[pivot], a[k]
        divisor = a[k][k]
        a[k] = [v/divisor for v in a[k]]
        for i in range(n):
            if i != k:
                multiplier = a[i][k]
                a[i] = [x - multiplier*y for x, y in zip(a[i], a[k])]
    return [row[n:] for row in a]


class Ring:
    def __init__(self, modulus):
        assert modulus[-1] == 1
        self.modulus = tuple(modulus)
        self.n = len(modulus) - 1
        self.zero = self.reduce((0,))
        self.one = self.reduce((1,))
        self.x = self.reduce((0, 1))

    def reduce(self, coefficients):
        a = list(coefficients) + [0]*max(0, self.n - len(coefficients))
        for i in range(len(a) - 1, self.n - 1, -1):
            leading = a[i]
            for j, coefficient in enumerate(self.modulus):
                a[i - self.n + j] -= leading*coefficient
        assert all(v == 0 for v in a[self.n:])
        return tuple(a[:self.n])

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def scale(self, a, integer):
        return tuple(integer*x for x in a)

    def mul(self, a, b):
        c = [0]*(2*self.n - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i+j] += x*y
        return self.reduce(c)

    def pow(self, a, exponent):
        out = self.one
        while exponent:
            if exponent & 1:
                out = self.mul(out, a)
            a = self.mul(a, a)
            exponent //= 2
        return out

    def multiplication_matrix(self, a):
        columns = [self.mul(a, self.reduce((0,)*i + (1,))) for i in range(self.n)]
        return list(zip(*columns))

    def norm(self, a):
        return determinant(self.multiplication_matrix(a))

    def divide_integrally(self, a, b):
        inverse = inverse_matrix(self.multiplication_matrix(b))
        answer = tuple(sum(x*y for x, y in zip(row, a)) for row in inverse)
        assert all(x.denominator == 1 for x in answer)
        answer = tuple(int(x) for x in answer)
        assert self.mul(answer, b) == a
        return answer


L = Ring((1, 8, 8, -6, -6, 1, 1))
M = Ring((1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1))
EMBEDDINGS = (1, 5, 11, 13, 17, 19)
z = M.x
t_in_m = M.add(z, M.pow(z, 41))
relative_basis = [M.pow(t_in_m, j) for j in range(6)]
relative_basis += [M.mul(z, v) for v in relative_basis]
change_of_basis = list(zip(*relative_basis))
assert abs(determinant(change_of_basis)) == 1
inverse_change = inverse_matrix(change_of_basis)
assert all(value.denominator == 1 for row in inverse_change for value in row)


def relative_coordinates(a):
    values = tuple(int(sum(x*y for x, y in zip(row, a))) for row in inverse_change)
    reconstructed = M.zero
    for coefficient, basis_vector in zip(values, relative_basis):
        reconstructed = M.add(reconstructed, M.scale(basis_vector, coefficient))
    assert reconstructed == a
    return values[:6], values[6:]


def polynomial_at(coefficients, value):
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer*value + coefficient
    return answer


def interval_value(coefficients, interval):
    lo, hi = Q(0), Q(0)
    for coefficient in reversed(coefficients):
        products = (lo*interval[0], lo*interval[1], hi*interval[0], hi*interval[1])
        lo, hi = min(products) + coefficient, max(products) + coefficient
    return lo, hi


def root_intervals():
    # Six disjoint sign-changing rational brackets exhaust a degree-six f.
    starts = ((197, 198), (146, 147), (-15, -14), (-74, -73), (-166, -165), (-192, -191))
    intervals = []
    for a, b in starts:
        lo, hi = Q(a, 100), Q(b, 100)
        flo, fhi = polynomial_at(L.modulus, lo), polynomial_at(L.modulus, hi)
        assert flo*fhi < 0
        for _ in range(36):
            middle = (lo + hi)/2
            fm = polynomial_at(L.modulus, middle)
            assert fm != 0
            if flo*fm < 0:
                hi, fhi = middle, fm
            else:
                lo, flo = middle, fm
        intervals.append((lo, hi))
    assert all(intervals[i][0] > intervals[i+1][1] for i in range(5))
    return intervals


INTERVALS = root_intervals()


def sign_vector(a):
    result = []
    for interval in INTERVALS:
        lo, hi = interval_value(a, interval)
        assert lo > 0 or hi < 0
        result.append(1 if lo > 0 else -1)
    return tuple(result)


def pair_mul(x, y, p):
    return ((x[0]*y[0]) % p, (x[0]*y[1] + x[1]*y[0]) % p)


def pair_pow(x, exponent, p):
    out = (1, 0)
    for _ in range(exponent):
        out = pair_mul(out, x, p)
    return out


def finite_local_data(p, root, exponent):
    if exponent == 1:
        roots = [a for a in range(p) if (a*a-root*a+1) % p == 0]
        assert len(roots) == 2 and roots[0] != roots[1]
        a = roots[0]
        epsilon = -a*a % p
        order = next(d for d in range(1, p) if pow(epsilon, d, p) == 1)
        assert order == p-1
        assert pow(epsilon, 22, p) == a*a % p
        return {"root_z": a, "epsilon_ratio": epsilon, "quotient_order": p-1,
                "epsilon_order": order}
    assert (p, root, exponent) == (7, 1, 2)
    roots = [(a, b) for a, b in product(range(7), repeat=2)
             if (a*a-a+1) % 7 == 0 and ((2*a-1)*b-a) % 7 == 0]
    assert roots == [(3, 2), (5, 6)]
    a = roots[0]
    square = pair_mul(a, a, 7)
    epsilon = tuple(-x % 7 for x in square)
    order = next(d for d in range(1, 43) if pair_pow(epsilon, d, 7) == (1, 0))
    assert order == 42 and pair_pow(epsilon, 22, 7) == square
    return {"root_z": a, "epsilon_ratio": epsilon, "quotient_order": 42,
            "epsilon_order": order}


def serialize(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def ring_polynomial(ring, coefficients, value):
    result = ring.zero
    for coefficient in reversed(coefficients):
        result = ring.add(ring.mul(result, value), ring.scale(ring.one, coefficient))
    return result


def prime_norms_below_105():
    result = {}
    for p in range(2, 105):
        if any(p % q == 0 for q in range(2, p)):
            continue
        if p == 3:
            degree = 6  # order of 3 modulo 7
            assert next(d for d in range(1, 7) if pow(3, d, 7) == 1) == 6
        elif p == 7:
            degree = 1  # order of 7 modulo 3
        else:
            degree = next(d for d in range(1, 13) if pow(p, d, 21) == 1)
        if p**degree < 105:
            result[p] = p**degree
    assert result == {2: 64, 7: 7, 43: 43}
    return result


def polynomial_product(a, b):
    result = [0]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return tuple(result)


def polynomial_remainder_mod(a, b, p):
    a = [x % p for x in a]
    assert b[-1] % p == 1
    for i in range(len(a)-1, len(b)-2, -1):
        leading = a[i]
        for j, coefficient in enumerate(b):
            a[i-len(b)+1+j] = (a[i-len(b)+1+j] - leading*coefficient) % p
    return tuple(a[:len(b)-1])


def verify_real_small_prime_factors():
    for degree in (1, 2, 3):
        for coefficients in product(range(2), repeat=degree):
            assert any(polynomial_remainder_mod(L.modulus, coefficients+(1,), 2))
    factors = {3: ((1,1,-1,1), (1,1,-1,1)),
               5: ((-2,-1,0,1), (2,0,1,1))}
    for p, pair in factors.items():
        assert all((a-b) % p == 0 for a,b in zip(L.modulus, polynomial_product(*pair)))
        for factor in pair:
            assert all(polynomial_at(factor, a) % p != 0 for a in range(p))
    power = (1,)
    for _ in range(6):
        power = polynomial_product(power, (-1,1))
    assert all((a-b) % 7 == 0 for a,b in zip(L.modulus, power))


def main(output_directory):
    pi = L.reduce((5, -2, -4, 1, 1))
    gamma7 = L.reduce((-1, 1))
    positive = L.reduce((4, 0, -1))
    assert L.norm(pi) == 43
    assert L.norm(gamma7) == 7 and L.norm(positive) == 1
    assert sign_vector(positive) == (1,)*6
    assert polynomial_at(L.modulus, 10) % 41 == 0
    assert polynomial_at(positive, 10) % 41 == 27
    assert pow(27, 20, 41) == 40  # exact positive nonsquare witness
    assert Q(factorial(6)**2*453789, 6**12) < Q(52, 5)**2
    assert Q(factorial(12)*453789, 12**12)*Q(200, 157)**6 < 105
    assert ring_polynomial(M, L.modulus, t_in_m) == M.zero
    derivative_l = tuple(i*c for i,c in enumerate(L.modulus))[1:]
    derivative_m = tuple(i*c for i,c in enumerate(M.modulus))[1:]
    assert -L.norm(L.reduce(derivative_l)) == 453789
    assert M.norm(M.reduce(derivative_m)) == 453789**2
    epsilon = M.add(M.one, M.scale(M.pow(z, 2), -1))
    epsilon_bar = M.add(M.one, M.scale(M.pow(z, 40), -1))
    assert M.norm(epsilon) == 1
    assert epsilon == M.mul(M.pow(z, 23), epsilon_bar)
    verify_real_small_prime_factors()

    alpha0 = M.add(M.add(M.one, z), M.pow(z, 4))
    all_data, d3, maximal = [], [], []
    for k in range(42):
        a, b = relative_coordinates(M.mul(M.pow(z, k), alpha0))
        c = L.add(L.scale(a, 2), L.mul(L.x, b))
        norm = L.add(L.add(L.mul(a, a), L.mul(L.x, L.mul(a, b))), L.mul(b, b))
        assert norm == pi
        nb = abs(L.norm(b))
        signs = sign_vector(b)
        phi = [r if s > 0 else 42-r for r, s in zip(EMBEDDINGS, signs)]
        if sum(r % 3 == 1 for r in phi) == 3:
            d3.append(k)
        if nb == 1:
            maximal.append(k)
        all_data.append((a, b, c, nb))
    assert len({data[2] for data in all_data}) == 42
    assert d3 == [0,1,4,8,11,12,15,16,17,18,19,20,21,22,25,29,32,33,36,37,38,39,40,41]
    assert maximal == [1,7,11,13,16,17,19,20,22,28,32,34,37,38,40,41]
    selected = [k for k in d3 if k not in maximal]
    assert selected == [0,4,8,12,15,18,21,25,29,33,36,39]

    expected = {k: (7, 1, 1) for k in (0,15,18,21,36,39)}
    expected.update({k: (7, 1, 2) for k in (12,33)})
    expected.update({k: (43, 20, 1) for k in (4,25)})
    expected.update({k: (43, 31, 1) for k in (8,29)})
    for k in selected:
        p, root, exponent = expected[k]
        _, b, _, nb = all_data[k]
        assert nb == p**exponent
        roots = [r for r in range(p) if polynomial_at(L.modulus, r) % p == 0]
        vanishing = [r for r in roots if polynomial_at(b, r) % p == 0]
        assert vanishing == [root]
        if p == 7:
            quotient = L.divide_integrally(b, L.pow(gamma7, exponent))
            assert abs(L.norm(quotient)) == 1
    assert polynomial_at(L.modulus, 1) == 7
    assert polynomial_at(tuple(i*c for i,c in enumerate(L.modulus))[1:], 1) == -7

    rho = M.pow(z, 2)
    generators, expected_signs = [], []
    for a in (5, 11, 13, 17):
        geometric_sum = M.zero
        for j in range(a):
            geometric_sum = M.add(geometric_sum, M.pow(rho, j))
        xi = M.mul(M.pow(rho, ((1-a)//2) % 21), geometric_sum)
        coefficients, zero = relative_coordinates(xi)
        assert zero == L.zero and abs(L.norm(coefficients)) == 1
        signs = tuple((-1)**((a*k)//21) for k in EMBEDDINGS)
        assert sign_vector(coefficients) == signs
        generators.append(coefficients)
        expected_signs.append(signs)
    generators.extend((L.scale(L.one, -1), positive))
    inverses = [L.divide_integrally(L.one, a) for a in generators]
    units, inverse_units = [], []
    for mask in range(64):
        value, inverse = L.one, L.one
        for bit, (g, gi) in enumerate(zip(generators, inverses)):
            if mask & (1 << bit):
                value, inverse = L.mul(value, g), L.mul(inverse, gi)
        assert L.mul(value, inverse) == L.one
        units.append(value)
        inverse_units.append(inverse)
    assert len({sign_vector(v) for v in units[:32]}) == 32
    assert len(set(units)) == 64

    locals_ = {f"{p}:{root}:{exponent}": finite_local_data(p, root, exponent)
               for p, root, exponent in sorted(set(expected.values()))}
    manifest, counts = [], {}
    for k in selected:
        p, root, exponent = expected[k]
        a, b, c, _ = all_data[k]
        for j in range(exponent + 1):
            if j == 0:
                template = (a, L.scale(b, -1), b, L.add(a, L.mul(L.x, b)))
            elif j == exponent:
                template = (L.zero, L.scale(pi, -1), L.one, c)
            else:
                assert (p, root, exponent, j) == (7, 1, 2, 1)
                template = (a, L.scale(L.mul(gamma7, b), -1),
                            L.divide_integrally(b, gamma7), L.add(a, L.mul(L.x, b)))
            for mask in range(32 if j == 0 else 64):
                aa, bb, cc, dd = template
                matrix = (aa, L.mul(units[mask], bb),
                          L.mul(inverse_units[mask], cc), dd)
                aa, bb, cc, dd = matrix
                assert L.add(aa, dd) == c
                assert L.add(L.mul(aa, dd), L.scale(L.mul(bb, cc), -1)) == pi
                manifest.append({"trace_index": k, "conductor": [p, root, exponent],
                                 "multiplier": "OM" if j == 0 else [p, root, j],
                                 "orientation_mask": mask, "unit": units[mask],
                                 "unit_signs": sign_vector(units[mask]),
                                 "matrix_row_major_L_coefficients": matrix})
        counts[k] = 32 + 64*exponent
    assert len(manifest) == sum(counts.values()) == 1280
    assert sum(1+2*expected[k][2] for k in selected) == 40

    # A norm-64 generator for either prime above 2 in M.
    prime2_generator = M.one
    for exponent in (6, 12, 24):
        prime2_generator = M.add(prime2_generator, M.pow(z, exponent))
    assert M.norm(prime2_generator) == 64
    prime7_generator = M.add(M.add(M.one, z), M.scale(M.pow(z, 3), -1))
    assert M.norm(prime7_generator) == 7
    assert M.norm(alpha0) == 43
    small_prime_norms = prime_norms_below_105()

    raw_manifest = "".join(serialize(row) for row in manifest)
    expected_digest = "fc87260544e8d8418e7bc6fe3ca80ef7d9715a82fcdaea9b2a9716016e00004e"
    assert sha256(raw_manifest.encode()).hexdigest() == expected_digest
    retained_manifest = Path(__file__).with_name("d3_matrix_manifest_v60.jsonl")
    assert retained_manifest.read_bytes() == raw_manifest.encode()
    summary = {"status": "PASS_EXACT_D3_FOUNDATIONS_AND_1280_CLASS_REPRESENTATIVES",
               "scope": "1280 trace-labelled SL2 matrix classes; no geometric point count",
               "geometric_search_performed": False,
               "arithmetic_d3_indices": d3, "selected_nonmaximal_indices": selected,
               "conductors": {k: expected[k] for k in selected}, "counts_per_trace": counts,
               "matrix_count": len(manifest), "positive_sector_count": 40,
               "local_unit_quotients": locals_,
               "sign_generator_coefficients": generators[:5],
               "sign_generator_signs": [sign_vector(a) for a in generators[:5]],
               "positive_nonsquare_coefficients": positive,
               "prime2_generator_in_z_basis": prime2_generator,
               "prime7_generator_in_z_basis": prime7_generator,
               "prime43_generator_in_z_basis": alpha0,
               "prime_ideal_norms_below_minkowski_bound": small_prime_norms,
               "discriminants": {"L": 453789, "M": 453789**2},
               "implementation_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
               "manifest_sha256": sha256(raw_manifest.encode()).hexdigest(),
               "root_brackets": [[str(a),str(b)] for a,b in INTERVALS]}
    if output_directory:
        output_directory.mkdir(parents=True, exist_ok=True)
        (output_directory/"matrix_manifest_2026_09_18.jsonl").write_text(raw_manifest)
        (output_directory/"foundation_output_2026_09_18.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path)
    arguments = parser.parse_args()
    main(arguments.output_directory)
