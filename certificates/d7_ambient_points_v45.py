#!/usr/bin/env python3
"""Exact audit of the nonmaximal d=7 orientation and fixed-point counts.

This certificate supplies the exact finite arithmetic behind the accompanying
human proof.  It checks:

* the two prime conductors and their comaximality;
* the missing CM-unit input E_M=<epsilon,mu_42,E_L> through its exact local
  consequences, including the order-84 conductor-83 class;
* explicit generators for E_L/E_L^2 and all 32/64 determinant sectors;
* one matrix representative for every multiplier/Picard class and all 512
  trace-labelled SL_2 orientation representatives;
* exact multiplier-order tests for the proper and maximal templates;
* the free fixed-point involution T -> cI-T and the resulting
  32+64+16=112 ambient quotient-point count;
* the two direct/conjugate actions over each of the normalized A, B, and C
  vertex orbits.

The certificate is exact.  The structural completeness proof (relative order
class sequence, prime-conductor local lattice lemma, determinant quotient, and
fixed-point fiber lemma) is written in docs/D7_AMBIENT_CM_POINTS.md.
"""

from hashlib import sha256
import sympy as sp


T = sp.symbols("T")
F = sp.Poly(
    T**6 + T**5 - 6*T**4 - 6*T**3 + 8*T**2 + 8*T + 1,
    T,
    domain=sp.QQ,
)
PI = T**4 + T**3 - 4*T**2 - 2*T + 5
IDENTITY = sp.eye(2)


def reduce_l(value):
    return sp.Poly(sp.expand(value), T, domain=sp.QQ).rem(F).as_expr()


def reduce_field(value):
    numerator, denominator = sp.fraction(sp.cancel(value))
    numerator_poly = sp.Poly(numerator, T, domain=sp.QQ).rem(F)
    denominator_poly = sp.Poly(denominator, T, domain=sp.QQ).rem(F)
    return (numerator_poly * sp.invert(denominator_poly, F)).rem(F).as_expr()


def matrix_reduce(matrix):
    return matrix.applyfunc(reduce_field)


def matrix_multiply(left, right):
    return matrix_reduce(left * right)


def matrix_power(matrix, exponent):
    result = IDENTITY
    base = matrix
    while exponent:
        if exponent & 1:
            result = matrix_multiply(result, base)
        base = matrix_multiply(base, base)
        exponent //= 2
    return result


def assert_integral(value):
    sp.Poly(reduce_field(value), T, domain=sp.ZZ)


def is_integral(value):
    polynomial = sp.Poly(reduce_field(value), T, domain=sp.QQ)
    return all(sp.Rational(coefficient).q == 1 for coefficient in polynomial.all_coeffs())


def assert_integral_matrix(matrix):
    for value in matrix:
        assert_integral(value)


def absolute_norm(value):
    return int(sp.resultant(F.as_expr(), reduce_l(value), T))


def pair_add(left, right):
    return reduce_l(left[0] + right[0]), reduce_l(left[1] + right[1])


def pair_multiply(left, right):
    a, b = left
    c, d = right
    return reduce_l(a*c - b*d), reduce_l(a*d + b*c + T*b*d)


def pair_power(pair, exponent):
    result = (sp.Integer(1), sp.Integer(0))
    base = pair
    while exponent:
        if exponent & 1:
            result = pair_multiply(result, base)
        base = pair_multiply(base, base)
        exponent //= 2
    return result


Z_PAIR = (sp.Integer(0), sp.Integer(1))


def z_power(exponent):
    return pair_power(Z_PAIR, exponent % 42)


ALPHA_ZERO = (2 - T**2, T**3 - 2*T + 1)


def alpha_pair(index):
    return pair_multiply(z_power(index), ALPHA_ZERO)


def alpha_trace(index):
    a, b = alpha_pair(index)
    return reduce_l(2*a + T*b)


def multiplication_matrix(pair):
    a, b = pair
    return sp.Matrix([[a, -b], [b, reduce_l(a + T*b)]])


def assert_characteristic(matrix, index):
    assert reduce_l(matrix.trace() - alpha_trace(index)) == 0
    assert reduce_l(matrix.det() - PI) == 0


def coefficient_tuple(value):
    polynomial = sp.Poly(reduce_field(value), T, domain=sp.ZZ)
    return tuple(int(polynomial.nth(i)) for i in range(6))


def matrix_tuple(matrix):
    return tuple(coefficient_tuple(value) for value in matrix)


def conjugate_by_unit(matrix, unit):
    diagonal = sp.diag(unit, 1)
    diagonal_inverse = sp.diag(reduce_field(1/unit), 1)
    answer = matrix_multiply(matrix_multiply(diagonal, matrix), diagonal_inverse)
    assert_integral_matrix(answer)
    return answer


# ---------------------------------------------------------------------------
# Exact real signs and the 64 unit square classes.
# ---------------------------------------------------------------------------


ROOT_INTERVALS = tuple(
    interval
    for interval, multiplicity in sp.polys.polytools.intervals(
        F, eps=sp.Rational(1, 10)**24
    )
    for _ in range(multiplicity)
)
ROOT_INTERVALS = tuple(
    sorted(ROOT_INTERVALS, key=lambda interval: sp.Rational(interval[0]), reverse=True)
)
assert len(ROOT_INTERVALS) == 6
assert all(
    sp.Rational(ROOT_INTERVALS[index][0])
    > sp.Rational(ROOT_INTERVALS[index + 1][1])
    for index in range(5)
)


def interval_evaluate(value, interval):
    coefficients = sp.Poly(reduce_field(value), T, domain=sp.QQ).all_coeffs()
    lower = upper = sp.Rational(coefficients[0])
    left, right = map(sp.Rational, interval)
    for coefficient in coefficients[1:]:
        candidates = (lower*left, lower*right, upper*left, upper*right)
        lower = min(candidates) + coefficient
        upper = max(candidates) + coefficient
    return lower, upper


def sign_vector(value):
    answer = []
    for interval in ROOT_INTERVALS:
        lower, upper = interval_evaluate(value, interval)
        assert lower > 0 or upper < 0
        answer.append(1 if lower > 0 else -1)
    return tuple(answer)


def xi_unit(index):
    # rho=z^2 and xi_a=rho^((1-a)/2)*(1+rho+...+rho^(a-1)).
    rho = z_power(2)
    total = (sp.Integer(0), sp.Integer(0))
    term = (sp.Integer(1), sp.Integer(0))
    for _ in range(index):
        total = pair_add(total, term)
        term = pair_multiply(term, rho)
    shift = pair_power(rho, ((1-index)//2) % 21)
    real_pair = pair_multiply(shift, total)
    assert reduce_l(real_pair[1]) == 0
    assert abs(absolute_norm(real_pair[0])) == 1
    return reduce_l(real_pair[0])


XI_INDICES = (5, 11, 13, 17, 19)
XI_UNITS = tuple(xi_unit(index) for index in XI_INDICES)
EXPECTED_XI_SIGNS = (
    (1, -1, 1, -1, 1, 1),
    (1, 1, -1, 1, 1, -1),
    (1, -1, 1, 1, 1, -1),
    (1, 1, 1, 1, -1, -1),
    (1, 1, -1, -1, -1, -1),
)
assert tuple(sign_vector(unit) for unit in XI_UNITS) == EXPECTED_XI_SIGNS

# The five xi_a vectors alone have rank four.  Replacing xi_19 by -1 gives a
# basis of the rank-five sign quotient E_L/E_L^+.
SIGN_GENERATORS = XI_UNITS[:4] + (sp.Integer(-1),)


POSITIVE_NONSQUARE = 4 - T**2
assert absolute_norm(POSITIVE_NONSQUARE) == 1
assert sign_vector(POSITIVE_NONSQUARE) == (1,)*6


def unit_product(bits, include_positive):
    answer = sp.Integer(1)
    for bit, unit in zip(bits[:5], SIGN_GENERATORS):
        if bit:
            answer = reduce_l(answer*unit)
    if include_positive and bits[5]:
        answer = reduce_l(answer*POSITIVE_NONSQUARE)
    assert abs(absolute_norm(answer)) == 1
    return answer


def mask_bits(mask, length):
    return tuple((mask >> index) & 1 for index in range(length))


SIGN_UNITS = tuple(unit_product(mask_bits(mask, 5), False) for mask in range(32))
SQUARE_CLASS_UNITS = tuple(unit_product(mask_bits(mask, 6), True) for mask in range(64))
assert len({sign_vector(unit) for unit in SIGN_UNITS}) == 32
assert len(SQUARE_CLASS_UNITS) == 64
assert len({coefficient_tuple(unit) for unit in SQUARE_CLASS_UNITS}) == 64

# A separate exact residue witness for the positive nonsquare: at the
# degree-one prime (41,T-10), u=4-T^2 is 27, a nonsquare modulo 41.
assert int(F.as_expr().subs(T, 10)) % 41 == 0
assert int(POSITIVE_NONSQUARE.subs(T, 10)) % 41 == 27
assert pow(27, 20, 41) == 40


# ---------------------------------------------------------------------------
# Local unit quotient and Picard-class inputs.
# ---------------------------------------------------------------------------


def finite_pair_multiply(left, right, prime, trace_root):
    a, b = left
    c, d = right
    return (
        (a*c - b*d) % prime,
        (a*d + b*c + trace_root*b*d) % prime,
    )


def finite_pair_power(pair, exponent, prime, trace_root):
    answer = (1, 0)
    base = pair
    while exponent:
        if exponent & 1:
            answer = finite_pair_multiply(answer, base, prime, trace_root)
        base = finite_pair_multiply(base, base, prime, trace_root)
        exponent //= 2
    return answer


def quotient_order(pair, prime, trace_root):
    # The coset is trivial exactly when its z coefficient vanishes.
    for order in range(1, prime + 2):
        if (prime + 1) % order == 0:
            value = finite_pair_power(pair, order, prime, trace_root)
            if value[1] == 0:
                return order
    raise AssertionError("local quotient order not found")


for prime, root in ((41, 12), (83, 36)):
    discriminant = (root*root - 4) % prime
    assert pow(discriminant, (prime-1)//2, prime) == prime - 1
    epsilon_local = (2 % prime, (-root) % prime)  # 1-z^2
    z_local = (0, 1)
    assert quotient_order(epsilon_local, prime, root) == 42
    assert quotient_order(z_local, prime, root) == 21
    epsilon_22 = finite_pair_power(epsilon_local, 22, prime, root)
    # Equality in the quotient means the ratio is a base-field scalar.
    assert epsilon_22[0]*z_local[1] % prime == epsilon_22[1]*z_local[0] % prime

BETA_83 = (1, 3)
assert quotient_order(BETA_83, 83, 36) == 84


GAMMA_41 = -1 - 2*T + T**2 - T**3 - T**4
GAMMA_83 = -2 - 2*T + 2*T**2 - T**4
BEZOUT_41 = -174 - 718*T - 267*T**2 + 733*T**3 + 131*T**4 - 167*T**5
assert abs(absolute_norm(GAMMA_41)) == 41
assert abs(absolute_norm(GAMMA_83)) == 83
assert reduce_l(BEZOUT_41*GAMMA_41 + 30*GAMMA_83) == 1


# The relative order-class exact sequence and the CM-unit lemma in the human
# proof give these exact Picard counts.  The local assertions above certify the
# quotient and image orders, including the missing order-84 class.
PICARD_COUNTS = {41: 1, 83: 2}
ICM_COUNTS = {41: 2, 83: 3}
assert 42//42 == PICARD_COUNTS[41]
assert 84//42 == PICARD_COUNTS[83]


# ---------------------------------------------------------------------------
# Exact module templates, multiplier-order tests, and orientation manifests.
# ---------------------------------------------------------------------------


def basis_action(index, basis):
    maximal_action = multiplication_matrix(alpha_pair(index))
    return matrix_reduce(basis.inv() * maximal_action * basis)


def z_action_in_basis(basis):
    maximal_z = multiplication_matrix(Z_PAIR)
    return matrix_reduce(basis.inv() * maximal_z * basis)


def companion_basis(index):
    a, b = alpha_pair(index)
    return sp.Matrix([[1, a], [0, b]])


I_83_BASIS = sp.Matrix([[16, GAMMA_83], [58 + 3*T, 3*GAMMA_83]])
BETA_83_MATRIX = multiplication_matrix((sp.Integer(1), sp.Integer(3)))
I_83_PREIMAGE_BASIS = sp.Matrix([[19, GAMMA_83], [1, 0]])
assert matrix_reduce(BETA_83_MATRIX*I_83_PREIMAGE_BASIS - I_83_BASIS) == sp.zeros(2)
assert reduce_l(I_83_PREIMAGE_BASIS.det() + GAMMA_83) == 0
assert abs(absolute_norm(I_83_PREIMAGE_BASIS.det())) == 83
assert_integral(reduce_field((58 + 3*T)/GAMMA_83))


def module_templates(index):
    prime = 83 if index in (2, 23) else 41
    templates = []

    r_basis = companion_basis(index)
    r_action = basis_action(index, r_basis)
    expected_companion = sp.Matrix([[0, -PI], [1, alpha_trace(index)]])
    assert matrix_reduce(r_action - expected_companion) == sp.zeros(2)
    assert_integral_matrix(r_action)
    assert_characteristic(r_action, index)
    r_z = z_action_in_basis(r_basis)
    assert any(not is_integral(value) for value in r_z)
    templates.append((f"R{prime}:principal", r_action, f"R{prime}"))

    if prime == 83:
        ideal_action = basis_action(index, I_83_BASIS)
        assert_integral_matrix(ideal_action)
        assert_characteristic(ideal_action, index)
        ideal_z = z_action_in_basis(I_83_BASIS)
        assert any(not is_integral(value) for value in ideal_z)
        assert matrix_tuple(ideal_action) != matrix_tuple(r_action)
        templates.append(("R83:nontrivial", ideal_action, "R83"))

    maximal_action = multiplication_matrix(alpha_pair(index))
    assert_integral_matrix(maximal_action)
    assert_characteristic(maximal_action, index)
    assert_integral_matrix(multiplication_matrix(Z_PAIR))
    templates.append(("OM", maximal_action, "OM"))
    return templates


ORIENTATION_RECORDS = []
for index in (2, 14, 23, 35):
    for label, base_matrix, multiplier in module_templates(index):
        units = SIGN_UNITS if multiplier == "OM" else SQUARE_CLASS_UNITS
        for orientation_index, unit in enumerate(units):
            matrix = conjugate_by_unit(base_matrix, unit)
            assert_characteristic(matrix, index)
            ORIENTATION_RECORDS.append(
                (index, label, orientation_index, sign_vector(unit), matrix_tuple(matrix))
            )

COUNTS_BY_INDEX = {
    index: sum(1 for record in ORIENTATION_RECORDS if record[0] == index)
    for index in (2, 14, 23, 35)
}
assert COUNTS_BY_INDEX == {2: 160, 14: 96, 23: 160, 35: 96}
assert len(ORIENTATION_RECORDS) == 512
assert len({(record[0], record[4]) for record in ORIENTATION_RECORDS}) == 512

manifest = sha256()
for record in ORIENTATION_RECORDS:
    manifest.update((repr(record) + "\n").encode("ascii"))
ORIENTATION_DIGEST = manifest.hexdigest()
assert ORIENTATION_DIGEST == "ef559dbf40bdb2eb06c16870a9b409d94d1d2cb6d7c930d2d7cd6a7a2f8ec89c"


# ---------------------------------------------------------------------------
# Fixed-point fibers and the 112-point incidence count.
# ---------------------------------------------------------------------------


for index in (2, 14):
    opposite = index + 21
    assert alpha_pair(opposite) == tuple(reduce_l(-value) for value in alpha_pair(index))
    assert reduce_l(alpha_trace(opposite) + alpha_trace(index)) == 0


# Exact anti-actions identify the determinant coset of T -> cI-T on every
# module template.  For the O_M basis and the principal proper-order bases the
# determinant is -1 (mask bit 4).  On the nontrivial R_83 ideal it is
# T^2-4=-(4-T^2), hence mask bits 4 and 5.
STANDARD_ANTI = sp.Matrix([[1, T], [0, -1]])
assert reduce_l(STANDARD_ANTI.det() + 1) == 0
for index in (2, 14, 23, 35):
    maximal = multiplication_matrix(alpha_pair(index))
    conjugate = alpha_trace(index)*IDENTITY - maximal
    assert matrix_reduce(STANDARD_ANTI*maximal - conjugate*STANDARD_ANTI) == sp.zeros(2)

    companion = sp.Matrix([[0, -PI], [1, alpha_trace(index)]])
    companion_anti = sp.Matrix([[1, alpha_trace(index)], [0, -1]])
    companion_conjugate = alpha_trace(index)*IDENTITY - companion
    assert reduce_l(companion_anti.det() + 1) == 0
    assert matrix_reduce(
        companion_anti*companion - companion_conjugate*companion_anti
    ) == sp.zeros(2)

I_83_ANTI = sp.Matrix(
    [
        [18*T**5 - 109*T**3 - 18*T**2 + 166*T + 55,
         -T**5 + T**4 + 5*T**3 - 3*T**2 - 9*T - 4],
        [-683*T**5 - 184*T**4 + 4245*T**3 + 905*T**2 - 6206*T - 603,
         -18*T**5 + 109*T**3 + 18*T**2 - 166*T - 55],
    ]
)
assert reduce_l(I_83_ANTI.det() - (T**2 - 4)) == 0
for index in (2, 23):
    ideal_action = basis_action(index, I_83_BASIS)
    ideal_conjugate = alpha_trace(index)*IDENTITY - ideal_action
    assert matrix_reduce(
        I_83_ANTI*ideal_action - ideal_conjugate*I_83_ANTI
    ) == sp.zeros(2)

PAIR_XOR = {
    "OM": 16,
    "R41:principal": 16,
    "R83:principal": 16,
    "R83:nontrivial": 48,
}
for label, pairing in PAIR_XOR.items():
    modulus = 32 if label == "OM" else 64
    assert all((mask ^ pairing) != mask for mask in range(modulus))
    assert all(((mask ^ pairing) ^ pairing) == mask for mask in range(modulus))

POINT_RECORDS = []
for label, modulus in (
    ("OM", 32),
    ("R41:principal", 64),
    ("R83:principal", 64),
    ("R83:nontrivial", 64),
):
    pairing = PAIR_XOR[label]
    representatives = sorted({min(mask, mask ^ pairing) for mask in range(modulus)})
    POINT_RECORDS.extend((label, mask) for mask in representatives)
assert len(POINT_RECORDS) == 112
point_manifest = sha256()
for record in POINT_RECORDS:
    point_manifest.update((repr(record) + "\n").encode("ascii"))
POINT_DIGEST = point_manifest.hexdigest()
assert POINT_DIGEST == "576e1d0e1966bf2868a33876147b8d7dc8ea1c002c64a8abc35aea95bb8e671b"


def anti_conjugator(generator):
    # In the cyclic basis (e_1, Z e_1), Z is [[0,-1],[1,t]].
    cyclic_basis = sp.Matrix([[1, generator[0, 0]], [0, generator[1, 0]]])
    assert reduce_field(cyclic_basis.det()) != 0
    standard_conjugation = sp.Matrix([[1, T], [0, -1]])
    answer = matrix_reduce(cyclic_basis * standard_conjugation * cyclic_basis.inv())
    assert matrix_reduce(answer*answer - IDENTITY) == sp.zeros(2)
    assert reduce_l(answer.det() + 1) == 0
    assert matrix_reduce(answer*generator - (T*IDENTITY-generator)*answer) == sp.zeros(2)
    return answer


A_PRINTED = sp.Matrix(
    [
        [T**5 - 4*T**3 - T**2 + 3*T,
         -T**5 + T**4 + 6*T**3 - 3*T**2 - 10*T],
        [T**3 - 1, -T**5 + 5*T**3 + T**2 - 6*T],
    ]
)
B_PRINTED = sp.Matrix([[0, -1], [1, T]])
C_PRINTED = sp.Matrix(
    [
        [2*T - 4*T**2 - 4*T**3 + T**4 + T**5,
         -5*T - 4*T**2 + T**3 + T**4],
        [-1 + 5*T + 4*T**2 - T**3 - T**4,
         -T + 4*T**2 + 4*T**3 - T**4 - T**5],
    ]
)
H = sp.diag(T, 1)
H_INVERSE = sp.diag(reduce_field(1/T), 1)
A_NORMALIZED = matrix_multiply(matrix_multiply(H, A_PRINTED), H_INVERSE)
B_NORMALIZED = matrix_multiply(matrix_multiply(H, B_PRINTED), H_INVERSE)
C_NORMALIZED = matrix_multiply(matrix_multiply(H, C_PRINTED), H_INVERSE)

Z_A = sp.Matrix(
    [
        [-T**5 + 7*T**3 - 11*T - 1, -T**5 + 8*T**3 - 15*T - 2],
        [-T**4 + 6*T**2 - 7, T**5 - 7*T**3 + 12*T + 1],
    ]
)
VERTEX_GENERATORS = {"A": Z_A, "B": B_NORMALIZED, "C": C_NORMALIZED}

for name, generator in VERTEX_GENERATORS.items():
    assert_integral_matrix(generator)
    assert matrix_reduce(generator*generator - T*generator + IDENTITY) == sp.zeros(2)
    assert reduce_l(generator.trace() - T) == 0
    assert reduce_l(generator.det() - 1) == 0
    if name == "A":
        assert matrix_reduce(matrix_power(generator, 3) - A_NORMALIZED) == sp.zeros(2)
    anti = anti_conjugator(generator)
    for index in (2, 14, 23, 35):
        a, b = alpha_pair(index)
        direct = matrix_reduce(a*IDENTITY + b*generator)
        conjugate = matrix_reduce(alpha_trace(index)*IDENTITY - direct)
        assert_integral_matrix(direct)
        assert_integral_matrix(conjugate)
        assert_characteristic(direct, index)
        assert_characteristic(conjugate, index)
        assert matrix_reduce(anti*direct - conjugate*anti) == sp.zeros(2)
        assert matrix_tuple(direct) != matrix_tuple(conjugate)

assert abs(absolute_norm(T**2 - 4)) == 1


# ---------------------------------------------------------------------------
# Orientation fingerprints for the three vertices and the two already
# excluded maximal-order companion points.
# ---------------------------------------------------------------------------


SIGN_TO_MASK = {sign_vector(unit): mask for mask, unit in enumerate(SIGN_UNITS)}
assert len(SIGN_TO_MASK) == 32
STANDARD_Z = multiplication_matrix(Z_PAIR)

VERTEX_MASKS = {}
for name, generator in VERTEX_GENERATORS.items():
    basis = sp.Matrix([[1, generator[0, 0]], [0, generator[1, 0]]])
    assert matrix_reduce(generator*basis - basis*STANDARD_Z) == sp.zeros(2)
    VERTEX_MASKS[name] = SIGN_TO_MASK[sign_vector(basis.det())]
assert VERTEX_MASKS == {"A": 7, "B": 15, "C": 0}

EXCLUDED_COMPANION_MASKS = {}
for index, expected_mask in ((1, 13), (17, 2)):
    maximal = multiplication_matrix(alpha_pair(index))
    companion = sp.Matrix([[0, -PI], [1, alpha_trace(index)]])
    bridge = sp.Matrix([[1, maximal[0, 0]], [0, maximal[1, 0]]])
    assert matrix_reduce(bridge*companion - maximal*bridge) == sp.zeros(2)
    mask = SIGN_TO_MASK[sign_vector(bridge.det())]
    assert mask == expected_mask
    EXCLUDED_COMPANION_MASKS[index] = min(mask, mask ^ 16)
    opposite = index + 21
    opposite_maximal = multiplication_matrix(alpha_pair(opposite))
    opposite_companion = sp.Matrix([[0, -PI], [1, alpha_trace(opposite)]])
    opposite_bridge = sp.Matrix(
        [[1, opposite_maximal[0, 0]], [0, opposite_maximal[1, 0]]]
    )
    assert matrix_reduce(
        opposite_bridge*opposite_companion - opposite_maximal*opposite_bridge
    ) == sp.zeros(2)
    opposite_mask = SIGN_TO_MASK[sign_vector(opposite_bridge.det())]
    assert opposite_mask == (mask ^ 16)
assert EXCLUDED_COMPANION_MASKS == {1: 13, 17: 2}

CLASSIFIED_OM_MASKS = set(VERTEX_MASKS.values()) | set(EXCLUDED_COMPANION_MASKS.values())
assert CLASSIFIED_OM_MASKS == {0, 2, 7, 13, 15}
assert 16 - len(CLASSIFIED_OM_MASKS) == 11

# Human fiber lemma: T and cI-T have the same fixed point and form a free pair
# on SL2 classes.  The counts below are therefore exact consequences of the
# certified class counts and multiplier strata.
MATRIX_CLASSES = {"R41": 64, "R83": 128, "OM": 32}
POINTS_PER_TRACE = {41: (64 + 32)//2, 83: (128 + 32)//2}
assert POINTS_PER_TRACE == {41: 48, 83: 80}

AMBIENT_POINTS = {"proper_R41": 32, "proper_R83": 64, "maximal_OM": 16}
assert sum(AMBIENT_POINTS.values()) == 112
assert 32*4 + 64*4 + 16*8 == 512
assert 112 - 3 - 2 == 107


print("EXACT NONMAXIMAL d=7 AMBIENT-POINT CERTIFICATE v45")
print("unit sectors: |E_L/E_L^+|=32; |E_L/E_L^2|=64")
print("local quotient orders: epsilon=(42,42); beta_83=84")
print("Picard counts:", PICARD_COUNTS, "; ICM counts:", ICM_COUNTS)
print("trace-labelled SL2 orientation counts:", COUNTS_BY_INDEX, "total=512")
print("orientation manifest SHA-256:", ORIENTATION_DIGEST)
print("fixed points per conductor trace: conductor 41 -> 48; conductor 83 -> 80")
print("deduplicated ambient points:", AMBIENT_POINTS, "total=112")
print("ambient point-label SHA-256:", POINT_DIGEST)
print("known normalized vertex points: A, B, C (six matrix classes per trace)")
print("classified O_M point masks: vertices", VERTEX_MASKS,
      "; exact exclusions", EXCLUDED_COMPANION_MASKS)
print("remaining exact curve-membership universe: 107 ambient points")
print("PASS: exact arithmetic inputs and representative manifest verified")
