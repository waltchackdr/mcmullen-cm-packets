#!/usr/bin/env python3
"""Independent exact reconstruction of the nonmaximal d=7 ambient packet.

This verifier uses only the standard-library quotient-field and Sturm core in
``quadrilateral_k17_stdlib_v41.py``.  It does not import the primary v45
certificate.  It reconstructs the four orders, the proper module templates,
all 512 determinant-orientation representatives, the fixed-point root
involutions, and the resulting 112 ambient Hilbert-modular CM points.

The script deliberately does *not* claim that the 107 points not already
classified by earlier vertex/exclusion certificates lie on, or off, the
McMullen curve.  That is the next geometric orbit problem.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256

from quadrilateral_k17_stdlib_v41 import (
    EMBEDDINGS,
    IDENTITY,
    ONE,
    T,
    ZERO,
    LElement,
    Matrix,
    all_entries_integral,
    column_intertwiner,
    companion,
    conjugate_by,
    elt,
    embedded_roots,
    field_discriminant,
    madd,
    mdet,
    minverse,
    mmul,
    mscale,
    msub,
    mtrace,
    matrix,
    sign_vector,
)


Pair = tuple[LElement, LElement]
ZERO_MATRIX = matrix(0, 0, 0, 0)
INDICES = (2, 14, 23, 35)
PI = T**4 + T**3 - 4 * T**2 - 2 * T + 5
BASE_DISCRIMINANT = T**2 - 4


def pair_multiply(left: Pair, right: Pair) -> Pair:
    """Multiply A+B*z pairs using z^2-t*z+1=0."""
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c + T * b * d


def z_power(exponent: int) -> Pair:
    answer = (ONE, ZERO)
    base = (ZERO, ONE)
    for _ in range(exponent % 42):
        answer = pair_multiply(answer, base)
    return answer


def alpha_pair(index: int) -> Pair:
    z4 = z_power(4)
    alpha0 = (ONE + z4[0], ONE + z4[1])
    return pair_multiply(z_power(index), alpha0)


def pair_trace(value: Pair) -> LElement:
    a, b = value
    return 2 * a + T * b


def pair_norm(value: Pair) -> LElement:
    a, b = value
    return a * a + T * a * b + b * b


def multiplication_matrix(value: Pair) -> Matrix:
    a, b = value
    return matrix(a, -b, b, a + T * b)


def diagonal(value: LElement) -> Matrix:
    return matrix(value, 0, 0, 1)


def twist(value: Matrix, unit: LElement) -> Matrix:
    return conjugate_by(diagonal(unit), value)


def f2_rank(rows: tuple[tuple[int, ...], ...]) -> int:
    values = [sum(bit << j for j, bit in enumerate(row)) for row in rows]
    rank = 0
    for bit in range(6):
        pivot = next((i for i in range(rank, len(values)) if (values[i] >> bit) & 1), None)
        if pivot is None:
            continue
        values[rank], values[pivot] = values[pivot], values[rank]
        for i in range(len(values)):
            if i != rank and ((values[i] >> bit) & 1):
                values[i] ^= values[rank]
        rank += 1
    return rank


def sign_bits(signs: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(0 if sign > 0 else 1 for sign in signs)


def finite_multiply(p: int, trace: int, left: tuple[int, int],
                    right: tuple[int, int]) -> tuple[int, int]:
    a, b = left
    c, d = right
    return ((a * c - b * d) % p, (a * d + b * c + trace * b * d) % p)


def finite_power(p: int, trace: int, value: tuple[int, int],
                 exponent: int) -> tuple[int, int]:
    answer = (1, 0)
    base = value
    n = exponent
    while n:
        if n & 1:
            answer = finite_multiply(p, trace, answer, base)
        base = finite_multiply(p, trace, base, base)
        n >>= 1
    return answer


def local_quotient_order(p: int, trace: int, value: tuple[int, int],
                         expected: int) -> tuple[int, int]:
    """Return value^(p-1), after certifying its exact quotient order."""
    quotient_value = finite_power(p, trace, value, p - 1)
    assert finite_power(p, trace, quotient_value, expected) == (1, 0)
    for prime in (2, 3, 7):
        if expected % prime == 0:
            assert finite_power(p, trace, quotient_value, expected // prime) != (1, 0)
    return quotient_value


def polynomial_modulus_value(x: int, p: int) -> int:
    # f in ascending coefficient order.
    coefficients = (1, 8, 8, -6, -6, 1, 1)
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * x + coefficient) % p
    return value


def element_modulus_value(value: LElement, x: int, p: int) -> int:
    answer = 0
    for coefficient in reversed(value.coefficients):
        assert coefficient.denominator % p != 0
        residue = coefficient.numerator * pow(coefficient.denominator, -1, p)
        answer = (answer * x + residue) % p
    return answer


def coefficient_string(value: LElement) -> str:
    # Independent format: high-to-low coefficients, explicit numerator/denominator.
    return ".".join(
        f"{coefficient.numerator}/{coefficient.denominator}"
        for coefficient in reversed(value.coefficients)
    )


def manifest_record(tag: str, index: int, mask: int, value: Matrix) -> bytes:
    # Column-major matrix order intentionally differs from the primary path.
    entries = (value[0][0], value[1][0], value[0][1], value[1][1])
    text = f"{tag}:{index:02d}:{mask:02x}:" + "/".join(
        coefficient_string(entry) for entry in entries
    ) + "\n"
    return text.encode("ascii")


def root_conjugate(value: Matrix) -> Matrix:
    return msub(mscale(mtrace(value), IDENTITY), value)


def exact_fixed_invariants(value: Matrix) -> tuple[LElement, LElement]:
    """Return x and y^2 for the upper fixed point of a determinant-pi matrix."""
    a, c, d = value[0][0], value[1][0], value[1][1]
    return (a - d) / (2 * c), (4 * mdet(value) - mtrace(value) ** 2) / (4 * c**2)


def main() -> None:
    roots = embedded_roots()
    assert field_discriminant() == 453789
    assert PI.norm() == 43
    assert sign_vector(PI, roots) == (1,) * 6
    assert BASE_DISCRIMINANT.norm() == 1
    assert sign_vector(BASE_DISCRIMINANT, roots) == (-1,) * 6

    # Four real cyclotomic units, -1, and the missing totally-positive
    # square-class direction.  The first five give all 32 sign classes.
    xi5 = T**4 - 3 * T**2 + 1
    xi11 = T**3 + T**2 - 2 * T - 1
    xi13 = T**5 - 4 * T**3 + 3 * T + 1
    xi17 = T**3 - 2 * T
    minus_one = -ONE
    u = 4 - T**2
    sign_units = (xi5, xi11, xi13, xi17, minus_one)
    unit_basis = sign_units + (u,)
    expected_signs = (
        (1, -1, 1, -1, 1, 1),
        (1, 1, -1, 1, 1, -1),
        (1, -1, 1, 1, 1, -1),
        (1, 1, 1, 1, -1, -1),
        (-1, -1, -1, -1, -1, -1),
    )
    assert tuple(sign_vector(value, roots) for value in sign_units) == expected_signs
    assert f2_rank(tuple(sign_bits(value) for value in expected_signs)) == 5
    assert all(abs(value.norm()) == 1 for value in unit_basis)
    assert sign_vector(u, roots) == (1,) * 6
    # A global square remains a square at every odd residue prime.  At the
    # degree-one prime (41,t-10), u reduces to the nonsquare 27.
    assert polynomial_modulus_value(10, 41) == 0
    assert (4 - 10 * 10) % 41 == 27
    assert pow(27, 20, 41) == 40

    deltas: list[LElement] = []
    for mask in range(64):
        value = ONE
        for bit, generator in enumerate(unit_basis):
            if (mask >> bit) & 1:
                value *= generator
        assert value.integral_power_basis() and abs(value.norm()) == 1
        deltas.append(value)
    assert len(set(deltas)) == 64
    assert len({sign_vector(value, roots) for value in deltas[:32]}) == 32

    # Exact local unit computations.  Both conductor primes are degree one
    # and inert in M/L.  epsilon=1-z^2=(2,-t) has quotient order 42.
    for p, residue_trace in ((41, 12), (83, 36)):
        assert polynomial_modulus_value(residue_trace, p) == 0
        assert pow((-7) % p, (p - 1) // 2, p) == p - 1
        epsilon = (2, (-residue_trace) % p)
        local_quotient_order(p, residue_trace, epsilon, 42)
        z_class = local_quotient_order(p, residue_trace, (0, 1), 21)
        epsilon_class = finite_power(p, residue_trace, epsilon, p - 1)
        assert z_class == finite_power(p, residue_trace, epsilon_class, 22)

    # beta supplies the nontrivial conductor-83 proper ideal class.
    beta_class = local_quotient_order(83, 36, (1, 3), 84)
    assert beta_class == (57, 38)
    assert finite_power(83, 36, beta_class, 42) == (82, 0)

    picard_counts = {41: 1, 83: 2}
    icm_counts = {41: 2, 83: 3}
    oriented_counts = {41: 96, 83: 160}
    assert picard_counts[41] == (41 + 1) // 42
    assert picard_counts[83] == (83 + 1) // 42
    assert icm_counts == {p: picard_counts[p] + 1 for p in picard_counts}
    assert oriented_counts == {
        41: 64 * picard_counts[41] + 32,
        83: 64 * picard_counts[83] + 32,
    }

    # Module templates in the independently reconstructed (1,z) basis.
    z_matrix = matrix(0, -1, 1, T)
    gamma83 = elt(-2, -2, 2, 0, -1)
    v_matrix = matrix(16, gamma83, elt(58, 3), 3 * gamma83)
    assert all_entries_integral(v_matrix)
    assert mdet(v_matrix).norm() == -30551387
    beta_matrix = multiplication_matrix((ONE, LElement(3)))
    beta_coordinates = mmul(minverse(beta_matrix), v_matrix)
    assert beta_coordinates == matrix(19, gamma83, 1, 0)
    assert all_entries_integral(beta_coordinates)
    # Thus V lies in beta*O_M.  Its relative index there is the conductor
    # prime: det(beta_coordinates)=-gamma83 has absolute norm 83.
    assert mdet(beta_coordinates) == -gamma83
    assert abs(mdet(beta_coordinates).norm()) == 83

    bases: dict[tuple[int, str], Matrix] = {}
    anti: dict[tuple[int, str], tuple[Matrix, int]] = {}
    order_discriminant_norms: dict[int, int] = {}

    for index in INDICES:
        alpha = alpha_pair(index)
        maximal_matrix = multiplication_matrix(alpha)
        trace = pair_trace(alpha)
        assert pair_norm(alpha) == PI
        assert (mtrace(maximal_matrix), mdet(maximal_matrix)) == (trace, PI)
        assert all_entries_integral(maximal_matrix)
        discriminant = trace**2 - 4 * PI
        assert discriminant == alpha[1] ** 2 * BASE_DISCRIMINANT
        order_discriminant_norms[index] = abs(int(discriminant.norm()))

        regular_matrix = companion(trace, PI)
        assert all_entries_integral(regular_matrix)
        assert (mtrace(regular_matrix), mdet(regular_matrix)) == (trace, PI)

        bases[index, "O"] = maximal_matrix
        bases[index, "R"] = regular_matrix

        # z is integral on the O_M lattice but not on the regular R lattice.
        q_regular = column_intertwiner(maximal_matrix)
        assert mmul(q_regular, regular_matrix) == mmul(maximal_matrix, q_regular)
        z_on_regular = conjugate_by(minverse(q_regular), z_matrix)
        assert not all_entries_integral(z_on_regular)
        if index in (2, 23):
            # The displayed V columns also lie in R_k.  Combined with the
            # exact conductor-83 index above, this is the local intersection
            # construction V = beta*O_M intersect R_k.
            assert all_entries_integral(mmul(minverse(q_regular), v_matrix))

        j_z = matrix(1, T, 0, -1)
        j_regular = matrix(1, trace, 0, -1)
        assert mdet(j_z) == -ONE and mdet(j_regular) == -ONE
        assert mmul(j_z, maximal_matrix) == mmul(root_conjugate(maximal_matrix), j_z)
        assert mmul(j_regular, regular_matrix) == mmul(root_conjugate(regular_matrix), j_regular)
        anti[index, "O"] = (j_z, 16)
        anti[index, "R"] = (j_regular, 16)

        if index in (2, 23):
            ideal_matrix = conjugate_by(minverse(v_matrix), maximal_matrix)
            assert all_entries_integral(ideal_matrix)
            assert (mtrace(ideal_matrix), mdet(ideal_matrix)) == (trace, PI)
            bases[index, "I"] = ideal_matrix

            # The multiplier cannot be O_M: z fails to preserve this lattice.
            z_on_ideal = conjugate_by(minverse(v_matrix), z_matrix)
            assert not all_entries_integral(z_on_ideal)
            assert z_on_ideal[1][0].coefficients[0] == Fraction(71454, 83)

            p = elt(55, 166, -18, -109, 0, 18)
            q = elt(-4, -9, -3, 5, 1, -1)
            r = elt(-603, -6206, 905, 4245, -184, -683)
            p_ideal = matrix(p, q, r, -p)
            assert all_entries_integral(p_ideal)
            assert mdet(p_ideal) == T**2 - 4 == -u
            assert mmul(p_ideal, ideal_matrix) == mmul(root_conjugate(ideal_matrix), p_ideal)
            anti[index, "I"] = (p_ideal, 48)

    assert order_discriminant_norms == {2: 83**2, 14: 41**2, 23: 83**2, 35: 41**2}
    for p, residue_trace, relevant_indices in (
        (83, 36, (2, 23)),
        (41, 12, (14, 35)),
    ):
        roots_mod_p = tuple(x for x in range(p) if polynomial_modulus_value(x, p) == 0)
        assert len(roots_mod_p) == 6
        for index in relevant_indices:
            discriminant = mtrace(bases[index, "O"]) ** 2 - 4 * PI
            vanishing_roots = tuple(
                x for x in roots_mod_p if element_modulus_value(discriminant, x, p) == 0
            )
            assert vanishing_roots == (residue_trace,)

    # The two traces for each conductor differ by alpha -> -alpha.  This is
    # the exact cross-trace point-set identification used in the 112 count.
    assert alpha_pair(23) == (-alpha_pair(2)[0], -alpha_pair(2)[1])
    assert alpha_pair(35) == (-alpha_pair(14)[0], -alpha_pair(14)[1])
    assert bases[23, "O"] == mscale(-ONE, bases[2, "O"])
    assert bases[35, "O"] == mscale(-ONE, bases[14, "O"])
    assert bases[23, "I"] == mscale(-ONE, bases[2, "I"])
    assert bases[23, "R"] == twist(mscale(-ONE, bases[2, "R"]), -ONE)
    assert bases[35, "R"] == twist(mscale(-ONE, bases[14, "R"]), -ONE)
    common_o_invariants = exact_fixed_invariants(bases[2, "O"])
    assert all(exact_fixed_invariants(bases[index, "O"]) == common_o_invariants for index in INDICES)

    # Complete deterministic orientation manifest and free root involutions.
    digest = sha256()
    record_count = 0
    per_index: dict[int, int] = {index: 0 for index in INDICES}
    invariant_positive_checks = 0
    for tag in ("I", "O", "R"):
        for index in INDICES:
            if (index, tag) not in bases:
                continue
            base = bases[index, tag]
            masks = range(32) if tag == "O" else range(64)
            anti_matrix, involution_mask = anti[index, tag]
            for mask in masks:
                value = twist(base, deltas[mask])
                assert all_entries_integral(value)
                assert (mtrace(value), mdet(value)) == (mtrace(base), PI)
                digest.update(manifest_record(tag, index, mask, value))
                record_count += 1
                per_index[index] += 1

                _, height_squared = exact_fixed_invariants(value)
                assert sign_vector(height_squared, roots) == (1,) * 6
                invariant_positive_checks += 6

                partner_mask = mask ^ involution_mask
                partner = twist(base, deltas[partner_mask])
                # When the I-mask already contains u, canonical xor removes
                # u rather than multiplying by it.  Scaling the anti-map by
                # u^-1 changes its determinant by u^-2 and gives the exact
                # canonical determinant, without changing its relation.
                effective_anti = anti_matrix
                if tag == "I" and ((mask >> 5) & 1):
                    effective_anti = mscale(u.inverse(), anti_matrix)
                # S sends the same fixed point, encoded by cI-value, to its
                # selected partner representative and has determinant one.
                s_matrix = mmul(
                    mmul(diagonal(deltas[partner_mask]), minverse(effective_anti)),
                    minverse(diagonal(deltas[mask])),
                )
                assert all_entries_integral(s_matrix) and mdet(s_matrix) == ONE
                assert conjugate_by(s_matrix, root_conjugate(value)) == partner

    assert record_count == 512
    assert per_index == {2: 160, 14: 96, 23: 160, 35: 96}
    assert invariant_positive_checks == 512 * 6
    manifest_digest = digest.hexdigest()
    expected_digest = "5212fa42df26b1ebdbad04bb4f8f2edba473d7a1ce733bbd23ebfd94fca7d82b"
    assert manifest_digest == expected_digest

    # Quotient-point count.  The root involutions above are free: an
    # orientation-preserving real Mobius map fixing the upper point cannot
    # exchange the two CM eigenvalues.  k and k+21 give the same point sets;
    # the O_M stratum is common to both conductor packets.
    ambient_points = {
        "O": 32 // 2,
        "R41": 64 // 2,
        "R83": 64 // 2,
        "I83": 64 // 2,
    }
    assert ambient_points == {"O": 16, "R41": 32, "R83": 32, "I83": 32}
    assert sum(ambient_points.values()) == 112
    point_digest_builder = sha256()
    point_record_count = 0
    for tag, modulus, pairing in (
        ("I83", 64, 48),
        ("O", 32, 16),
        ("R41", 64, 16),
        ("R83", 64, 16),
    ):
        representatives = sorted({min(mask, mask ^ pairing) for mask in range(modulus)})
        for mask in representatives:
            point_digest_builder.update(f"{tag}>{mask:02x}\n".encode("ascii"))
            point_record_count += 1
    point_digest = point_digest_builder.hexdigest()
    assert point_record_count == 112
    assert point_digest == "67a7f316ab84a67b34409db7e91e0af1598ecb6583e3634855c4fa7067ccb428"

    # Known O_M point masks.  Masks are reduced modulo xor 16 by choosing
    # the representative with the -1 bit zero.
    sign_to_mask = {sign_vector(deltas[mask], roots): mask for mask in range(32)}
    b1 = alpha_pair(1)[1]
    b17 = alpha_pair(17)[1]
    determinant_c = elt(13, 12, -7, -7, 1, 1)
    determinant_a = elt(-7, 0, 6, 0, -1)
    determinant_b = T.inverse()
    known_masks = {
        "C": sign_to_mask[sign_vector(determinant_c, roots)],
        "excluded_17_38": sign_to_mask[sign_vector(b17, roots)],
        "A": sign_to_mask[sign_vector(determinant_a, roots)],
        "excluded_1_22": sign_to_mask[sign_vector(b1, roots)],
        "B": sign_to_mask[sign_vector(determinant_b, roots)],
    }
    assert known_masks == {
        "C": 0,
        "excluded_17_38": 2,
        "A": 7,
        "excluded_1_22": 13,
        "B": 15,
    }
    assert alpha_pair(38)[1] == -b17
    assert alpha_pair(22)[1] == -b1
    assert sign_vector(determinant_a, roots) == (1, 1, -1, -1, 1, 1)
    assert sign_vector(determinant_b, roots) == (1, 1, -1, -1, -1, -1)
    assert sign_vector(determinant_c, roots) == (1, 1, 1, 1, 1, 1)

    classified_o_masks = set(known_masks.values())
    unresolved_o_masks = set(range(16)) - classified_o_masks
    unresolved_proper_points = ambient_points["R41"] + ambient_points["R83"] + ambient_points["I83"]
    assert len(unresolved_o_masks) == 11
    assert unresolved_proper_points == 96
    assert len(unresolved_o_masks) + unresolved_proper_points == 107

    print("INDEPENDENT d=7 AMBIENT-POINT STURM CERTIFICATE v45")
    print("unit sign rank=5; u positive nonsquare witness=(41,t-10,u=27)")
    print("local quotient orders: epsilon=42 at p=41,83; beta=84 at p=83")
    print("orientation counts:", per_index, "total=512")
    print("manifest digest (independent column-major format):", manifest_digest)
    print("free root involutions: O/R xor16; I83 xor48")
    print("ambient quotient points:", ambient_points, "total=112")
    print("ambient point digest (independent label format):", point_digest)
    print("known O masks:", known_masks)
    print("unresolved O masks:", tuple(sorted(unresolved_o_masks)))
    print("unresolved geometry: 11 O + 96 proper-order = 107")
    print("PASS: independent exact d=7 ambient arithmetic and point count")


if __name__ == "__main__":
    main()
