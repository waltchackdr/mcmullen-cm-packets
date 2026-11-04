#!/usr/bin/env python3
"""Exact twelve-class CM-type and Hodge reconstruction.

This extends ``quadrilateral_k17_stdlib_v41``.  The only non-stdlib import is
that independent module's fresh Fraction quotient ring, Sturm isolator, and
rational interval evaluator.  No repository certificate or experiment is
imported or executed.

Every decisive CM-type sign is certified by an interval excluding zero.
Hodge eigenlines are then enumerated by the CM cocharacter criterion, and
divisor products by disjoint unions of degree-two eigenline index sets.
"""

from __future__ import annotations

from itertools import combinations
from math import gcd

from quadrilateral_k17_stdlib_v41 import (
    EMBEDDINGS,
    ONE,
    T,
    embedded_roots,
    interval_eval,
    power_trace,
    sign_vector,
    sine_quotient,
)


MAXIMAL_33_INDICES = (1, 11, 16, 17, 19, 20, 22, 32, 37, 38, 40, 41)
GALOIS_GROUP = tuple(a for a in range(42) if gcd(a, 42) == 1)


def cm_type_for_index(k: int, roots: dict) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return the six certified signs and the selected upper-line CM type.

    For alpha_k=zeta^k+zeta^(k+1)+zeta^(k+4), write
    alpha_k=x_k+y_k*zeta.  At an embedding r in EMBEDDINGS,

        Im(alpha_k(zeta -> zeta^r)) = sigma_r(y_k) sin(r*pi/21).

    The sine factor is positive for all six printed r's, so the exact sign of
    sigma_r(y_k) selects r versus -r for the upper companion eigenline.
    """
    yk = sine_quotient(k) + sine_quotient(k + 1) + sine_quotient(k + 4)
    signs = sign_vector(yk, roots)
    phi = tuple(r if s > 0 else (-r) % 42 for r, s in zip(EMBEDDINGS, signs))
    return signs, phi


def translate_type(phi: frozenset[int], multiplier: int) -> frozenset[int]:
    return frozenset((multiplier * exponent) % 42 for exponent in phi)


def type_orbit(phi: frozenset[int]) -> frozenset[frozenset[int]]:
    return frozenset(translate_type(phi, multiplier) for multiplier in GALOIS_GROUP)


def type_stabilizer(phi: frozenset[int]) -> tuple[int, ...]:
    return tuple(multiplier for multiplier in GALOIS_GROUP if translate_type(phi, multiplier) == phi)


def hodge_eigenlines(phi: frozenset[int], codimension: int) -> set[frozenset[int]]:
    """Enumerate wedge eigenlines satisfying the CM cocharacter criterion."""
    conjugate_types = type_orbit(phi)
    answer: set[frozenset[int]] = set()
    for subset_tuple in combinations(GALOIS_GROUP, 2 * codimension):
        subset = frozenset(subset_tuple)
        if all(len(subset & conjugate_phi) == codimension for conjugate_phi in conjugate_types):
            answer.add(subset)
    return answer


def divisor_generated_eigenlines(phi: frozenset[int], codimension: int) -> set[frozenset[int]]:
    divisor_lines = hodge_eigenlines(phi, 1)
    answer: set[frozenset[int]] = set()
    for factors in combinations(divisor_lines, codimension):
        union = frozenset().union(*factors)
        if len(union) == 2 * codimension:
            answer.add(union)
    return answer


def main() -> None:
    roots = embedded_roots()
    pi = T**4 + T**3 - 4 * T**2 - 2 * T + 5
    base_discriminant = T**2 - 4
    assert sign_vector(base_discriminant, roots) == (-1,) * 6
    assert base_discriminant.norm() == 1

    sign_vectors: dict[int, tuple[int, ...]] = {}
    cm_types: dict[int, frozenset[int]] = {}
    certified_signs = 0

    for k in MAXIMAL_33_INDICES:
        ck = power_trace(k) + power_trace(k + 1) + power_trace(k + 4)
        yk = sine_quotient(k) + sine_quotient(k + 1) + sine_quotient(k + 4)
        delta_k = ck**2 - 4 * pi

        # This identity is also the order-generation check: y_k is an
        # integral norm-one unit, so alpha_k=x_k+y_k*zeta generates O_M over
        # O_L exactly when zeta does.
        assert yk.integral_power_basis()
        assert yk.norm() == 1
        assert delta_k == yk**2 * base_discriminant
        assert delta_k.norm() == 1
        assert sign_vector(delta_k, roots) == (-1,) * 6

        signs, phi_tuple = cm_type_for_index(k, roots)
        # sign_vector raises unless every rational evaluation interval lies
        # strictly on one side of zero.
        certified_signs += len(signs)
        sign_vectors[k] = signs
        phi = frozenset(phi_tuple)
        assert len(phi) == 6
        assert all(((exponent + 21) % 42) not in phi for exponent in phi)
        cm_types[k] = phi

    assert certified_signs == len(MAXIMAL_33_INDICES) * len(EMBEDDINGS)

    blocks_by_orbit: dict[frozenset[frozenset[int]], list[int]] = {}
    for k, phi in cm_types.items():
        blocks_by_orbit.setdefault(type_orbit(phi), []).append(k)
    blocks = tuple(
        tuple(sorted(indices))
        for indices in sorted(blocks_by_orbit.values(), key=lambda values: (len(values), values))
    )

    block_data: dict[tuple[int, ...], dict[str, object]] = {}
    for block in blocks:
        representative = cm_types[block[0]]
        stabilizer = type_stabilizer(representative)
        hodge_sets = tuple(hodge_eigenlines(representative, p) for p in range(1, 7))
        divisor_sets = tuple(divisor_generated_eigenlines(representative, p) for p in range(1, 7))
        counts = tuple(len(lines) for lines in hodge_sets)
        divisor_counts = tuple(len(lines) for lines in divisor_sets)
        exceptional_middle = hodge_sets[2] - divisor_sets[2]
        block_data[block] = {
            "stabilizer": stabilizer,
            "counts": counts,
            "divisor_counts": divisor_counts,
            "exceptional_middle": exceptional_middle,
        }

    imprimitive_block = next(block for block, data in block_data.items() if data["stabilizer"] == (1, 29))
    primitive_block = next(block for block, data in block_data.items() if data["stabilizer"] == (1,))
    primitive_exceptional = block_data[primitive_block]["exceptional_middle"]
    expected_weil_lines = {
        frozenset(exponent for exponent in GALOIS_GROUP if exponent % 3 == 1),
        frozenset(exponent for exponent in GALOIS_GROUP if exponent % 3 == 2),
    }
    assert primitive_exceptional == expected_weil_lines

    # Expected-output comparisons are deliberately delayed until every sign,
    # type orbit, and Hodge subset has been reconstructed from definitions.
    expected_blocks = (
        (1, 20, 22, 41),
        (11, 16, 17, 19, 32, 37, 38, 40),
    )
    assert blocks == expected_blocks
    assert imprimitive_block == expected_blocks[0]
    assert primitive_block == expected_blocks[1]
    assert block_data[imprimitive_block]["counts"] == (12, 51, 88, 51, 12, 1)
    assert block_data[imprimitive_block]["divisor_counts"] == (12, 51, 88, 51, 12, 1)
    assert block_data[primitive_block]["counts"] == (6, 15, 22, 15, 6, 1)
    assert block_data[primitive_block]["divisor_counts"] == (6, 15, 20, 15, 6, 1)
    assert len(primitive_exceptional) == 2

    print("INDEPENDENT MAXIMAL-(3,3) CM/HODGE STURM CHECK")
    print(f"certified CM-type signs: {certified_signs}/{certified_signs} (all rational intervals exclude zero)")
    print("per-index certified sign vectors and upper-line CM types:")
    for k in MAXIMAL_33_INDICES:
        phi = tuple(r if s > 0 else (-r) % 42 for r, s in zip(EMBEDDINGS, sign_vectors[k]))
        print(f"  k={k:02d}: signs={sign_vectors[k]} Phi={phi} stabilizer={type_stabilizer(cm_types[k])}")
    print("Galois blocks and Hodge data:")
    for block in blocks:
        data = block_data[block]
        print(
            f"  block={block} stabilizer={data['stabilizer']} "
            f"Hodge={data['counts']} divisor-products={data['divisor_counts']}"
        )
    exceptional_display = tuple(sorted(tuple(sorted(line)) for line in primitive_exceptional))
    print("primitive middle-degree exceptional eigenlines:", exceptional_display)
    print("primitive middle-degree deficit: 22 - 20 = 2")
    print("PASS: C14 decisive CM-type inputs are exact interval-certified for all twelve indices")


if __name__ == "__main__":
    main()
