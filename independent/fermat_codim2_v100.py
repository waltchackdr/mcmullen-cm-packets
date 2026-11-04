"""Independent exact Fermat B x D codimension-two verification.

No primary implementation is imported. Enumerate complementary sums of pairs
of residue-triple signatures, rather than four-subsets of CM weight columns.
Construct partitions from opposite pairs or one Galois orbit representative;
no recursive partition search is used. Stdlib only; preparation, not a release.
"""
import collections
import hashlib
import itertools
import json
import math
from pathlib import Path

M = 42
UNITS = tuple(a for a in range(1, M) if math.gcd(a, M) == 1)
TRIPLES = {"B": (10, 12, 20), "D": (1, 14, 27)}
EXPECTED_TYPES = {"B": (1, 5, 11, 13, 17, 19),
                  "D": (1, 5, 11, 13, 19, 25)}
REP_SUPPORT = (("B", 1), ("B", 19), ("D", 31), ("D", 41))
REP_BLOCKS = ((2, 10, 31, 41), (12, 15, 18, 39), (14, 28), (20, 22))


def residues(triple, multiplier=1, modulus=M):
    return tuple(multiplier * a % modulus for a in triple)


def cm_type(triple, modulus=M):
    return tuple(g for g in UNITS if sum(residues(triple, g, modulus)) == modulus)


def balanced(block):
    return (len(block) % 2 == 0 and all(0 < x < M for x in block)
            and all(sum(residues(block, g)) == M * len(block) // 2 for g in UNITS))


def acceptable_block(block):
    if len(block) not in (2, 4, 6) or not balanced(block):
        return False
    return len(block) != 6 or any(sum(t) % M == 0 for t in itertools.combinations(block, 3))


def canonical_type(phi):
    return min(tuple(sorted(g * a % M for a in phi)) for g in UNITS)


def triple_census(modulus, orbit_representatives):
    """Unordered triples then permutation multiplicities, unlike primary loops."""
    counts = [0] * len(orbit_representatives)
    exact_types = {factor: [] for factor in TRIPLES}
    for triple in itertools.combinations_with_replacement(range(1, modulus), 3):
        if sum(triple) % modulus:
            continue
        phi = cm_type(triple, modulus)
        for factor in exact_types:
            if phi == EXPECTED_TYPES[factor]:
                exact_types[factor].append(list(triple))
        if math.gcd(modulus, *triple) == 1:
            multiplicity = math.factorial(3)
            for count in collections.Counter(triple).values():
                multiplicity //= math.factorial(count)
            counts[orbit_representatives.index(canonical_type(phi))] += multiplicity
    return {"primitive_ordered_triples": sum(counts), "counts_by_census_type": counts}, exact_types


def check(primary_path=None):
    types = {factor: cm_type(triple) for factor, triple in TRIPLES.items()}
    assert types == EXPECTED_TYPES
    type_stabilizers = {f: [g for g in UNITS if tuple(sorted(g*a % M for a in phi)) == phi]
                        for f, phi in types.items()}
    assert all(s == [1] for s in type_stabilizers.values())
    assert canonical_type(types["B"]) != canonical_type(types["D"])
    labels = tuple((factor, g) for factor in TRIPLES for g in UNITS)
    actual_triples = {label: residues(TRIPLES[label[0]], label[1]) for label in labels}
    assert all(len(set(residues(triple, g) for g in UNITS)) == 12 for triple in TRIPLES.values())
    signatures = {label: tuple(sum(residues(actual_triples[label], h)) for h in UNITS) for label in labels}
    pairs_by_signature = collections.defaultdict(list)
    for a, b in itertools.combinations(labels, 2):
        signature = tuple(x + y for x, y in zip(signatures[a], signatures[b]))
        pairs_by_signature[signature].append(frozenset((a, b)))
    invariant_supports = set()
    for signature, pair_list in pairs_by_signature.items():
        complement = tuple(6 * M - x for x in signature)
        for p in pair_list:
            for q in pairs_by_signature.get(complement, ()):
                if p.isdisjoint(q):
                    invariant_supports.add(tuple(sorted(p | q)))
    divisor_pairs = pairs_by_signature[(3 * M,) * len(UNITS)]
    expected_pairs = {frozenset(((f, g), (f, M - g))) for f in TRIPLES for g in UNITS}
    assert set(divisor_pairs) == expected_pairs
    divisor_supports = {tuple(sorted(p | q)) for p, q in itertools.combinations(divisor_pairs, 2) if p.isdisjoint(q)}
    assert len(divisor_pairs) == 12 and len(divisor_supports) == 66
    exceptional = invariant_supports - divisor_supports
    orbit = {tuple(sorted((f, g * h % M) for f, g in REP_SUPPORT)): h for h in UNITS}
    assert len(orbit) == 12 and set(orbit) == exceptional
    assert len(invariant_supports) == 78
    certificates = []
    for support in sorted(invariant_supports):
        gamma = sorted(a for label in support for a in actual_triples[label])
        if support in orbit:
            h = orbit[support]
            blocks = [list(residues(block, h)) for block in REP_BLOCKS]
        else:
            blocks = []
            for f, g in support:
                if g < M // 2:
                    assert (f, M-g) in support
                    blocks.extend([[a, M-a] for a in actual_triples[f, g]])
        assert sorted(a for block in blocks for a in block) == gamma
        assert all(len(block) in (2, 4) and balanced(block) for block in blocks)
        certificates.append({"support": support, "blocks": blocks})
    conjugate_reps = tuple(g for g in UNITS if g < M // 2)
    type_reps = sorted({canonical_type(tuple((M - g if mask & (1 << i) else g)
                                            for i, g in enumerate(conjugate_reps)))
                        for mask in range(64)})
    censuses, presentations = {}, {}
    for modulus in (21, 42):
        censuses[str(modulus)], presentations[str(modulus)] = triple_census(modulus, type_reps)
    result = {
        "status": "PASS_PREPARATION_INDEPENDENT_EXACT_CHECK",
        "types": types, "type_stabilizers": type_stabilizers,
        "effective_conductors": {f: M // math.gcd(M, *t) for f, t in TRIPLES.items()},
        "pair_domain": math.comb(24, 2), "four_support_domain": math.comb(24, 4),
        "invariant_support_count": len(invariant_supports), "divisor_pair_count": len(divisor_pairs),
        "divisor_product_count": len(divisor_supports), "exceptional_orbit_size": len(orbit),
        "primitive_triple_census": censuses,
        "B_exact_type_triples": presentations["42"]["B"],
        "D_exact_type_triples": presentations["42"]["D"],
        "constructive_certificates": certificates,
    }
    if primary_path is not None:
        primary = json.loads(Path(primary_path).read_text())["exact_output_payload"]
        assert primary["primitive_triple_census"] == censuses
        for f in TRIPLES:
            assert primary[f + "_exact_type_triples"] == presentations["42"][f]
            assert primary[f + "_type"] == list(types[f])
        supplied = primary["selected_presentation_certificates"]
        assert len(supplied) == 78
        seen = set()
        for cert in supplied:
            support = tuple(sorted(tuple(x) for x in cert["support"]))
            assert support not in seen
            seen.add(support)
            gamma = sorted(a for label in support for a in actual_triples[label])
            assert cert["gamma"] == gamma
            assert sorted(a for block in cert["blocks"] for a in block) == gamma
            assert all(acceptable_block(block) for block in cert["blocks"])
        assert seen == invariant_supports
        result["primary_certificates_all_validated"] = True
        result["primary_sha256"] = hashlib.sha256(Path(primary_path).read_bytes()).hexdigest()
    return result


if __name__ == "__main__":
    import argparse
    import resource
    resource.setrlimit(resource.RLIMIT_CPU, (18, 19))
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary", type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.primary), sort_keys=True, separators=(",", ":")))
