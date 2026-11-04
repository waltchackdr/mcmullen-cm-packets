#!/usr/bin/env python3
"""Independent one-copy primitive-pair exterior atoms via pair differences.

Each factor uses 3^6 signed conjugate-pair differences. Hash-joining
opposite synchronized CRT weights gives every difference relation.
Difference zero has degree polynomial 1+u^2; either nonzero difference
has degree polynomial u. Atoms containing an opposite pair are exactly
that divisor pair; all other atoms are the inclusion-minimal signed
supports among the joined difference relations. This is not an integer
Hilbert-basis calculation for arbitrary powers or an algebraicity proof.
"""
from collections import defaultdict, Counter
from hashlib import sha256
from itertools import product
from math import gcd, comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "experiments/fixed_field_cm_census_v98/independent.json"
SOURCE_SHA = "8d2a5d643e3fca75f06270d9fa5d102d7d7711f1dfebd1dc7184c75445a1d4fe"
CRT = tuple(product(range(2), range(6)))
UNITS = tuple(a for a in range(42) if gcd(a, 42) == 1)
EXP = {g: next(a for a in UNITS if a % 3 == (-1)**g[0] % 3
               and a % 7 == pow(3, g[1], 7)) for g in CRT}
PAIR_REPRESENTATIVES = tuple(EXP[0, j] for j in range(6))


def difference_records(phi):
    columns = [tuple(1 if EXP[(-h[0] % 2, (j-h[1]) % 6)] in phi else -1
                     for h in CRT) for j in range(6)]
    result = defaultdict(list)
    for d in product((-1, 0, 1), repeat=6):
        total = tuple(sum(d[j]*columns[j][h] for j in range(6)) for h in range(12))
        result[total].append(d)
    assert sum(map(len, result.values())) == 729
    return result


def selected_support(vector):
    return frozenset(("B" if k < 6 else "D",
                      PAIR_REPRESENTATIVES[k % 6] if sign > 0
                      else -PAIR_REPRESENTATIVES[k % 6] % 42)
                     for k, sign in enumerate(vector) if sign)


def check():
    raw = SOURCE.read_bytes()
    assert sha256(raw).hexdigest() == SOURCE_SHA
    types = json.loads(raw)["orbits"]
    left = difference_records(set(types[0]["representative"]))
    right = difference_records(set(types[2]["representative"]))
    relations = [a+b for weight, aa in left.items()
                 for a in aa for b in right.get(tuple(-x for x in weight), ())]
    assert len(set(relations)) == len(relations)
    dimensions = [0]*25
    signed_supports = set()
    for relation in relations:
        occupied = sum(x != 0 for x in relation)
        zero_pairs = 12-occupied
        for k in range(zero_pairs+1):
            dimensions[occupied+2*k] += comb(zero_pairs, k)
        if occupied:
            signed_supports.add(selected_support(relation))
    assert not any(dimensions[1::2]) and dimensions == dimensions[::-1]
    assert sum(dimensions) == 11856
    # Signed support inclusion, not mere unsigned direction inclusion:
    # a subrelation with an unavailable opposite sign cannot factor a wedge.
    nonpair_atoms = [s for s in signed_supports if not any(t < s for t in signed_supports)]
    pair_atoms = [frozenset(((factor, a), (factor, -a % 42)))
                  for factor in ("B", "D") for a in PAIR_REPRESENTATIVES]
    atoms = pair_atoms+nonpair_atoms
    assert len(set(atoms)) == len(atoms)
    counts = Counter(len(s)//2 for s in atoms)
    assert counts == {1: 12, 2: 12, 3: 64}
    ordered = sorted(atoms, key=lambda s: (len(s), sorted(s)))
    degree6 = [s for s in atoms if len(s) == 6]
    pure = [s for s in degree6 if len({label for label, _ in s}) == 1]
    assert len(pure) == 4 and len(degree6)-len(pure) == 60
    assert all(any(atom <= support for atom in nonpair_atoms) for support in signed_supports)
    assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    return {"status": "PASS_INDEPENDENT_ONE_COPY_ATOMS_NO_ALGEBRAICITY_OR_ALL_POWERS_CLAIM",
            "input_sha256": {str(SOURCE.relative_to(ROOT)): SOURCE_SHA},
            "method": "Synchronized additive CRT weights;729 differences per factor; signed-support inclusion and degree polynomials",
            "pair_representatives_mod42": list(PAIR_REPRESENTATIVES),
            "joined_difference_relation_count_including_zero": len(relations),
            "invariant_dimensions": {str(k): v for k, v in enumerate(dimensions[::2])},
            "invariant_support_count": sum(dimensions),
            "atom_counts_by_codimension": {str(k): v for k, v in sorted(counts.items())},
            "atoms": [{"support": [list(x) for x in sorted(s)]} for s in ordered],
            "degree6_pure_atoms": [[list(x) for x in sorted(s)] for s in sorted(pure, key=lambda s: sorted(s))],
            "degree6_mixed_count": 60,
            "atom_counts_by_bidegree": {f"{a},{b}": v for (a, b), v in sorted(Counter(
                (sum(f == "B" for f, _ in s), sum(f == "D" for f, _ in s)) for s in atoms).items())},
            "scope": "Only exterior-support atoms on one copy of each of the two exceptional primitive sixfold types. Neither rational algebraicity nor arbitrary-multiplicity semigroup generation is asserted."}


if __name__ == "__main__":
    import resource
    import signal
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    signal.alarm(30)
    print(json.dumps(check(), sort_keys=True, indent=2))
