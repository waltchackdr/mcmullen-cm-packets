#!/usr/bin/env python3
"""Independent CRT six-plus-six hash join, with signed-multiset minimality.

Global completeness requires the separate signed-graph proof of the box and
l1 bounds. This executable neither solves pivot coordinates nor imports the
primary enumeration or its result.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "experiments/fixed_field_cm_census_v98/independent.json"
SOURCE_SHA = "8d2a5d643e3fca75f06270d9fa5d102d7d7711f1dfebd1dc7184c75445a1d4fe"
GROUP = tuple(product(range(2), range(6)))
UNITS = tuple(x for x in range(42) if gcd(x, 42) == 1)
CRT = {h: next(x for x in UNITS if x % 3 == pow(2, h[0], 3)
              and x % 7 == pow(3, h[1], 7)) for h in GROUP}
LOG = {x: h for h, x in CRT.items()}
REPRESENTATIVES = (1, 5, 11, 13, 17, 19)


def digest(value):
    return sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def column_weights(cm_type):
    return [tuple(1 if CRT[((LOG[a][0]-h[0]) % 2,
                           (LOG[a][1]-h[1]) % 6)] in cm_type else -1
                  for h in GROUP) for a in REPRESENTATIVES]


def half_table(columns):
    table = defaultdict(list)
    for vector in product(range(-2, 3), repeat=6):
        weight = tuple(sum(vector[i]*columns[i][r] for i in range(6))
                       for r in range(12))
        table[weight].append(vector)
    assert sum(map(len, table.values())) == 15625
    return table


def multiset(vector):
    # Initial segments encode multiplicity; inclusion is exactly conformality.
    return frozenset((i, 1 if x > 0 else -1, j)
                     for i, x in enumerate(vector) for j in range(abs(x)))


def check():
    raw = (ROOT / SOURCE).read_bytes()
    assert sha256(raw).hexdigest() == SOURCE_SHA
    orbits = json.loads(raw)["orbits"]
    columns = [column_weights(set(orbits[i]["representative"])) for i in (0, 2)]
    left, right = [half_table(c) for c in columns]
    candidates = set()
    for weight, vectors in left.items():
        for a in vectors:
            for b in right.get(tuple(-x for x in weight), ()):
                z = a + b
                if 0 < sum(map(abs, z)) <= 12:
                    candidates.add(z)
    assert all(tuple(-x for x in z) in candidates for z in candidates)
    supports = {z: multiset(z) for z in candidates}
    # Compare against ALL candidate supports, not just already discovered atoms.
    atoms = {z for z, s in supports.items()
             if not any(t < s for t in supports.values())}
    ordered_atoms = sorted(atoms, key=lambda z: (sum(map(abs, z)), z))
    witnesses = []
    for z in sorted(candidates):
        residual = z
        summands = []
        while any(residual):
            s = multiset(residual)
            g = next(g for g in ordered_atoms if supports[g] <= s)
            summands.append(g)
            residual = tuple(x-y for x, y in zip(residual, g))
        assert tuple(sum(g[i] for g in summands) for i in range(12)) == z
        assert sum(sum(map(abs, g)) for g in summands) == sum(map(abs, z))
        witnesses.append([z, summands])
    canonical = lambda values: sorted(z for z in values if next(x for x in z if x) > 0)
    graver = canonical(atoms)
    nonzero = canonical(candidates)
    full_columns = columns[0]+columns[1]
    assert all(all(sum(z[j]*full_columns[j][r] for j in range(12)) == 0
                   for r in range(12)) for z in candidates)
    return {
        "status": "EXACT_HASH_JOIN_AND_SIGNED_MULTISET_CHECK_NO_ALGEBRAICITY",
        "input_sha256": {SOURCE: SOURCE_SHA},
        "column_order": [[factor, a] for factor in (0, 2) for a in REPRESENTATIVES],
        "crt_row_order": [list(h) for h in GROUP],
        "crt_full_weight_columns": full_columns,
        "half_tuples_each": 15625,
        "candidate_count_up_to_sign": len(nonzero),
        "candidate_vectors_up_to_sign": nonzero,
        "candidate_digest": digest(nonzero),
        "graver_count_up_to_sign": len(graver),
        "graver_vectors_up_to_sign": graver,
        "graver_digest": digest(graver),
        "counts_by_degree": dict(sorted(Counter(sum(map(abs, z)) for z in graver).items())),
        "counts_by_maximum_coefficient": dict(sorted(Counter(max(map(abs, z)) for z in graver).items())),
        "fully_decomposed_signed_candidates": len(witnesses),
        "full_conformal_decomposition_digest": digest(witnesses),
        "scope": "Graver completeness uses the independently reviewed signed-graph bound l1<=12 and coefficient<=2; no extrapolation or algebraicity.",
    }


if __name__ == "__main__":
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_CPU, (20, 22))
    signal.alarm(30)
    assert __debug__
    print(json.dumps(check(), indent=2, sort_keys=True))
