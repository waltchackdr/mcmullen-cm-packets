#!/usr/bin/env python3
"""Independent CRT/permutation-determinant audit of the B x J7 weights.

Uses no SymPy and reads no primary implementation or result. Complete
signed-unit enumeration and all minors of a normalized configuration
support a separate human conformal-circuit generation argument.
"""
from collections import Counter
from hashlib import sha256
from itertools import permutations, combinations, product
from math import gcd
from pathlib import Path
import json

SOURCE = Path(__file__).resolve().parents[1] / "experiments/fixed_field_cm_census_v98/independent.json"
SOURCE_SHA = "8d2a5d643e3fca75f06270d9fa5d102d7d7711f1dfebd1dc7184c75445a1d4fe"
G = tuple(a for a in range(42) if gcd(a, 42) == 1)
CRT = tuple(product(range(2), range(6)))
EXP = {g: next(a for a in G if a % 3 == (-1)**g[0] % 3
               and a % 7 == pow(3, g[1], 7)) for g in CRT}
COORD = {a: g for g, a in EXP.items()}


def determinant(matrix):
    n = len(matrix)
    total = 0
    for p in permutations(range(n)):
        value = (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        for i in range(n):
            value *= matrix[i][p[i]]
        total += value
    return total


def replace_column(matrix, column, values):
    result = [r[:] for r in matrix]
    for i, value in enumerate(values):
        result[i][column] = value
    return result


def support(v):
    return frozenset(i for i, x in enumerate(v) if x)


def check():
    raw = SOURCE.read_bytes()
    assert sha256(raw).hexdigest() == SOURCE_SHA
    types = json.loads(raw)["orbits"]
    columns = []
    multiplicities = {}
    for name, index in (("B", 0), ("J7", 4)):
        phi = set(types[index]["representative"])
        full = []
        for a in G:
            ga = COORD[a]
            values = []
            for h in CRT:
                gh = ((ga[0]-h[0]) % 2, (ga[1]-h[1]) % 6)
                values.append(1 if EXP[gh] in phi else -1)
            # Primary ga convention is obtained by common row inversion.
            order = [CRT.index(COORD[pow(g, -1, 42)]) for g in G]
            column = tuple(values[k] for k in order)
            full.append(column)
            existing = [tuple(c["weights"]) for c in columns]
            if column not in existing and tuple(-x for x in column) not in existing:
                columns.append({"factor": name, "embedding_representative": a,
                                "weights": list(column)})
        multiplicities[name] = sorted(Counter(full).values())
    assert len(columns) == 9
    full_matrix = [[c["weights"][h] for c in columns] for h in range(12)]
    W = full_matrix[:5]
    D = [r[:5] for r in W]
    basis_det = determinant(D)
    assert basis_det == -16
    C = [[0]*9 for _ in range(5)]
    for k in range(9):
        for q in range(5):
            v = determinant(replace_column(D, q, [r[k] for r in W]))
            assert v % basis_det == 0
            C[q][k] = v // basis_det
    assert all(sum(D[i][q]*C[q][k] for q in range(5)) == W[i][k]
               for i in range(5) for k in range(9))
    # Verify every discarded row is in the first five rows' rational span.
    Dt = [list(r) for r in zip(*D)]
    for row in full_matrix:
        numerators = [determinant(replace_column(Dt, q, row[:5])) for q in range(5)]
        assert all(sum(numerators[q]*W[q][k] for q in range(5)) == basis_det*row[k]
                   for k in range(9))
    minors = [{"columns": list(ix),
               "determinant": determinant([[row[i] for i in ix] for row in W])}
              for ix in combinations(range(9), 5)]
    assert {abs(v["determinant"]) for v in minors} == {0, 16}
    histogram = Counter()
    for k in range(1, 6):
        for rows in combinations(range(5), k):
            for cols in combinations(range(9), k):
                histogram[determinant([[C[i][j] for j in cols] for i in rows])] += 1
    assert set(histogram) <= {-1, 0, 1} and sum(histogram.values()) == 2001
    relations = []
    for v in product((-1, 0, 1), repeat=9):
        if next((x for x in v if x), 0) != 1:
            continue
        if all(sum(row[i]*v[i] for i in range(9)) == 0 for row in full_matrix):
            relations.append(v)
    circuits = sorted(v for v in relations if not any(support(w) < support(v) for w in relations))
    assert len(circuits) == 15
    assert Counter(len(support(v)) for v in circuits) == {4: 9, 6: 6}
    assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    return {"status": "PASS_EXACT_UNIMODULAR_ARITHMETIC_NOT_ALGEBRAICITY",
            "input_sha256": SOURCE_SHA, "columns": columns,
            "factor_signed_weight_multiplicities": multiplicities,
            "rank": 5, "reduced_matrix": W,
            "maximal_minors": minors, "basis_determinant": basis_det,
            "normalized_matrix": C,
            "normalized_square_minor_histogram": {str(k): v for k, v in sorted(histogram.items())},
            "signed_unit_relations_up_to_sign": len(relations),
            "circuits_up_to_sign": [list(v) for v in circuits],
            "circuit_counts_by_degree": {str(k): v for k, v in sorted(Counter(len(support(v)) for v in circuits).items())},
            "qualification": "Raw W has maximal minors 0 or +/-16. Normalized D^-1 W is totally unimodular; no claim raw W is totally unimodular. Circuit completeness uses the separate human TU-to-unit-circuits proof.",
            "algebraicity_proved": False}


if __name__ == "__main__":
    import resource
    import signal
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    signal.alarm(30)
    print(json.dumps(check(), sort_keys=True, indent=2))
