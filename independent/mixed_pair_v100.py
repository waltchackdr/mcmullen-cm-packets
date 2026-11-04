#!/usr/bin/env python3
"""Independent mixed-pair counts from CRT characters and pair differences.

For each of six conjugate embedding pairs the selection difference is
-1, 0 or 1. Difference zero has exterior-degree polynomial 1+u^2;
either nonzero difference has polynomial u. This uses 729 difference
vectors per factor, not the primary producer's 4096 subset-bitset loop.
The same additive CRT Galois row acts on both factors of each product.
No primary producer or output is imported or read by this implementation.
These finite rational Hodge dimensions assert no mixed-cycle algebraicity.
"""
from collections import defaultdict
from hashlib import sha256
from itertools import product, combinations
from math import comb, gcd
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "experiments/fixed_field_cm_census_v98/independent.json"
SOURCE_SHA = "8d2a5d643e3fca75f06270d9fa5d102d7d7711f1dfebd1dc7184c75445a1d4fe"
GROUP = tuple(product(range(2), range(6)))
UNITS = tuple(a for a in range(42) if gcd(a, 42) == 1)
EXP = {g: next(a for a in UNITS if a % 3 == (-1)**g[0] % 3
               and a % 7 == pow(3, g[1], 7)) for g in GROUP}
ZERO = (0,) * 12


def minus(g, h):
    return ((g[0]-h[0]) % 2, (g[1]-h[1]) % 6)


def weights(phi):
    return {EXP[g]: tuple(1 if EXP[minus(g, h)] in phi else -1
                         for h in GROUP) for g in GROUP}


def difference_polynomial(phi):
    w = weights(phi)
    positive = [w[EXP[0, j]] for j in range(6)]
    assert all(w[EXP[1, (j+3) % 6]] == tuple(-x for x in positive[j])
               for j in range(6))
    answer = defaultdict(lambda: [0] * 13)
    for d in product((-1, 0, 1), repeat=6):
        z = d.count(0)
        total = tuple(sum(d[j]*positive[j][i] for j in range(6))
                      for i in range(12))
        for k in range(z+1):
            answer[total][6-z+2*k] += comb(z, k)
    assert sum(sum(v) for v in answer.values()) == 4096
    return answer, w


def convolution(a, b):
    out = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def polynomial_power(a, n):
    out = [1]
    for _ in range(n):
        out = convolution(out, a)
    return out


def self_control(rank):
    if rank == 5:
        out = polynomial_power([1, 4, 1], 6)
        for d in (1, 2):
            term = polynomial_power([comb(2, k+d)*comb(2, k)
                                     for k in range(3-d)], 6)
            for k, v in enumerate(term):
                out[3*d+k] += 2*v
        return out
    multiplicity = {6: 1, 3: 2, 1: 6}[rank]
    return polynomial_power([comb(2*multiplicity, k)**2
                             for k in range(2*multiplicity+1)], rank)


def check():
    raw = SOURCE.read_bytes()
    assert sha256(raw).hexdigest() == SOURCE_SHA
    records = json.loads(raw)["orbits"]
    assert len(records) == 8
    data = [difference_polynomial(set(x["representative"])) for x in records]
    pairs = []
    for i in range(8):
        for j in range(i, 8):
            left, lw = data[i]
            right, rw = data[j]
            bidegrees = [[0]*13 for _ in range(13)]
            for total, polynomial in left.items():
                opposite = right.get(tuple(-x for x in total))
                if opposite is None:
                    continue
                for a, v in enumerate(polynomial):
                    if not v:
                        continue
                    for b, w in enumerate(opposite):
                        bidegrees[a][b] += v*w
            separate = [[left[ZERO][a]*right[ZERO][b] for b in range(13)]
                        for a in range(13)]
            dims = [0]*25
            for a in range(13):
                for b in range(13):
                    dims[a+b] += bidegrees[a][b]
            assert not any(dims[1::2]) and dims == dims[::-1]
            sep = convolution(left[ZERO], right[ZERO])[::2]
            result = {"pair": [i, j], "hodge_dimensions_by_codimension": dims[::2],
                      "separate_factor_algebra_dimensions": sep,
                      "additional_mixed_dimensions": [a-b for a, b in zip(dims[::2], sep)],
                      "mixed_counts_by_bidegree": {
                          f"{a},{b}": bidegrees[a][b]-separate[a][b]
                          for a in range(13) for b in range(13)
                          if bidegrees[a][b] != separate[a][b]}}
            if i == j:
                assert dims[::2] == self_control(records[i]["hodge_group_rank"])
                result["self_n2_closed_formula_matches"] = True
            if (i, j) in ((0, 4), (0, 2)):
                # Direct enumeration of all degree-two and degree-four supports
                # on 24 separately indexed lines, with no difference grouping.
                columns = list(lw.values())+list(rw.values())
                low = {}
                for degree in (2, 4):
                    count = sum(all(sum(columns[k][h] for k in subset) == 0
                                    for h in range(12))
                                for subset in combinations(range(24), degree))
                    assert count == dims[degree]
                    low[str(degree)] = count
                result["direct_small_degree_support_counts"] = low
            pairs.append(result)
    assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    return {"status": "PASS_INDEPENDENT_FINITE_MIXED_PAIR_COUNTS_NO_ALGEBRAICITY_PROOF",
            "input_sha256": {str(SOURCE.relative_to(ROOT)): SOURCE_SHA},
            "method": "Additive CRT synchronized signed weights and 3^6 conjugate-pair difference vectors with binomial exterior-degree polynomials",
            "pairs": pairs, "mixed_algebraicity_proved": False}


if __name__ == "__main__":
    import resource
    import signal
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    signal.alarm(30)
    print(json.dumps(check(), sort_keys=True, indent=2))
