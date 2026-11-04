#!/usr/bin/env python3
"""Exact reconstruction of the reviewed full 42-trace arithmetic census.

All local rings, ideals and anti-actions are reconstructed in the retained
standard-library quadratic path. Sources are pinned before and after execution;
all five mathematical artifacts must match their retained bytes. No geometry
search is performed. Human proof: docs/FULL_PACKET_ARITHMETIC_V62.md.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import cache
from hashlib import sha256
import json
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError("Assertions must be enabled")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "certificates/full_packet_v62"
ENTRY_SHA256 = sha256(Path(__file__).read_bytes()).hexdigest()
SOURCE_PINS = {'independent/d3_completeness_quadratic_v60.py': '7b34e8c676a0045b8c02d1ae065b395003bb48191d1cddf8d085a458db4f0363', 'independent/packet_order_coverage_stdlib_v58.py': 'ed328d11d87d2cc9e93322bddb2b3207b1edf8c5d8dc66dea2026f0a07ea07c3', 'independent/quadrilateral_k17_stdlib_v41.py': 'c5bd45f334a028820a071b5c758acfb8e48afc1e8454b18a3174dd04ca58a38c'}
ARTIFACT_PINS = {'all_42_matrix_classes.jsonl': '08492a608ce13b99b6e9470241379a1a01746e78539b126aeb785a15af9b9761', 'ambient_points.json': 'c0e84da6c266c087cbce65a6e8aca8e3f30350ce56b5a1089c4ae4978a381f36', 'local_quotients.json': '9ef01be9df80366e398b75a7aad39b4705ea06c9a84ff884d6b4a08b69f05613', 'module_templates.json': '09ccff217511a67e946567135f9ffb74da020843621d0895ecce8746d25c245d', 'root_bridges.jsonl': '13dcf16cd44e5359659b307dc9f0eb45d8d0864b8dda9c3aabb4291bb620ece6'}


def verify_sources():
    assert sha256(Path(__file__).read_bytes()).hexdigest() == ENTRY_SHA256, "entrypoint changed during execution"
    for relative, expected in SOURCE_PINS.items():
        assert sha256((ROOT/relative).read_bytes()).hexdigest() == expected, relative


verify_sources()
for name, expected in ARTIFACT_PINS.items():
    assert sha256((DATA/name).read_bytes()).hexdigest() == expected, name
sys.path.insert(0, str(ROOT / "independent"))
import packet_order_coverage_stdlib_v58 as coverage
import d3_completeness_quadratic_v60 as foundation

F = coverage.field
for module, relative in ((coverage, "independent/packet_order_coverage_stdlib_v58.py"),
                         (foundation, "independent/d3_completeness_quadratic_v60.py"),
                         (F, "independent/quadrilateral_k17_stdlib_v41.py")):
    assert Path(module.__file__).resolve() == (ROOT/relative).resolve(), relative
T, ONE, ZERO, PI = F.T, F.ONE, F.ZERO, coverage.PI
IDENTITY = F.IDENTITY
U = 4-T*T


@cache
def inverse(a):
    return a.inverse()


def ints(a):
    return coverage.coefficients(a)


def rational(a):
    return [str(c) for c in a.coefficients]


def flat(a):
    assert F.all_entries_integral(a)
    return [ints(c) for row in a for c in row]


def pair_mul(x, y):
    return foundation.multiply(x, y)


def action_pair(pair):
    a, b = pair
    return F.matrix(a, -b, b, a+T*b)


def other_root(matrix):
    a, b = matrix[0]
    c, d = matrix[1]
    return F.matrix(d, -b, -c, a)


def local_ring(p, r):
    def multiply(a, b):
        x, y = a
        v, w = b
        return ((x*v-y*w) % p, (x*w+y*v+r*y*w) % p)

    def norm(a):
        x, y = a
        return (x*x+r*x*y+y*y) % p

    def normalize(a):
        assert norm(a)
        return (1, a[1]*pow(a[0], -1, p) % p) if a[0] else (0, 1)

    def powers(a):
        answer, value = [], (1, 0)
        while value not in answer:
            answer.append(value)
            value = normalize(multiply(value, a))
        assert value == (1, 0)
        return answer

    quotient = {(1, s) for s in range(p) if norm((1, s))} | {(0, 1)}
    epsilon = normalize((2, -r % p))
    image = powers(epsilon)
    zz = powers((0, 1))
    assert image[22 % len(image)] == (0, 1)
    roots = [s for s in range(p) if (s*s-r*s+1) % p == 0]
    assert len(roots) in (0, 2)
    assert len(quotient) == p-1 if roots else len(quotient) == p+1
    primitive = next((1, s) for s in range(p)
                     if norm((1, s)) and len(powers((1, s))) == len(quotient))
    assert set(powers(primitive)) == quotient
    data = dict(prime=p, trace_root=r, splitting="split" if roots else "inert",
                z_roots=roots, quotient_order=len(quotient),
                epsilon_projective=list(epsilon), epsilon_order=len(image),
                z_order=len(zz), primitive_beta=list(primitive),
                picard_order=len(quotient)//len(image))
    return data, multiply, normalize, powers, image


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    arithmetic = coverage.reconstruct()
    powers = coverage.z_powers()
    alpha = [tuple(sum((powers[k+j][i] for j in (0, 1, 4)), ZERO) for i in (0, 1))
             for k in range(42)]
    old_raw, _ = foundation.reconstruct()
    units = {}
    for line in old_raw.splitlines():
        row = json.loads(line)
        units[row["orientation_mask"]] = F.LElement(row["unit"])
    assert set(units) == set(range(64))
    unit_inverse = {m: inverse(u) for m, u in units.items()}
    conductors = {row["trace_index"]: tuple(row["conductor_key"])
                  for row in arithmetic["traces"]}
    strata = {
        "OM": ((1, 0, 0), ONE, 1, (0,)),
        "R7": ((7, 1, 1), alpha[0][1], 0, (0,)),
        "R7sq": ((7, 1, 2), alpha[12][1], 12, (0,)),
        "R41_12": ((41, 12, 1), alpha[14][1], 14, (0,)),
        "R41_22": ((41, 22, 1), alpha[5][1], 5, (0,)),
        "R43_20": ((43, 20, 1), alpha[4][1], 4, (0,)),
        "R43_31": ((43, 31, 1), alpha[8][1], 8, (0,)),
        "R83": ((83, 36, 1), alpha[2][1], 2, (0, 3)),
        "R127": ((127, 94, 1), alpha[10][1], 10, (0, 9, 100)),
    }

    local = {}
    for name, ((p, r, exponent), g, k, slopes) in strata.items():
        if name in ("OM", "R7sq"):
            continue
        data, mul, normalize, finite_powers, image = local_ring(p, r)
        local[name] = data
        assert data["picard_order"] == len(slopes)
        assert data["epsilon_order"] in (6, 42)
        cosets = [{normalize(mul((1, s), value)) for value in image} for s in slopes]
        assert all(len(coset) == len(image) for coset in cosets)
        assert sum(map(len, cosets)) == len(set().union(*cosets)) == data["quotient_order"]
        data["module_slopes"] = list(slopes)
        data["coset_representatives_are_disjoint_and_exhaustive"] = True
        if name == "R127":
            assert data["splitting"] == "split" and data["z_roots"] == [40, 54]
            assert len(finite_powers((1, 9))) == 126
            assert normalize((1+9*r, -9 % p)) == (1, 100)
            assert cosets[1] != cosets[2]
            data["conjugation_on_picard_classes"] = [0, 2, 1]
        if name == "R83":
            assert len(finite_powers((1, 3))) == 84
            beta_square = normalize(mul((1, 3), (1, 3)))
            assert image.index(beta_square) == 29
            data["epsilon_power_equal_to_beta_square"] = 29
            data["conjugation_on_picard_classes"] = [0, 1]
    assert local["R41_22"]["splitting"] == "inert"
    assert local["R41_22"]["quotient_order"] == local["R41_22"]["epsilon_order"] == 42

    # The square conductor has already been reviewed; recheck its primitive
    # ratio explicitly in F7[v]/v^2, t=1+v and z=3+2v.
    dual_epsilon, value, order = (5, 2), (1, 0), 0
    while True:
        value = (value[0]*5 % 7, (value[0]*2+value[1]*5) % 7)
        order += 1
        if value == (1, 0):
            break
        assert order < 42
    assert order == 42
    local["R7sq"] = dict(prime=7, trace_root=1, exponent=2,
                          quotient_order=42, epsilon_order=42, picard_order=1,
                          residue_ring="F7[v]/(v^2)", epsilon_ratio=list(dual_epsilon))

    @cache
    def action(k, name, j, mask):
        _, g, _, slopes = strata[name]
        s = slopes[j]
        a, b = alpha[k]
        n = 1+s*T+s*s
        delta = units[mask]
        result = F.matrix(a-s*b, -delta*g*b,
                          b*inverse(g)*n*unit_inverse[mask], a+(T+s)*b)
        assert F.all_entries_integral(result)
        assert F.mdet(result) == PI
        assert result[0][0]+result[1][1] == 2*a+T*b
        return result

    def basis(name, j):
        _, g, _, slopes = strata[name]
        return F.matrix(1, 0, slopes[j], g)

    def partner(name, j, mask):
        if name == "R127" and j:
            return 3-j, mask^16
        return j, mask^(48 if name == "R83" and j else 16)

    def canonical(name, j, mask):
        return min((j, mask), partner(name, j, mask))

    def label(name, j, mask):
        return f"{name}:pic{j}:{mask:02d}"

    epsilon = (ONE-powers[2][0], -powers[2][1])
    eta83 = pair_mul(epsilon, powers[28])
    assert foundation.relative_norm(eta83) == U
    standard_conjugation = F.matrix(1, T, 0, -1)
    anti, module_records = {}, []
    for name, (key, g, k, slopes) in strata.items():
        for j, s in enumerate(slopes):
            bmat = basis(name, j)
            z_action = F.mmul(F.mmul(F.minverse(bmat), action_pair((ZERO, ONE))), bmat)
            assert F.all_entries_integral(z_action) == (name == "OM")
            target_j = 3-j if name == "R127" and j else j
            eta = eta83 if name == "R83" and j else (ONE, ZERO)
            c = F.mmul(F.mmul(F.mmul(F.minverse(basis(name, target_j)), action_pair(eta)),
                              standard_conjugation), bmat)
            assert F.all_entries_integral(c)
            assert F.mdet(c) == (-U if name == "R83" and j else -ONE)
            assert F.all_entries_integral(F.minverse(c))
            assert F.mmul(c, action(k, name, j, 0)) == F.mmul(other_root(action(k, name, target_j, 0)), c)
            anti[name, j] = c
            module_records.append(dict(stratum=name, conductor=list(key), picard_class=j,
                slope=s, basis_columns_in_one_z=flat(bmat), index_norm=int(abs(g.norm())),
                representative_trace=k, norm_pi_action=flat(action(k, name, j, 0)),
                conjugate_picard_class=target_j, anti_intertwiner=flat(c),
                anti_determinant=ints(F.mdet(c))))

    @cache
    def root_bridge(name, j, mask):
        next_j, next_mask = partner(name, j, mask)
        bridge = F.mmul(F.mmul(F.matrix(units[next_mask], 0, 0, 1), anti[name, j]),
                       F.matrix(unit_inverse[mask], 0, 0, 1))
        if name == "R83" and j and not mask & 32:
            bridge = F.mscale(inverse(U), bridge)
        assert F.all_entries_integral(bridge) and F.mdet(bridge) == ONE
        return next_j, next_mask, bridge

    points = []
    for name, (key, g, k, slopes) in strata.items():
        for j, s in enumerate(slopes):
            for mask in range(32 if name == "OM" else 64):
                if (j, mask) != canonical(name, j, mask):
                    continue
                q = g*units[mask]
                n = 1+s*T+s*s
                x2 = -q*(T+2*s)*inverse(n)
                h4 = q*q*U*inverse(n)**2
                mod2 = q*q*inverse(n)
                assert x2*x2+h4 == 4*mod2
                m = action(k, name, j, mask)
                assert (m[0][0]-m[1][1])*inverse(m[1][0]) == x2
                points.append(dict(point_label=label(name, j, mask), stratum=name,
                    conductor=list(key), picard_class=j, orientation_mask=mask,
                    trace_for_matrix=k, norm_pi_matrix=flat(m),
                    root_partner=list(partner(name, j, mask)),
                    twice_real_part=rational(x2), four_times_height_squared=rational(h4),
                    modulus_squared=rational(mod2)))
    point_labels = {p["point_label"] for p in points}
    assert len(points) == len(point_labels) == 368

    full_manifest, bridges = [], []
    incidences = Counter()
    for k in range(42):
        cp, cr, ce = conductors[k]
        for name, ((p, r, e), _, _, slopes) in strata.items():
            if name != "OM" and not (cp == p and cr == r and ce >= e):
                continue
            for j in range(len(slopes)):
                for mask in range(32 if name == "OM" else 64):
                    matrix = action(k, name, j, mask)
                    canonical_j, canonical_mask = canonical(name, j, mask)
                    destination = label(name, canonical_j, canonical_mask)
                    assert destination in point_labels
                    next_j, next_mask, bridge = root_bridge(name, j, mask)
                    assert F.mmul(bridge, matrix) == F.mmul(other_root(action(k, name, next_j, next_mask)), bridge)
                    row = dict(trace_index=k, stratum=name, picard_class=j,
                               orientation_mask=mask, matrix_row_major=flat(matrix),
                               canonical_point_label=destination,
                               root_partner=[next_j, next_mask])
                    full_manifest.append(row)
                    incidences[destination] += 1
                    if k == strata[name][2]:
                        bridges.append(dict(stratum=name, picard_class=j, orientation_mask=mask,
                            target_picard_class=next_j, target_orientation_mask=next_mask,
                            reference_trace=k, SL2_root_bridge=flat(bridge)))
    assert len(full_manifest) == 3520
    assert len(bridges) == 736
    matrix_labels = {(x["trace_index"], x["stratum"], x["picard_class"], x["orientation_mask"])
                     for x in full_manifest}
    assert len(matrix_labels) == 3520
    point_counts = dict(Counter(p["stratum"] for p in points))
    assert point_counts == dict(OM=16, R7=32, R7sq=32, R41_12=32, R41_22=32,
                               R43_20=32, R43_31=32, R83=64, R127=96)
    expected_incidence = dict(OM=84, R7=28, R7sq=4, R41_12=4, R41_22=4,
                              R43_20=4, R43_31=4, R83=4, R127=4)
    assert all(incidences[p["point_label"]] == expected_incidence[p["stratum"]] for p in points)
    per_trace = dict(Counter(row["trace_index"] for row in full_manifest))
    table = []
    for group in arithmetic["embedded_order_groups"]:
        indices = group["trace_indices"]
        counts = {per_trace[k] for k in indices}
        assert len(counts) == 1
        table.append(dict(conductor=group["conductor_key"], trace_indices=indices,
                          matrix_classes_per_trace=counts.pop(),
                          total_matrix_classes=sum(per_trace[k] for k in indices)))

    artifacts = {
        "local_quotients.json": (json.dumps(local, indent=2, sort_keys=True)+"\n").encode(),
        "module_templates.json": (json.dumps(module_records, indent=2, sort_keys=True)+"\n").encode(),
        "all_42_matrix_classes.jsonl": "".join(encode(row) for row in full_manifest).encode(),
        "ambient_points.json": (json.dumps(points, indent=2, sort_keys=True)+"\n").encode(),
        "root_bridges.jsonl": "".join(encode(row) for row in bridges).encode(),
    }
    summary = dict(status="PASS_EXACT_FULL_PACKET_ARITHMETIC_3520_CLASSES_368_POINTS",
        trace_count=42, embedded_orders=9, matrix_class_count=3520, ambient_point_count=368,
        point_counts_by_multiplier=point_counts, incidences_per_point=expected_incidence,
        table_by_generated_order=table, root_bridge_count=len(bridges),
        scope="Arithmetic classes and ambient SL2 quotient points; no curve membership or geometric CM claims",
        artifact_sha256={name: sha256(raw).hexdigest() for name, raw in artifacts.items()},
        source_sha256={str(Path(module.__file__).relative_to(ROOT)): sha256(Path(module.__file__).read_bytes()).hexdigest()
                       for module in (coverage, foundation, F)},
        implementation_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    assert set(artifacts) == set(ARTIFACT_PINS)
    for name, raw in artifacts.items():
        assert sha256(raw).hexdigest() == ARTIFACT_PINS[name], name
        assert (DATA/name).read_bytes() == raw, name
    verify_sources()
    print(json.dumps(summary, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
