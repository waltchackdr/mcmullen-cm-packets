#!/usr/bin/env python3
"""Exact d3 point representatives and 1,280 SL2 transport witnesses.

Registered finite representative/transport check; no curve searches. Uses the
pinned foundation ring implementation, so this is not an independent replay
of its arithmetic foundations. Reconstructed artifacts must match retained data.
"""
from collections import Counter
from functools import cache
from hashlib import sha256
from pathlib import Path
import json
import resource
import runpy
import time

if not __debug__:
    raise RuntimeError("Assertions must be enabled")

START_CPU = time.process_time()
resource.setrlimit(resource.RLIMIT_CPU, (25, 25))
ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "d3_points_v60"
PINS = {
    "certificates/d3_completeness_v60.py":
        "7341c4dddc68e06e3288cda31372e6ad3a7617074a01cb2292526eba17baaf17",
    "certificates/d3_points_v60/foundation_inputs.json":
        "378559344333a813c4b94e721bafbf44581232a489a0cd16f24751a15569adf1",
    "certificates/d3_matrix_manifest_v60.jsonl":
        "fc87260544e8d8418e7bc6fe3ca80ef7d9715a82fcdaea9b2a9716016e00004e",
    "certificates/d3_points_v60/trace_order_inputs.json":
        "2aef609b6ce72dfb8b6a9017b31674995f0a3cb2cd3ce80aaf7b763e3edfdcd2",
    "certificates/d7_ambient_points_v45.py":
        "34cb35db91b953162b595fb61e031f551e6df0a0d7b5192dc15a6c00c39274a2",
}
for filename, expected in PINS.items():
    assert sha256((ROOT/filename).read_bytes()).hexdigest() == expected, filename
ARITH = runpy.run_path(str(ROOT/next(iter(PINS))))
L, M = ARITH["L"], ARITH["M"]
L.divide_integrally = cache(L.divide_integrally)
L.norm = cache(L.norm)
relative = ARITH["relative_coordinates"]
sign_vector = cache(ARITH["sign_vector"])
ZERO, ONE, T = L.zero, L.one, L.x
NEG = lambda a: L.scale(a, -1)
SUB = lambda a, b: L.add(a, NEG(b))
I = (ONE, ZERO, ZERO, ONE)


def mm(a, b):
    return (L.add(L.mul(a[0], b[0]), L.mul(a[1], b[2])),
            L.add(L.mul(a[0], b[1]), L.mul(a[1], b[3])),
            L.add(L.mul(a[2], b[0]), L.mul(a[3], b[2])),
            L.add(L.mul(a[2], b[1]), L.mul(a[3], b[3])))


def det(a):
    return SUB(L.mul(a[0], a[3]), L.mul(a[1], a[2]))


def inverse_sl2(a):
    assert det(a) == ONE
    return (a[3], NEG(a[1]), NEG(a[2]), a[0])


def conjugate(a, b):
    return mm(mm(a, b), inverse_sl2(a))


def other_root(a):
    return (a[3], NEG(a[1]), NEG(a[2]), a[0])


def diag(a, b):
    return (a, ZERO, ZERO, b)


foundation = json.loads((ROOT/list(PINS)[1]).read_text())
raw_input = (ROOT/list(PINS)[2]).read_text().splitlines()
input_rows = [json.loads(line) for line in raw_input]
coverage = json.loads((ROOT/list(PINS)[3]).read_text())
generators = [tuple(a) for a in foundation["sign_generator_coefficients"]]
generators.append(tuple(foundation["positive_nonsquare_coefficients"]))
assert generators[4] == NEG(ONE)
U = generators[5]
assert U == SUB(L.scale(ONE, 4), L.mul(T, T))
assert sign_vector(U) == (1,)*6
generator_inverses = [L.divide_integrally(ONE, a) for a in generators]
units, inverse_units = [], []
for mask in range(64):
    unit, inverse = ONE, ONE
    for bit in range(6):
        if mask & (1 << bit):
            unit = L.mul(unit, generators[bit])
            inverse = L.mul(inverse, generator_inverses[bit])
    assert L.mul(unit, inverse) == ONE
    units.append(unit)
    inverse_units.append(inverse)
assert len(set(units)) == 64
assert len({sign_vector(a) for a in units[:32]}) == 32
for mask in range(64):
    assert units[mask ^ 16] == NEG(units[mask])

z = M.x
alpha0 = M.add(M.add(M.one, z), M.pow(z, 4))
alpha = [relative(M.mul(M.pow(z, k), alpha0)) for k in range(42)]
PI = L.reduce((5, -2, -4, 1, 1))
GAMMA = SUB(T, ONE)
STRATA = {
    "OM": (None, ONE, 0),
    "R7": ((7, 1, 1), alpha[0][1], 0),
    "R7sq": ((7, 1, 2), alpha[12][1], 12),
    "R43_20": ((43, 20, 1), alpha[4][1], 4),
    "R43_31": ((43, 31, 1), alpha[8][1], 8),
}
by_conductor = {conductor: name for name, (conductor, _, _) in STRATA.items()}
for group in coverage["embedded_order_groups"]:
    key = tuple(group["conductor_key"]) if group["conductor_key"] else None
    if key not in by_conductor or key is None:
        continue
    reference_a, reference_b = alpha[group["representative"]]
    assert reference_b == STRATA[by_conductor[key]][1]
    for witness in group["affine_unit_witnesses"]:
        a, b = alpha[witness["trace_index"]]
        e, inv = tuple(witness["unit_ratio"]), tuple(witness["inverse_unit"])
        assert L.mul(e, inv) == ONE
        assert b == L.mul(e, reference_b)
        assert SUB(a, L.mul(e, reference_a)) == tuple(witness["affine_shift"])


def action(k, g, mask):
    a, b = alpha[k]
    q = L.mul(g, units[mask])
    result = (a, NEG(L.mul(q, b)),
              L.mul(L.divide_integrally(b, g), inverse_units[mask]),
              L.add(a, L.mul(T, b)))
    assert det(result) == PI
    return result


def invariants(q):
    square = L.mul(q, q)
    return {"twice_real_part": NEG(L.mul(q, T)),
            "four_times_height_squared": L.mul(square, U),
            "modulus_squared": square}


def label(name, mask):
    return f"{name}:{mask:02d}"


points = []
for name, (conductor, g, k) in STRATA.items():
    masks = range(16) if name == "OM" else list(range(16))+list(range(32, 48))
    for mask in masks:
        q = L.mul(g, units[mask])
        matrix = action(k, g, mask)
        # A unique upper fixed tuple: x=-qt/2, h^2=q^2(4-t^2)/4.
        inv = invariants(q)
        translation = (ONE, L.mul(q, T), ZERO, ONE)
        opposite = action(k, g, mask ^ 16)
        assert conjugate(translation, other_root(matrix)) == opposite
        assert inv["twice_real_part"] != invariants(NEG(q))["twice_real_part"]
        points.append({"point_label": label(name, mask), "stratum": name,
                       "multiplier_conductor": conductor, "canonical_mask": mask,
                       "root_pair_masks": [mask, mask ^ 16], "g": g, "q": q,
                       "exact_tuple": inv, "trace_index_for_matrix": k,
                       "norm_pi_matrix_row_major": matrix,
                       "root_translation_row_major": translation,
                       "same_existing_v45_label": label(name, mask) if name == "OM" else None})
assert len(points) == len({tuple(p["q"]) for p in points}) == 144
assert len({tuple(p["exact_tuple"]["twice_real_part"]) for p in points}) == 144
points_by_label = {p["point_label"]: p for p in points}

# e = delta_a * v^2 for each nontrivial old-to-canonical basis determinant.
V15 = L.reduce((1, 5, 0, -5, 0, 1))
V18 = L.reduce((3, 13, -1, -11, 0, 2))
VINTER = L.reduce((1, 18, 3, -16, -1, 3))
square_witnesses = {0: (0, ONE), 15: (33, V15), 18: (39, V18),
                    21: (16, ONE), 36: (49, V15), 39: (55, V18),
                    4: (0, ONE), 25: (16, ONE),
                    8: (0, ONE), 29: (16, ONE),
                    12: (0, ONE), 33: (16, ONE)}
unit_witness_records, bridges = {}, []
incidences, literal_count, seen = Counter(), 0, set()
for row_number, row in enumerate(input_rows, 1):
    k, mask = row["trace_index"], row["orientation_mask"]
    key = None if row["multiplier"] == "OM" else tuple(row["multiplier"])
    name = by_conductor[key]
    _, g, _ = STRATA[name]
    old = tuple(tuple(a) for a in row["matrix_row_major_L_coefficients"])
    assert tuple(row["unit"]) == units[mask]
    assert tuple(row["unit_signs"]) == sign_vector(units[mask])
    assert (k, name, mask) not in seen
    seen.add((k, name, mask))
    if name == "OM":
        assert mask < 32
        Q, e, amask, v = I, ONE, 0, ONE
        assert old == action(k, ONE, mask)
        witness_name = "OM_identity"
    elif tuple(row["conductor"]) == key:
        a, b = alpha[k]
        e = L.divide_integrally(b, g)
        amask, v = square_witnesses[k]
        Q = (ONE, a, ZERO, e)
        trace = L.add(L.scale(a, 2), L.mul(T, b))
        expected_old = (ZERO, NEG(L.mul(units[mask], PI)), inverse_units[mask], trace)
        assert old == expected_old
        witness_name = f"regular_trace_{k}"
    else:
        assert k in (12, 33) and name == "R7"
        e = L.divide_integrally(GAMMA, g)
        amask, v = 7, VINTER
        Q = diag(ONE, e)
        assert old == action(k, GAMMA, mask)
        witness_name = "intermediate_gamma_to_B0"
    assert e == L.mul(units[amask], L.mul(v, v))
    assert abs(L.norm(e)) == abs(L.norm(v)) == 1
    unit_witness_records[witness_name] = {
        "basis_matrix_Q": Q, "determinant_unit_e": e,
        "unit_square_mask": amask, "square_root_v": v,
        "identity": "e = delta[unit_square_mask] * square_root_v^2"}
    target_mask = mask ^ amask
    h = ONE
    for bit in range(6):
        if (amask & (1 << bit)) and not (mask & (1 << bit)):
            h = L.mul(h, generators[bit])
    scalar = L.divide_integrally(ONE, L.mul(v, h))
    P0 = mm(mm(diag(units[target_mask], ONE), Q), diag(inverse_units[mask], ONE))
    P = tuple(L.mul(scalar, value) for value in P0)
    assert det(P) == ONE
    assert conjugate(P, old) == action(k, g, target_mask)
    canonical_mask = target_mask & ~16
    root_swapped = bool(target_mask & 16)
    if root_swapped:
        q = L.mul(g, units[target_mask])
        root_translation = (ONE, L.mul(q, T), ZERO, ONE)
        G = mm(root_translation, P)
    else:
        G = P
    target = action(k, g, canonical_mask)
    assert conjugate(G, old) == (other_root(target) if root_swapped else target)
    assert G[2] == ZERO and det(G) == ONE
    # Exact tuple transport, independently of inferring it from conjugacy.
    # 2x=(a-d)/c; 4h^2=(4pi-trace^2)/c^2 for the old action.
    # These quotients are integral for our templates, even when c is not a unit.
    old_x2 = L.divide_integrally(SUB(old[0], old[3]), old[2])
    old_trace = L.add(old[0], old[3])
    old_h4 = L.divide_integrally(SUB(L.scale(PI, 4), L.mul(old_trace, old_trace)), L.mul(old[2], old[2]))
    # G=[[a,b],[0,d]], ad=1: 2x'=a^2*2x+2ab, 4h'^2=a^4*4h^2.
    aa = L.mul(G[0], G[0])
    target_label = label(name, canonical_mask)
    target_inv = points_by_label[target_label]["exact_tuple"]
    assert L.add(L.mul(aa, old_x2), L.scale(L.mul(G[0], G[1]), 2)) == target_inv["twice_real_part"]
    assert L.mul(L.mul(aa, aa), old_h4) == target_inv["four_times_height_squared"]
    literal_equal = (old_x2 == target_inv["twice_real_part"] and
                     old_h4 == target_inv["four_times_height_squared"])
    literal_count += literal_equal
    incidences[target_label] += 1
    bridges.append({"input_line": row_number, "trace_index": k,
                    "old_orientation_mask": mask, "stratum": name,
                    "basis_unit_witness": witness_name,
                    "transported_mask_before_root_pair": target_mask,
                    "canonical_point_label": target_label,
                    "old_chosen_tuple_equals_canonical_tuple": literal_equal,
                    "old_twice_real_part": old_x2,
                    "old_four_times_height_squared": old_h4,
                    "SL2_transport_row_major": G,
                    "conjugates_to_other_root_of_canonical_action": root_swapped})

assert len(bridges) == len(seen) == 1280
for point in points:
    expected_incidence = 24 if point["stratum"] == "OM" else 16 if point["stratum"] == "R7" else 4
    assert incidences[point["point_label"]] == expected_incidence
point_counts = dict(Counter(p["stratum"] for p in points))
assert point_counts == {"OM": 16, "R7": 32, "R7sq": 32, "R43_20": 32, "R43_31": 32}


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))+"\n"


point_bytes = (json.dumps(points, indent=2, sort_keys=True)+"\n").encode()
bridge_bytes = "".join(encoded(row) for row in bridges).encode()
witness_bytes = (json.dumps(unit_witness_records, indent=2, sort_keys=True)+"\n").encode()
artifacts = {"canonical_points_v59.json": point_bytes,
             "matrix_to_point_bridges_v59.jsonl": bridge_bytes,
             "basis_unit_witnesses_v59.json": witness_bytes}
for filename, data in artifacts.items():
    assert (OUT/filename).read_bytes() == data, filename
summary = {
    "status": "PASS_EXACT_144_D3_POINTS_AND_1280_SL2_BRIDGES",
    "scope": "144 SL2(O_L) quotient points using C37 and the reviewed two-root/multiplier proof; no curve decisions",
    "curve_search_performed": False,
    "imports_foundation_arithmetic_not_an_independent_foundation_check": True,
    "canonical_chosen_tuples_pairwise_distinct": True,
    "counts_by_multiplier": point_counts,
    "canonical_point_count": 144,
    "old_trace_labelled_matrix_count": len(bridges),
    "old_representatives_literally_equal_to_canonical_tuple": literal_count,
    "old_representatives_transported_by_verified_SL2_matrix": len(bridges),
    "exact_existing_d7_population_overlap": {"stratum": "OM", "count": 16},
    "source_sha256": PINS,
    "implementation_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "artifact_sha256": {filename: sha256(data).hexdigest() for filename, data in artifacts.items()},
    "ordering_of_real_embeddings": [1, 5, 11, 13, 17, 19],
    "invariants": {"twice_real_part": "2x", "four_times_height_squared": "4h^2", "modulus_squared": "|w|^2"},
    "cpu_seconds": round(time.process_time()-START_CPU, 6),
}
assert summary["cpu_seconds"] < 25
for filename, expected in PINS.items():
    assert sha256((ROOT/filename).read_bytes()).hexdigest() == expected, filename
print(json.dumps(summary, indent=2, sort_keys=True))
