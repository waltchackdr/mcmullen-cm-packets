#!/usr/bin/env python3
"""Independently verify the d3 ambient-point representatives and SL2 bridges.

Definitions and old matrix classes come from the retained independent
quadratic model.  JSON files provide only proposed points and witnesses.
The two-root and intrinsic-multiplier arguments remain the reviewed human
proof; no curve membership, field-automorphism quotient or geometric CM
signature is certified here.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import cache
from hashlib import sha256
import json
from pathlib import Path

import d3_completeness_quadratic_v60 as foundation

if not __debug__:
    raise RuntimeError("Assertions must be enabled; do not use -O/-OO")

F = foundation.field
L, T, ONE, ZERO = F.LElement, F.T, F.ONE, F.ZERO
PI = foundation.PI
COEFF = foundation.coefficients
DEFAULT_DATA = Path(__file__).resolve().parents[1] / "certificates/d3_points_v60"


@cache
def inverse(value):
    return value.inverse()


def matrix(coefficients):
    assert len(coefficients) == 4 and all(len(row) == 6 for row in coefficients)
    answer = F.matrix(*(L(row) for row in coefficients))
    assert F.all_entries_integral(answer)
    return answer


def flat(value):
    return [COEFF(entry) for row in value for entry in row]


def other_root(value):
    return F.matrix(value[1][1], -value[0][1], -value[1][0], value[0][0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-directory", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()

    # No proposed point/bridge data have been opened at this stage.
    old_bytes, arithmetic = foundation.reconstruct()
    old_rows = [json.loads(line) for line in old_bytes.splitlines()]
    alpha, units = {}, {}
    for row in old_rows:
        k, mask = row["trace_index"], row["orientation_mask"]
        units[mask] = L(row["unit"])
        if row["multiplier"] == "OM" and mask == 0:
            entries = matrix(row["matrix_row_major_L_coefficients"])
            alpha[k] = (entries[0][0], entries[1][0])
    assert set(units) == set(range(64))
    strata = {"OM": (None, ONE, 0),
              "R7": ([7, 1, 1], alpha[0][1], 0),
              "R7sq": ([7, 1, 2], alpha[12][1], 12),
              "R43_20": ([43, 20, 1], alpha[4][1], 4),
              "R43_31": ([43, 31, 1], alpha[8][1], 8)}

    def action(k, g, mask):
        a, b = alpha[k]
        q = g*units[mask]
        value = F.matrix(a, -q*b, b*inverse(q), a+T*b)
        assert F.mdet(value) == PI and F.all_entries_integral(value)
        return value

    predicted = {}
    for name, (conductor, g, k) in strata.items():
        masks = range(16) if name == "OM" else tuple(range(16))+tuple(range(32, 48))
        for mask in masks:
            q = g*units[mask]
            label = f"{name}:{mask:02d}"
            translate = F.matrix(1, q*T, 0, 1)
            target = action(k, g, mask)
            assert F.mdet(translate) == ONE
            assert F.mmul(translate, other_root(target)) == F.mmul(action(k, g, mask^16), translate)
            predicted[label] = dict(point_label=label, stratum=name,
                multiplier_conductor=conductor, canonical_mask=mask,
                root_pair_masks=[mask, mask^16], g=COEFF(g), q=COEFF(q),
                exact_tuple=dict(twice_real_part=COEFF(-q*T),
                                 four_times_height_squared=COEFF(q*q*(4-T*T)),
                                 modulus_squared=COEFF(q*q)),
                trace_index_for_matrix=k, norm_pi_matrix_row_major=flat(target),
                root_translation_row_major=flat(translate),
                same_existing_v45_label=label if name == "OM" else None)
    assert len(predicted) == 144
    assert len({tuple(p["q"]) for p in predicted.values()}) == 144

    paths = {name: args.data_directory/name for name in
             ("canonical_points_v59.json", "basis_unit_witnesses_v59.json",
              "matrix_to_point_bridges_v59.jsonl")}
    raw = {name: path.read_bytes() for name, path in paths.items()}
    points = json.loads(raw["canonical_points_v59.json"])
    assert len(points) == len(predicted)
    assert {p["point_label"] for p in points} == set(predicted)
    for point in points:
        assert point == predicted[point["point_label"]]
    witnesses = json.loads(raw["basis_unit_witnesses_v59.json"])
    bridges = [json.loads(line) for line in raw["matrix_to_point_bridges_v59.jsonl"].splitlines()]
    assert len(bridges) == len(old_rows) == 1280

    unit_data = {}
    for name, witness in witnesses.items():
        q = matrix(witness["basis_matrix_Q"])
        e, v = L(witness["determinant_unit_e"]), L(witness["square_root_v"])
        amask = witness["unit_square_mask"]
        assert amask in units
        assert e == units[amask]*v*v == F.mdet(q)
        assert e.integral_power_basis() and inverse(e).integral_power_basis()
        assert v.integral_power_basis() and inverse(v).integral_power_basis()
        unit_data[name] = q, e, amask

    by_conductor = {tuple(key) if key else None: name for name, (key, _, _) in strata.items()}
    used_witnesses, incidences, literal_count, seen = set(), Counter(), 0, set()
    for line, (old_row, bridge) in enumerate(zip(old_rows, bridges), 1):
        k, mask = old_row["trace_index"], old_row["orientation_mask"]
        key = None if old_row["multiplier"] == "OM" else tuple(old_row["multiplier"])
        name = by_conductor[key]
        _, g, _ = strata[name]
        assert (k, name, mask) not in seen
        seen.add((k, name, mask))
        assert bridge["input_line"] == line and bridge["trace_index"] == k
        assert bridge["old_orientation_mask"] == mask and bridge["stratum"] == name
        a, b = alpha[k]
        if name == "OM":
            witness_name = "OM_identity"
            expected_q = F.IDENTITY
        elif tuple(old_row["conductor"]) == key:
            witness_name = f"regular_trace_{k}"
            expected_q = F.matrix(1, a, 0, b*inverse(g))
        else:
            assert name == "R7" and k in (12, 33)
            witness_name = "intermediate_gamma_to_B0"
            expected_q = F.matrix(1, 0, 0, (T-1)*inverse(g))
        assert bridge["basis_unit_witness"] == witness_name
        used_witnesses.add(witness_name)
        q, e, amask = unit_data[witness_name]
        assert q == expected_q
        target_mask = mask ^ amask
        canonical_mask = target_mask & ~16
        target_label = f"{name}:{canonical_mask:02d}"
        root_swap = bool(target_mask & 16)
        assert bridge["transported_mask_before_root_pair"] == target_mask
        assert bridge["canonical_point_label"] == target_label
        assert bridge["conjugates_to_other_root_of_canonical_action"] == root_swap

        old = matrix(old_row["matrix_row_major_L_coefficients"])
        transport = matrix(bridge["SL2_transport_row_major"])
        assert F.mdet(transport) == ONE
        # Intertwining avoids computing an inverse or constructing the
        # primary algorithm's square correction. The supplied witness must
        # satisfy the exact identity in a separate field implementation.
        target = action(k, g, canonical_mask)
        target_root = other_root(target) if root_swap else target
        assert F.mmul(transport, old) == F.mmul(target_root, transport)
        assert transport[1][0] == ZERO
        x2 = (old[0][0]-old[1][1])*inverse(old[1][0])
        trace = old[0][0]+old[1][1]
        h4 = (4*PI-trace*trace)*inverse(old[1][0])**2
        assert COEFF(x2) == bridge["old_twice_real_part"]
        assert COEFF(h4) == bridge["old_four_times_height_squared"]
        destination = predicted[target_label]["exact_tuple"]
        aa, bb = transport[0]
        assert aa**2*x2+2*aa*bb == L(destination["twice_real_part"])
        assert aa**4*h4 == L(destination["four_times_height_squared"])
        literal = x2 == L(destination["twice_real_part"]) and h4 == L(destination["four_times_height_squared"])
        assert bridge["old_chosen_tuple_equals_canonical_tuple"] == literal
        literal_count += literal
        incidences[target_label] += 1
    assert used_witnesses == set(witnesses)
    for label, point in predicted.items():
        expected_count = 24 if point["stratum"] == "OM" else 16 if point["stratum"] == "R7" else 4
        assert incidences[label] == expected_count
    assert literal_count == 200
    print(json.dumps(dict(status="PASS_INDEPENDENT_D3_POINT_BRIDGES",
        point_count=144, matrix_classes=1280, literal_tuple_equalities=literal_count,
        counts_by_multiplier=dict(Counter(p["stratum"] for p in predicted.values())),
        matrix_manifest_sha256=sha256(old_bytes).hexdigest(),
        artifact_sha256={name: sha256(value).hexdigest() for name, value in raw.items()},
        implementation_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        foundation_sha256=sha256(Path(foundation.__file__).read_bytes()).hexdigest(),
        scope="Exact canonical representatives, unit-square witnesses and 1280 integral determinant-one bridges; quotient distinctness uses the reviewed intrinsic-multiplier/two-root proof; no curve geometry"),
        sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
