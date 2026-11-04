#!/usr/bin/env python3
"""Portable independent candidate certificate for 32 proper-(41,t-12) points.

Pending replay/review/registration. Uses the local hash-checked separate
standard-library field/Sturm/trace enumerator; no primary arithmetic or
experiment import. Only explicit witnesses and expected output constants
are shared with the primary path. No build-directory dependency.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"certificates"))
import d7_r41_expected_v63 as frozen


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def vectors_digest(vectors):
    h = sha256()
    for vector in sorted(vectors):
        h.update((",".join(map(str, vector))+"\n").encode("ascii"))
    return h.hexdigest()


def candidate_digest(candidates):
    return vectors_digest(tuple(c for entry in gamma for c in entry) for gamma in candidates)


def selected_masks(value):
    if value == "all":
        return frozen.MASKS
    try:
        values = tuple(int(item) for item in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use all or comma-separated proper-R41 masks") from exc
    if not values or len(values) != len(set(values)) or not set(values) <= set(frozen.MASKS):
        raise argparse.ArgumentTypeError(f"Select distinct masks from {frozen.MASKS}")
    return tuple(sorted(values))


def sources_and_domain():
    assert frozen.MASKS == tuple(range(16))+tuple(range(32, 48))
    root_pairs = [(m, m ^ 16) for m in frozen.MASKS]
    assert sorted(v for pair in root_pairs for v in pair) == list(range(64))
    batch = sha256()
    for mask in frozen.MASKS:
        expected = frozen.EXPECTED[mask]
        assert expected["survivor_count"] == 0
        assert sum(n for _, n in expected["polygon_histogram"]) == expected["candidate_count"]
        assert expected["candidate_count"] == expected["nonzero_c_count"]+expected["zero_c_count"]
        assert expected["zero_c_audit"][-1] == expected["zero_c_count"]
        batch.update(f'{mask:02d}:{expected["classification_digest"]}\n'.encode("ascii"))
    assert batch.hexdigest() == frozen.EXPECTED_BATCH_DIGEST
    for data in (frozen.BALANCING_EXPONENTS, frozen.EXPECTED_BOUNDS, frozen.EXPECTED):
        assert tuple(sorted(data)) == frozen.MASKS
    expected_sources = frozen.FOUNDATION_SHA256["sturm"]
    for relative, expected in expected_sources.items():
        assert digest(ROOT/relative) == expected, relative
    assert Path(frozen.__file__).resolve() == ROOT/"certificates/d7_r41_expected_v63.py"
    return {**expected_sources,
            "certificates/d7_r41_expected_v63.py": digest(Path(frozen.__file__)),
            "independent/d7_r41_sturm_v63.py": digest(Path(__file__))}


def reconstruct(retained, mask):
    field = retained.ambient
    t, one = field.T, field.ONE
    # Independent quadratic-pair arithmetic reconstructs alpha_14.
    a, b = field.alpha_pair(14)
    assert field.pair_norm((a, b)) == field.PI
    trace = field.pair_trace((a, b))
    assert trace == t**5-7*t**3-t*t+11*t+2
    assert abs(b.norm()) == 41
    assert field.element_modulus_value(b, 12, 41) == 0
    assert field.element_modulus_value(b, 22, 41) != 0
    base = field.matrix(0, -field.PI, 1, trace)
    basis = field.matrix(1, a, 0, b)
    assert field.mmul(basis, base) == field.mmul(field.multiplication_matrix((a, b)), basis)
    z_action = field.matrix(0, -1, 1, t)
    assert not field.all_entries_integral(field.conjugate_by(field.minverse(basis), z_action))
    xi = (t**4-3*t*t+1, t**3+t*t-2*t-1,
          t**5-4*t**3+3*t+1, t**3-2*t, t)
    orientation = xi[:4]+(-one, 4-t*t)
    delta, unit = one, one
    for bit, generator in enumerate(orientation):
        if mask & (1 << bit):
            delta *= generator
    for generator, exponent in zip(xi, frozen.BALANCING_EXPONENTS[mask]):
        unit *= generator**exponent
    for value in (delta, unit):
        assert value.integral_power_basis() and value.inverse().integral_power_basis()
    raw = field.twist(base, delta)
    opposite = field.twist(base, -delta)
    bridge = field.matrix(-1, -delta*trace, 0, -1)
    assert field.mdet(bridge) == one and field.all_entries_integral(bridge)
    conjugate_action = field.matrix(trace-raw[0][0], -raw[0][1], -raw[1][0], trace-raw[1][1])
    assert field.mmul(bridge, conjugate_action) == field.mmul(opposite, bridge)
    balancing = field.matrix(unit, 0, 0, unit.inverse())
    assert field.mdet(balancing) == one
    assert field.all_entries_integral(balancing) and field.all_entries_integral(field.minverse(balancing))
    balanced = field.conjugate_by(balancing, raw)
    for matrix in (raw, opposite, balanced):
        assert field.all_entries_integral(matrix)
        assert (field.mtrace(matrix), field.mdet(matrix)) == (trace, field.PI)
    center, height = field.exact_fixed_invariants(balanced)
    modulus = center**2+height
    eta = delta*unit**2
    assert center == -eta*trace/2
    assert height == eta**2*(4*field.PI-trace**2)/4
    assert modulus == eta**2*field.PI
    roots = field.embedded_roots()
    assert field.sign_vector(height, roots) == field.sign_vector(modulus, roots) == (1,)*6
    point = retained.Point(balanced, center, height, modulus)
    bounds = retained.derive_bounds(point)
    assert bounds == {key: tuple(map(Q, values)) for key, values in frozen.EXPECTED_BOUNDS[mask].items()}
    return point, bounds


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O/-OO")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--masks", type=selected_masks, default=frozen.MASKS)
    parser.add_argument("--check-sources", action="store_true", help="No arithmetic import or enumeration")
    parser.add_argument("--bounds-only", action="store_true", help="Reconstruct points/bounds only; no exclusion pass")
    args = parser.parse_args()
    sources = sources_and_domain()
    if args.check_sources:
        print("PASS_SOURCE_CHECK_ONLY " + json.dumps(sources, sort_keys=True))
        return
    # The retained independent module imports only the legacy constant table
    # from certificates/, never primary arithmetic or search code.
    import d7_r83_principal_sturm_v55 as retained
    for relative in frozen.FOUNDATION_SHA256["sturm"]:
        module = sys.modules[Path(relative).stem]
        assert Path(module.__file__).resolve() == (ROOT/relative).resolve()
    completed = {}
    for mask in args.masks:
        point, bounds = reconstruct(retained, mask)
        if args.bounds_only:
            print(f"BOUNDS_ONLY R41:principal:{mask:02d}", flush=True)
            continue
        print(f"START R41:principal:{mask:02d} entry_lists", flush=True)
        a_values, c_values = retained.entry_lists(point, bounds)
        bounded = retained.enumerate_candidates(mask, point, a_values, c_values)
        candidates = bounded.nonzero+bounded.zero
        histogram, classification_digest, survivors = retained.polygon_result(point, candidates)
        assert not survivors
        observed = {
            "a_count": len(a_values), "a_digest": vectors_digest(a_values),
            "c_count": len(c_values), "c_digest": vectors_digest(c_values),
            "bottom_rows": bounded.bottom_rows,
            "nonzero_c_count": len(bounded.nonzero), "zero_c_count": len(bounded.zero),
            "candidate_count": len(candidates), "candidate_digest": candidate_digest(candidates),
            "classification_digest": classification_digest,
            "polygon_histogram": [[list(key), n] for key, n in histogram],
            "zero_c_audit": [bounded.zero_bottom_rows, bounded.zero_d_rows,
                bounded.integral_inverse_rows, bounded.a_admissible_rows, len(bounded.zero)],
            "survivor_count": len(survivors), "survivor_digest": candidate_digest(survivors),
        }
        assert observed == frozen.EXPECTED[mask], (mask, observed)
        completed[mask] = observed
        print("POINT_PASS " + json.dumps({"mask": mask, **observed}, sort_keys=True), flush=True)
    for relative, expected in sources.items():
        assert digest(ROOT/relative) == expected, ("Source changed", relative)
    if args.bounds_only:
        print("PASS_BOUNDS_ONLY_NO_ORBIT_DECISION")
        return
    if args.masks == frozen.MASKS:
        batch = sha256(''.join(f'{m:02d}:{completed[m]["classification_digest"]}\n' for m in args.masks).encode("ascii")).hexdigest()
        assert batch == frozen.EXPECTED_BATCH_DIGEST
        print(f"BATCH_PASS thirty_two_proper_R41_points digest={batch}")
    else:
        print("SUBSET_PASS masks="+",".join(map(str, args.masks)))
    print("SOURCE_SHA256 " + json.dumps(sources, sort_keys=True))


if __name__ == "__main__":
    main()
