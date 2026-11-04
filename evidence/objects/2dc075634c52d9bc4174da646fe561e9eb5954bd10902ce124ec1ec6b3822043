#!/usr/bin/env python3
"""Portable primary certificate candidate for the 32 proper-(41,t-12) points.

Pending replay/review/registration: source creation is not theorem promotion.
Uses only local, hash-checked retained primary foundations. The enumeration
is the reviewed v59 primary path, with frozen expectations and no build-tree
or experiment import. See docs/D7_R41_COMPLETION_V63.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
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


def source_manifest():
    expected = frozen.FOUNDATION_SHA256["primary"]
    assert Path(frozen.__file__).resolve() == ROOT/"certificates/d7_r41_expected_v63.py"
    for relative, value in expected.items():
        assert digest(ROOT/relative) == value, relative
    return {**expected,
            "certificates/d7_r41_expected_v63.py": digest(Path(frozen.__file__)),
            "certificates/d7_r41_v63.py": digest(Path(__file__))}


def masks_argument(value):
    if value == "all":
        return frozen.MASKS
    try:
        values = tuple(int(item) for item in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Use all or comma-separated proper-R41 masks") from exc
    if not values or len(set(values)) != len(values) or any(m not in frozen.MASKS for m in values):
        raise argparse.ArgumentTypeError(f"Select distinct masks from {frozen.MASKS}")
    return tuple(sorted(values))


def check_domain():
    assert frozen.MASKS == tuple(range(16))+tuple(range(32, 48))
    fibers = {frozenset((m, m ^ 16)) for m in frozen.MASKS}
    assert len(fibers) == 32 and set().union(*fibers) == set(range(64))
    for mapping in (frozen.EXPECTED, frozen.EXPECTED_BOUNDS, frozen.BALANCING_EXPONENTS):
        assert set(mapping) == set(frozen.MASKS)
    batch = sha256()
    for mask in frozen.MASKS:
        expected = frozen.EXPECTED[mask]
        assert expected["survivor_count"] == 0
        assert sum(n for _, n in expected["polygon_histogram"]) == expected["candidate_count"]
        assert expected["candidate_count"] == expected["nonzero_c_count"]+expected["zero_c_count"]
        assert expected["zero_c_audit"][-1] == expected["zero_c_count"]
        batch.update(f'{mask:02d}:{expected["classification_digest"]}\n'.encode("ascii"))
    assert batch.hexdigest() == frozen.EXPECTED_BATCH_DIGEST


def reconstruct(primary, mask):
    a, sp = primary.ambient, primary.sp
    t = a.T
    aa, bb = a.alpha_pair(14)
    trace = a.alpha_trace(14)
    assert a.reduce_l(trace-(t**5-7*t**3-t*t+11*t+2)) == 0
    assert a.reduce_l(aa*aa+t*aa*bb+bb*bb-a.PI) == 0
    # The norm and residue fix the embedded conductor (41,t-12).
    assert abs(a.absolute_norm(bb)) == 41
    assert int(bb.subs(t, 12)) % 41 == 0
    assert int(bb.subs(t, 22)) % 41 != 0
    assert ("R41:principal", mask) in a.POINT_RECORDS
    base = sp.Matrix([[0, -a.PI], [1, trace]])
    companion_basis = sp.Matrix([[1, aa], [0, bb]])
    maximal_action = a.multiplication_matrix((aa, bb))
    assert a.matrix_reduce(companion_basis*base-maximal_action*companion_basis) == sp.zeros(2)
    z_action = sp.Matrix([[0, -1], [1, t]])
    z_in_basis = a.matrix_reduce(companion_basis.inv()*z_action*companion_basis)
    assert not all(a.is_integral(v) for v in z_in_basis)
    delta = a.SQUARE_CLASS_UNITS[mask]
    raw = a.conjugate_by_unit(base, delta)
    opposite = a.conjugate_by_unit(base, -delta)
    # Compare the conjugate root action, not raw, with mask m xor16.
    root_bridge = sp.Matrix([[-1, -delta*trace], [0, -1]])
    assert a.reduce_l(root_bridge.det()-1) == 0
    a.assert_integral_matrix(root_bridge)
    assert a.matrix_reduce(root_bridge*(trace*sp.eye(2)-raw)-opposite*root_bridge) == sp.zeros(2)
    unit = sp.Integer(1)
    for generator, exponent in zip(a.XI_UNITS, frozen.BALANCING_EXPONENTS[mask]):
        unit = a.reduce_l(unit*primary.field_power(generator, exponent))
    inverse = a.reduce_field(1/unit)
    a.assert_integral(unit)
    a.assert_integral(inverse)
    balancing = sp.diag(unit, inverse)
    assert a.reduce_l(balancing.det()-1) == 0
    balanced = a.matrix_multiply(a.matrix_multiply(balancing, raw), sp.diag(inverse, unit))
    for matrix in (raw, opposite, balanced):
        a.assert_integral_matrix(matrix)
        a.assert_characteristic(matrix, 14)
    c = balanced[1, 0]
    x = a.reduce_field((balanced[0, 0]-balanced[1, 1])/(2*c))
    height = a.reduce_field((4*a.PI-trace*trace)/(4*c*c))
    modulus = a.reduce_l(x*x+height)
    eta = a.reduce_l(delta*unit*unit)
    assert a.reduce_l(x+eta*trace/2) == 0
    assert a.reduce_l(height-eta*eta*(4*a.PI-trace*trace)/4) == 0
    assert a.reduce_l(modulus-eta*eta*a.PI) == 0
    assert a.sign_vector(height) == a.sign_vector(modulus) == (1,)*6
    center, height, modulus = map(primary.to_core, (x, height, modulus))
    bounds = primary.derive_bounds(height, modulus)
    assert bounds == tuple(tuple(map(Q, frozen.EXPECTED_BOUNDS[mask][k]))
                           for k in ("DEN", "NUM", "A", "B", "C", "D"))
    primary.enumeration.install_point(height, bounds)
    return center, height, bounds


class Progress:
    def __init__(self, mask):
        self.mask, self.stage = mask, None

    def report(self, stage, **counts):
        completed, total = counts.get("completed", 0), counts.get("total", 0)
        if stage != self.stage or completed == total or completed % 250 == 0:
            self.stage = stage
            print(f"PROGRESS R41:principal:{self.mask:02d} {stage} {json.dumps(counts, sort_keys=True)}", flush=True)


# Same primary arithmetic/enumeration path as the reviewed v59 driver.
def enumerate_candidates(primary, x, a_values, c_values, progress):
    core = primary.core
    den_squares = core.square_bounds(core.DEN_BOUND)
    num_squares = core.square_bounds(core.NUM_BOUND)
    den_radius, num_radius = sum(den_squares, Q(0)), sum(num_squares, Q(0))
    bottom_rows = []
    x_boxes = tuple(core.scaled_field_box(x, i) for i in range(6))
    y2_boxes = tuple(core.scaled_field_box(core.Y2, i) for i in range(6))
    a_boxes = {a: tuple(core.scaled_eval_integral(a, i) for i in range(6))
               for a in a_values}

    def bounded(v, w, boxes, squares, exact):
        unresolved = False
        for i in range(6):
            real = core.scaled_add(core.scaled_mul(boxes[i], x_boxes[i]),
                                   core.scaled_eval_integral(w, i))
            sq = core.scaled_add(core.scaled_mul(real, real),
                core.scaled_mul(core.scaled_mul(boxes[i], boxes[i]), y2_boxes[i]))
            decision = core.scaled_le(sq, squares[i])
            if decision is False:
                return False
            unresolved |= decision is None
        if not unresolved:
            return True
        value = exact(core.as_field(v), core.as_field(w), x)
        return all(core.sign_at(core.sub(core.scale(bound, core.ONE), value), i) >= 0
                   for i, bound in enumerate(squares))

    def numerator(a, b):
        return bounded(a, b, a_boxes[a], num_squares, core.numerator_square)

    for n, c in enumerate(c_values, 1):
        cq = core.as_field(c)
        boxes = tuple(core.scaled_eval_integral(c, i) for i in range(6))
        values = core.enumerate_trace_ellipsoid(core.neg(core.mul(cq, x)), den_radius,
            lambda d: bounded(c, d, boxes, den_squares, core.denominator_square))
        bottom_rows.extend((c, d) for d in values)
        progress.report("bottom_rows", completed=n, total=len(c_values),
                        bottom_rows=len(bottom_rows))
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != core.ZERO]
    zero_rows = [(c, d) for c, d in bottom_rows if c == core.ZERO]
    candidates = []
    for n, (c, d) in enumerate(nonzero_rows, 1):
        affine = core.affine_b_data(c, d)
        for a in a_values:
            b = core.affine_integral_b_fast(a, affine)
            if b is None or not numerator(a, b):
                continue
            assert core.sub(core.mul(core.as_field(a), core.as_field(d)),
                            core.mul(core.as_field(b), core.as_field(c))) == core.ONE
            candidates.append((a, b, c, d))
        progress.report("determinant_nonzero_c", completed=n,
                        total=len(nonzero_rows), candidates=len(candidates))
    nonzero_count = len(candidates)
    zero_audit = Counter()
    for n, (c, d) in enumerate(zero_rows, 1):
        if d == core.ZERO:
            zero_audit["d_zero"] += 1
        else:
            dq = core.as_field(d)
            aq = core.inverse(dq)
            if not core.integral(aq):
                zero_audit["nonintegral_inverse"] += 1
            else:
                a = core.integral_tuple(aq)
                if not primary.within_box(a, core.A_BOUND):
                    zero_audit["a_outside_box"] += 1
                else:
                    assert a in a_boxes
                    zero_audit["unit_rows"] += 1
                    bs = core.enumerate_trace_ellipsoid(core.neg(core.mul(aq, x)),
                                                        num_radius, lambda b: numerator(a, b))
                    assert core.mul(aq, dq) == core.ONE
                    candidates.extend((a, b, c, d) for b in bs)
        progress.report("determinant_zero_c", completed=n, total=len(zero_rows),
                        candidates=len(candidates)-nonzero_count)
    assert len(candidates) == len(set(candidates))
    return candidates, dict(bottom_rows=len(bottom_rows), nonzero_candidates=nonzero_count,
        zero_candidates=len(candidates)-nonzero_count, zero_rows=len(zero_rows),
        zero_row_audit=dict(zero_audit))


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O/-OO")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--masks", type=masks_argument, default=frozen.MASKS)
    parser.add_argument("--check-sources", action="store_true", help="No arithmetic import or enumeration")
    parser.add_argument("--bounds-only", action="store_true", help="Reconstruct points/bounds only; no exclusion pass")
    args = parser.parse_args()
    check_domain()
    sources = source_manifest()
    if args.check_sources:
        print("PASS_SOURCE_CHECK_ONLY " + json.dumps(sources, sort_keys=True))
        return
    import d7_r83_principal_v55 as primary
    for relative in frozen.FOUNDATION_SHA256["primary"]:
        module = sys.modules[Path(relative).stem]
        assert Path(module.__file__).resolve() == (ROOT/relative).resolve()
    completed = {}
    for mask in args.masks:
        center, height, bounds = reconstruct(primary, mask)
        if args.bounds_only:
            print(f"BOUNDS_ONLY R41:principal:{mask:02d}", flush=True)
            continue
        print(f"START R41:principal:{mask:02d} entry_lists", flush=True)
        a_values, c_values = primary.entry_lists(height, bounds)
        candidates, counts = enumerate_candidates(primary, center, a_values, c_values, Progress(mask))
        histogram, classification_digest = primary.polygon_result(center, candidates)
        audit = counts["zero_row_audit"]
        observed = {
            "a_count": len(a_values), "a_digest": vectors_digest(a_values),
            "c_count": len(c_values), "c_digest": vectors_digest(c_values),
            "bottom_rows": counts["bottom_rows"],
            "nonzero_c_count": counts["nonzero_candidates"],
            "zero_c_count": counts["zero_candidates"],
            "candidate_count": len(candidates), "candidate_digest": candidate_digest(candidates),
            "classification_digest": classification_digest,
            "polygon_histogram": [[list(key), n] for key, n in histogram],
            "zero_c_audit": [counts["zero_rows"], audit.get("d_zero", 0),
                counts["zero_rows"]-audit.get("d_zero", 0)-audit.get("nonintegral_inverse", 0),
                audit.get("unit_rows", 0), counts["zero_candidates"]],
            "survivor_count": 0, "survivor_digest": candidate_digest(()),
        }
        # polygon_result asserts that every candidate fails a necessary side.
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
