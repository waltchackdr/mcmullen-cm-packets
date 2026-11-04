#!/usr/bin/env python3
"""Portable primary replay for R127 Pic1 mask zero and its principal control.

This integrates reviewed v73 entry units, v74 DEN enumeration, and v79
primitive/sign geometry through the retained primary arithmetic lineage.
Every full run recomputes all finite entry, determinant, and polygon stages,
then checks the exact expected counts and eight canonical artifact hashes.
A preparation-only run makes no geometric decision. Candidate status does
not itself promote a theorem; paired frozen replay and review remain gates.
No runtime source or data is loaded from discovery/build directories.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from math import ceil
import json
import os
from pathlib import Path
import resource
import sys
import time
from tempfile import TemporaryDirectory

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
# trace index, Picard class, basis slope, conductor, inversion, unit exponents.
# These specify witnesses only; their field coefficients are derived below.
PILOTS = {
    "R127:pic0:00": (10, 0, 0, (127, 94, 1), False, (0, 0, 0, 0, 0)),
    "R43_20:pic0:00": (4, 0, 0, (43, 20, 1), False, (0, 0, 0, 0, 0)),
    "R7sq:pic0:00": (12, 0, 0, (7, 1, 2), False, (0, 0, 0, 0, 0)),
    "R41_22:pic0:00": (5, 0, 0, (41, 22, 1), False, (0, 0, 0, 0, 0)),
    "R43_31:pic0:00": (8, 0, 0, (43, 31, 1), False, (1, -1, 0, 0, 1)),
    "R7:pic0:00": (0, 0, 0, (7, 1, 1), False, (0, 0, 0, 0, 0)),
    "R127:pic1:00": (10, 1, 9, (127, 94, 1), True, (1, 0, 0, 0, 0)),
}
BOUND_KEYS = ("DEN", "NUM", "A", "B", "C", "D")
INPUT_RELATIVE = "certificates/broader_principal_v75/pilot_inputs.json"
WITNESSES_RELATIVE = "certificates/broader_principal_v75/unit_witnesses.json"
ENTRY_WITNESSES_RELATIVE = "certificates/r127_pic1_v80/entry_unit_witnesses.json"
EXPECTED_RELATIVE = "certificates/r127_pic1_v80/expected.json"
ZERO = (0,) * 6
SUPPORTED_POINTS = ("R127:pic0:00", "R127:pic1:00")
CANONICAL_ARTIFACTS = (
    "preparation.json", "A-entries.json", "C-entries.json", "bottom-rows.json",
    "candidates-nonzero.json", "candidates-zero.json", "classification.json", "survivors.json",
)
COUNT_KEYS = (
    "a_count", "c_count", "bottom_rows", "nonzero_candidates", "zero_candidates",
    "candidate_count", "survivor_count", "classification_digest", "c_zero_audit",
)

# Literal dictionary: the supervisor reads it without executing this file.
# The six mathematical modules are the complete local primary import closure.
SOURCE_PINS = {
    "certificates/d7_ambient_points_v45.py": "34cb35db91b953162b595fb61e031f551e6df0a0d7b5192dc15a6c00c39274a2",
    "certificates/d7_r83_mask14_v46.py": "70e9c7114adea2b319c3e430046ae3e917717142e7a78e4d43a043950f1147f3",
    "certificates/d7_r83_masks14_30_v47.py": "e78ba5e793134f28d8113d8129809421a5a7754c439afcd64ba0385543799941",
    "certificates/d7_r83_principal_expected_v55.py": "0fd36a67b004e6a250188aa8731a207456b5e49faa955005d6b12cd36665bea7",
    "certificates/d7_r83_principal_v55.py": "a937b69037b5d53f40a600aefb177bf89766d3529b22f73529b2b861ec9238cc",
    "certificates/residual_orbit_k17_k38_v42.py": "5186066f21a2f840ffc71936a54372d7ad51eb6cc919bf7812bb1bec258cd59d",
    "certificates/full_packet_v62/ambient_points.json": "c0e84da6c266c087cbce65a6e8aca8e3f30350ce56b5a1089c4ae4978a381f36",
    "certificates/full_packet_v62/module_templates.json": "09ccff217511a67e946567135f9ffb74da020843621d0895ecce8746d25c245d",
    "certificates/broader_principal_v75/pilot_inputs.json": "0eb217f70238f79308f589cc15a6de7c317e6962a712bd5598b2b6aaa638c386",
    "certificates/broader_principal_v75/unit_witnesses.json": "040ddafb120f5144ecc31ccb3e532121c9f8d588e58bd3483522659e52c43977",
    "certificates/r127_pic1_v80/entry_unit_witnesses.json": "9899ce81e27f2aadbf8bf41ef5ef435e0e92ec9ad5e099378326e07f220f37f5",
    "certificates/r127_pic1_v80/expected.json": "19a23cd3a229dc8c71825166f9e83d13f09d6ece389d3859df49755e9e317cfc"
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()



def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")
    temporary.replace(path)
    return digest(path)



def vector(value):
    return [str(Q(coefficient)) for coefficient in value]



def gamma_record(value):
    return [vector(entry) for entry in value]



def canonical_rows(values):
    return [gamma_record(row) for row in sorted(values)]



def check_sources():
    for relative, expected in SOURCE_PINS.items():
        assert digest(ROOT / relative) == expected, relative



def load_primary():
    check_sources()
    sys.path.insert(0, str(ROOT / "certificates"))
    import d7_r83_principal_v55 as primary

    for relative in SOURCE_PINS:
        if relative.endswith(".py"):
            module = sys.modules[Path(relative).stem]
            assert Path(module.__file__).resolve() == (ROOT / relative).resolve(), relative
    return primary



class Progress:
    def __init__(self, output, point_label):
        assert point_label in PILOTS
        self.output = output
        self.point_label = point_label
        self.started = time.monotonic()
        self.cpu_started = time.process_time()
        self.last = -10.0
        self.stage = None

    def report(self, stage, force=False, **counts):
        now = time.monotonic()
        if not force and stage == self.stage and now - self.last < 5:
            return
        self.stage, self.last = stage, now
        record = {
            "point_label": self.point_label, "backend": "primary", "pid": os.getpid(),
            "stage": stage, "updated_utc": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": now - self.started,
            "cpu_seconds": time.process_time() - self.cpu_started, **counts,
        }
        atomic_json(self.output / "progress.json", record)
        print(json.dumps(record, sort_keys=True), flush=True)



def primary_unit(primary, exponents):
    """Build a real unit from the fixed primary polynomial generators."""
    ambient = primary.ambient
    t = ambient.T
    generators = (t**4 - 3*t**2 + 1, t**3 + t**2 - 2*t - 1,
                  t**5 - 4*t**3 + 3*t + 1, t**3 - 2*t, t)
    assert len(exponents) == 5 and all(type(exponent) is int for exponent in exponents)
    answer = primary.sp.Integer(1)
    for generator, retained, exponent in zip(generators, ambient.XI_UNITS, exponents):
        assert ambient.reduce_l(generator - retained) == 0
        answer = ambient.reduce_l(answer * primary.field_power(generator, exponent))
    inverse = ambient.reduce_field(1 / answer)
    ambient.assert_integral(answer)
    ambient.assert_integral(inverse)
    assert ambient.reduce_l(answer * inverse) == 1
    assert abs(ambient.absolute_norm(answer)) == abs(ambient.absolute_norm(inverse)) == 1
    return answer, inverse



def reconstruct(primary, point_label):
    """Derive the common-basis action and its selected SL2 transport."""
    if point_label not in PILOTS:
        raise ValueError("Only the seven reviewed mask-zero v65 pilots are supported")
    a, sp, core = primary.ambient, primary.sp, primary.core
    trace_index, picard, slope, conductor, inversion_first, exponents = PILOTS[point_label]
    stratum = point_label.split(":", 1)[0]
    prime, residue, exponent = conductor
    t = a.T
    aa, bb = a.alpha_pair(trace_index)
    trace = a.alpha_trace(trace_index)
    assert a.reduce_l(aa * aa + t * aa * bb + bb * bb - a.PI) == 0
    assert abs(a.absolute_norm(bb)) == prime**exponent
    assert int(bb.subs(t, residue)) % prime == 0

    # The columns are 1+s*z and B_k*z. For Pic1, s=9 rather than 0.
    basis = sp.Matrix([[1, 0], [slope, bb]])
    maximal_action = a.multiplication_matrix((aa, bb))
    raw = a.matrix_reduce(basis.inv() * maximal_action * basis)
    n = a.reduce_l(1 + slope*t + slope*slope)
    assert a.matrix_reduce(basis * raw - maximal_action * basis) == sp.zeros(2)
    explicit = sp.Matrix([[aa - slope*bb, -bb*bb], [n, aa + (t + slope)*bb]])
    assert a.matrix_reduce(raw - explicit) == sp.zeros(2)
    a.assert_integral_matrix(raw)
    a.assert_characteristic(raw, trace_index)
    z_in_basis = a.matrix_reduce(basis.inv() * sp.Matrix([[0, -1], [1, t]]) * basis)
    assert not all(a.is_integral(value) for value in z_in_basis)
    assert a.SQUARE_CLASS_UNITS[0] == 1

    def invariants(matrix):
        c = matrix[1, 0]
        center = a.reduce_field((matrix[0, 0] - matrix[1, 1]) / (2*c))
        height = a.reduce_field((4*a.PI - trace*trace) / (4*c*c))
        modulus = a.reduce_l(center*center + height)
        assert a.sign_vector(height) == a.sign_vector(modulus) == (1,) * 6
        return center, height, modulus

    raw_x, raw_y, raw_mod = invariants(raw)
    assert a.reduce_field(raw_x + bb*(t + 2*slope)/(2*n)) == 0
    assert a.reduce_field(raw_y - bb*bb*(4 - t*t)/(4*n*n)) == 0
    assert a.reduce_field(raw_mod - bb*bb/n) == 0
    unit, inverse_unit = primary_unit(primary, exponents)
    inversion = sp.Matrix([[0, -1], [1, 0]]) if inversion_first else sp.eye(2)
    conjugator = a.matrix_multiply(sp.diag(unit, inverse_unit), inversion)
    conjugator_inverse = a.matrix_reduce(conjugator.inv())
    a.assert_integral_matrix(conjugator)
    a.assert_integral_matrix(conjugator_inverse)
    assert a.reduce_l(conjugator.det()) == 1
    assert a.matrix_multiply(conjugator, conjugator_inverse) == sp.eye(2)
    matrix = a.matrix_multiply(a.matrix_multiply(conjugator, raw), conjugator_inverse)
    a.assert_integral_matrix(matrix)
    a.assert_characteristic(matrix, trace_index)
    center_sym, height_sym, modulus_sym = invariants(matrix)

    # Independently check the direct Möbius formulas on the raw point.
    ga, gb, gc, gd = tuple(conjugator)
    denominator = a.reduce_l((gc*raw_x + gd)**2 + gc*gc*raw_y)
    assert a.sign_vector(denominator) == (1,) * 6
    mobius_x = a.reduce_field(((ga*raw_x + gb)*(gc*raw_x + gd) + ga*gc*raw_y)/denominator)
    mobius_y = a.reduce_field(raw_y/(denominator*denominator))
    mobius_mod = a.reduce_field(((ga*raw_x + gb)**2 + ga*ga*raw_y)/denominator)
    assert a.reduce_l(center_sym - mobius_x) == 0
    assert a.reduce_l(height_sym - mobius_y) == 0
    assert a.reduce_l(modulus_sym - mobius_mod) == 0
    center, height, modulus = map(primary.to_core, (center_sym, height_sym, modulus_sym))
    raw_center, raw_height, raw_modulus = map(primary.to_core, (raw_x, raw_y, raw_mod))
    matrix_coefficients = tuple(primary.to_core(entry) for entry in matrix)
    raw_coefficients = tuple(primary.to_core(entry) for entry in raw)
    conjugator_coefficients = tuple(primary.to_core(entry) for entry in conjugator)
    bounds = primary.derive_bounds(height, modulus)
    raw_bounds = primary.derive_bounds(raw_height, raw_modulus)
    assert core.EMBEDDINGS == (1, 5, 11, 13, 17, 19)

    # Prepared data and promoted manifests are comparisons, not definitions.
    prepared = json.loads((ROOT / INPUT_RELATIVE).read_text())
    assert prepared["ordered_embeddings"] == list(core.EMBEDDINGS)
    assert prepared["matrix_entry_order"] == ["a", "b", "c", "d"]
    assert {row["point_label"] for row in prepared["pilots"]} == set(PILOTS)
    pilots = [row for row in prepared["pilots"] if row["point_label"] == point_label]
    assert len(pilots) == 1
    pilot = pilots[0]
    assert pilot["trace_index"] == trace_index
    assert pilot["conductor"] == list(conductor)
    assert pilot["stratum"] == stratum
    assert pilot["picard_class"] == picard and pilot["orientation_mask"] == 0
    selected = pilot["selected"]
    assert selected["inversion_first"] is inversion_first
    assert selected["unit_exponents"] == list(exponents)
    assert tuple(map(Q, selected["unit"])) == primary.to_core(unit)
    assert tuple(map(Q, selected["unit_inverse"])) == primary.to_core(inverse_unit)
    assert tuple(tuple(map(Q, value)) for value in selected["conjugator"]) == conjugator_coefficients
    for point, expected_matrix, xx, yy, mm, expected_bounds in (
            (pilot["raw"], raw_coefficients, raw_center, raw_height, raw_modulus, raw_bounds),
            (selected["point"], matrix_coefficients, center, height, modulus, bounds)):
        assert tuple(tuple(map(Q, value)) for value in point["matrix"]) == expected_matrix
        for key, actual in (("center", xx), ("height_squared", yy), ("modulus_squared", mm)):
            assert tuple(map(Q, point[key])) == actual, key
        assert tuple(tuple(map(Q, point["bounds"][key])) for key in BOUND_KEYS) == expected_bounds

    manifest = json.loads((ROOT / "certificates/full_packet_v62/ambient_points.json").read_text())
    matches = [row for row in manifest if row["point_label"] == point_label]
    assert len(matches) == 1
    record = matches[0]
    assert record["conductor"] == list(conductor)
    assert record["picard_class"] == picard and record["orientation_mask"] == 0
    assert record["trace_for_matrix"] == trace_index
    assert record["root_partner"] == [0 if picard == 0 else 2, 16]
    assert tuple(tuple(map(Q, value)) for value in record["norm_pi_matrix"]) == raw_coefficients
    assert tuple(Q(value) / 2 for value in record["twice_real_part"]) == raw_center
    assert tuple(Q(value) / 4 for value in record["four_times_height_squared"]) == raw_height
    assert tuple(map(Q, record["modulus_squared"])) == raw_modulus
    templates = json.loads((ROOT / "certificates/full_packet_v62/module_templates.json").read_text())
    matches = [row for row in templates if row["stratum"] == stratum and row["picard_class"] == picard]
    assert len(matches) == 1
    template = matches[0]
    assert template["representative_trace"] == trace_index and template["slope"] == slope
    assert template["index_norm"] == prime**exponent and template["conductor"] == list(conductor)
    assert tuple(tuple(map(Q, value)) for value in template["basis_columns_in_one_z"]) == tuple(primary.to_core(entry) for entry in basis)
    assert tuple(tuple(map(Q, value)) for value in template["norm_pi_action"]) == raw_coefficients

    # The primary core has only these point-dependent globals. Center is an
    # explicit argument throughout; its sole lru cache stores field inverses.
    primary.enumeration.install_point(height, bounds)
    assert core.Y2 == height
    assert tuple(getattr(core, key + "_BOUND") for key in BOUND_KEYS) == bounds
    preparation = {
        "point_label": point_label, "trace_index": trace_index,
        "matrix": gamma_record(matrix_coefficients), "center": vector(center),
        "height_squared": vector(height), "modulus_squared": vector(modulus),
        "bounds": {key: vector(values) for key, values in zip(BOUND_KEYS, bounds)},
        "conjugator": gamma_record(conjugator_coefficients),
    }
    return center, height, bounds, preparation



def entry_lists(primary, bounds, progress, output, artifacts, entry_scalings):
    """Enumerate A/C via integral units, testing and retaining original v."""
    lists = []
    core = primary.core
    for label, coordinate_bounds in (("A", bounds[2]), ("C", bounds[4])):
        progress.report("entries_" + label, force=True)
        scaling = entry_scalings[label]
        # q=0 in the shifted-unit lemma: w=u*v, center zero. The
        # original box predicate receives v=u^-1*w, never w.
        values = inverse_mapped_ellipsoid(
            core, core.ZERO, scaling["gross_radius"], scaling["inverse_rows"],
            lambda value: primary.within_box(value, coordinate_bounds),
        )
        assert len(values) == len(set(values))
        assert {core.neg(value) for value in values} == set(values)
        values = sorted(values)
        name = label + "-entries.json"
        artifacts[name] = atomic_json(output / name, [vector(value) for value in values])
        lists.append(values)
        progress.report("entries_" + label, force=True, count=len(values))
    return lists



def entry_unit_scaling(primary, entry_bounds, preparation, bound_kind):
    """Derive and verify a DEN, NUM, A or C unit in primary arithmetic."""
    if bound_kind not in ("DEN", "NUM", "A", "C"):
        raise ValueError("Only the reviewed DEN, NUM, A and C witnesses are supported")
    witness_relative = ENTRY_WITNESSES_RELATIVE if bound_kind in ("A", "C") else WITNESSES_RELATIVE
    ambient, sp, core = primary.ambient, primary.sp, primary.core
    point_label = preparation["point_label"]
    assert point_label in PILOTS
    bundle = json.loads((ROOT / witness_relative).read_text())
    assert set(bundle["pilots"]) == set(PILOTS)
    witness = bundle["pilots"][point_label]
    assert witness["point_label"] == point_label
    assert witness["ordered_embeddings"] == list(core.EMBEDDINGS)
    assert witness["unit_generator_order"] == ["xi5", "xi11", "xi13", "xi17", "xi19=t"]
    selected = witness["witnesses"][bound_kind]["selected"]
    exponents = tuple(selected["unit_exponents"])
    t = ambient.T
    unit_sym, inverse_sym = primary_unit(primary, exponents)
    unit, inverse = map(primary.to_core, (unit_sym, inverse_sym))

    # Integer column j is multiplication applied to t^j, independently
    # reconstructed in the primary power-basis arithmetic.
    def multiplication_columns(value):
        result = []
        for j in range(6):
            product = ambient.reduce_l(value * t**j)
            ambient.assert_integral(product)
            coefficients = primary.to_core(product)
            assert all(coefficient.denominator == 1 for coefficient in coefficients)
            result.append(tuple(map(int, coefficients)))
        return tuple(result)

    columns = multiplication_columns(unit_sym)
    inverse_columns = multiplication_columns(inverse_sym)
    forward = sp.Matrix.hstack(*(sp.Matrix(column) for column in columns))
    backward = sp.Matrix.hstack(*(sp.Matrix(column) for column in inverse_columns))
    assert forward * backward == backward * forward == sp.eye(6)
    assert abs(forward.det()) == abs(backward.det()) == 1
    inverse_rows = tuple(tuple(int(backward[i, j]) for j in range(6)) for i in range(6))

    transformed_bounds = []
    for i, bound in enumerate(entry_bounds):
        assert bound > 0
        lo, hi = core.interval_eval(unit, core.ROOT_INTERVALS[i])
        enclosing = Q(ceil(1000000 * max(abs(lo), abs(hi)) * bound), 1000000)
        assert enclosing > 0
        # Separate exact sign check certifies each absolute-value bound,
        # including either sign of the real unit embedding.
        difference = core.sub(core.scale(enclosing * enclosing, core.ONE),
                              core.scale(bound * bound, core.mul(unit, unit)))
        assert core.sign_at(difference, i) >= 0
        transformed_bounds.append(enclosing)
    transformed_bounds = tuple(transformed_bounds)
    gross_radius = sum((bound * bound for bound in transformed_bounds), Q(0))

    for key in ("matrix", "center", "height_squared", "modulus_squared"):
        assert witness["point"][key] == preparation[key]
    assert tuple(map(Q, selected["unit"])) == unit
    assert tuple(map(Q, selected["unit_inverse"])) == inverse
    assert tuple(tuple(row) for row in selected["multiplication_matrix_columns"]) == columns
    assert tuple(tuple(row) for row in selected["inverse_matrix_columns"]) == inverse_columns
    assert tuple(map(Q, witness["witnesses"][bound_kind]["original_bounds"])) == entry_bounds
    assert tuple(map(Q, selected["transformed_bounds"])) == transformed_bounds
    assert Q(selected["transformed_gross_radius"]) == gross_radius
    assert Q(selected["original_gross_radius"]) == sum((bound * bound for bound in entry_bounds), Q(0))
    audit = {
        "scope": "ENTRY_ENUMERATION_COORDINATE_CHANGE_ONLY_NO_GEOMETRIC_CONJUGATION",
        "point_label": point_label, "unit": vector(unit), "unit_inverse": vector(inverse),
        "unit_exponents": list(exponents),
        "multiplication_matrix_columns": [list(column) for column in columns],
        "inverse_matrix_columns": [list(column) for column in inverse_columns],
        "original_bounds": vector(entry_bounds), "transformed_bounds": vector(transformed_bounds),
        "transformed_gross_radius": str(gross_radius),
        "witnesses_sha256": SOURCE_PINS[witness_relative],
    }
    return {"unit": unit, "inverse_rows": inverse_rows, "gross_radius": gross_radius, "audit": audit}



def num_unit_scaling(primary, num_bounds, preparation):
    """Preserve the v71 NUM verification and return contract."""
    return entry_unit_scaling(primary, num_bounds, preparation, "NUM")



def den_unit_scaling(primary, den_bounds, preparation):
    """Verify the existing DEN witness by the separate primary code path."""
    return entry_unit_scaling(primary, den_bounds, preparation, "DEN")



def inverse_mapped_ellipsoid(core, center, radius, inverse_rows, original_predicate):
    """Enumerate integral w, then test and return its original field entry."""
    if radius < 0:
        return []

    def inverse_map(w):
        return tuple(sum(row[j] * w[j] for j in range(6)) for row in inverse_rows)

    # Radius zero is intentional: the retained recursion includes equality.
    values = core.enumerate_trace_ellipsoid(
        center, radius, lambda w: original_predicate(inverse_map(w)))
    result = [inverse_map(w) for w in values]
    assert len(result) == len(set(result))
    return result



class PrimaryLattices:
    """Exact integral power-basis multiplication and column-lattice HNF."""
    def __init__(self, core, sp, hnf):
        self.core, self.sp, self.hnf = core, sp, hnf
        self.basis = [tuple(Q(int(i == j)) for i in range(6)) for j in range(6)]
        self.matrices = {}
        self.ideals = {}

    def multiplication_matrix(self, value):
        value = tuple(value)
        assert len(value) == 6 and all(type(v) is int for v in value)
        if value not in self.matrices:
            columns = [self.core.mul(self.core.as_field(value), basis) for basis in self.basis]
            assert all(q.denominator == 1 for column in columns for q in column)
            matrix = self.sp.Matrix.hstack(*(self.sp.Matrix([int(q) for q in column]) for column in columns))
            assert tuple(matrix[:, 0]) == value
            self.matrices[value] = matrix
        return self.matrices[value]

    def column_lattice(self, matrix):
        assert matrix.rows == 6 and all(value.is_Integer for value in matrix)
        h = self.hnf(matrix)
        rank = h.cols
        assert h.rows == 6 and rank in (0, 6)
        if rank == 0:
            assert matrix.is_zero_matrix
            return {"rank": 0, "index": None, "hnf": [[] for _ in range(6)]}
        for i in range(6):
            assert h[i, i] > 0
            assert all(h[i, j] == 0 for j in range(i))
            assert all(0 <= h[i, j] < h[i, i] for j in range(i + 1, 6))
        index = abs(int(h.det()))
        assert index >= 1
        # Independently check every input column is an integral combination
        # of the resulting upper-triangular HNF columns (no rational inverse).
        for j in range(matrix.cols):
            coefficients = [0] * 6
            for i in range(5, -1, -1):
                residual = int(matrix[i, j]) - sum(int(h[i, k])*coefficients[k] for k in range(i + 1, 6))
                diagonal = int(h[i, i])
                assert residual % diagonal == 0
                coefficients[i] = residual // diagonal
        return {"rank": 6, "index": str(index),
                "hnf": [[int(h[i, j]) for j in range(6)] for i in range(6)]}

    def row_ideal(self, row):
        c, d = row
        data = self.column_lattice(self.multiplication_matrix(c).row_join(self.multiplication_matrix(d)))
        assert (data["index"] is None) == (row == (ZERO, ZERO))
        return data

    def principal_ideal(self, c):
        assert c != ZERO
        key = min(c, tuple(-value for value in c))
        if key not in self.ideals:
            self.ideals[key] = self.column_lattice(self.multiplication_matrix(key))
        return self.ideals[key]



def verify_affine_data(lattices, row, affine):
    """Certify ad-bc=1 for every integral b returned by the retained affine map.

    For b=(R*a+offset)/D, M_c*R=D*M_d and M_c*offset=-D*e0
    imply D*c*b=D*(a*d-1). Since D>0, this proves the determinant
    identity for all a, including a outside the finite test list.
    """
    c,d = row
    assert c != ZERO
    assert isinstance(affine,tuple) and len(affine) == 3
    denominator,rows,offset = affine
    assert type(denominator) is int and denominator > 0
    assert len(rows) == 6 and all(len(r) == 6 for r in rows) and len(offset) == 6
    assert all(type(v) is int for r in rows for v in r)
    assert all(type(v) is int for v in offset)
    sp = lattices.sp
    mc,md = lattices.multiplication_matrix(c),lattices.multiplication_matrix(d)
    assert mc*sp.Matrix(rows) == denominator*md
    assert mc*sp.Matrix(offset) == sp.Matrix([-denominator,0,0,0,0,0])
    return {"denominator":denominator,"rows":[list(r) for r in rows],"offset":list(offset),
            "verified_identities":["M_c*R=D*M_d","M_c*offset=-D*e0"],
            "conclusion":"Every integral b=(R*a+offset)/D satisfies a*d-b*c=1"}



def negate_matrix(value):
    return tuple(tuple(-coefficient for coefficient in entry) for entry in value)



def append_checkpoint(path, record):
    with path.open("a") as stream:
        stream.write(json.dumps({"status": "INCOMPLETE_GEOMETRY_ROW_CHECKPOINT", **record},
                                sort_keys=True, separators=(",", ":")) + "\n")
        stream.flush()



class TimedRows:
    """Observe the unchanged c=0 loop without replacing its arithmetic body."""
    def __init__(self, rows, journal):
        self.rows, self.journal, self.completed = rows, journal, 0
        self.cpu_seconds = 0.0

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        assert self.completed == 0
        for index, row in enumerate(self.rows):
            record = {"stage": "zero_c", "row_index": index, "bottom_row": gamma_record(row)}
            append_checkpoint(self.journal, {**record, "event": "ROW_START"})
            started = time.process_time()
            yield row
            elapsed = time.process_time() - started
            self.cpu_seconds += elapsed
            self.completed += 1
            append_checkpoint(self.journal, {**record, "event": "ROW_COMPLETE", "cpu_seconds": elapsed})



def bottom_rows_only(primary, x, c_values, progress, output, artifacts, den_scaling):
    """Extract the unchanged v72/v73 DEN stage, with per-C partial journals."""
    core = primary.core
    den_squares = core.square_bounds(core.DEN_BOUND)
    den_radius = sum(den_squares, Q(0))
    bottom_rows = []
    original_records, scaled_den_records = [], []
    original_signs, scaled_den_signs = Counter(), Counter()
    x_boxes = tuple(core.scaled_field_box(x, i) for i in range(6))
    y2_boxes = tuple(core.scaled_field_box(core.Y2, i) for i in range(6))

    def bounded(v, w, boxes, squares, exact):
        unresolved = False
        for i in range(6):
            real = core.scaled_add(core.scaled_mul(boxes[i], x_boxes[i]), core.scaled_eval_integral(w, i))
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

    stage_started = time.process_time()
    total_visits = 0
    journal_path = output / "partial-bottom-rows.jsonl"
    progress.report("bottom_rows", force=True, completed=0, total=len(c_values))
    with journal_path.open("x") as journal:
        for n, c in enumerate(c_values, 1):
            cq = core.as_field(c)
            boxes = tuple(core.scaled_eval_integral(c, i) for i in range(6))
            imaginary_trace = core.trace(core.mul(core.mul(cq, cq), core.Y2))
            assert imaginary_trace >= 0
            original_radius = den_radius - imaginary_trace
            original_sign = "negative" if original_radius < 0 else "zero" if original_radius == 0 else "positive"
            original_signs[original_sign] += 1
            original_records.append({"entry": vector(c), "old_radius": str(den_radius),
                                     "imaginary_trace": str(imaginary_trace),
                                     "radius": str(original_radius), "sign": original_sign})
            scaled_c = core.mul(den_scaling["unit"], cq)
            scaled_center = core.neg(core.mul(scaled_c, x))
            weighted_trace = core.trace(core.mul(core.mul(scaled_c, scaled_c), core.Y2))
            assert weighted_trace >= 0
            scaled_radius = den_scaling["gross_radius"] - weighted_trace
            sign = "negative" if scaled_radius < 0 else "zero" if scaled_radius == 0 else "positive"
            scaled_den_signs[sign] += 1
            radius_record = {
                "c": vector(c), "coordinate": "w=u_DEN*d", "scaled_c": vector(scaled_c),
                "transformed_center": vector(scaled_center),
                "raw_original_gross_radius": str(den_radius), "raw_original_radius": str(original_radius),
                "transformed_gross_radius": str(den_scaling["gross_radius"]),
                "fixed_term": "Tr((u_DEN*c)^2*h^2)", "weighted_trace": str(weighted_trace),
                "radius": str(scaled_radius), "sign": sign,
            }
            scaled_den_records.append(radius_record)
            c_started = time.process_time()
            last_checkpoint = time.monotonic()
            visits = 0
            accepted = []

            def checkpoint(event):
                record = {
                    "status": "INCOMPLETE_STAGE_CHECKPOINT", "event": event,
                    "point_label": progress.point_label, "index": n, "total": len(c_values),
                    **radius_record, "visited_leaves": visits,
                    "retained_d_count": len(accepted), "retained_d": [vector(d) for d in sorted(accepted)],
                    "c_cpu_seconds": time.process_time() - c_started,
                    "stage_cpu_seconds": time.process_time() - stage_started,
                }
                journal.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
                journal.flush()

            def original_denominator(d):
                nonlocal visits, last_checkpoint
                visits += 1
                keep = bounded(c, d, boxes, den_squares, core.denominator_square)
                if keep:
                    accepted.append(d)
                if visits % 4096 == 0 and time.monotonic() - last_checkpoint >= 5:
                    checkpoint("C_PROGRESS")
                    progress.report("bottom_rows", force=True, completed=n - 1, total=len(c_values),
                                    active_c_index=n, active_c=vector(c), active_radius=str(scaled_radius),
                                    active_visits=visits, active_retained=len(accepted),
                                    bottom_rows=len(bottom_rows), visited_leaves=total_visits + visits)
                    last_checkpoint = time.monotonic()
                return keep

            checkpoint("C_START")
            values = inverse_mapped_ellipsoid(
                core, scaled_center, scaled_radius, den_scaling["inverse_rows"], original_denominator)
            assert values == accepted
            total_visits += visits
            bottom_rows.extend((c, d) for d in values)
            checkpoint("C_COMPLETE")
            progress.report("bottom_rows", completed=n, total=len(c_values), bottom_rows=len(bottom_rows),
                            visited_leaves=total_visits)
    assert len(bottom_rows) == len(set(bottom_rows))
    artifacts["bottom-rows.json"] = atomic_json(output / "bottom-rows.json", canonical_rows(bottom_rows))
    artifacts["den-unit-radii.json"] = atomic_json(output / "den-unit-radii.json", {
        "scope": "DEN_ENTRY_ENUMERATION_COORDINATE_CHANGE_ONLY_ORIGINAL_PREDICATE_UNCHANGED",
        "unit": vector(den_scaling["unit"]),
        "summary": {sign: scaled_den_signs[sign] for sign in ("negative", "zero", "positive")},
        "records": scaled_den_records,
    })
    artifacts["den-original-radii.json"] = atomic_json(output / "den-original-radii.json", {
        "scope": "ORIGINAL_DEN_TRACE_RADIUS_COMPATIBILITY_AUDIT_ONLY",
        "summary": {sign: original_signs[sign] for sign in ("negative", "zero", "positive")},
        "records": original_records,
    })
    artifacts["partial-bottom-rows.jsonl"] = digest(journal_path)
    return sorted(bottom_rows), {"bottom_rows": len(bottom_rows), "visited_leaves": total_visits,
            "bottom_stage_cpu_seconds": time.process_time() - stage_started}



def enumerate_candidates(primary, x, inputs, progress, output, artifacts, num_scaling):
    """Fresh complete domain; exact primitive sign rows and original predicates."""
    from sympy.matrices.normalforms import hermite_normal_form
    core = primary.core
    a_values, c_values, bottom_rows = inputs["a_values"], inputs["c_values"], inputs["bottoms"]
    den_squares = core.square_bounds(core.DEN_BOUND)
    num_squares = core.square_bounds(core.NUM_BOUND)
    num_radius = sum(num_squares, Q(0))
    radius_records = {"b_zero_c": []}
    radius_signs = {"b_zero_c": Counter()}
    timings = {name: 0.0 for name in ("box_cache_cpu", "bottom_validation_cpu", "primitive_check_cpu",
                                      "affine_setup_cpu", "affine_certificate_cpu", "nonzero_scan_cpu",
                                      "retained_determinant_cpu", "zero_c_cpu")}
    started = time.process_time()
    x_boxes = tuple(core.scaled_field_box(x, i) for i in range(6))
    y2_boxes = tuple(core.scaled_field_box(core.Y2, i) for i in range(6))
    a_boxes = {a: tuple(core.scaled_eval_integral(a, i) for i in range(6)) for a in a_values}
    timings["box_cache_cpu"] = time.process_time() - started

    def bounded(v, w, boxes, squares, exact):
        unresolved = False
        for i in range(6):
            real = core.scaled_add(core.scaled_mul(boxes[i], x_boxes[i]), core.scaled_eval_integral(w, i))
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


    def reduced_radius(branch, entry, field_entry, old_radius, **context):
        imaginary_trace = core.trace(core.mul(core.mul(field_entry, field_entry), core.Y2))
        assert imaginary_trace >= 0
        radius = old_radius - imaginary_trace
        sign = "negative" if radius < 0 else "zero" if radius == 0 else "positive"
        radius_signs[branch][sign] += 1
        radius_records[branch].append({"entry": vector(entry), "old_radius": str(old_radius),
                                       "imaginary_trace": str(imaginary_trace),
                                       "radius": str(radius), "sign": sign, **context})
        return radius


    started = time.process_time()
    c_boxes = {c: tuple(core.scaled_eval_integral(c, i) for i in range(6)) for c in {row[0] for row in bottom_rows}}
    assert all(bounded(c, d, c_boxes[c], den_squares, core.denominator_square) for c, d in bottom_rows)
    timings["bottom_validation_cpu"] = time.process_time() - started
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != core.ZERO]
    representatives = sorted({min(row, negate_matrix(row)) for row in nonzero_rows})
    assert len(representatives) * 2 == len(nonzero_rows)
    journal = output / "geometry-rows-partial.jsonl"
    journal.touch(exist_ok=False)
    zero_rows = TimedRows([(c, d) for c, d in bottom_rows if c == core.ZERO], journal)
    lattices = PrimaryLattices(core, primary.sp, hermite_normal_form)
    nonzero, row_audits = [], []
    fresh_indices = {}
    audit_counts = Counter()
    bottom_set, a_set = set(bottom_rows), set(a_values)
    progress.report("determinant_nonzero_c", force=True, completed=0, total=len(representatives))
    for n, (c, d) in enumerate(representatives, 1):
        row_start = time.process_time()
        context = {"stage": "nonzero_c", "row_index": n - 1, "bottom_row": gamma_record((c, d))}
        append_checkpoint(journal, {**context, "event": "ROW_START"})
        row_times = {}
        started = time.process_time()
        cinv = core.inverse(core.as_field(c))
        unit_c = core.integral(cinv)
        if unit_c:
            assert core.mul(core.as_field(c), cinv) == core.ONE
            index = "1"
            primitive_witness = {"unit_c_inverse": vector(cinv)}
        else:
            ideal = lattices.row_ideal((c, d))
            index = ideal["index"]
            primitive_witness = {"row_ideal_hnf": ideal["hnf"]}
        assert index is not None and int(index) >= 1
        fresh_indices[(c, d)] = fresh_indices[negate_matrix((c, d))] = index
        row_times["primitive_check_cpu"] = time.process_time() - started
        affine_certificate = None
        accepted = []
        integral_matches = numerator_rejections = 0
        if index == "1":
            audit_counts["primitive_sign_rows"] += 1
            audit_counts["unit_c_sign_rows"] += unit_c
            started = time.process_time()
            affine = core.affine_b_data(c, d)
            row_times["affine_setup_cpu"] = time.process_time() - started
            started = time.process_time()
            affine_certificate = verify_affine_data(lattices, (c, d), affine)
            row_times["affine_certificate_cpu"] = time.process_time() - started
            started = time.process_time()
            for a in a_values:
                b = core.affine_integral_b_fast(a, affine)
                if b is None:
                    continue
                integral_matches += 1
                if not numerator(a, b):
                    numerator_rejections += 1
                    continue
                accepted.append((a, b, c, d))
            row_times["nonzero_scan_cpu"] = time.process_time() - started
            audit_counts["a_probes"] += len(a_values)
            started = time.process_time()
            for gamma in accepted:
                opposite = negate_matrix(gamma)
                assert gamma != opposite
                for restored in (gamma, opposite):
                    aa, bb, cc, dd = restored
                    assert aa in a_set and (cc, dd) in bottom_set
                    assert numerator(aa, bb)
                    assert core.sub(core.mul(core.as_field(aa), core.as_field(dd)),
                                    core.mul(core.as_field(bb), core.as_field(cc))) == core.ONE
                    nonzero.append(restored)
            row_times["retained_determinant_cpu"] = time.process_time() - started
        else:
            assert int(index) > 1
            audit_counts["nonprimitive_sign_rows"] += 1
        audit_counts["integral_matches"] += integral_matches
        audit_counts["numerator_rejections"] += numerator_rejections
        audit_counts["bounded_sign_candidates"] += len(accepted)
        for key, elapsed in row_times.items():
            timings[key] += elapsed
        row_record = {**context, "event": "ROW_COMPLETE", "unit_c": unit_c,
                      "ideal_index": index, "primitive_witness": primitive_witness,
                      "affine_certificate": affine_certificate, "integral_matches": integral_matches,
                      "numerator_rejections": numerator_rejections, "bounded_sign_candidates": len(accepted),
                      "restored_candidates": 2 * len(accepted), "timing_seconds": row_times,
                      "cpu_seconds": time.process_time() - row_start}
        row_audits.append(row_record)
        append_checkpoint(journal, row_record)
        progress.report("determinant_nonzero_c", completed=n, total=len(representatives), candidates=len(nonzero))
    nonzero_set = set(nonzero)
    assert len(nonzero) == len(nonzero_set)
    assert {negate_matrix(gamma) for gamma in nonzero_set} == nonzero_set
    assert len(nonzero) == 2 * audit_counts["bounded_sign_candidates"]
    assert audit_counts["primitive_sign_rows"] + audit_counts["nonprimitive_sign_rows"] == len(representatives)
    assert audit_counts["integral_matches"] == audit_counts["numerator_rejections"] + audit_counts["bounded_sign_candidates"]
    artifacts["candidates-nonzero.json"] = atomic_json(output / "candidates-nonzero.json", canonical_rows(nonzero))

    # Record every c=0 ideal index afresh; no saved primitive data is read.
    started = time.process_time()
    for c, d in zero_rows.rows:
        if d == core.ZERO:
            fresh_indices[(c, d)] = None
        else:
            dinv = core.inverse(core.as_field(d))
            if core.integral(dinv):
                assert core.mul(core.as_field(d), dinv) == core.ONE
                fresh_indices[(c, d)] = "1"
            else:
                fresh_indices[(c, d)] = lattices.row_ideal((c, d))["index"]
    timings["zero_primitive_check_cpu"] = time.process_time() - started
    assert set(fresh_indices) == bottom_set
    assert all(fresh_indices[row] == fresh_indices[negate_matrix(row)] for row in bottom_rows)
    artifacts["row-indices.json"] = atomic_json(output / "row-indices.json", [gamma_record(row) + [fresh_indices[row]] for row in sorted(bottom_rows)])
    artifacts["primitive-rows.json"] = atomic_json(output / "primitive-rows.json", canonical_rows([row for row in bottom_rows if fresh_indices[row] == "1"]))

    zero, zero_audit = [], Counter()
    progress.report("determinant_zero_c", force=True, completed=0, total=len(zero_rows))
    for n, (c, d) in enumerate(zero_rows, 1):
        if d == core.ZERO:
            zero_audit["d_zero"] += 1
        else:
            dq = core.as_field(d)
            aq = core.inverse(dq)
            if not core.integral(aq):
                zero_audit["nonintegral_inverse"] += 1
            else:
                zero_audit["integral_inverse"] += 1
                a = core.integral_tuple(aq)
                if not primary.within_box(a, core.A_BOUND):
                    zero_audit["a_outside_box"] += 1
                else:
                    assert a in a_boxes
                    zero_audit["unit_rows"] += 1
                    scaled_a = core.mul(num_scaling["unit"], aq)
                    scaled_center = core.neg(core.mul(scaled_a, x))
                    radius = reduced_radius(
                        "b_zero_c", a, scaled_a, num_scaling["gross_radius"], d=vector(d),
                        coordinate="w=u*b", raw_original_gross_radius=str(num_radius),
                        transformed_gross_radius=str(num_scaling["gross_radius"]),
                        transformed_center=vector(scaled_center), scaled_a=vector(scaled_a),
                        fixed_term="Tr((u*a)^2*h^2)")
                    bs = inverse_mapped_ellipsoid(
                        core, scaled_center, radius, num_scaling["inverse_rows"], lambda b: numerator(a, b))
                    assert core.mul(aq, dq) == core.ONE
                    zero.extend((a, b, c, d) for b in bs)
        progress.report("determinant_zero_c", completed=n, total=len(zero_rows), candidates=len(zero))
    assert zero_audit["d_zero"] == 1
    assert zero_audit["d_zero"] + zero_audit["nonintegral_inverse"] + zero_audit["integral_inverse"] == len(zero_rows)
    assert zero_audit["a_outside_box"] + zero_audit["unit_rows"] == zero_audit["integral_inverse"]
    assert len(zero) == len(set(zero)) and not (set(zero) & set(nonzero))
    artifacts["candidates-zero.json"] = atomic_json(output / "candidates-zero.json", canonical_rows(zero))
    assert zero_rows.completed == len(zero_rows)
    timings["zero_c_cpu"] = zero_rows.cpu_seconds
    started = time.process_time()
    zero_set = set(zero)
    assert {negate_matrix(gamma) for gamma in zero_set} == zero_set
    for aa, bb, cc, dd in zero:
        assert numerator(aa, bb)
        assert core.sub(core.mul(core.as_field(aa), core.as_field(dd)),
                        core.mul(core.as_field(bb), core.as_field(cc))) == core.ONE
    timings["retained_determinant_cpu"] += time.process_time() - started
    radius_summary = {branch: {sign: counts[sign] for sign in ("negative", "zero", "positive")}
                      for branch, counts in radius_signs.items()}
    artifacts["trace-radii.json"] = atomic_json(output / "trace-radii.json", {
        "scope": "C_ZERO_RADIUS_REPLAY_WITH_SEPARATE_FRESH_DEN_RADIUS_AUDIT",
        "summary": radius_summary, "records": radius_records,
    })
    artifacts["geometry-row-audit.json"] = atomic_json(output / "geometry-row-audit.json", {
        "scope": "EXACT_PRIMITIVITY_AND_FULL_A_SCAN_NO_RESIDUE_BUCKETS",
        "counts": dict(audit_counts), "rows": row_audits, "timing_seconds": timings,
        "original_nonzero_rows": len(nonzero_rows), "nonzero_sign_representatives": len(representatives),
        "full_candidate_sign_closure": True,
    })
    artifacts["geometry-rows-partial.jsonl"] = digest(journal)
    audit = [len(zero_rows), zero_audit["d_zero"], zero_audit["integral_inverse"], zero_audit["unit_rows"], len(zero)]
    return sorted(nonzero + zero), {
        "bottom_rows": len(bottom_rows), "nonzero_candidates": len(nonzero),
        "zero_candidates": len(zero), "candidate_count": len(nonzero) + len(zero),
        "c_zero_audit": audit, "trace_radius_summary": radius_summary,
        "nonzero_audit": dict(audit_counts), "timing_seconds": timings,
    }



def polygon_result(primary, center, candidates, progress, output, artifacts):
    """Classify all candidates without assuming the survivor set is empty."""
    exclusions, survivors, histogram = [], [], Counter()
    classification_digest = sha256()
    checkpoint = output / "survivors-partial.jsonl"
    checkpoint.touch(exist_ok=False)
    progress.report("polygon", force=True, completed=0, total=len(candidates))
    for n, gamma in enumerate(candidates, 1):
        inside, failure = primary.core.in_all_coordinate_polygons(gamma, center)
        if inside:
            assert failure is None
            survivors.append(gamma)
            with checkpoint.open("a") as stream:
                stream.write(json.dumps({"status": "INCOMPLETE_STAGE_SURVIVOR_CHECKPOINT",
                                         "gamma": gamma_record(gamma)},
                                        sort_keys=True, separators=(",", ":")) + "\n")
                stream.flush()
                os.fsync(stream.fileno())
        else:
            assert failure is not None
            side, embedding = failure
            histogram[failure] += 1
            exclusions.append({"gamma": gamma_record(gamma), "failure": [side, embedding]})
            flat = tuple(coefficient for entry in gamma for coefficient in entry)
            classification_digest.update((",".join(map(str, flat)) + f";{side},{embedding}\n").encode("ascii"))
        progress.report("polygon", completed=n, total=len(candidates), survivors=len(survivors))
    assert sum(histogram.values()) + len(survivors) == len(candidates)
    classification = {
        "histogram": [{"side_index": side, "embedding_index": embedding, "count": count}
                      for (side, embedding), count in sorted(histogram.items())],
        "excluded": exclusions, "survivors": canonical_rows(survivors),
    }
    artifacts["classification.json"] = atomic_json(output / "classification.json", classification)
    artifacts["survivors.json"] = atomic_json(output / "survivors.json", canonical_rows(survivors))
    artifacts["survivors-partial.jsonl"] = digest(checkpoint)
    return {"survivor_count": len(survivors), "classification_digest": classification_digest.hexdigest(),
            "polygon_histogram": classification["histogram"]}



def load_expected():
    """Read checking data only; it never supplies any enumerated object."""
    check_sources()
    expected = json.loads((ROOT / EXPECTED_RELATIVE).read_text())
    assert expected["schema_version"] == 1
    assert expected["scope"] == "R127_PIC1_MASK_ZERO_WITH_PRINCIPAL_CONTROL"
    assert set(expected["points"]) == set(SUPPORTED_POINTS)
    assert set(SUPPORTED_POINTS) < set(PILOTS)
    for point in expected["points"].values():
        assert set(point["counts"]) == set(COUNT_KEYS)
        assert set(point["canonical_artifact_sha256"]) == set(CANONICAL_ARTIFACTS)
        assert point["counts"]["survivor_count"] == 0
    return expected



def compare_expected(output, result, expected_point):
    """Check freshly computed counts and bytes after the full search ends."""
    observed_counts = {key: result[key] for key in COUNT_KEYS}
    assert observed_counts == expected_point["counts"], (observed_counts, expected_point["counts"])
    assert observed_counts["survivor_count"] == 0
    observed_hashes = {name: digest(output / name) for name in CANONICAL_ARTIFACTS}
    assert observed_hashes == expected_point["canonical_artifact_sha256"]
    assert all(result["artifact_sha256"][name] == value for name, value in observed_hashes.items())
    return observed_hashes



def finite_cpu_limits():
    """Keep inherited finite limits; cap an otherwise unlimited process."""
    soft, hard = resource.getrlimit(resource.RLIMIT_CPU)
    if hard in (-1, resource.RLIM_INFINITY):
        soft = 600 if soft in (-1, resource.RLIM_INFINITY) else min(soft, 600)
        resource.setrlimit(resource.RLIMIT_CPU, (soft, 601))
    return resource.getrlimit(resource.RLIMIT_CPU)



def validate_generated_domain(a_values, c_values, bottom_rows):
    """Check fresh original-coordinate domains before applying sign reduction."""
    for values in (a_values, c_values):
        assert values == sorted(values) and len(values) == len(set(values))
        assert all(len(value) == 6 and all(type(q) is int for q in value) for value in values)
        assert ZERO in values
        assert {tuple(-q for q in value) for value in values} == set(values)
    assert bottom_rows == sorted(bottom_rows) and len(bottom_rows) == len(set(bottom_rows))
    assert (ZERO, ZERO) in bottom_rows
    c_set = set(c_values)
    assert all(c in c_set and len(d) == 6 and all(type(q) is int for q in d) for c, d in bottom_rows)
    assert {negate_matrix(row) for row in bottom_rows} == set(bottom_rows)


def run_point(point_label, output, prepare_only):
    if point_label not in SUPPORTED_POINTS:
        raise ValueError("This candidate covers only the two named R127 mask-zero points")
    cpu_limits = finite_cpu_limits()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    progress = Progress(output, point_label)
    progress.report("initializing", force=True)
    try:
        driver_sha256 = digest(Path(__file__))
        executable = Path(sys.executable).resolve()
        executable_sha256 = digest(executable)
        started = time.process_time()
        primary = load_primary()
        stage_times = {"arithmetic_load_cpu": time.process_time() - started}
        progress.report("reconstruction", force=True)
        started = time.process_time()
        center, height, bounds, preparation = reconstruct(primary, point_label)
        stage_times["reconstruction_cpu"] = time.process_time() - started
        artifacts = {"preparation.json": atomic_json(output / "preparation.json", preparation)}
        started = time.process_time()
        num_scaling = num_unit_scaling(primary, bounds[1], preparation)
        artifacts["num-unit-scaling.json"] = atomic_json(output / "num-unit-scaling.json", num_scaling["audit"])
        den_scaling = den_unit_scaling(primary, bounds[0], preparation)
        artifacts["den-unit-scaling.json"] = atomic_json(output / "den-unit-scaling.json", den_scaling["audit"])
        entry_scalings = {}
        for kind, index in (("A", 2), ("C", 4)):
            entry_scalings[kind] = entry_unit_scaling(primary, bounds[index], preparation, kind)
            name = kind + "-unit-scaling.json"
            artifacts[name] = atomic_json(output / name, entry_scalings[kind]["audit"])
        stage_times["unit_verification_cpu"] = time.process_time() - started
        result = {
            "point_label": point_label, "backend": "primary", "trace_index": preparation["trace_index"],
            "status": "PREPARATION_ONLY_NO_GEOMETRY", "scope": "R127_PIC1_MASK_ZERO_WITH_PRINCIPAL_CONTROL",
            "promotion_status": "CANDIDATE_PENDING_PAIRED_REPLAY_AND_REVIEW",
            "assertions_enabled": __debug__, "source_sha256": SOURCE_PINS,
            "driver_sha256": driver_sha256, "artifact_sha256": artifacts,
            "expected_sha256": SOURCE_PINS[EXPECTED_RELATIVE],
            "runtime": {"python_executable": str(executable), "python_sha256": executable_sha256,
                        "python_version": sys.version, "sympy_version": primary.sp.__version__,
                        "cpu_limits_seconds": list(cpu_limits)},
        }
        if not prepare_only:
            started = time.process_time()
            a_values, c_values = entry_lists(primary, bounds, progress, output, artifacts, entry_scalings)
            stage_times["entry_enumeration_cpu"] = time.process_time() - started
            result.update(a_count=len(a_values), c_count=len(c_values))
            bottom_rows, bottom_counts = bottom_rows_only(primary, center, c_values, progress, output, artifacts, den_scaling)
            stage_times["bottom_enumeration_cpu"] = bottom_counts["bottom_stage_cpu_seconds"]
            result["den_visited_leaves"] = bottom_counts["visited_leaves"]
            validate_generated_domain(a_values, c_values, bottom_rows)
            fresh_inputs = {"a_values": a_values, "c_values": c_values, "bottoms": bottom_rows}
            candidates, counts = enumerate_candidates(primary, center, fresh_inputs, progress, output, artifacts, num_scaling)
            result.update(counts)
            started = time.process_time()
            result.update(polygon_result(primary, center, candidates, progress, output, artifacts))
            stage_times["polygon_cpu"] = time.process_time() - started
            # No expected count/hash influences any preceding stage or branch.
            expected_point = load_expected()["points"][point_label]
            result["canonical_artifact_sha256"] = compare_expected(output, result, expected_point)
            result["status"] = "POINT_PASS"
        result.setdefault("timing_seconds", {}).update(stage_times)
        check_sources()
        assert digest(Path(__file__)) == driver_sha256
        assert digest(executable) == executable_sha256
        result.update(sources_unchanged=True, executable_unchanged=True,
                      wall_seconds=time.monotonic() - progress.started,
                      cpu_seconds=time.process_time() - progress.cpu_started)
        atomic_json(output / "result.json", result)
        progress.report("completed", force=True, result_status=result["status"],
                        survivor_count=result.get("survivor_count"))
        print(result["status"] + " " + json.dumps(result, sort_keys=True, separators=(",", ":")), flush=True)
        return result
    except BaseException as exc:
        progress.report("failed", force=True, error=repr(exc))
        raise


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not use -O/-OO")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point", choices=SUPPORTED_POINTS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--check-sources", action="store_true")
    args = parser.parse_args()
    check_sources()
    if args.check_sources:
        print(json.dumps({
            "status": "SOURCES_CHECKED_NO_GEOMETRY", "supported_points": list(SUPPORTED_POINTS),
            "source_sha256": SOURCE_PINS, "driver_sha256": digest(Path(__file__)),
            "scope": "R127_PIC1_MASK_ZERO_WITH_PRINCIPAL_CONTROL",
        }, sort_keys=True, separators=(",", ":")), flush=True)
        return
    if args.point is None:
        parser.error("--point is required unless --check-sources is selected")
    if args.output is None:
        with TemporaryDirectory(prefix="r127-pic1-v80-") as temporary:
            run_point(args.point, Path(temporary) / "artifacts", args.prepare_only)
    else:
        run_point(args.point, args.output, args.prepare_only)


if __name__ == "__main__":
    main()
