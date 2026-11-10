#!/usr/bin/env python3
"""Portable independent R127 Pic1 exclusion candidate with principal control.

Fresh point reconstruction, A/C unit entry enumeration, DEN bottom rows,
independent integer primitive indices, certified affine maps, both determinant
branches and the original exact polygon classifier are recomputed. Expected
counts/hashes are checked only after completion. No experiment/build file is
read or imported; no cached finite list defines the search. This is the retained
standard-library/Sturm arithmetic lineage, not a third independent proof.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import resource
import sys
import time
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PILOTS = {
    "R7:pic0:00": ("R7", 0, 0, 0, (7, 1, 1)),
    "R7sq:pic0:00": ("R7sq", 0, 12, 0, (7, 1, 2)),
    "R41_22:pic0:00": ("R41_22", 0, 5, 0, (41, 22, 1)),
    "R43_20:pic0:00": ("R43_20", 0, 4, 0, (43, 20, 1)),
    "R43_31:pic0:00": ("R43_31", 0, 8, 0, (43, 31, 1)),
    "R127:pic0:00": ("R127", 0, 10, 0, (127, 94, 1)),
    "R127:pic1:00": ("R127", 1, 10, 9, (127, 94, 1)),
}
ALLOWED_POINTS = ("R127:pic0:00", "R127:pic1:00")
EXPECTED_PATH = "certificates/r127_pic1_v80/expected.json"
CANONICAL_ARTIFACTS = ("preparation.json", "A-entries.json", "C-entries.json",
                       "bottom-rows.json", "candidates-nonzero.json", "candidates-zero.json",
                       "classification.json", "survivors.json")
EMBEDDINGS = (1, 5, 11, 13, 17, 19)
ZERO = (0,) * 6
NUM_WITNESS = "certificates/broader_principal_v75/unit_witnesses.json"
ENTRY_WITNESS = "certificates/r127_pic1_v80/entry_unit_witnesses.json"
ONE = (1,0,0,0,0,0)
SOURCE_PINS = {
    "certificates/broader_principal_v75/pilot_inputs.json": "0eb217f70238f79308f589cc15a6de7c317e6962a712bd5598b2b6aaa638c386",
    "certificates/broader_principal_v75/unit_witnesses.json": "040ddafb120f5144ecc31ccb3e532121c9f8d588e58bd3483522659e52c43977",
    "certificates/d7_r83_principal_expected_v55.py": "0fd36a67b004e6a250188aa8731a207456b5e49faa955005d6b12cd36665bea7",
    "certificates/full_packet_v62/ambient_points.json": "c0e84da6c266c087cbce65a6e8aca8e3f30350ce56b5a1089c4ae4978a381f36",
    "certificates/full_packet_v62/module_templates.json": "09ccff217511a67e946567135f9ffb74da020843621d0895ecce8746d25c245d",
    "certificates/r127_pic1_v80/entry_unit_witnesses.json": "9899ce81e27f2aadbf8bf41ef5ef435e0e92ec9ad5e099378326e07f220f37f5",
    "certificates/r127_pic1_v80/expected.json": "19a23cd3a229dc8c71825166f9e83d13f09d6ece389d3859df49755e9e317cfc",
    "independent/d7_ambient_points_sturm_v45.py": "3eb68f57e6d0d88b2e548f605742d434c8714b25eb1404b3c27d32415b61f42a",
    "independent/d7_r83_principal_sturm_v55.py": "894339e82389ff651bf25b04693821c4dfeedc427dec2fa4fd8993f195bb663f",
    "independent/quadrilateral_k17_stdlib_v41.py": "c5bd45f334a028820a071b5c758acfb8e48afc1e8454b18a3174dd04ca58a38c",
    "independent/residual_orbit_k17_k38_sturm_v42.py": "12a0d875e2715fe2c2ae123d2dc5202f4d7ce4e024f114ad163126654219fc13"
}


def file_digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_pins():
    for relative, expected in SOURCE_PINS.items():
        assert file_digest(ROOT / relative) == expected, relative


def coefficients(value):
    return [str(Q(q)) for q in value.coefficients]


def vector(value):
    assert len(value) == 6
    return [str(Q(q)) for q in value]


def gamma(value):
    return [vector(entry) for entry in value]


def matrix(value):
    return [coefficients(entry) for row in value for entry in row]


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def write_json(directory, name, value):
    with (directory/name).open("x") as stream:
        checked_text_write(stream, canonical_json(value))


class Progress:
    def __init__(self,directory,label):
        self.label=label
        self.stream=(directory/"progress.jsonl").open("x")
        self.stage_cpu={}

    def emit(self,event,**data):
        record={"event":event,"point_label":self.label,"utc":datetime.now(timezone.utc).isoformat(),
                "cpu_seconds":time.process_time(),**data}
        line=canonical_json(record)
        checked_text_write(self.stream,line)
        self.stream.flush()
        print(line,end="",flush=True)

    @contextmanager
    def stage(self,name):
        start=time.process_time()
        self.emit("START_STAGE",stage=name)
        try:
            yield
        except BaseException:
            self.emit("INTERRUPTED_STAGE",stage=name,stage_cpu_seconds=time.process_time()-start)
            raise
        self.stage_cpu[name]=time.process_time()-start
        self.emit("COMPLETE_STAGE",stage=name,stage_cpu_seconds=self.stage_cpu[name])


def load_foundation():
    check_pins()
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(ROOT / "independent"))
    import d7_r83_principal_sturm_v55 as retained
    for relative in SOURCE_PINS:
        if relative.endswith(".py"):
            module = sys.modules[Path(relative).stem]
            assert Path(module.__file__).resolve() == ROOT / relative
    assert not any(name == "sympy" or name.startswith("sympy.") for name in sys.modules)
    assert retained.ambient.EMBEDDINGS == retained.orbit.EMBEDDINGS == EMBEDDINGS
    return retained


def cyclotomic_units(field):
    t = field.T
    return (t**4-3*t*t+1, t**3+t*t-2*t-1,
            t**5-4*t**3+3*t+1, t**3-2*t, t)


def unit_from_exponents(field, powers):
    assert len(powers) == 5 and all(type(power) is int for power in powers)
    result = field.ONE
    for generator, power in zip(cyclotomic_units(field), powers):
        assert generator.integral_power_basis() and generator.inverse().integral_power_basis()
        assert abs(generator.norm()) == 1
        result *= generator**power
    assert result.integral_power_basis() and result.inverse().integral_power_basis()
    assert abs(result.norm()) == 1
    return result


def reconstruct(retained, point_label):
    """Derive raw module and selected SL2 tuple before comparing approved data."""
    assert point_label in PILOTS
    stratum, picard, index, slope, conductor = PILOTS[point_label]
    field, orbit = retained.ambient, retained.orbit
    t, one = field.T, field.ONE
    alpha = field.alpha_pair(index)
    aa, g = alpha
    assert field.pair_norm(alpha) == field.PI
    assert abs(g.norm()) == conductor[0]**conductor[2]
    assert field.element_modulus_value(g, conductor[1], conductor[0]) == 0
    # In (1+s*z,g*z), a pair r+p*z has coordinates (r,(p-s*r)/g).
    # Pic1 is nonprincipal as a CM-order ideal; this remains an O_L-free basis.
    basis = ((one, field.LElement(slope)), (field.ZERO, g))
    images = [field.pair_multiply(alpha, entry) for entry in basis]
    columns = [(real, (zpart-slope*real)/g) for real, zpart in images]
    raw = field.matrix(columns[0][0], columns[1][0], columns[0][1], columns[1][1])
    module_basis = field.matrix(1, 0, slope, g)
    assert field.mmul(module_basis, raw) == field.mmul(field.multiplication_matrix(alpha), module_basis)
    assert raw == field.matrix(aa-slope*g, -g*g, 1+slope*t+slope*slope, aa+(t+slope)*g)
    assert field.all_entries_integral(raw)
    assert field.mdet(raw) == field.PI and field.mtrace(raw) == field.pair_trace(alpha)
    raw_x, raw_y2 = field.exact_fixed_invariants(raw)
    n = 1+slope*t+slope*slope
    raw_modulus = raw_x*raw_x+raw_y2
    assert raw_x == -g*(t+2*slope)/(2*n)
    assert raw_y2 == g*g*(4-t*t)/(4*n*n) and raw_modulus == g*g/n
    assert orbit.signs(raw_y2) == orbit.signs(raw_modulus) == (1,)*6

    # Compare the derived RAW basis and action with the promoted arithmetic census.
    manifest = json.loads((ROOT / "certificates/full_packet_v62/ambient_points.json").read_text())
    matches = [row for row in manifest if row["point_label"] == point_label]
    assert len(matches) == 1
    row = matches[0]
    assert (row["stratum"], row["picard_class"], row["orientation_mask"], row["trace_for_matrix"]) == (stratum, picard, 0, index)
    assert row["conductor"] == list(conductor)
    assert row["root_partner"] == [2 if picard == 1 else 0, 16]
    assert matrix(raw) == [[str(q) for q in entry] for entry in row["norm_pi_matrix"]]
    assert coefficients(2*raw_x) == row["twice_real_part"]
    assert coefficients(4*raw_y2) == row["four_times_height_squared"]
    assert coefficients(raw_modulus) == row["modulus_squared"]
    templates = json.loads((ROOT / "certificates/full_packet_v62/module_templates.json").read_text())
    matches = [row for row in templates if row["stratum"] == stratum and row["picard_class"] == picard]
    assert len(matches) == 1
    template = matches[0]
    assert template["slope"] == slope and template["representative_trace"] == index
    assert matrix(module_basis) == [[str(q) for q in entry] for entry in template["basis_columns_in_one_z"]]
    assert matrix(raw) == [[str(q) for q in entry] for entry in template["norm_pi_action"]]

    approved = json.loads((ROOT / "certificates/broader_principal_v75/pilot_inputs.json").read_text())
    assert len(approved["pilots"]) == 7 and {row["point_label"] for row in approved["pilots"]} == set(PILOTS)
    row = next(row for row in approved["pilots"] if row["point_label"] == point_label)
    assert (row["stratum"], row["picard_class"], row["orientation_mask"], row["trace_index"], row["conductor"]) == (stratum, picard, 0, index, list(conductor))
    selected = row["selected"]
    assert type(selected["inversion_first"]) is bool
    v = unit_from_exponents(field, selected["unit_exponents"])
    diagonal = field.matrix(v, 0, 0, v.inverse())
    inversion = field.matrix(0, -1, 1, 0)
    conjugator = field.mmul(diagonal, inversion) if selected["inversion_first"] else diagonal
    assert field.all_entries_integral(conjugator) and field.all_entries_integral(field.minverse(conjugator))
    assert field.mdet(conjugator) == one
    action = field.conjugate_by(conjugator, raw)
    assert field.all_entries_integral(action)
    assert field.mtrace(action) == field.mtrace(raw) and field.mdet(action) == field.PI
    center, height_squared = field.exact_fixed_invariants(action)
    modulus_squared = center*center+height_squared
    # Independent direct Mobius transport of the raw complex tuple.
    a, b, c, d = (entry for matrix_row in conjugator for entry in matrix_row)
    den = (c*raw_x+d)**2+c*c*raw_y2
    assert orbit.signs(den) == (1,)*6
    assert center == ((a*raw_x+b)*(c*raw_x+d)+a*c*raw_y2)/den
    assert height_squared == raw_y2/den**2
    assert modulus_squared == ((a*raw_x+b)**2+a*a*raw_y2)/den
    assert orbit.signs(height_squared) == orbit.signs(modulus_squared) == (1,)*6
    point = retained.Point(action, center, height_squared, modulus_squared)
    bounds = retained.derive_bounds(point)
    preparation = {
        "point_label": point_label, "trace_index": index, "matrix": matrix(action),
        "center": coefficients(center), "height_squared": coefficients(height_squared),
        "modulus_squared": coefficients(modulus_squared),
        "bounds": {key: list(map(str, values)) for key, values in bounds.items()},
        "conjugator": matrix(conjugator),
    }
    assert selected["unit"] == coefficients(v) and selected["unit_inverse"] == coefficients(v.inverse())
    assert [[str(q) for q in entry] for entry in selected["conjugator"]] == matrix(conjugator)
    assert matrix(raw) == [[str(q) for q in entry] for entry in row["raw"]["matrix"]]
    assert coefficients(raw_x) == row["raw"]["center"]
    assert coefficients(raw_y2) == row["raw"]["height_squared"]
    assert coefficients(raw_modulus) == row["raw"]["modulus_squared"]
    for key in ("center", "height_squared", "modulus_squared", "bounds"):
        assert preparation[key] == selected["point"][key], key
    assert preparation["matrix"] == [[str(q) for q in entry] for entry in selected["point"]["matrix"]]
    retained.install_point_data(point, bounds)
    assert orbit.Y2 == point.height_squared
    for name in ("DEN", "NUM", "A", "B", "C", "D"):
        assert getattr(orbit, name+"_BOUND") == bounds[name]
    assert orbit.DEN_RADIUS == orbit.sum_of_squares(bounds["DEN"])
    assert orbit.NUM_RADIUS == orbit.sum_of_squares(bounds["NUM"])
    assert orbit.Y2_SCALED == tuple(orbit.scaled_interval(orbit.interval_eval(point.height_squared, root)) for root in orbit.ROOTS)
    return point, bounds, preparation


def inverse_map(value, columns):
    return tuple(sum(columns[j][i]*value[j] for j in range(6)) for i in range(6))


def certify_entry_unit(retained, point, bounds, point_label, kind):
    assert point_label in PILOTS
    assert kind in ("DEN", "NUM", "A", "C")
    field, orbit = retained.ambient, retained.orbit
    witness_path = ENTRY_WITNESS if kind in ("A", "C") else NUM_WITNESS
    witness_data = json.loads((ROOT / witness_path).read_text())
    assert set(witness_data["pilots"]) == set(PILOTS)
    record = witness_data["pilots"][point_label]
    assert record["point_label"] == point_label and record["ordered_embeddings"] == list(EMBEDDINGS)
    assert record["unit_generator_order"] == ["xi5", "xi11", "xi13", "xi17", "xi19=t"]
    expected_point = {"matrix": matrix(point.matrix), "center": coefficients(point.center),
                      "height_squared": coefficients(point.height_squared),
                      "modulus_squared": coefficients(point.modulus_squared)}
    for key, expected in expected_point.items():
        assert record["point"][key] == expected
    selected = record["witnesses"][kind]["selected"]
    t = field.T
    u = unit_from_exponents(field, selected["unit_exponents"])
    inverse = u.inverse()
    assert coefficients(u) == selected["unit"] and coefficients(inverse) == selected["unit_inverse"]
    assert u*inverse == field.ONE
    assert u.integral_power_basis() and inverse.integral_power_basis()
    assert abs(u.norm()) == abs(inverse.norm()) == 1

    def integral_coeff(value):
        assert value.integral_power_basis()
        return tuple(int(q) for q in value.coefficients)

    columns = [integral_coeff(u*t**j) for j in range(6)]
    inverse_columns = [integral_coeff(inverse*t**j) for j in range(6)]
    assert columns == [tuple(row) for row in selected["multiplication_matrix_columns"]]
    assert inverse_columns == [tuple(row) for row in selected["inverse_matrix_columns"]]
    for j in range(6):
        expected = tuple(int(i == j) for i in range(6))
        assert inverse_map(columns[j], inverse_columns) == expected
        assert inverse_map(inverse_columns[j], columns) == expected
    transformed = tuple(map(Q, selected["transformed_bounds"]))
    assert len(transformed) == 6 and all(bound > 0 for bound in transformed)
    for bound, original, root in zip(transformed, bounds[kind], orbit.ROOTS):
        box = orbit.interval_eval(u, root)
        assert bound >= max(abs(box.lo), abs(box.hi))*original
    gross = sum((bound*bound for bound in transformed), Q(0))
    assert gross == Q(selected["transformed_gross_radius"])
    diagnostic = {
        "scope": {"NUM": "C_ZERO_B_ENUMERATION_ONLY_NO_POINT_TRANSFORMATION",
                  "DEN": "BOTTOM_D_ENUMERATION_ONLY_NO_POINT_TRANSFORMATION",
                  "A": "A_ENTRY_ENUMERATION_ONLY_NO_POINT_TRANSFORMATION",
                  "C": "C_ENTRY_ENUMERATION_ONLY_NO_POINT_TRANSFORMATION"}[kind],
        "witness_path": witness_path, "witness_sha256": SOURCE_PINS[witness_path],
        "unit_exponents": selected["unit_exponents"], "unit": coefficients(u),
        "unit_inverse": coefficients(inverse), "multiplication_matrix_columns": columns,
        "inverse_matrix_columns": inverse_columns, "transformed_bounds": list(map(str, transformed)),
        "gross_radius": str(gross), "original_bounds": list(map(str, bounds[kind])),
    }
    return (u, inverse_columns, gross), diagnostic


def certify_num_unit(retained, point, bounds, point_label):
    return certify_entry_unit(retained, point, bounds, point_label, "NUM")


def certify_den_unit(retained, point, bounds, point_label):
    return certify_entry_unit(retained, point, bounds, point_label, "DEN")


def preconditioned_ball(orbit, center, radius, inverse_columns, predicate, visited):
    if radius < 0:
        return []

    def accept(w):
        visited[0] += 1
        return predicate(inverse_map(w, inverse_columns))

    values = orbit.enumerate_trace_ball(center, radius, accept)
    return [inverse_map(w, inverse_columns) for w in values]


def entries(retained, bounds, directory, progress, entry_units):
    orbit = retained.orbit
    result, entry_audit = {}, {}
    for name in ("C", "A"):
        with progress.stage(f"entries-{name}"):
            visited = [0]
            values = entry_values(retained, bounds, name, entry_units[name], visited)
            values.sort()
            assert len(values) == len(set(values))
            assert {tuple(-q for q in value) for value in values} == set(values)
            write_json(directory, f"{name}-entries.json", [vector(value) for value in values])
            entry_audit[name] = {"gross_radius": str(entry_units[name][2]),
                                 "lattice_points_visited": visited[0], "retained_entries": len(values)}
            progress.emit("ENTRY_COUNT", entry=name, count=len(values), lattice_points_visited=visited[0])
            result[name] = values
    write_json(directory, "entry-unit-audit.json", entry_audit)
    return result["A"], result["C"]


def enumerate_bottom_rows(retained,point,c_values,directory,progress,den_unit):
    """Fresh complete DEN-unit enumeration, retaining the original predicate."""
    orbit, x = retained.orbit, point.center
    x_boxes = tuple(orbit.scaled_interval(orbit.interval_eval(x, root)) for root in orbit.ROOTS)
    bottom_rows, bottom_radius_audit, den_unit_radius_audit = [], [], []
    den_u, den_inverse_columns, den_gross = den_unit
    with progress.stage("bottom-rows"):
        for index, c in enumerate(c_values, 1):
            ce = orbit.vector_element(c)
            # Preserve the v71 unweighted radius audit as a compatibility
            # diagnostic; enumeration uses the separately recorded DEN unit.
            fixed_trace = (ce**2 * point.height_squared).trace()
            assert fixed_trace >= 0
            remaining_radius = orbit.DEN_RADIUS - fixed_trace
            bottom_radius_audit.append({"c": vector(c), "fixed_trace": str(fixed_trace),
                                        "remaining_radius": str(remaining_radius)})
            c_boxes = tuple(orbit.scaled_vector_box(c, i) for i in range(6))
            center = tuple((-den_u*ce*x).coefficients)
            weighted_trace = ((den_u*ce)**2 * point.height_squared).trace()
            assert weighted_trace >= 0
            weighted_radius = den_gross - weighted_trace
            visited = [0]
            values = preconditioned_ball(
                orbit, center, weighted_radius, den_inverse_columns,
                lambda d, c=c, c_boxes=c_boxes: orbit.bottom_bounded(c, d, x, c_boxes, x_boxes),
                visited)
            assert len(values) == len(set(values))
            for d in values:
                assert all(type(q) is int for q in d)
                assert orbit.bottom_bounded(c, d, x, c_boxes, x_boxes)
                bottom_rows.append((c, d))
            den_unit_radius_audit.append({"c": vector(c), "center": list(map(str, center)),
                                          "fixed_trace": str(weighted_trace),
                                          "remaining_radius": str(weighted_radius),
                                          "lattice_points_visited": visited[0], "retained_d": len(values)})
            if index % 5 == 0 or index == len(c_values):
                progress.emit("BOTTOM_PROGRESS", c_done=index, c_total=len(c_values), bottom_rows=len(bottom_rows))
        bottom_rows.sort()
        assert len(bottom_rows) == len(set(bottom_rows))
        write_json(directory, "bottom-rows.json", [[vector(c), vector(d)] for c, d in bottom_rows])
        write_json(directory, "bottom-radius-audit.json", bottom_radius_audit)
        write_json(directory, "den-unit-radius-audit.json", den_unit_radius_audit)

    assert (ZERO,ZERO) in bottom_rows
    assert {negate_row(row) for row in bottom_rows}==set(bottom_rows)
    return bottom_rows



def classify(retained, point, values, directory, progress):
    survivors, excluded = [], []
    histogram = Counter()
    # A hard cap may interrupt the polygon stage. Preserve every positive
    # candidate immediately; this partial journal never certifies completeness.
    with progress.stage("polygon"), (directory / "survivors-partial.jsonl").open("x") as journal:
        for index, value in enumerate(sorted(values), 1):
            failure = retained.orbit.first_polygon_exclusion(value, point.center)
            if failure is None:
                survivors.append(value)
                journal.write(canonical_json({"status": "INCOMPLETE_STAGE_SURVIVOR_CHECKPOINT",
                                               "gamma": gamma(value)}))
                journal.flush()
            else:
                excluded.append((value, failure))
                histogram[failure] += 1
            if index % 100 == 0 or index == len(values):
                progress.emit("POLYGON_PROGRESS", candidates_done=index,
                              candidate_count=len(values), survivors=len(survivors))
        assert sum(histogram.values()) + len(survivors) == len(values)
        write_json(directory, "survivors.json", [gamma(value) for value in survivors])
        write_json(directory, "classification.json", {
            "histogram": [{"side_index": side, "embedding_index": embedding, "count": count}
                          for (side, embedding), count in sorted(histogram.items())],
            "excluded": [{"gamma": gamma(value), "failure": list(failure)} for value, failure in excluded],
            "survivors": [gamma(value) for value in survivors],
        })
        return {"survivor_count": len(survivors),
                "classification_digest": retained.orbit.canonical_digest(excluded)}


def check_expected(point_label, result, directory):
    """Compare only after independent reconstruction and complete enumeration."""
    expected = json.loads((ROOT / EXPECTED_PATH).read_text())
    assert expected["scope"] == "R127_PIC1_MASK_ZERO_WITH_PRINCIPAL_CONTROL"
    assert set(expected["points"]) == set(ALLOWED_POINTS)
    wanted = expected["points"][point_label]
    keys = {"a_count", "c_count", "bottom_rows", "nonzero_candidates", "zero_candidates",
            "c_zero_audit", "candidate_count", "classification_digest", "survivor_count"}
    assert set(wanted["counts"]) == keys
    observed = {key: result[key] for key in keys}
    assert observed == wanted["counts"], (point_label, observed, wanted["counts"])
    assert observed["candidate_count"] == observed["nonzero_candidates"] + observed["zero_candidates"]
    assert observed["c_zero_audit"][-1] == observed["zero_candidates"]
    assert observed["survivor_count"] == 0
    assert set(wanted["canonical_artifact_sha256"]) == set(CANONICAL_ARTIFACTS)
    actual = {name: file_digest(directory / name) for name in CANONICAL_ARTIFACTS}
    assert actual == wanted["canonical_artifact_sha256"], (point_label, actual)
    return actual


def cpu_limits():
    """Preserve an inherited finite limit, otherwise impose a bounded default."""
    soft, hard = resource.getrlimit(resource.RLIMIT_CPU)
    if hard == resource.RLIM_INFINITY:
        # Never weaken an existing finite soft limit supplied by a supervisor.
        soft = 600 if soft == resource.RLIM_INFINITY else min(soft, 600)
        resource.setrlimit(resource.RLIMIT_CPU, (soft, 601))
    limits = resource.getrlimit(resource.RLIMIT_CPU)
    if limits[0] <= 0 or limits[1] == resource.RLIM_INFINITY or limits[1] <= 0:
        raise SystemExit("Arithmetic requires a finite positive CPU limit")
    return limits




# Local port of the reviewed independent v76 integer lattice algorithm.
def multiplication_basis(field):
    """Derive integer column matrices of multiplication by 1,t,...,t^5."""
    basis = tuple(field.T**j for j in range(6))
    tables = []
    for value in basis:
        columns = []
        for basis_vector in basis:
            product = value*basis_vector
            assert product.integral_power_basis()
            columns.append(tuple(int(q) for q in product.coefficients))
        tables.append(tuple(columns))
    assert tables[0] == tuple(tuple(int(i == j) for i in range(6)) for j in range(6))
    return tuple(tables)

def multiplication_columns(value, tables):
    assert len(value) == 6 and all(type(q) is int for q in value)
    nonzero = [(k, q) for k, q in enumerate(value) if q]
    return tuple(tuple(sum(q*tables[k][j][i] for k, q in nonzero) for i in range(6)) for j in range(6))

def apply_columns(columns, value):
    assert len(columns) == len(value)
    return tuple(sum(column[i]*q for column, q in zip(columns, value)) for i in range(len(columns[0])))

def bezout(a, b):
    """Euclidean algorithm with signed Bezout coefficients and positive gcd."""
    assert type(a) is int and type(b) is int and (a or b)
    r0, r1, s0, s1, t0, t1 = abs(a), abs(b), 1, 0, 0, 1
    while r1:
        q, r = divmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, s0-q*s1
        t0, t1 = t1, t0-q*t1
    s, t = s0 if a >= 0 else -s0, t0 if b >= 0 else -t0
    assert r0 > 0 and s*a+t*b == r0
    return r0, s, t

def column_hnf(columns, with_transform=False):
    """Integer column HNF via independent unimodular Bezout elimination.

    Process rows from bottom to top. Pending columns vanish below the current
    row. A Bezout 2-column operation has determinant one and replaces the pivot
    by gcd(p,q), clearing its partner. Column swaps/signs and quotient additions
    are also unimodular. Thus the generator lattice is invariant. Once every
    row is processed, discarded columns are zero and retained pivots give a
    basis. Positive pivots and quotient reductions make upper column HNF.
    Euclid terminates by strict remainder decrease; the outer loops are finite.
    """
    assert columns and all(len(v) == len(columns[0]) for v in columns)
    assert all(type(q) is int for v in columns for q in v)
    n, m = len(columns[0]), len(columns)
    values = [list(v) for v in columns]
    transform = [[int(i == j) for i in range(m)] for j in range(m)] if with_transform else None
    edge = m-1
    pivot_rows = []
    for row in range(n-1, -1, -1):
        chosen = next((j for j in range(edge+1) if values[j][row]), None)
        if chosen is None:
            continue
        values[edge], values[chosen] = values[chosen], values[edge]
        if with_transform:
            transform[edge], transform[chosen] = transform[chosen], transform[edge]
        for j in range(edge):
            q = values[j][row]
            if not q:
                continue
            p = values[edge][row]
            g, s, t = bezout(p, q)
            old_p, old_q = values[edge], values[j]
            values[edge] = [s*a+t*b for a, b in zip(old_p, old_q)]
            values[j] = [-(q//g)*a+(p//g)*b for a, b in zip(old_p, old_q)]
            if with_transform:
                old_p, old_q = transform[edge], transform[j]
                transform[edge] = [s*a+t*b for a, b in zip(old_p, old_q)]
                transform[j] = [-(q//g)*a+(p//g)*b for a, b in zip(old_p, old_q)]
            assert values[edge][row] == g and values[j][row] == 0
        if values[edge][row] < 0:
            values[edge] = [-v for v in values[edge]]
            if with_transform:
                transform[edge] = [-v for v in transform[edge]]
        pivot = values[edge][row]
        assert pivot > 0
        for j in range(edge+1, m):
            quotient = values[j][row]//pivot
            values[j] = [a-quotient*b for a, b in zip(values[j], values[edge])]
            if with_transform:
                transform[j] = [a-quotient*b for a, b in zip(transform[j], transform[edge])]
            assert 0 <= values[j][row] < pivot
        assert all(values[j][row] == 0 for j in range(edge))
        assert all(values[j][i] == 0 for j in range(edge+1) for i in range(row+1, n))
        pivot_rows.append(row)
        edge -= 1
    assert all(not any(v) for v in values[:edge+1])
    retained_columns = values[edge+1:]
    rank = len(retained_columns)
    hnf = [[retained_columns[j][i] for j in range(rank)] for i in range(n)]
    index = None
    if rank == n:
        index = 1
        for i in range(n):
            assert hnf[i][i] > 0
            assert all(hnf[i][j] == 0 for j in range(i))
            assert all(0 <= hnf[i][j] < hnf[i][i] for j in range(i+1, n))
            index *= hnf[i][i]
    if with_transform:
        assert [apply_columns(columns, v) for v in transform] == [tuple(v) for v in values]
    return {"rank": rank, "index": index, "hnf": hnf,
            "transform_columns": transform, "zero_columns": edge+1}

def checked_text_write(stream, text):
    """Fail explicitly before exceeding the 128 MiB per-file budget; never truncate."""
    assert text.isascii()  # JSON uses ensure_ascii=True; lengths are exact bytes.
    soft, _hard = resource.getrlimit(resource.RLIMIT_FSIZE)
    limit = 128*1024*1024
    if soft != resource.RLIM_INFINITY:
        limit = min(limit, soft)
    required = stream.tell()+len(text)
    if required > limit:
        raise RuntimeError(f"Artifact {stream.name} requires {required} bytes, exceeding file limit {limit}; output is incomplete, never truncated")
    stream.write(text)

def row_json(row):
    return [vector(v) for v in row]

def negate_row(row):
    return tuple(tuple(-q for q in v) for v in row)

def negative_gamma(value):
    return tuple(tuple(-q for q in entry) for entry in value)

def certify_affine_identity(lattice, tables, c, d, affine, mc=None):
    """Prove c*b=a*d-1 for every integral output of the original affine map.

    Its formula is b=(R*a+offset)/D. Integer basis identities
    M_c*R=D*M_d and M_c*offset=-D*e0 certify the formula for ALL a,
    before any success enumeration, without per-success field multiplication.
    """
    assert c != ZERO
    denominator = affine.denominator
    assert type(denominator) is int and denominator > 0
    assert len(affine.rows) == len(affine.offsets) == 6
    assert all(len(row) == 6 for row in affine.rows)
    assert all(type(q) is int for row in affine.rows for q in row)
    assert all(type(q) is int for q in affine.offsets)
    mc = mc if mc is not None else lattice["multiplication_columns"](c, tables)
    md = lattice["multiplication_columns"](d, tables)
    for j in range(6):
        r_column = tuple(affine.rows[i][j] for i in range(6))
        actual = lattice["apply_columns"](mc, r_column)
        assert actual == tuple(denominator*q for q in md[j]), ("Affine coefficient identity failed", j)
    actual_offset = lattice["apply_columns"](mc, affine.offsets)
    assert actual_offset == tuple(-denominator*q for q in ONE), "Affine offset identity failed"
    return {"status": "EXACT_AFFINE_IDENTITY_VERIFIED_FOR_ALL_A",
            "denominator": str(denominator), "matrix_rows": [list(row) for row in affine.rows],
            "offsets": list(affine.offsets)}

def make_boxes(orbit,point,a_values):
    x_boxes=tuple(orbit.scaled_interval(orbit.interval_eval(point.center,root)) for root in orbit.ROOTS)
    a_boxes={a:tuple(orbit.scaled_vector_box(a,i) for i in range(6)) for a in a_values}
    return a_boxes,x_boxes

def entry_values(retained, bounds, name, unit_data, visited):
    """Enumerate q=0 unit coordinates, returning original A or C entries."""
    assert name in ("A", "C")
    _unit, inverse_columns, gross = unit_data
    return preconditioned_ball(
        retained.orbit, (Q(0),)*6, gross, inverse_columns,
        lambda value: retained.orbit.field_box_membership(value, bounds[name]), visited)

def nonzero_candidates(retained,point,a_values,bottom_rows,a_boxes,x_boxes,directory,progress):
    lattice={"multiplication_basis":multiplication_basis,"multiplication_columns":multiplication_columns,
             "apply_columns":apply_columns,"column_hnf":column_hnf}
    indices={}
    assert a_values==sorted(set(a_values))
    assert {tuple(-q for q in a) for a in a_values}==set(a_values)
    assert bottom_rows==sorted(set(bottom_rows))
    assert {negate_row(row) for row in bottom_rows}==set(bottom_rows)
    orbit=retained.orbit
    tables=lattice["multiplication_basis"](retained.ambient)
    selected=sorted({min(row,negate_row(row)) for row in bottom_rows if row[0]!=ZERO})
    assert len(selected)*2==sum(c!=ZERO for c,d in bottom_rows)
    cache={}
    half=[]
    diagnostics=[]
    totals={"selected_nonzero_rows":len(selected),"primitive_nonzero_rows":0,"nonprimitive_nonzero_rows":0,
            "unit_c_rows":0,"a_probes":0,"integral_matches":0,"numerator_rejections":0}
    timing={"c_ideal_setup":0.0,"row_primitive_index":0.0,"affine_setup_and_certificate":0.0,
            "a_scan_with_original_numerator":0.0,"accepted_determinant_verification":0.0}
    with progress.stage("determinant-c-nonzero"),(directory/"partial-nonzero-rows.jsonl").open("x") as journal:
        def checkpoint(event,**record):
            checked_text_write(journal,canonical_json({"event":event,**record}))
            journal.flush()
        for ordinal,(c,d) in enumerate(selected,1):
            row_start=time.process_time()
            checkpoint("START_ROW",ordinal=ordinal,bottom_row=row_json((c,d)))
            started=time.process_time()
            if c not in cache:
                mc=lattice["multiplication_columns"](c,tables)
                c_ideal=lattice["column_hnf"](mc)
                assert c_ideal["rank"]==6 and c_ideal["index"]>=1
                cache[c]=(mc,c_ideal["index"]==1,c_ideal["index"])
            mc,unit_c,c_index=cache[c]
            timing["c_ideal_setup"]+=time.process_time()-started
            started=time.process_time()
            if unit_c:
                actual_index=1
            else:
                reduced=lattice["column_hnf"](mc+lattice["multiplication_columns"](d,tables))
                assert reduced["rank"]==6
                actual_index=reduced["index"]
            assert type(actual_index) is int and actual_index>=1
            indices[(c,d)]=indices[negate_row((c,d))]=actual_index
            timing["row_primitive_index"]+=time.process_time()-started
            diagnostic={"ordinal":ordinal,"bottom_row":row_json((c,d)),"ideal_index":str(actual_index),
                        "c_ideal_index":str(c_index),"unit_c":unit_c,"a_probes":0,
                        "integral_matches":0,"numerator_rejections":0,"half_candidates":0}
            if actual_index!=1:
                totals["nonprimitive_nonzero_rows"]+=1
                diagnostic.update(status="NONPRIMITIVE_SKIPPED_BEFORE_A_SCAN",row_cpu_seconds=time.process_time()-row_start)
                diagnostics.append(diagnostic)
                checkpoint("COMPLETE_ROW",**diagnostic,accepted=[])
                progress.emit("NONZERO_ROW_COMPLETE",**diagnostic)
                continue
            totals["primitive_nonzero_rows"]+=1
            totals["unit_c_rows"]+=int(unit_c)
            started=time.process_time()
            affine=orbit.determinant_one_affine_map(c,d)
            certificate=certify_affine_identity(lattice,tables,c,d,affine,mc)
            # The certified offset is -c^{-1}; denominator one exactly detects
            # the integral inverse. Unit rows need no quotient or Bezout work.
            assert (affine.denominator==1)==unit_c
            timing["affine_setup_and_certificate"]+=time.process_time()-started
            diagnostic["affine_identity"]=certificate
            accepted=[]
            chunk=[]
            scan_start=time.process_time()
            determinant_cpu=0.0
            for a_index,a in enumerate(a_values):
                b=affine.evaluate(a)
                if b is not None:
                    diagnostic["integral_matches"]+=1
                    if orbit.top_bounded(a,b,point.center,a_boxes[a],x_boxes):
                        value=(a,b,c,d)
                        verify_start=time.process_time()
                        ae,be,ce,de=map(orbit.vector_element,value)
                        assert ae*de-be*ce==orbit.ONE
                        determinant_cpu+=time.process_time()-verify_start
                        accepted.append(value)
                        chunk.append(value)
                    else:
                        diagnostic["numerator_rejections"]+=1
                diagnostic["a_probes"]=a_index+1
                if (a_index+1)%4096==0:
                    checkpoint("ROW_PROGRESS",ordinal=ordinal,bottom_row=row_json((c,d)),
                               a_done=a_index+1,integral_matches=diagnostic["integral_matches"],
                               half_candidates=len(accepted),accepted_since_checkpoint=[gamma(v) for v in chunk],
                               row_cpu_seconds=time.process_time()-row_start)
                    chunk=[]
            scan_cpu=time.process_time()-scan_start
            timing["a_scan_with_original_numerator"]+=scan_cpu-determinant_cpu
            timing["accepted_determinant_verification"]+=determinant_cpu
            half.extend(accepted)
            diagnostic.update(status="COMPLETE_PRIMITIVE_ROW",half_candidates=len(accepted),
                              row_cpu_seconds=time.process_time()-row_start)
            for key in ("a_probes","integral_matches","numerator_rejections"):
                totals[key]+=diagnostic[key]
            diagnostics.append(diagnostic)
            checkpoint("COMPLETE_ROW",**diagnostic,accepted_since_checkpoint=[gamma(v) for v in chunk])
            progress.emit("NONZERO_ROW_COMPLETE",ordinal=ordinal,total=len(selected),
                          integral_matches=diagnostic["integral_matches"],half_candidates=len(accepted),unit_c=unit_c)
        assert len(half)==len(set(half))
        values=sorted(half+[negative_gamma(value) for value in half])
        assert len(values)==len(set(values))==2*len(half)
        assert {negative_gamma(value) for value in values}==set(values)
        assert all((value[2],value[3]) in indices for value in values)
        write_json(directory,"candidates-nonzero.json",[gamma(value) for value in values])
        write_json(directory,"nonzero-row-audit.json",diagnostics)
        write_json(directory,"nonzero-row-indices.json",[row_json(row)+[str(indices[row])] for row in sorted(indices)])
    totals.update(half_nonzero_candidates=len(half),nonzero_candidates=len(values),timing_seconds=timing)
    return values,totals

def zero_candidates(retained,point,a_values,bottom_rows,a_boxes,x_boxes,nonzero,directory,progress,num_unit):
    """Retain the complete original v73 c=0 NUM-unit mathematics and audit."""
    orbit,x=retained.orbit,point.center
    zero, zero_rows, zero_d, integral_inverse, a_admissible = [], 0, 0, 0, 0
    zero_radius_audit, unit_radius_audit = [], []
    u, inverse_columns, gross = num_unit
    a_set = set(a_values)
    with progress.stage("determinant-c-zero"), (directory/"partial-zero-rows.jsonl").open("x") as journal:
        def checkpoint(event, **record):
            checked_text_write(journal,canonical_json({"event":event,**record}))
            journal.flush()
        for c, d in bottom_rows:
            if c != ZERO:
                continue
            zero_rows += 1
            row_start=time.process_time()
            checkpoint("START_ZERO_ROW",ordinal=zero_rows,bottom_row=row_json((c,d)))
            if d == ZERO:
                zero_d += 1
                checkpoint("COMPLETE_ZERO_ROW",ordinal=zero_rows,reason="ZERO_D",candidates=[])
                continue
            ae = orbit.vector_element(d).inverse()
            if not ae.integral_power_basis():
                checkpoint("COMPLETE_ZERO_ROW",ordinal=zero_rows,reason="NONINTEGRAL_INVERSE",candidates=[])
                continue
            integral_inverse += 1
            a = orbit.integral_vector(ae)
            if not orbit.field_box_membership(a, orbit.A_BOUND):
                checkpoint("COMPLETE_ZERO_ROW",ordinal=zero_rows,reason="A_OUTSIDE_BOX",candidates=[])
                continue
            assert a in a_set
            a_admissible += 1
            # Preserve the original v68 radius audit, while the new weighted
            # radius used for enumeration is recorded in a separate artifact.
            fixed_trace = (ae**2 * point.height_squared).trace()
            assert fixed_trace >= 0
            remaining_radius = orbit.NUM_RADIUS - fixed_trace
            zero_radius_audit.append({"a": vector(a), "d": vector(d),
                                      "fixed_trace": str(fixed_trace),
                                      "remaining_radius": str(remaining_radius)})
            center = tuple((-u*ae*x).coefficients)
            weighted_trace = ((u*ae)**2 * point.height_squared).trace()
            assert weighted_trace >= 0
            weighted_radius = gross - weighted_trace
            visited = [0]
            values = preconditioned_ball(
                orbit, center, weighted_radius, inverse_columns,
                lambda b, a=a: orbit.top_bounded(a, b, x, a_boxes[a], x_boxes), visited)
            assert len(values) == len(set(values))
            for b in values:
                assert all(type(q) is int for q in b)
                assert orbit.top_bounded(a, b, x, a_boxes[a], x_boxes)
                value = (a, b, c, d)
                aa, bb, cc, dd = map(orbit.vector_element, value)
                assert cc == orbit.ZERO and aa*dd - bb*cc == orbit.ONE
                zero.append(value)
            checkpoint("COMPLETE_ZERO_ROW",ordinal=zero_rows,reason="ENUMERATED_NUM_BALL",
                       candidates=[gamma((a,b,c,d)) for b in values],lattice_points_visited=visited[0],
                       row_cpu_seconds=time.process_time()-row_start)
            unit_radius_audit.append({"a": vector(a), "d": vector(d),
                                      "center": list(map(str, center)),
                                      "fixed_trace": str(weighted_trace),
                                      "remaining_radius": str(weighted_radius),
                                      "lattice_points_visited": visited[0], "retained_b": len(values)})
            progress.emit("ZERO_PROGRESS", bottom_rows_done=zero_rows, candidates=len(zero),
                          lattice_points_visited=visited[0])
        zero.sort()
        assert zero_d == 1
        assert len(zero) == len(set(zero)) and set(zero).isdisjoint(nonzero)
        write_json(directory, "candidates-zero.json", [gamma(value) for value in zero])
        write_json(directory, "zero-radius-audit.json", zero_radius_audit)
        write_json(directory, "unit-radius-audit.json", unit_radius_audit)
    return zero, {
        "bottom_rows": len(bottom_rows), "zero_candidates": len(zero),
        "c_zero_audit": [zero_rows, zero_d, integral_inverse, a_admissible, len(zero)],
    }

def main():
    if not __debug__:
        raise SystemExit("Assertions must be enabled")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point", choices=ALLOWED_POINTS)
    parser.add_argument("--output", type=Path,
                        help="Fresh artifact directory; default is a retained temporary directory")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--prepare-only", action="store_true")
    modes.add_argument("--check-sources", action="store_true", help="Hash checks only; no arithmetic import")
    args = parser.parse_args()
    if not args.check_sources and args.point is None:
        parser.error("--point is required except with --check-sources")
    source_hash = file_digest(Path(__file__))
    check_pins()
    if args.check_sources:
        print("PASS_SOURCE_CHECK_ONLY " + canonical_json({
            "status": "PASS_SOURCE_CHECK_ONLY", "allowed_points": list(ALLOWED_POINTS),
            "driver_sha256": source_hash, "source_sha256": SOURCE_PINS}).rstrip(), flush=True)
        return
    limits = cpu_limits()
    started = time.monotonic()
    if args.output is None:
        directory = Path(tempfile.mkdtemp(prefix="r127-pic1-sturm-v80-"))
    else:
        directory = args.output.resolve()
        directory.mkdir(parents=True, exist_ok=False)
    progress = Progress(directory, args.point)
    with progress.stage("reconstruct-and-bound"):
        retained = load_foundation()
        point, bounds, preparation = reconstruct(retained, args.point)
        num_unit, num_diagnostic = certify_num_unit(retained, point, bounds, args.point)
        den_unit, den_diagnostic = certify_den_unit(retained, point, bounds, args.point)
        write_json(directory, "preparation.json", preparation)
        write_json(directory, "unit-preconditioning.json", num_diagnostic)
        write_json(directory, "den-unit-preconditioning.json", den_diagnostic)
        entry_units,entry_diagnostics={},{}
        for kind in ("A","C"):
            entry_units[kind],entry_diagnostics[kind]=certify_entry_unit(retained,point,bounds,args.point,kind)
        write_json(directory,"entry-unit-preconditioning.json",entry_diagnostics)
    result = {
        "schema_version": 1, "point_label": args.point, "implementation": "stdlib-sturm",
        "status": "PREPARATION_ONLY_NO_GEOMETRY",
        "scope": "R127_PIC1_MASK_ZERO_WITH_PRINCIPAL_CONTROL",
        "ordered_embeddings": list(EMBEDDINGS), "coefficient_order": "ascending powers 1,t,...,t^5",
        "matrix_entry_order": ["a", "b", "c", "d"],
        "driver_sha256": source_hash, "source_sha256": SOURCE_PINS,
        "python": sys.version, "cpu_limits": list(limits), "artifact_directory": str(directory),
    }
    if not args.prepare_only:
        a_values, c_values = entries(retained, bounds, directory, progress, entry_units)
        result.update(a_count=len(a_values), c_count=len(c_values))
        bottom_rows=enumerate_bottom_rows(retained,point,c_values,directory,progress,den_unit)
        with progress.stage("original-coordinate-boxes"):
            a_boxes,x_boxes=make_boxes(retained.orbit,point,a_values)
        nonzero,nonzero_audit=nonzero_candidates(retained,point,a_values,bottom_rows,a_boxes,x_boxes,directory,progress)
        write_json(directory,"nonzero-summary.json",nonzero_audit)
        zero,counts=zero_candidates(retained,point,a_values,bottom_rows,a_boxes,x_boxes,nonzero,directory,progress,num_unit)
        assert {negative_gamma(value) for value in zero}==set(zero)
        result.update(nonzero_candidates=len(nonzero))
        result.update(counts)
        result.update(candidate_count=len(nonzero)+len(zero))
        result.update(classify(retained, point, nonzero + zero, directory, progress))
        result["canonical_artifact_sha256"] = check_expected(args.point, result, directory)
        result["status"] = "POINT_PASS"
    check_pins()
    assert file_digest(Path(__file__)) == source_hash
    assert not any(name == "sympy" or name.startswith("sympy.") for name in sys.modules)
    result.update(sources_unchanged=True, cpu_seconds=time.process_time(),
                  wall_seconds=time.monotonic()-started,stage_cpu_seconds=progress.stage_cpu)
    progress.emit("COMPLETE_PILOT", status=result["status"])
    progress.stream.close()
    result["artifact_sha256"] = {path.name: file_digest(path) for path in sorted(directory.iterdir())}
    write_json(directory, "result.json", result)
    print(result["status"] + " " + canonical_json(result).rstrip(), flush=True)


if __name__ == "__main__":
    main()
