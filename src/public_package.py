"""Render and verify the derived reader edition; never execute certificates."""
import collections
import hashlib
import json
from pathlib import Path

from reader_links import check_links


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def cell(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def public_schema_errors(data, index):
    """Fail closed if private-record fields or omitted-object aliases reappear."""
    errors = []
    if set(data) != {"schema_version", "scope", "parent_register_sha256", "editorial_overrides_sha256", "entries"} or data.get("schema_version") != "public-reader-1":
        errors.append("Unexpected public register schema")
    keys = {"id", "title", "status", "qualification", "paper_locations", "related_claims", "no_standalone_statement_reason", "direct_evidence"}
    for entry in data.get("entries", []):
        if set(entry) != keys:
            errors.append("Unexpected claim fields")
        for item in entry.get("direct_evidence", []):
            if set(item) != {"sha256", "availability", "locator"}:
                errors.append("Unexpected reference fields")
        for item in entry.get("paper_locations", []):
            if set(item) != {"label", "line", "path", "relation", "source_sha256", "statement_id"}:
                errors.append("Unexpected paper-location fields")
    if set(index) != {"scope", "objects", "historical_missing_versions"}:
        errors.append("Unexpected evidence-index fields")
    for row in index.get("objects", []):
        expected = {"sha256", "bytes", "status"}
        if row.get("status") == "INCLUDED":
            expected.add("logical_paths")
        elif row.get("status") != "OMITTED_FROM_EDITION":
            errors.append("Unknown evidence availability")
        if set(row) != expected:
            errors.append("Unexpected object fields or omitted aliases")
    missing = index.get("historical_missing_versions", [])
    if len(missing) != 2 or {x.get("logical_path") for x in missing} != {"paper/main.tex", "docs/CLAIM_LEDGER.md"} or any(set(x) != {"logical_path", "sha256"} for x in missing):
        errors.append("Unexpected missing-version metadata")
    return errors


def render_claims(data):
    rows = ["# Claim register", "", "A derived reader view with stable identifiers. Read the [claims guide](claims.md)",
            "and [limitations](limitations.md) for scope. The [paper-order index](claim-reading-order.md)",
            "uses the same records. This register is navigation, not a new proof review.", "",
            "| ID | Claim | Disposition |", "| --- | --- | --- |"]
    reading = ["# Claims in paper reading order", "", "Generated from the same [register](claim-register.md). Relations distinguish",
               "statements, inputs, consequences and refutations; a location is not a certificate.", ""]
    for entry in data["entries"]:
        rows.append(f"| [{entry['id']}](#{entry['id'].lower()}) | {cell(entry['title'])} | {entry['status']} |")
    for entry in data["entries"]:
        rows += ["", "## " + entry["id"], "", entry["title"], "", "**Disposition:** " + entry["status"], "", entry["qualification"], ""]
        for location in entry["paper_locations"]:
            rows.append(f"- [{location['label']}](../{location['path']}), line {location['line']}; {location['relation']}.")
        if not entry["paper_locations"]:
            rows.append(entry["no_standalone_statement_reason"])
        if entry["related_claims"]:
            rows += ["", "Related claims: " + ", ".join(f"[{x}](#{x.lower()})" for x in entry["related_claims"]) + "."]
        rows += ["", "Exact reference digests and inclusion status are recorded in [the data](../data/claims/register.json)."]
    for key, title in [("packet", "Paper I"), ("quotients", "Paper II"), ("hodge", "Hodge note"), ("supplement", "Supplement")]:
        locations = collections.defaultdict(list)
        for entry in data["entries"]:
            for item in entry["paper_locations"]:
                if "/" + key + "/" in item["path"]:
                    locations[(item["path"], item["line"], item["label"])].append((entry["id"], item["relation"]))
        reading += ["", "## " + title, "", "| Paper location | Claim and relationship |", "| --- | --- |"]
        for (path, line, label), claims in sorted(locations.items()):
            links = ", ".join(f"[{name}](claim-register.md#{name.lower()}) ({relation})" for name, relation in claims)
            reading.append(f"| [{label}](../{path}) · line {line} | {links} |")
    reading += ["", "## Entries without a standalone current statement", ""]
    for entry in data["entries"]:
        if not entry["paper_locations"]:
            reading.append(f"- [{entry['id']}](claim-register.md#{entry['id'].lower()}): {entry['no_standalone_statement_reason']}")
    reading += ["", "The [statement index](../papers/statement-index.csv) also includes supporting",
                "blocks without a separate claim ID. No additional IDs are inferred."]
    return {"docs/claim-register.md": "\n".join(rows) + "\n", "docs/claim-reading-order.md": "\n".join(reading) + "\n"}


def verify(root):
    root = Path(root).resolve()
    errors = []
    manifest = json.loads((root / "manifest.json").read_text())
    names = set()
    for item in manifest["files"]:
        name = item["path"]
        path = root / name
        if name in names or not path.resolve().is_relative_to(root) or path.is_symlink():
            errors.append("Unsafe or duplicate path: " + name)
            continue
        names.add(name)
        if not path.is_file() or digest(path.read_bytes()) != item["sha256"] or path.stat().st_size != item["bytes"]:
            errors.append("File identity differs: " + name)
    # A normal checkout has root Git metadata. It is outside the release
    # manifest; reject indirection and exclude only this exact directory.
    git = root / ".git"
    if git.is_symlink() or (git.exists() and not git.is_dir()):
        errors.append("Unsupported Git metadata indirection")
    actual = {p.relative_to(root).as_posix() for p in root.rglob("*")
              if p.is_file() and p.relative_to(root).parts[0] != ".git"}
    if actual != names | {"manifest.json"}:
        errors.append("Inventory differs from manifest")
    # This one manifested workflow is intentional public configuration.
    # Other hidden files remain forbidden, including other workflows.
    allowed_hidden_files = {".github/workflows/verify.yml"}
    if any(x.startswith(".") or x == "__pycache__" for name in actual - allowed_hidden_files for x in Path(name).parts) or any(name.endswith(".pyc") for name in actual):
        errors.append("Hidden configuration or history in distribution")
    data = json.loads((root / "data/claims/register.json").read_text())
    index = json.loads((root / "evidence/index.json").read_text())
    errors.extend(public_schema_errors(data, index))
    expected = [f"C{x:02}" for x in range(1, 103)] + [f"O{x:02}" for x in range(1, 5)] + [f"P{x:02}" for x in range(1, 13)]
    if [x["id"] for x in data["entries"]] != expected:
        errors.append("Claim roster differs")
    allowed_statuses = {"RECORDED_PROOF_WITH_PREMISES", "ADOPTED_SCOPED", "CONDITIONAL_EXECUTION", "REFUTED_AS_STATED", "SUPERSEDED_PARTIAL_RESULT", "SUPERSEDED_TRACKING", "HISTORICAL_PREPARATION"}
    availability = collections.Counter()
    for entry in data["entries"]:
        if entry["status"] not in allowed_statuses or not entry["qualification"]:
            errors.append("Missing status/qualification: " + entry["id"])
        if any(x not in expected for x in entry["related_claims"]):
            errors.append("Unknown related claim: " + entry["id"])
        for location in entry["paper_locations"]:
            source = root / location["path"]
            if not source.resolve().is_relative_to(root) or not source.is_file():
                errors.append("Missing paper source: " + entry["id"])
                continue
            raw = source.read_bytes()
            text = raw.decode()
            marker = "\\label{" + location["label"] + "}"
            if digest(raw) != location["source_sha256"] or text.count(marker) != 1 or text.count("\n", 0, text.find(marker)) + 1 != location["line"]:
                errors.append("Paper location differs: " + entry["id"])
        for reference in entry["direct_evidence"]:
            availability[reference["availability"]] += 1
            locator = reference.get("locator")
            if locator:
                source = root / locator
                if not source.resolve().is_relative_to(root) or not source.is_file() or digest(source.read_bytes()) != reference["sha256"]:
                    errors.append("Included reference differs: " + entry["id"])
            elif reference["availability"] != "OMITTED_FROM_EDITION":
                errors.append("Invalid omission state: " + entry["id"])
    for name, expected_text in render_claims(data).items():
        if (root / name).read_text() != expected_text:
            errors.append("Stale generated view: " + name)
    for row in index["objects"]:
        path = root / ("evidence/objects/" + row["sha256"])
        if row["status"] == "INCLUDED":
            if not path.is_file() or digest(path.read_bytes()) != row["sha256"]:
                errors.append("Object identity differs: " + row["sha256"])
        elif path.exists():
            errors.append("Omitted object present: " + row["sha256"])
    links = fragments = 0
    cache = {}
    for name in manifest["reader_guides"]:
        checked, anchors, issues = check_links(root / name, root, cache)
        links += checked
        fragments += anchors
        errors.extend(issues)
    return {"status": "PUBLIC_EDITION_INTEGRITY_CONSISTENT" if not errors else "FAIL", "files": len(names),
            "claims": len(data["entries"]), "paper_mappings": sum(len(x["paper_locations"]) for x in data["entries"]),
            "links": links, "markdown_anchors": fragments, "reference_occurrences": dict(availability), "errors": errors,
            "scope": "Included identities and derived navigation only. No proof, complete runtime closure, external endorsement or privacy guarantee."}
