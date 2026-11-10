#!/usr/bin/env python3
"""Check that the included certificate programs can find the files they pin.

This reads every Python program under certificates/ and independent/ without
running it. It reports two kinds of missing input:

1. a sibling module the program imports that is not in this edition;
2. a file the program pins by SHA256 that is not in this edition.

A pinned file is a string naming a .py, .json or .jsonl file that is paired
with a 64-digit hexadecimal string in the same dictionary. It is looked up by
content hash, so its location does not matter. A pin whose path lies in a
directory this edition does not ship usually records a source of a historical
run; it is listed for information and does not fail the check, so this check
cannot show that such a file is never read.

The check runs no certificate and establishes no mathematical claim. It
exists so that an edition which leaves out a file its own programs need fails
loudly. Standard library only.
"""

import ast
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
PROGRAM_DIRECTORIES = ("certificates", "independent")
DECLARED_REQUIREMENTS = {"sympy", "mpmath"}
SKIPPED_TREES = {".git", "evidence"}
HEX64 = re.compile(r"[0-9a-f]{64}")
PINNED_NAME = re.compile(r"\.(py|json|jsonl)$")


def main():
    programs = sorted(path for name in PROGRAM_DIRECTORIES
                      for path in (ROOT / name).rglob("*.py"))
    local_modules = {path.stem for path in programs}
    present = set()
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if path.is_file() and relative.parts[0] not in SKIPPED_TREES:
            present.add(sha256(path.read_bytes()).hexdigest())
    standard = set(sys.stdlib_module_names)

    missing_imports, missing_pins, outside, unparsed = [], [], [], []
    for path in programs:
        program = path.relative_to(ROOT).as_posix()
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError as error:
            unparsed.append({"program": program, "error": str(error)})
            continue
        seen = set()
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                modules = [node.module.split(".")[0]]
            for module in modules:
                if (module not in standard and module not in DECLARED_REQUIREMENTS
                        and module not in local_modules and ("import", module) not in seen):
                    seen.add(("import", module))
                    missing_imports.append({"program": program, "module": module})
            if not isinstance(node, ast.Dict):
                continue
            for key, value in zip(node.keys, node.values):
                if not (isinstance(key, ast.Constant) and isinstance(key.value, str)
                        and PINNED_NAME.search(key.value)
                        and isinstance(value, ast.Constant) and isinstance(value.value, str)
                        and HEX64.fullmatch(value.value)):
                    continue
                if value.value in present or (key.value, value.value) in seen:
                    continue
                seen.add((key.value, value.value))
                row = {"program": program, "name": key.value, "sha256": value.value}
                parent = Path(key.value).parent
                if parent != Path(".") and not (ROOT / parent).is_dir():
                    outside.append(row)
                else:
                    missing_pins.append(row)

    passed = not (missing_imports or missing_pins or unparsed)
    print(json.dumps({
        "status": "CERTIFICATE_INPUTS_PRESENT" if passed else "CERTIFICATE_INPUTS_MISSING",
        "programs": len(programs),
        "missing_imports": missing_imports,
        "missing_pinned_files": missing_pins,
        "unparsed_programs": unparsed,
        "pins_outside_this_edition": outside,
        "scope": "Static presence check only; no certificate is run.",
    }, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
