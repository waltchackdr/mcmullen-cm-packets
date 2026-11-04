#!/usr/bin/env python3
"""Build the four reading documents with a cached Tectonic; never run certificates."""
import argparse
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from repository_tools import ROOT, PAPERS, DOCUMENTS, sha256, auxiliary_numbers

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--tectonic", default="tectonic", help="Installed executable; no download is performed")
parser.add_argument("--cache", required=True, help="Existing populated Tectonic cache, copied into the build")
parser.add_argument("--destination", required=True, help="New output directory; never overwritten")
args = parser.parse_args()
exe = shutil.which(args.tectonic)
if not exe:
    parser.error("Tectonic is unavailable. Install it separately; see docs/reproducing-results.md.")
cache = Path(args.cache).expanduser().resolve()
dest = Path(args.destination).resolve()
if not cache.is_dir() or dest.exists():
    parser.error("An existing cache and a nonexistent output directory are required.")
if cache == dest or cache.is_relative_to(dest):
    parser.error("Cache source must not be inside the new output directory.")
start = time.monotonic()
dest.mkdir(parents=True)
stage = dest / "sources"
shutil.copytree(PAPERS, stage, ignore=shutil.ignore_patterns("pdf", "*.pdf"))
shutil.copytree(cache, dest / "cache")
env = {**os.environ, "TECTONIC_CACHE_DIR": str(dest / "cache")}
report = {"status": "INCOMPLETE", "scope": "Document build and numbering only; no mathematical certificate execution.",
          "compiler": str(Path(exe).resolve()), "compiler_sha256": sha256(exe),
          "source_sha256": {str(p.relative_to(PAPERS)): sha256(p) for p in sorted(PAPERS.rglob('*.tex'))},
          "passes": [], "documents": []}

def cpu_used():
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime

try:
    for iteration in range(3):
        numbers = {}
        for key, name in DOCUMENTS.items():
            remain = 27 - (time.monotonic() - start)
            cpu = int(20 - cpu_used())
            if remain < 1 or cpu < 1:
                raise RuntimeError("Local build bound reached; partial outputs retained. Do not split a large replay into local batches.")
            out = dest / "pdf" / name
            out.mkdir(parents=True, exist_ok=True)
            def limits():
                resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu))
            command = [exe, "--only-cached", "--keep-logs", "--keep-intermediates", "--outdir", str(out), name + ".tex"]
            run = subprocess.run(command, cwd=stage, env=env, capture_output=True, text=True, timeout=remain, preexec_fn=limits)
            (out / f"command-{iteration}.log").write_text(run.stdout + run.stderr)
            report["passes"].append({"document": name, "pass": iteration, "returncode": run.returncode})
            if run.returncode:
                raise RuntimeError(f"Build failed for {name}; inspect retained command log (cache may be incomplete).")
            current = auxiliary_numbers(out / (name + ".aux"))
            duplicate = set(numbers) & set(current)
            if duplicate:
                raise RuntimeError(f"Duplicate companion labels: {sorted(duplicate)}")
            numbers.update(current)
        catalog = "% Generated from native auxiliary files; do not invent result numbers.\n"
        catalog += "".join("\\expandafter\\def\\csname xnumber@" + k + "\\endcsname{" + v + "}\n" for k, v in sorted(numbers.items()))
        target = stage / "shared/external-numbers.tex"
        changed = target.read_text() != catalog
        target.write_text(catalog)
        if not changed:
            break
    else:
        raise RuntimeError("Companion number catalog failed to stabilize in three passes.")
    for name in DOCUMENTS.values():
        out = dest / "pdf" / name
        log = (out / (name + ".log")).read_text()
        undefined = len(re.findall(r"(?:Reference|Citation).*?undefined", log))
        overfull = log.count("Overfull")
        if undefined or overfull or "LaTeX Warning: Label(s) may have changed" in log:
            raise RuntimeError(f"Unresolved references or layout overflow in {name}")
        pdf = out / (name + ".pdf")
        report["documents"].append({"document": name, "pdf": str(pdf), "sha256": sha256(pdf), "undefined": undefined, "overfull": overfull})
    report["status"] = "READING_DOCUMENTS_BUILT"
    report["catalog_matches_reader_source"] = (stage / "shared/external-numbers.tex").read_bytes() == (PAPERS / "shared/external-numbers.tex").read_bytes()
except Exception as exc:
    report["error"] = str(exc)
finally:
    report["elapsed_seconds"] = round(time.monotonic() - start, 3)
    report["child_cpu_seconds"] = round(cpu_used(), 3)
    (dest / "build-report.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
sys.exit(report["status"] != "READING_DOCUMENTS_BUILT")
