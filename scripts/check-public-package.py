#!/usr/bin/env python3
"""Check a derived publication preview without executing mathematical code."""
import argparse
import json
from pathlib import Path
import resource
import signal
import sys

resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
signal.alarm(28)
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from public_package import verify

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("directory", nargs="?", default=str(Path(__file__).resolve().parents[1]))
args = parser.parse_args()
try:
    report = verify(args.directory)
except (OSError, ValueError, KeyError) as exc:
    report = {"status": "FAIL", "error_type": type(exc).__name__}
print(json.dumps(report, indent=2))
sys.exit(report["status"] == "FAIL")
