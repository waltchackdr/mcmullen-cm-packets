# Reader commands

Run `python3 scripts/check-public-package.py` from the snapshot root to check
included bytes, claim navigation and local links. It has finite CPU and wall
limits and performs no mathematical computation or remote action.

`python3 scripts/check-certificate-inputs.py` checks, without running any
certificate, that the files the included programs pin or import are present.

`python3 scripts/build-papers.py --help` describes the cached paper build.
See [reproduction instructions](../docs/reproducing-results.md) for dependencies
and the distinction between document compilation and mathematical replay.
