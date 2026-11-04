# Reproducing and inspecting the results

## First run

Download and unpack the complete edition, then open a terminal in the unpacked
folder containing `README.md`, `papers/` and `scripts/`. A single downloaded
Python or TeX file is not enough. The commands below use Python 3.9 or newer on
macOS or Linux; they need no account, credentials or paid service.

### Integrity and navigation

From the snapshot root, with Python 3.9 or newer on macOS or Linux:

```sh
python3 scripts/check-public-package.py
```

The bounded standard-library check validates the manifest, claim roster,
section anchors, paper labels and included reference bytes. It reports omitted
references separately. It executes no certificate and does not prove a theorem,
authenticate a historical execution or establish complete dependency closure.
The root `.git/` directory of a normal checkout is excluded from the content
manifest. No other hidden directory is exempted. Git history and configuration
require their own publication check; this command does not inspect them.

### Two exact mathematical examples

These included examples use only Python's standard library. No package
installation is needed. Run them from the same edition root:

```sh
python3 -B certificates/known_cm_comparison_v58.py
python3 -B independent/k11_vertex_bridge_stdlib_v60.py
```

Each prints a JSON report and exits with code zero on success. The expected
`status` values are, respectively, `PASS_EXACT_CM_COMPARISON_REGRESSION` and
`PASS_NEW_EXACT_K11_BRIDGE_STDLIB`. Both were run successfully from this edition.
The first checks supplied CM types and cyclotomic matrix identities; the second
checks the displayed arithmetic vertex bridge using rational arithmetic. Their
printed `scope` fields explain the limits: neither establishes the geometric
orbit claim or replaces its human proof. These are examples, not the full
packet campaign. Do not use Python's `-O` option: it disables assertions.

The `-B` option prevents bytecode files from being added to the checked edition.
To retain output, redirect it to a file outside the edition, for example
`python3 -B certificates/known_cm_comparison_v58.py > ../cm-comparison.json`.
A nonzero exit, traceback or different status is a failed run to investigate.

## Automatic checks on GitHub

The [workflow](../.github/workflows/verify.yml) runs the three commands above
on Python 3.12 and Ubuntu for pushes to `main` and pull requests. You can inspect
individual steps and their JSON output in the repository's
[Actions tab](https://github.com/waltchackdr/mcmullen-cm-packets/actions).
Maintainers can also start a run manually there.

The job has a five-minute limit, with a one-minute limit for each check. It uses
read-only repository permissions, no project secrets and no third-party Python
packages. It does not run the full historical campaign or build the papers.
A passing run confirms the stated package checks and two exact examples; it
does not certify all mathematical claims or constitute independent review.

The workflow itself is included in the manifest. It is the only permitted
hidden configuration file in the distribution. The checker does not update
hashes automatically; see [contributing](../CONTRIBUTING.md) for how proposed
changes and revised manifests are handled.

## Rebuild the papers

For an editor preview, open the complete `papers/` folder and choose the main
file listed in the [paper guide](../papers/README.md). Section files are inputs
to those documents. They carry `% !TeX root` hints for editors that support
them; other editors require choosing the main document explicitly.

Use an installed Tectonic executable and an existing populated cache:

```sh
python3 scripts/build-papers.py --tectonic /path/to/tectonic \
  --cache /path/to/existing/cache --destination build/reading-copy
```

Run this in a working copy. The helper uses cached resources and a new
destination, retains diagnostics and checks companion numbering. Compiler
subprocesses are bounded; inspect cache size before copying it. Rebuilt PDF
metadata may differ. Compare text, references, page counts and visual layout.

## Mathematical replay

Start with a precise statement in the papers, then follow its
[claim entry](claim-register.md) and [evidence availability](../evidence/README.md).
The included `requirements.txt` pins SymPy but does not specify the complete
environment for every historical certificate family. Some programs depend on
unshipped inputs or a separately assembled environment.

For a selected script that imports SymPy, create an environment outside the
edition and install the pinned dependency there:

```sh
python3 -m venv ../mcmullen-venv
../mcmullen-venv/bin/python -m pip install -r requirements.txt
```

This installation downloads SymPy and its dependencies. It is unnecessary for
the examples above and does not supply omitted historical inputs. Inspect the
chosen script's inputs and runtime requirements before executing it.

Define a finite workload, exact success criteria, source/runtime identity and
output retention before execution. Substantial enumeration belongs on suitable
compute resources with explicit time and cost limits. This edition supplies no
credentials and performs no automatic provisioning. Do not treat unrestricted
test discovery or a historical runner as a routine smoke check.

Preserve exact old versions. A missing input or changed source is a discrepancy
to investigate, not a reason to reset the expected hash. The
[limitations](limitations.md) distinguish missing historical versions from
deliberate distribution omissions.
