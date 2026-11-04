# Contributing

Corrections, questions and extensions are welcome. You can [open an issue](https://github.com/waltchackdr/mcmullen-cm-packets/issues)
to report an error or discuss an idea. Include the paper, theorem, claim ID or
file involved, and enough detail for someone else to follow your reasoning.
For a Python problem, include the command, Python version and error message.
Please leave credentials and personal information out of reports.

To propose an edit, fork the repository, make your changes on a branch and open
a pull request. A fork is your own copy; a pull request proposes changes to this
repository. The maintainer reviews proposals and decides what to merge.
Public access does not grant permission to edit this repository directly.

## Preparing a contribution

The [overview](docs/research-overview.md), [claims guide](docs/claims.md) and
[limitations](docs/limitations.md) explain the package and its scope.
Keep changes focused and describe what they improve.

For a mathematical change, identify the affected claims, hypotheses, evidence
and dependencies. Distinguish conjectures, computations and proofs, and retain
counterexamples and failed approaches when they explain a result. Preserve
existing theorem labels and historical certificate identifiers so references
remain usable. If a certificate needs correction, retain its original version
and document the replacement and affected claims.

Use descriptive filenames for new prose, four-space Python indentation and
UTF-8 text with LF line endings and a final newline.

## Checking changes

Follow the [reproduction guide](docs/reproducing-results.md) to run the package
check and two small exact Python examples. GitHub Actions runs these same
checks for pushes to `main` and pull requests. A passing run checks package
consistency and those examples; it is not a review of every proof.

The manifest records the files and hashes of a particular snapshot. Changes
will intentionally fail the integrity check until the manifest is updated.
You can submit a proposal with that expected failure: explain which files you
changed. The maintainer reviews the changes, preserves the previous snapshot
in Git history and updates the manifest before merging. Changes to mathematical
claims also need their corresponding evidence and review records updated;
regenerating a manifest alone does not validate them.

## Reuse and attribution

See [licensing](LICENSE.md) for MIT code/data and CC BY 4.0 writing terms, and
[third-party notices](THIRD-PARTY-NOTICES.md) for exclusions. Contribute material
you have the right to share under the applicable terms, and acknowledge any
sources. Citing third-party material does not itself grant permission to
redistribute it.
