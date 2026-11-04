# Papers

The four documents are [Paper I](packet-and-principal-models.tex),
[Paper II](quotients-and-descent.tex), the [Hodge note](hodge-classes.tex) and
the [computational supplement](arithmetic-and-computation.tex). Reading PDFs
are in [pdf/](pdf/); [shared/preamble.tex](shared/preamble.tex) supplies common
typography and theorem environments.

## Opening and editing sources

Open the complete `papers/` folder as the LaTeX project, then select one main
document below. Each main file contains its own `\documentclass` and
`\begin{document}`. Keep `shared/`, `references/` and `sections/` beside it.

| Section folder | Main document to compile |
| --- | --- |
| `sections/packet/` | [packet-and-principal-models.tex](packet-and-principal-models.tex) |
| `sections/quotients/` | [quotients-and-descent.tex](quotients-and-descent.tex) |
| `sections/hodge/` | [hodge-classes.tex](hodge-classes.tex) |
| `sections/supplement/` | [arithmetic-and-computation.tex](arithmetic-and-computation.tex) |

For example, `sections/packet/04-geometric-fibers.tex` is one section of Paper I,
so its preview is built through `packet-and-principal-models.tex`. Section and
bibliography files begin with a `% !TeX root` comment pointing to their parent.
Editors that do not recognize this convention need the main file selected
manually. A “no complete document” message when opening an isolated section
means the editor has not found that parent; it is not evidence of missing prose.
Shared macros serve all four papers and have no single main document.

Each document has its own abstract, contents, numbering and cited-only
bibliography. Companion references use I, II, H and S.
[statement-index.csv](statement-index.csv) records 140 statement locations;
the [paper-order claim index](../docs/claim-reading-order.md) connects them to
stable claim IDs. A location is not itself a proof certificate.

This edition changes publication-facing front matter, PDF metadata and editor
navigation relative to the retained reading edition. Navigation edits add root
comments and move the document-class declaration into each main file. The
mathematical text, references, labels and common macros are preserved. Exact
source transformations are recorded in [edition provenance](../edition-provenance.json).
The earlier split sources remain in
the private archive, not at their former paths in this distribution.

Use the [reproduction guide](../docs/reproducing-results.md) for builds and the
[overview](../docs/research-overview.md) for dependencies and reading order.
