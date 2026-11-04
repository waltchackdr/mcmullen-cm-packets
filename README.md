# CM packets, principal models and quotient descent

Research preprints on abelian varieties of dimension six with complex multiplication, their integral models and their fields of definition.

## Papers

| Document | PDF | Source |
| --- | --- | --- |
| I. A fixed CM packet and its principal models | [Read](papers/pdf/packet-and-principal-models.pdf) | [LaTeX](papers/packet-and-principal-models.tex) |
| II. Quotients, polarizations and descent | [Read](papers/pdf/quotients-and-descent.pdf) | [LaTeX](papers/quotients-and-descent.tex) |
| H. Hodge classes | [Read](papers/pdf/hodge-classes.pdf) | [LaTeX](papers/hodge-classes.tex) |
| S. Arithmetic and computation | [Read](papers/pdf/arithmetic-and-computation.pdf) | [LaTeX](papers/arithmetic-and-computation.tex) |

Start with Paper I for the packet and principal models, or Paper II for the
quotient and descent problem. The [research overview](docs/research-overview.md)
explains their relationship. Precise hypotheses, computational premises, review
status and open questions are collected in [scope and limitations](docs/limitations.md).

## Navigate the results

- [Claim register](docs/claim-register.md) and [paper-order index](docs/claim-reading-order.md)
- [Claims and evidence](docs/claims.md)
- [Prior work and contribution](docs/prior-work.md)
- [Reproduction instructions](docs/reproducing-results.md)
- [Certificate sources](certificates/README.md), [separate implementations](independent/README.md)
  and [data catalog](data/README.md)
- [Evidence selection](evidence/README.md) and [contributing](CONTRIBUTING.md)

## Run the examples

From the repository root, run the integrity check with Python 3.9 or newer
on macOS or Linux:

```sh
python3 scripts/check-public-package.py
```

It checks included files, claim navigation and source locations.

To run a small exact mathematical example without installing dependencies:

```sh
python3 -B certificates/known_cm_comparison_v58.py
```

Expect JSON with `status: PASS_EXACT_CM_COMPARISON_REGRESSION`. The
[quickstart](docs/reproducing-results.md#first-run) explains a second example,
output retention and the limits of each check. For section previews, use the
[main-document map](papers/README.md#opening-and-editing-sources).

## Reuse

Original code and computational data are licensed under MIT; original papers
and reader documentation are licensed under CC BY 4.0. See [licensing](LICENSE.md)
and [third-party notices](THIRD-PARTY-NOTICES.md), including PDF font terms.

## Feedback welcome

Please use [GitHub Issues](https://github.com/waltchackdr/mcmullen-cm-packets/issues)
for mathematical corrections, prior-work references, reproducibility problems
or suggestions that improve the exposition. Include the paper and theorem label
or claim ID, the precise issue, and a source or minimal example when possible.
For computational reports, include the command, software versions and relevant
output, with credentials and personal information removed.

The [review brief](docs/specialist-review.md) identifies focused proof and novelty
questions. Corrections will be recorded with the affected claims and supporting
evidence, preserving earlier versions in the repository history.
