# P030 — Reproducing the numerical tables and figures

This page concerns the **four printed numerical tables** and **two figures** of the
JMP manuscript. The exact-moment single-point continuum certificate uses the
separate `exact_m2_certificate/` programs; it has different assumptions from
older rounded-SVD numerical experiments.

## From a clean Linux checkout

Use Python 3.13. Install the figure dependencies in addition to the existing
scientific requirements:

```bash
python -m pip install -r requirements-figures.txt
python tools/verify_release.py
python tools/verify_paper_reproduction.py
```

The second test independently looks up published values in raw JSON, checks
historical bath-grid convergence and unresolved-point handling, generates four
LaTeX table-row fragments, generates both PDF figures, and deliberately corrupts
a source dataset to ensure that **a mismatched result is rejected**.

Manual regeneration, from any current working directory:

```bash
python REPRODUCIBILITY/code/make_tables.py --check
python REPRODUCIBILITY/code/make_figures.py
python REPRODUCIBILITY/code/make_figure_nonnormal.py
```

Outputs are written into `REPRODUCIBILITY/generated/tables/` and
`REPRODUCIBILITY/generated/figures/`; these build artifacts are untracked.
All generated table values are calculated from recorded data, then compared
cell-for-cell with `reference/published_table_rows.json` (the four numerical
tables in the submitted manuscript).

If the author-supplied LaTeX `numerics_section.tex` is available, reproduce its
**scientific numerical entries exactly**, without replacing any prose, captions,
or authorship/policy disclosure:

```bash
python REPRODUCIBILITY/code/make_tables.py --check \
  --manuscript-section /path/to/paper/numerics_section.tex \
  --rebuild-section /tmp/numerics_section_rebuilt.tex
sha256sum /path/to/paper/numerics_section.tex /tmp/numerics_section_rebuilt.tex
```

The paper's editorial text is not synthesized from numerical data; the script
changes **only** numeric rows inside the four recognized table labels.
With the exact manuscript-source revision used for the 23 November submission,
the two SHA-256 values are identical:
`aef659f5de5ed9ac5a7f40c0140dd277e07d621dfb1b9e891fc482b7d75c888d`.

## One-to-one provenance

| Manuscript label/output | Source data | Generation |
| --- | --- | --- |
| `tab:ep` | `data/testC_EP.json` at `c=2.2` for five couplings | `code/make_tables.py` |
| `tab:thresholdsummary` | `data/testA_m1.json`, `data/testB_poly_m2.json` | `code/make_tables.py` |
| `tab:m2hp` | `order2_high_precision/direct_m2_results.json`, rounded-SVD bath | `code/make_tables.py` |
| `tab:nonnormal` | `data/testN2_micro.json`, scalar-block initial state | `code/make_tables.py` |
| `fig_thresholds.pdf` | `data/testA_m1.json`, `data/testB_poly_m2.json` | `code/make_figures.py` |
| `fig_nonnormal.pdf` | `data/testN2_micro.json` | `code/make_figure_nonnormal.py` |

The plots reproduce the **same underlying curves, markers, threshold positions
and unresolved-point filtering**. Embedded PDF font metrics and vector paths
can differ across Matplotlib versions. No byte-identical PDF claim is made.

## Boundaries

1. The generated numerical-table rows are checked against an independently
   captured source reference; generating identical rows merely by copying a
   manuscript file is prohibited by the regression test.
2. The reference scientific tables do not infer a continuum bound from
   historical bath-grid agreement. Only the *separate* exactly orthogonal
   parameter point has the conditional certified continuum enclosure.
3. This workflow generates the numerical scientific content, not an entire
   LaTeX article or supplementary manuscript; human-authored exposition and
   cross-references are preserved when rebuilding the author's section.
