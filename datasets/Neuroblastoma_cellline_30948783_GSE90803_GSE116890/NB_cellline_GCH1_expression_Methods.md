# Methods

## Data source

Affymetrix HG-U133 Plus 2.0 (GPL570) MAS5.0 log2 microarray data (PMID
30948783), two GEO series:

- `GSE90803.rds`: paired ADRN/MES neuroblastoma cell lines (SH-SY5Y /
  SH-EP2, plus 3 patient-derived pairs: 691, 700, 717).
- `GSE116890.rds`: a NOTCH3-IC induction time-course (Day 0/1/7/14/21,
  ON = induced/MES-like vs. OFF = uninduced/ADRN-like).
- `probe_info.rds`: probeset-to-gene mapping (BH4 panel + control genes)
  shared by both series.

This notebook does not redo any upstream processing (normalization,
probe summarization) -- it reads the 3 pre-built `.rds` objects directly
and extracts per-gene expression via `extract_affy_genes()`
(`emt_code/bh4_pathway_shared_setup.R`), which looks up each gene's
probeset in the raw MAS5.0 matrix.

## Plotting

- **Extended 4a/4b**: `GCH1` expression by cell-line pair (GSE90803),
  split into a "3 pairs" panel (691/700/717, ADRN vs. MES) and a
  standalone SH-SY5Y/SH-EP2 panel, combined via `patchwork` with widths
  matched to their explicit x-axis position gaps so bar widths are
  visually identical across both panels.
- **Extended 4c**: `GCH1` expression over the NOTCH3-IC induction
  time-course (GSE116890), OFF vs. ON, with a shared Day-0 baseline point.

No hypothesis tests are run in this notebook -- both panels are
descriptive (raw/summarized expression values only).

## Saved figures

| Figure | File(s) in `figures/` |
| --- | --- |
| ExtData4a/4b | `Extended_4ab_barplot_GCH1_GeneExpression.pdf`/`.svg` |
| ExtData4c | `Extended_4c_lineplot_GCH1_GeneExpression.pdf`/`.svg` |

## Software

See `NB_cellline_GCH1_expression_sessionInfo.txt` (placeholder until the
notebook is re-knit in the project's RStudio Server container -- see
README's "Running the analyses").
