# Methods

## Data source

`park_bh4_seurat.rds` is a pre-processed human PDAC scRNA-seq Seurat
object (Park et al., PMID 38702773), already subset to malignant cells
via CopyKAT (so the whole object is the malignant compartment -- no
further cell-state filtering needed before computing the UCell
`maxRank`). This notebook does not redo any upstream processing.

## Cell-state classification (Raghavan 2021 rule)

Malignant cells are classified into `Classical`/`IC`/`Basal` states using
the `AddModuleScore()`-based rule from Raghavan et al. 2021 (PMID
34890551): scBasal, scClassical, and IC module scores are computed from
`34890551_Raghavan_state_signatures.xlsx`, then
`IC_Score = IC - |Classical - Basal|`; a sample is called `IC` if
`IC_Score >= 0.2`, else `Classical` if `Classical > Basal`, else `Basal`.

## BH4 signature scoring

The human BH4 gene panel (`GCH1`, `PTS`, `SPR`, `QDPR`, `DHFR`, `PCBD1`,
`PCBD2`) is scored with `UCell` in 3 variants (see
`emt_code/bh4_pathway_shared_setup.R`). `maxRank` is the median
`nFeature_RNA` across this object (already-malignant-only, so no
additional state filter is applied before taking the median).

## Statistics

Fig5f only plots primary-tumor samples (`Origin %in% c("Pm0", "Pm1")`).
Boxplots (`make_boxplot_state()`) run a Kruskal-Wallis test across cell
states, and where significant (p < 0.05), pairwise Wilcoxon tests with
Benjamini-Hochberg correction; only adjusted p < 0.05 comparisons are
bracket-annotated.

## Saved figures

| Figure | File(s) in `figures/` |
| --- | --- |
| Fig5f | `5f_Boxplot_byCellState_BH4_score.pdf`/`.svg` |

## Software

See `PDAC_human_park_BH4_scoring_sessionInfo.txt` (placeholder until the
notebook is re-knit in the project's RStudio Server container -- see
README's "Running the analyses").
