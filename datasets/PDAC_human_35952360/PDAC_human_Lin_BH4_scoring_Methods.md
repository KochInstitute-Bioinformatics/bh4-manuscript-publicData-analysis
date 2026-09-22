# Methods

## Data source

`Lin2020_Pancreas_malignant.rds` is a pre-processed human PDAC scRNA-seq
Seurat object (Lin et al. 2020, PMID 35952360 -- the same publication
that provides the mouse dataset in the companion `PDAC_mouse_35952360/`
directory), subset to malignant cells. This notebook does not redo any
upstream processing.

## Cell-state classification (Raghavan 2021 rule)

Identical classification rule to `PDAC_human_38702773/` (Park et al.):
scBasal/scClassical/IC module scores via `AddModuleScore()`, then
`IC_Score = IC - |Classical - Basal|`, called `IC` if `IC_Score >= 0.2`,
else `Classical` if `Classical > Basal`, else `Basal`. Samples flagged
`is_met = "yes"` (metastatic, `sample` starting with `MET`) are excluded
from the Extended 6d boxplot -- primary tumor only.

## Known Gaps (resolved during the `datasets/` reorg)

The source script
(`tmp/WH_BH4pathway_final/src/Analysis_script.Rmd`, Data3 section) read
the Raghavan state signatures from `state_signatures.xlsx` -- a file that
does not exist anywhere in the exported analysis bundle. Data2
(`PDAC_human_park_BH4_scoring.Rmd`, Park et al.) runs this identical
classification code against `34890551_Raghavan_state_signatures.xlsx`,
and since both datasets use the same Raghavan 2021 rule with the same
`ic_threshold = 0.2` and the same 3 signature columns, this was most
likely a copy-paste filename shortening rather than a genuinely different
file. This notebook was updated to read
`34890551_Raghavan_state_signatures.xlsx` (copied into this directory)
instead, per the repo owner's confirmation -- if this assumption turns
out to be wrong, the original (missing) `state_signatures.xlsx` reference
is preserved in git history (see the `reorg` branch's parent commit).

## BH4 signature scoring

Human BH4 gene panel scored with `UCell` in 3 variants (see
`emt_code/bh4_pathway_shared_setup.R`); `maxRank` is the median
`nFeature_RNA` across the object.

## Statistics

Boxplots (`make_boxplot_state()`) run a Kruskal-Wallis test across cell
states, and where significant (p < 0.05), pairwise Wilcoxon tests with
Benjamini-Hochberg correction; only adjusted p < 0.05 comparisons are
bracket-annotated.

## Saved figures

| Figure | File(s) in `figures/` |
| --- | --- |
| ExtData6d | `Extended_6d_Boxplot_byCellState_BH4_score.pdf`/`.svg` |

## Software

See `PDAC_human_Lin_BH4_scoring_sessionInfo.txt` (placeholder until the
notebook is re-knit in the project's RStudio Server container -- see
README's "Running the analyses").
