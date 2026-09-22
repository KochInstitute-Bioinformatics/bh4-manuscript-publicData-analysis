# Methods

## Data source

`TCGA_bulkRNAseq_skcm.rds` is a pre-processed TCGA SKCM (melanoma) bulk
RNA-seq object, one sample per patient. This notebook does not redo any
upstream processing -- it loads the object, normalizes, and scores.

## Subtype classification (Tsoi 2018 4-state trajectory)

Each patient is assigned one of 4 transcriptional subtypes
(`Melanocytic`, `Transitory`, `Neural crest-like`, `Undifferentiated`)
from Tsoi et al. 2018 (PMID 29657129), read from
`SKCM_Tsoi2018_suppTable_s4.xlsx`. Patient IDs are derived from TCGA
barcodes (first 3 hyphen-delimited fields); when a patient has multiple
rows in the supplementary table, the first is kept (`slice_head(n = 1)`).

## BH4 signature scoring

BH4 signatures are scored with ssGSEA (`GSVA::ssgseaParam(..., normalize
= TRUE)` + `gsva(param)` -- the modern param-object API; see Known Gaps)
on the `data` (log-normalized) assay layer -- ssGSEA is rank-based, so
the choice of layer doesn't affect the result as long as it preserves
gene-expression rank order within each sample.

## Known Gaps (fixed during in-container verification)

The source script's `gsva(expr, sigs, method = "ssgsea", ssgsea.norm =
TRUE, verbose = FALSE)` call used GSVA's old function-argument API,
removed in GSVA >= 2.x (the container has GSVA 2.4.8). This errored with
`unable to find an inherited method for function 'gsva' for signature
'param = "matrix"'`. The source script's own comment ("Legacy GSVA API
(matches your other scripts)") was incorrect -- `BRCA_TCGA` and
`BRCA_atlas_CellxGene_Chen` already use the modern
`ssgseaParam()`/`gsvaParam()` + `gsva(param)` pattern, which this
notebook now matches. The same bug affected
`Neuroblastoma_human_25150838_SEQC_GSE62564` and
`Neuroblastoma_human_40603900_Treehouse_GSE294351`'s ssGSEA calls; both
were fixed the same way.

## Statistics

The boxplot (`make_boxplot_state()`) runs a Kruskal-Wallis test across
the 4 Tsoi subtypes (abbreviated on the x-axis: Melano/Transit/NCL/
Undiff), and where significant (p < 0.05), pairwise Wilcoxon tests with
Benjamini-Hochberg correction; only adjusted p < 0.05 comparisons are
bracket-annotated.

## Saved figures

| Figure | File(s) in `figures/` |
| --- | --- |
| Fig5g | `5g_Boxplot_byCellState_BH4_score.pdf`/`.svg` |

## Software

See `SKCM_TCGA_BH4_scoring_sessionInfo.txt` (placeholder until the
notebook is re-knit in the project's RStudio Server container -- see
README's "Running the analyses").
