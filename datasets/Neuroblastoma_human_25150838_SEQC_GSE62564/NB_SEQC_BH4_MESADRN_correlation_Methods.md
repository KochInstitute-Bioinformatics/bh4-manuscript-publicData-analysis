# Methods

## Data source

`GSE62564.rds` is a pre-processed bundle (`list(gene_matrix =
<genes x samples>, clinical_metadata = <data.frame>)`) of the SEQC
neuroblastoma cohort (Zhang et al. 2015, PMID 25150838, GEO accession
GSE62564). This notebook does not redo any upstream processing.

## MES/ADRN and BH4 signature scoring

MES and ADRN signature gene lists come from van Groningen et al. 2017
(PMID 28650485, `28650485_Groningen_SuppTable2.xlsx`). MES, ADRN, and
all 3 BH4 signatures are scored together in one ssGSEA call
(`GSVA::ssgseaParam(..., normalize = TRUE)` + `gsva(param)` -- the
modern param-object API; see Known Gaps) on the gene matrix, so all
scores share the same sample background/normalization.
`MES_ADRN_score = -(ADRN - MES)` (i.e. `MES - ADRN`); this is the x-axis
of the correlation and scatter plot below.

## Known Gaps (fixed during in-container verification)

The source script's `gsva(matrix, sets, method = "ssgsea", ssgsea.norm =
TRUE, verbose = FALSE)` call used GSVA's old function-argument API,
removed in GSVA >= 2.x (the container has GSVA 2.4.8) -- errors with
`unable to find an inherited method for function 'gsva' for signature
'param = "matrix"'`. Fixed to use `ssgseaParam()` + `gsva(param)`,
matching `BRCA_TCGA`/`BRCA_atlas_CellxGene_Chen`. The same bug was found
and fixed the same way in `SKCM_TCGA` and the companion
`Neuroblastoma_human_40603900_Treehouse_GSE294351` notebook.

## Statistics

Spearman correlation (`cor.test(..., method = "spearman")`) between each
of the 3 BH4 signature scores and `MES_ADRN_score`, across all samples;
p-values are Benjamini-Hochberg corrected across the 3 tests
(`p_value_BH < 0.05` flags `significant`). Full results are written to
`BH4_signatures_spearman_correlation_stats_ssgsea.csv`.

## Saved figures

| Figure | File(s) |
| --- | --- |
| Fig5h | `figures/5h_scatterplot_BH4vsMESADRN_score.pdf`/`.svg` |
| (stats) | `BH4_signatures_spearman_correlation_stats_ssgsea.csv` |

## Software

See `NB_SEQC_BH4_MESADRN_correlation_sessionInfo.txt` (placeholder until
the notebook is re-knit in the project's RStudio Server container -- see
README's "Running the analyses").
