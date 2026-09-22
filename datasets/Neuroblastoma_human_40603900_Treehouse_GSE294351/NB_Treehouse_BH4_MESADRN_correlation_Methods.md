# Methods

## Data source

`Treehouse_NB.rds` is a pre-processed bundle (`list(gene_matrix =
<genes x samples>, clinical_metadata = <data.frame>)`) of a Treehouse
neuroblastoma cohort (PMID 40603900, GEO accession GSE294351). This
notebook does not redo any upstream processing.

## Methodology

Identical to the companion SEQC notebook
(`Neuroblastoma_human_25150838_SEQC_GSE62564/`) -- same MES/ADRN
signature source (van Groningen et al. 2017, PMID 28650485), same
combined ssGSEA call across MES/ADRN/BH4 signatures (`ssgseaParam()` +
`gsva(param)` -- see that notebook's Known Gaps entry for a GSVA-API
bug found and fixed here too during in-container verification), same
`MES_ADRN_score = -(ADRN - MES)` derivation, same Spearman correlation +
BH-corrected statistics. See that notebook's Methods.md for the full
write-up; only the cohort differs.

## Saved figures

| Figure | File(s) |
| --- | --- |
| Fig5i | `figures/5i_scatterplot_BH4vsMESADRN_score.pdf`/`.svg` |
| (stats) | `BH4_signatures_spearman_correlation_stats_ssgsea.csv` |

## Software

See `NB_Treehouse_BH4_MESADRN_correlation_sessionInfo.txt` (placeholder
until the notebook is re-knit in the project's RStudio Server container
-- see README's "Running the analyses").
