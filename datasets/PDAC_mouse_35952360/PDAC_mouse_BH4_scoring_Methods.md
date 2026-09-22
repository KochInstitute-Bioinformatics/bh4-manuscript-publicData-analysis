# Methods

## Data source

`pitter_bh4_seurat.rds` is a pre-processed mouse PDAC scRNA-seq Seurat
object (Pitter et al., PMID 35952360), already annotated with a
`cell_state` column (`Luminal cells`, `Basal cells`, `Mesenchymal cells`,
`Other`) and a fitted UMAP embedding. This notebook does not redo any
upstream processing (QC, integration, clustering) -- it starts from that
checkpoint.

## BH4 signature scoring

Mouse orthologs of the human BH4 gene panel (`Gch1`, `Pts`, `Spr`,
`Qdpr`, `Dhfr`, `Pcbd1`, `Pcbd2`) are scored with `UCell` in 3 variants
(de novo synthesis, recycling/regeneration, and the full composite
pathway -- see `emt_code/bh4_pathway_shared_setup.R`). `maxRank` is set
per-dataset to the median `nFeature_RNA` of the malignant compartment
(`cell_state != "Other"`), not a fixed default, so the ranking depth
scales with this dataset's own sequencing depth.

## Statistics

Boxplots (`make_boxplot_state()`, `emt_code/bh4_pathway_shared_setup.R`)
run a Kruskal-Wallis test across cell states, and where significant
(p < 0.05), pairwise Wilcoxon tests with Benjamini-Hochberg correction;
only comparisons with adjusted p < 0.05 are bracket-annotated.

## Known Gaps (resolved -- author-confirmed fix)

Extended 6a's `p1` was undefined when this notebook was first split out
of the source script (`tmp/WH_BH4pathway_final/src/Analysis_script.Rmd`,
Data1 section) -- `plot_grid(p1, p3, rel_widths = c(1.5, 6))` referenced
`p1`, but it was never defined anywhere in this dataset's code path.
Author (Evelyn) confirmed by email this was an accidental deletion and
supplied the original definition -- a narrow UMAP colored by
`cell_state` (`scale_color_manual(values = state_color, ...)`, Other
cells drawn first/behind via `order(cell_state == "Other", decreasing =
TRUE)`), titled "Autochthonous KPCT cells". Restored verbatim; see the
notebook's Extended 6a chunk.

## Saved figures

| Figure | File(s) in `figures/` |
| --- | --- |
| Fig5d | `5d_UMAP_BH4_score.pdf`/`.svg` |
| ExtData6a | `Extended_6a_UMAP_byCellState_GeneExpression.pdf`/`.svg` |
| Fig5e | `5e_Boxplot_byCellState_BH4_score.pdf`/`.svg` |
| ExtData6c | `Extended_6c_Boxplot_byCellState_GeneExpression.pdf`/`.svg` |

## Software

See `PDAC_mouse_BH4_scoring_sessionInfo.txt` (placeholder until the
notebook is re-knit in the project's RStudio Server container -- see
README's "Running the analyses").
