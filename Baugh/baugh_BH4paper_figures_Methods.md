# Methods

## Data source and processing pipeline

scRNA-seq data from Baugh et al. (8 samples: `NT2.5_V1`–`V4`, mouse
mammary tumor; `NT2.5_LM_V1`–`V4`, the lung-metastasis-derived condition)
were processed upstream by `baugh_BH4paper_prep.Rmd` and loaded here from
`singlet.integrated.tumor_chakrobortyMethod.rds`; this notebook does not
redo any of that processing. Briefly, in five stages:

1. **Per-sample import and doublet calling**: each sample imported
   independently (`Read10X`, `min.cells = 3`, `min.features = 200`),
   given a first-pass clustering, and doublet-annotated with
   `DoubletFinder` (expected doublet count from a Poisson per-cell
   duplication-rate model, `mu = 0.144` for standard 10x GEM chemistry,
   adjusted for the homotypic doublet proportion within each sample's own
   clustering).
2. **Merge, filter, integrate**: the 8 doublet-annotated samples merged,
   doublets dropped, filtered to `200 <= nFeature_RNA <= 8000` and
   `percent.mt <= 15`, log-normalized, 2000 HVGs, scaled (regressing out
   `nCount_RNA` and `percent.mt`), then CCA-integrated across all 8
   samples (`IntegrateLayers`, `CCAIntegration`) and clustered
   (resolution 0.4) on the integrated reduction.
3. **Cell-type annotation**: clusters manually assigned Tumor (clusters
   0, 1, 2, 4, 9, 10, identified via the source publication's tumor
   marker panel -- `Lcn2`, `Wfdc2`, `Cd24a`, `Cd276`, `Col9a1`, `Erbb2`)
   or a canonical non-tumor identity (Fibroblast, Monocyte, Neutrophil,
   Endothelial, T cell, Macrophage) from marker expression.
4. **Tumor subsetting and re-integration**: subset to the tumor clusters,
   re-split by sample, and re-integrated from scratch (a fresh CCA
   reduction, kept distinct from the whole-dataset integration), then
   re-clustered at three resolutions (0.1, 0.2, 0.4) in one
   `FindClusters()` call. Figure 1c below uses the resolution-0.4
   clustering.
5. **Single-cell EMT KS scoring** (see below), producing the
   `EMT_KS_Score` metadata column used in Figure 1a.

## EMT scoring method

**KS method** (Chakraborty/Bell et al.;
`emt_code/Chakraborty_et_al_KS_code_vCAW.R`), using the tumor-calibrated
`Epi.tum`/`Mes.tum` signature (`emt_code/EM_gene_signature_tumor_KS.xlsx`).
Because this dataset is mouse, the human signature genes are mapped to
1:1 mouse orthologs (`Mouse_Human_Orthologs_OnetoOne.xlsx`) before
scoring -- `KSScore()` reads its gene lists from a two-column xlsx file
path rather than as function arguments, so the mouse-mapped signature is
written to `emt_code/EM_gene_signature_tumor_mouse_KS.xlsx` for it to
read.

**Two different granularities of KS scoring feed the two scored figures**,
both from the same method and signature but not numerically the same
values:

- **Figure 1a** (feature plot) uses the **single-cell** KS score: a
  non-vectorized `KSScore()` call run per cell (tens of thousands of
  `ks.test()` triples) on the tumor object's log-normalized expression,
  computed once upstream in the prep notebook's Stage 5 and cached as
  `EMT_KS_Score` in `singlet.integrated.tumor_chakrobortyMethod.rds`
  (re-running it here would be slow; it also goes stale silently if the
  signature, ortholog mapping, or upstream object changes, since nothing
  in this notebook re-triggers that stage).
- **Figure 1b** (box plot) instead uses a **replicate-level pseudobulk**
  KS score, computed directly in this notebook: raw counts summed per
  `orig.ident` (8 pseudobulk replicates, one per sample), TMM-normalized
  (`edgeR::calcNormFactors()`) and converted to log-CPM
  (`edgeR::cpm(..., log = TRUE, prior.count = 1)`), then scored once per
  replicate with the same `KSScore()` and mouse-mapped tumor signature.
  This is recomputed fresh every run, so it isn't subject to the same
  staleness risk as Figure 1a's cached per-cell score.

## Gene sets (Figure 1c)

Three manually curated marker gene groups, per the manuscript's requested
panel:

- **EMT TF**: `Snai1`, `Zeb1`, `Vim`
- **Stemness**: `Cd44`, `Sox9`, `Cd24a`, `Itga6`, `Trp63`
- **Epithelial**: `Cdh1`, `Epcam`, `Krt7`, `Itgb4`, `Krt8`, `Krt14`

## Statistical comparisons

Only Figure 1b carries a significance test: a two-sided **Welch's
t-test** (`ggpubr::stat_compare_means(method = "t.test")` --
`t.test()`'s unequal-variance default), comparing the replicate-level
pseudobulk KS scores between the two `Model` conditions (`NT2.5` vs.
`NT2.5_LM`, n = 4 replicates each), drawn as a single bracket labeled
with the raw p-value (`label = "p.format"`). No multiple-comparison
adjustment is applied, since there is only one comparison to make (two
`Model` levels). Figures 1a and 1c are purely descriptive (a feature plot
and a mean-expression/percent-expressed dot plot); neither carries a
statistical test.

## Visualization

- **Figure 1a**: single-cell `EMT_KS_Score` UMAP feature plot
  (`reduction = "umap.cca.new"`), one panel per `Model` built separately
  and combined with `patchwork` (rather than Seurat's `split.by`, which
  would let each panel pick its own color limits) so both panels share
  one manually clamped RdBu diverging color scale (ColorBrewer, reversed
  so red = high/mesenchymal, blue = low/epithelial), limits fixed at
  `[-0.2, 0.2]` with out-of-range cells squished to the limit color
  (`scales::squish`) rather than left unplotted or off-scale.
- **Figure 1b**: box + jitter plot of the replicate-level pseudobulk KS
  score by `Model`, with the Welch's t-test bracket described above.
- **Figure 1c**: gene-level dot plot at the resolution-0.4 tumor
  clustering, one panel row per `Model`, dot color = per-gene z-scored
  mean expression (z-scored across cluster x `Model` groups), dot size =
  percent of cells expressing, genes grouped into the three marker
  categories above (dashed vertical separators between groups, a
  color-coded annotation bar below the main panel), clusters manually
  reordered top-to-bottom (`1, 2, 3, 8, 0, 7, 6, 9, 5, 4`) per reviewer
  feedback rather than left in numeric order.

## Saved figures

| File | Panel | Description |
| --- | --- | --- |
| `Fig1a_EMT_KS_Score_FeaturePlot_chakroborty_RdBu.svg` | Fig1a | Single-cell EMT KS score feature plot, by condition |
| `Fig1b_EMT_KS_Score_byCondition_boxplot_chakroborty.svg` | Fig1b | Pseudobulk-by-replicate EMT KS score box plot, by condition, with Welch's t-test |
| `Fig1c_EMT_Stem_Epi_Dotplot_byCluster_res0.4.svg` | Fig1c | EMT TF/Stemness/Epithelial marker dot plot, resolution-0.4 clusters x condition |

## Software

Analysis was performed in R (see `baugh_BH4paper_prep_sessionInfo.txt` and
`baugh_BH4paper_figures_sessionInfo.txt` for the full session/package
version listings for the prep and figures notebooks respectively),
primarily using `Seurat`, `DoubletFinder` (prep only), `edgeR`
(pseudobulk TMM normalization, figures only), `patchwork`, and `ggpubr`.
