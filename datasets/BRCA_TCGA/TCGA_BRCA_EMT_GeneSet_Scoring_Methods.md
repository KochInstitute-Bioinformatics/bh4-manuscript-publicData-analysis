# Methods

## Data source and BRCA subset construction

Pan-cancer expression data were prepared upstream by `TCGA_BRCA_prep.Rmd`
and loaded here from `TCGA_BRCA_seurat.rds`; this notebook does not redo
that assembly. Briefly: the prep notebook loads the Broad Firehose
`EBPlusPlusAdjustPANCAN` batch-corrected RSEM matrix
(`geneexpressionmatrix_uniqueIDs.txt`; continuous, upper-quartile-normalized
and empirical-Bayes batch-corrected across cancer types -- not raw counts),
replaces non-finite values with 0 without rounding (the EMT/gene-set scoring
methods below all expect continuous expression), builds one Seurat object
with one "cell" per tumor sample, attaches the matched clinical/metadata
table (`Updated_Metadata_fixed.txt`, joined on TCGA barcode), and subsets
the ~11,069-sample pan-cancer object down to the ~1,215 samples with
`type == "BRCA"`, saved as `TCGA_BRCA_seurat.rds`.

The GDC/Broad Firehose pages that originally hosted these two raw input
files are no longer live; provenance is instead tracked via the associated
publication (Liu, J., Lichtenberg, T., Hoadley, K.A. et al. An Integrated
TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival
Outcome Analytics. *Cell* 173(2), 400–416.e11 (2018).
<https://doi.org/10.1016/j.cell.2018.02.052>), documented in the prep
notebook's "Raw data" section.

No additional per-sample normalization (e.g. Seurat's `NormalizeData()`) is
layered on top of the Firehose matrix -- that would risk reintroducing the
composition bias the upstream UQ normalization was designed to remove. The
prep notebook checks this decision by confirming `nCount_RNA` is
reasonably stable across cancer types before subsetting to BRCA.

This directory only carries the code needed to reproduce the BRCA-specific
figures. An earlier pan-cancer companion notebook
(`TCGA_EMT_GeneSet_Scoring.Rmd`, referenced in this notebook's own
"Purpose" section as prior context for why BRCA-only rescoring isn't a
no-op for every method) originally scored all ~11,069 samples but is not
included in this repository; every score used in the figures below is
recomputed directly on the BRCA-only cohort in this notebook, from
`TCGA_BRCA_seurat.rds` alone.

## Expression normalization for scoring

`CreateSeuratObject()` only populates the `counts` layer. Rather than
route through `NormalizeData()` (whose total-count normalization is the
extra step decided against above), the `data` layer used for all scoring
below is built directly as `log2(counts + 1)`.

## Gene sets

Same 9-gene-set list as the companion notebooks
(`ChenAtlas_malignant_pseudobulk.Rmd`/`SUM149PT_EMT_scoring_plotting.Rmd`):
`BH4` (`DHFR`, `GCH1`, `PCBD1`, `PCBD2`, `PTS`, `QDPR`, `SPR`) and its
`BH4.synth` (`GCH1`, `PTS`, `SPR`) / `BH4.recyc` (`QDPR`, `DHFR`, `PCBD1`,
`PCBD2`) subsets, `Glutathione` (`SLC7A11`, `SLC3A2`, `GCLC`, `GCLM`,
`GSS`, `GSR`), `CoQ10` (13 genes: `COQ2`-`COQ9`, `HPDL`, `PDSS1`, `PDSS2`,
`STARD7`), `Thioredoxin` (`TXN`, `TXNRD1`, `TXN2`, `TXNRD2`),
`Peroxiredoxin` (`PRDX1`-`PRDX5`), and the `Epi.tum`/`Mes.tum` gene sets
from the tumor-calibrated KS EMT signature
(`emt_code/EM_gene_signature_tumor_KS.xlsx`). A separate `goi` list
(`GPX4`, `AIFM2`, `DHODH`, `PRDX6`) of individual genes is scored as raw
expression only (not as a gene set) for the full gene-set-score heatmap
below.

## EMT state scoring and tertile stratification

**KS method** (Chakraborty et al.; `emt_code/Chakraborty_et_al_KS_code_vCAW.R`),
computed via `KSScore()` -- a per-sample, two-sample Kolmogorov-Smirnov
test comparing each sample's own expression across the tumor-calibrated
`Mes.tum` genes against its own expression across the `Epi.tum` genes,
producing `EMT_KS_TumorSig`. Because this is a per-sample test that never
references other samples, BRCA-only and pan-cancer-context KS scores for
the same sample are expected to be numerically identical -- rescoring here
is for consistency/self-containment, not because the method requires it.

Samples are stratified into low/mid/high tertiles of `EMT_KS_TumorSig`
(`KS_Group`), computed immediately afterward on the BRCA-only distribution
(quantile breaks at 1/3 and 2/3 of the score) -- i.e. tertile boundaries
are BRCA-calibrated, not inherited from any pan-cancer scoring.

## Gene-set activity scoring

Pathway/gene-set activity was scored three ways on the BRCA-only `data`
layer described above, all computed directly in this notebook against the
9-gene-set list:

- **GSVA** (Hänzelmann et al., 2013), via `gsvaParam(..., kcdf = "none")`
  and `gsva()`.
- **ssGSEA** (Barbie et al., 2009), via `ssgseaParam(..., normalize =
  TRUE)` and `gsva()`, same input matrix and gene sets as GSVA.
- **UCell** (Andreatta & Carmona, 2021), via
  `UCell::AddModuleScore_UCell()`, `maxRank = 10000`.

Unlike the KS method and UCell (both rank/compare within a single sample
and are therefore unaffected by cohort composition), GSVA and ssGSEA
rank-transform each gene's expression *across all samples in the input
matrix* before scoring -- so their BRCA-only scores are expected to differ
from any pan-cancer-context scores for the same sample, reflecting each
sample's position relative to the BRCA cohort specifically rather than the
full pan-cancer background.

The re-scored object (`EMT_KS_TumorSig`, `KS_Group`, and all
GSVA/ssGSEA/UCell columns) is saved as `TCGA_BRCA_seurat_EMTscored.rds`.

## Statistical comparisons

- **Pearson correlation** (`ggpubr::stat_cor()`, default method), reported
  on three families of scatter plots: GSVA score vs. ssGSEA score for each
  of the 9 gene sets (one panel per set); BH4 GSVA score vs.
  `EMT_KS_TumorSig`; and raw log2-normalized expression of each of BH4's 7
  member genes vs. `EMT_KS_TumorSig` (one panel per gene).
- **All pairwise two-sided Wilcoxon rank-sum tests** among the 3
  `KS_Group` tertile levels (low-vs-mid, mid-vs-high, low-vs-high), with
  p-values adjusted for multiple comparisons across those 3 tests by the
  Benjamini-Hochberg (BH) procedure. Applied per gene-set/scoring-method
  combination (BH4 and, more generally, `Glutathione`/`CoQ10`/
  `Thioredoxin`/`Peroxiredoxin`, each scored by GSVA, ssGSEA, and UCell).
  Every pair is shown on the boxplots (not just the significant ones),
  labeled with its BH-adjusted p-value. No omnibus Kruskal-Wallis test is
  reported: with ~1,215 BRCA samples split into 3 tertiles by construction,
  it returns significant on essentially every panel regardless of effect
  size and doesn't gate the pairwise brackets, so it would add a second,
  less specific p-value rather than new information.

## Visualization

- **Boxplots**: BH4 (+ `BH4.synth`/`BH4.recyc`) scores by `KS_Group`
  (GSVA), shown together in one multi-panel figure without significance
  testing; and, separately, a single-panel-per-gene-set/method design
  (BH4, `Glutathione`, `CoQ10`, `Thioredoxin`, `Peroxiredoxin`, each for
  GSVA/ssGSEA/UCell) with a tight whisker-based y-axis (`boxplot.stats()`)
  and the BH-adjusted pairwise Wilcoxon brackets described above.
- **Scatter plots**: GSVA-vs-ssGSEA agreement per gene set; BH4 GSVA score
  vs. EMT KS score; each BH4 gene's expression vs. EMT KS score -- all with
  an OLS trend line and the Pearson correlation annotation described above.
- **Dot plot**: mean expression / percent-expressed summary of BH4's 7
  genes by `KS_Group`, colored by min-max normalized mean expression per
  gene (chosen over a per-gene z-score, which is unstable with only 3
  groups).
- **Heatmaps** (`ComplexHeatmap`): two variants, both with samples ordered
  by `EMT_KS_TumorSig` (no clustering) and rows z-scored --
  - Full heatmap: all 9 gene-set GSVA scores plus `goi`'s 4 genes' raw
    expression, with a row-side strip distinguishing gene-set-score rows
    from `goi`'s raw-expression rows.
  - Simplified heatmap: 7 gene-set scores only (`Epi.tum`/`Mes.tum` plus
    the 5 non-BH4-subset pathway sets; `BH4.synth`/`BH4.recyc` and `goi`
    dropped), run once per scoring method (GSVA, ssGSEA, UCell).

## Saved figures

Of all the panels generated, only the ssGSEA BH4-by-`KS_Group` single-panel
boxplot is saved to `figures/` (as SVG); every other panel (the combined
BH4-subsets boxplot, GSVA/UCell versions, scatter plots, dot plot, and both
heatmap variants) is rendered inline only:

| File | Gene set | Method |
| --- | --- | --- |
| `Fig5c_ssgsea_BH4_TCGA_BRCA.svg` | BH4 | ssGSEA |

## Software

Analysis was performed in R (see `TCGA_BRCA_prep_sessionInfo.txt` and
`TCGA_BRCA_EMT_GeneSet_Scoring_sessionInfo.txt` for the full
session/package version listings for the prep and scoring/plotting
notebooks respectively), primarily using `tidyverse`, `Seurat`, `GSVA`,
`UCell`, `ComplexHeatmap`, and `ggpubr`. Scoring methodology mirrors the
companion `ChenAtlas_malignant_pseudobulk.Rmd`/
`SUM149PT_EMT_scoring_plotting.Rmd` analyses, with grouping here based on
a data-driven `EMT_KS_TumorSig` tertile split (`KS_Group`) rather than
ground-truth experimental labels, since this cohort (unlike SUM149PT) has
no annotated EMT state.
