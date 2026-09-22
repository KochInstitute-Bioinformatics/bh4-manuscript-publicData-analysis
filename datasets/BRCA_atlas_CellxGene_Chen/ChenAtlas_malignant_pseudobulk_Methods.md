# Methods

## Data source and pseudobulk construction

Donor-level pseudobulk expression for the malignant-cell compartment of
this scRNA-seq atlas was prepared upstream by
`ChenAtlas_malignant_pseudobulk_prep.Rmd` and loaded here from
`Atlas_malignant_pseudobulk_donor.rds`; this notebook does not redo that
aggregation. Briefly: the prep notebook downloads the source h5ad,
resolves genes to symbols, subsets to malignant cells, averages each
malignant cell's log-normalized expression per donor
(`AverageExpression()`, since this atlas carries no raw counts to sum),
and attaches donor-level clinical metadata (batch, disease, tissue, sex,
ethnicity, development stage, grade) and malignant cell counts per donor
before saving the object this notebook loads. No score columns are
computed upstream -- all EMT/gene-set scoring (below) happens in this
notebook.

## EMT state scoring and tertile stratification

**KS method** (Chakraborty et al.; `emt_code/Chakraborty_et_al_KS_code_vCAW.R`),
computed in this notebook via `KSScore()` -- a per-sample, two-sample
Kolmogorov-Smirnov test comparing each donor's own expression across the
tumor-calibrated `Mes.tum` genes against its own expression across the
`Epi.tum` genes (`emt_code/EM_gene_signature_tumor_KS.xlsx`), producing
`EMT_KS_TumorSig`. MLR is not used anywhere in this notebook.

Donors are stratified into low/mid/high tertiles of `EMT_KS_TumorSig`
(`KS_Group`), computed immediately afterward in this notebook (quantile
breaks at 1/3 and 2/3 of the score).

## Gene sets

Same 9-gene-set list as the companion bulk-data notebooks
(`TCGA_BRCA_EMT_GeneSet_Scoring.Rmd`/`SUM149PT_EMT_scoring_plotting.Rmd`):
`BH4` (`DHFR`, `GCH1`, `PCBD1`, `PCBD2`, `PTS`, `QDPR`, `SPR`) and its
`BH4.synth` (`GCH1`, `PTS`, `SPR`) / `BH4.recyc` (`QDPR`, `DHFR`, `PCBD1`,
`PCBD2`) subsets, `Glutathione` (`SLC7A11`, `SLC3A2`, `GCLC`, `GCLM`,
`GSS`, `GSR`), `CoQ10` (13 genes: `COQ2`-`COQ9`, `HPDL`, `PDSS1`, `PDSS2`,
`STARD7`), `Thioredoxin` (`TXN`, `TXNRD1`, `TXN2`, `TXNRD2`),
`Peroxiredoxin` (`PRDX1`-`PRDX5`), and the `Epi.tum`/`Mes.tum` gene sets
from the tumor-calibrated KS EMT signature. These gene lists are built
directly in this notebook, since `AddModuleScore_UCell()` needs the
actual gene lists, not just resulting score columns.

## Gene-set activity scoring

Pathway/gene-set activity was scored three ways on `pseudobulk`'s
log-normalized `data` layer, all computed directly in this notebook:

- **GSVA** (Hänzelmann et al., 2013), via `gsvaParam(expr_matrix_pb,
  pathways_to_score)` and `gsva()` -- no `kcdf` override (unlike the bulk
  RNA-seq notebooks).
- **ssGSEA** (Barbie et al., 2009), via `ssgseaParam(..., normalize =
  TRUE)` and `gsva()`, same `expr_matrix_pb`/`pathways_to_score` as GSVA.
- **UCell** (Andreatta & Carmona, 2021), via
  `UCell::AddModuleScore_UCell()`, `maxRank = 10000` (matching the
  bulk-data notebooks, well above UCell's single-cell-tuned default of
  1,500). UCell's per-cell-ranking approach could be judged inapplicable
  at the single-cell level, but its rank-based score is still
  well-defined per-sample at the pseudobulk level -- the same reasoning
  that applies to GSVA/ssGSEA -- so it's scored here alongside them.

None of these three scoring steps write back to
`Atlas_malignant_pseudobulk_donor.rds`: all scoring in this notebook
(EMT KS scoring included) is in-memory only, for the current run's
figures.

## Statistical comparisons

For each gene set/scoring-method combination, all pairwise two-sided
Wilcoxon rank-sum tests were run among the 3 `KS_Group` levels
(low-vs-mid, mid-vs-high, low-vs-high), with p-values adjusted for
multiple comparisons across those 3 tests by the Benjamini-Hochberg (BH)
procedure. Every pair is shown on the boxplots (not just the significant
ones), labeled with its BH-adjusted p-value. No omnibus Kruskal-Wallis
test is reported, for consistency with the companion bulk-data notebooks
-- it doesn't gate the pairwise brackets, and the pairwise comparisons are
the more specific answer to which tertiles differ.

## Visualization

- **Boxplots**: one panel per gene set, box plot by `KS_Group`, y-axis
  limited to a tight whisker-based range (`boxplot.stats()`), with the
  BH-adjusted pairwise Wilcoxon brackets described above. Run once per
  scoring method (ssGSEA, GSVA, UCell) across all 9 scored gene sets.
- **Heatmap** (`ComplexHeatmap`): ssGSEA scores for 7 gene sets
  (`Epi.tum`, `Mes.tum`, `BH4`, `Glutathione`, `CoQ10`, `Thioredoxin`,
  `Peroxiredoxin` -- `BH4.synth`/`BH4.recyc` excluded as subsets of `BH4`'s
  own genes), donors ordered by `EMT_KS_TumorSig` (no clustering),
  z-scored per row, with `EMT_KS_TumorSig` as a continuous column color
  bar.

## Saved figures

Of all the boxplots generated, only the 5 ssGSEA panels below are saved to
`figures/` (as SVG); every other panel (GSVA, UCell, and the remaining
ssGSEA gene sets) is rendered inline only, and nothing else in this
notebook writes to disk:

| File | Gene set | Method |
| --- | --- | --- |
| `Fig5b_ssgsea_BH4_Chen.svg` | BH4 | ssGSEA |
| `ExtData5a_ssgsea_Glutathione_Chen.svg` | Glutathione | ssGSEA |
| `ExtData5b_ssgsea_CoQ10_Chen.svg` | CoQ10 | ssGSEA |
| `ExtData5c_ssgsea_Thioredoxin_Chen.svg` | Thioredoxin | ssGSEA |
| `ExtData5d_ssgsea_Peroxiredoxin_Chen.svg` | Peroxiredoxin | ssGSEA |

## Software

Analysis was performed in R (see
`ChenAtlas_pseudobulk_figures_sessionInfo.txt` for the full
session/package version listing), primarily using `tidyverse`, `GSVA`,
`UCell`/`Seurat`, `ComplexHeatmap`, and `ggpubr`.
