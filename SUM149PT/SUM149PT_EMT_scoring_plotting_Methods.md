# Methods

## Data source and sample selection

This dataset is from:

Brown, M.S., Abdollahi, B., Wilkins, O.M. et al. Phenotypic heterogeneity
driven by plasticity of the intermediate EMT state governs disease
progression and metastasis in breast cancer. *Sci. Adv.* 8, eabj8002 (2022).
<https://doi.org/10.1126/sciadv.abj8002>

RNA-seq expression data for the SUM149PT EMT-spectrum cell line panel were
obtained from GEO accession
[GSE172609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172609)
(`GSE172609_norm_counts_TPM_GRCh38.p13_NCBI.tsv`
for TPM-normalized counts, `GSE172609_raw_counts_GRCh38.p13_NCBI.tsv` for raw
counts), aligned to GRCh38.p13, together with the corresponding NCBI gene
annotation table (`Human.GRCh38.p13.annot.tsv`). GSE172609 is the RNA-seq
SubSeries of the study's SuperSeries,
[GSE172613](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172613).
Sample-to-condition mapping (GEO accession to experimental `EMT_State`
label) was read from `Sample_to_Condition.xlsx`.

The dataset comprises 21 samples spanning 7 experimentally induced EMT states
(`Parental`, `E`, `EM1`, `EM2`, `EM3`, `M1`, `M2`; 3 replicates each). The 3
`Parental` ("P") replicates were excluded prior to all downstream scoring,
differential expression, and plotting, leaving 18 samples across the 6
remaining states (`E` -> `EM1` -> `EM2` -> `EM3` -> `M1` -> `M2`), reflecting
the induced epithelial-to-mesenchymal spectrum.

## Expression matrix preparation

Gene IDs were annotated with gene symbol and gene type (`GeneID`, `Symbol`,
`GeneType`) via a left join against the NCBI annotation table. For genes
represented by multiple GeneIDs mapping to the same gene symbol, a single
representative "max-per-gene" (MPG) GeneID was retained per symbol -- the
GeneID with the highest average TPM across all samples (ties broken by table
order). Both the TPM and raw-count matrices were then filtered to
protein-coding genes only (`GeneType == "protein-coding"`) with one row per
gene symbol (`MPG == "Yes"`), using the same MPG selection for both.

The TPM matrix was log2-transformed (`log2(TPM + 1)`) to put values on a
scale appropriate for the KS/GSVA/ssGSEA scoring methods below. The raw-count
matrix was kept separately, unnormalized, as input to DESeq2 differential
expression and to UCell scoring (below); DESeq2 performs its own
size-factor normalization internally, and UCell ranks genes within each
sample regardless of the input scale.

## Gene sets

The following gene sets were used for pathway/gene-set activity scoring:

- **BH4** (7 genes: `DHFR`, `GCH1`, `PCBD1`, `PCBD2`, `PTS`, `QDPR`, `SPR`) --
  tetrahydrobiopterin synthesis and recycling
  - **BH4.synth** (`GCH1`, `PTS`, `SPR`) -- de novo synthesis subset
  - **BH4.recyc** (`QDPR`, `DHFR`, `PCBD1`, `PCBD2`) -- recycling/regeneration subset
- **Glutathione** (`SLC7A11`, `SLC3A2`, `GCLC`, `GCLM`, `GSS`, `GSR`) --
  cystine uptake and glutathione synthesis/recycling
- **CoQ10** (13 genes: `COQ2`-`COQ9`, `HPDL`, `PDSS1`, `PDSS2`, `STARD7`) --
  CoQ/ubiquinol biosynthesis
- **Thioredoxin** (`TXN`, `TXNRD1`, `TXN2`, `TXNRD2`)
- **Peroxiredoxin** (`PRDX1`-`PRDX5`, canonical peroxiredoxins)
- **Epi.tum** / **Mes.tum** -- the epithelial and mesenchymal gene sets from
  the tumor-calibrated KS EMT signature (`emt_code/EM_gene_signature_tumor_KS.xlsx`),
  used both as scoring inputs to the KS method (below) and as standalone
  gene sets for GSVA/ssGSEA/UCell scoring. The tumor-calibrated signature was
  used (rather than the cell-line-calibrated alternative also available in
  `emt_code/`) to match the signature used for TCGA BRCA tumor samples in
  the companion `TCGA_BRCA_EMT_GeneSet_Scoring.Rmd` analysis.

## Differential expression (DESeq2 + apeglm)

DESeq2 differential expression was run on the raw, protein-coding,
`Parental`-excluded count matrix (design `~ EMT_State`, `E` set as the
reference level), comparing each of the other 5 states (`EM1`, `EM2`, `EM3`,
`M1`, `M2`) to `E`. Log2 fold changes were shrunk with `apeglm`
(`lfcShrink(..., type = "apeglm")`) for each of the 5 `state vs. E`
coefficients. The resulting per-gene, per-state log2 fold changes and
BH-adjusted p-values were assembled into a `BH4 x EMT_State` matrix and
visualized as a heatmap, with cells flagged with an asterisk where a gene's
fold change is both "large" (`|log2FC| >= 1`) and statistically significant
(BH-adjusted p <= 0.05). This heatmap is saved as
`figures/ExtData5f_heatmap_bh4genes_lfc_sig.svg`.

## EMT state scoring

**KS method** (Chakraborty et al.; `emt_code/Chakraborty_et_al_KS_code_vCAW.R`).
For each sample, a two-sample Kolmogorov-Smirnov test compares that sample's
own expression values across the Mes.tum genes against its own expression
values across the Epi.tum genes. The signed KS statistic (positive =
mesenchymal-shifted, negative = epithelial-shifted) is taken as the sample's
EMT score, with the sign/magnitude resolved from one-sided KS tests
(`greater` in each direction) at p < 0.05; if neither one-sided test reaches
significance, the larger of the two one-sided statistics is used, signed
accordingly. This is a per-sample test that does not reference other samples
in the cohort, so scores are directly comparable across studies/cohorts.

## Gene-set activity scoring

Pathway/gene-set activity was scored three ways against the 9-gene-set list
above (`BH4`, `BH4.synth`, `BH4.recyc`, `Glutathione`, `CoQ10`, `Thioredoxin`,
`Peroxiredoxin`, `Epi.tum`, `Mes.tum`):

- **GSVA** (Hänzelmann et al., 2013), on the log2-transformed expression
  matrix, via `gsvaParam(..., kcdf = "none")` (Gaussian-kernel density
  estimation disabled, appropriate for this already-continuous log-TPM data)
  and `gsva()`.
- **ssGSEA** (Barbie et al., 2009), on the same log2-transformed matrix, via
  `ssgseaParam(..., normalize = TRUE)` and `gsva()`.
- **UCell** (Andreatta & Carmona, 2021), via `UCell::AddModuleScore_UCell()`,
  run on a `Seurat` object built from the raw, protein-coding,
  `Parental`-excluded count matrix (the same one used for DESeq2 above) --
  UCell ranks genes per sample, so its scores are unaffected by whether the
  input is raw or log-transformed. `maxRank` was set to 10,000, well above
  UCell's single-cell-tuned default of 1,500, since bulk RNA-seq samples are
  far denser than single cells (most protein-coding genes are detected at
  some level) and a low `maxRank` would truncate the per-sample ranking
  before genes further down the list -- some of which may still be pathway
  members -- are ever ranked.

GSVA and ssGSEA rank-transform each gene's expression across all samples in
the input matrix before scoring, so those two methods' scores reflect each
sample's position relative to the other 17 samples in this cohort. UCell, by
contrast, ranks genes within each sample only (like the KS method above), so
its scores do not depend on which other samples are included in the cohort.

## Statistical comparisons

Two complementary approaches were used to test whether the BH4 gene-set
score (GSVA, ssGSEA, or UCell) differs between the baseline `E` state and
each of the other 5 `EMT_State` levels (`EM1`, `EM2`, `EM3`, `M1`, `M2`):

- **Pairwise Welch's t-tests**, BH-adjusted across the 5 `state vs. E`
  comparisons. Welch's t-test (rather than the nonparametric Wilcoxon
  rank-sum test) was used because, at n=3 replicates per group, the exact
  Wilcoxon test's null distribution is too coarse to reach significance at
  all: with only `choose(6,3) = 20` possible ways to split 6 values into two
  groups of 3, the most extreme possible rank separation yields a two-sided
  p-value of `2/20 = 0.1`, making p < 0.05 mathematically unreachable
  regardless of effect size. A parametric test on the actual score
  magnitudes does not have this discreteness floor, at the cost of assuming
  approximate normality of the per-group scores (a reasonable assumption for
  continuous GSVA/ssGSEA/UCell scores, though a stronger one than the
  Wilcoxon test's rank-based approach).
- **One-way ANOVA + Dunnett's test**, an alternative to the pairwise
  Welch's t-tests above: an ordinary one-way ANOVA (`aov()`) followed by
  Dunnett's multiple-comparisons test (`multcomp::glht()`,
  `mcp(Group = "Dunnett")`), with `E` set as the reference level. Dunnett's
  test is built for exactly this "every group vs. one reference" design --
  it accounts for the correlation among comparisons that share a common
  reference group, so its adjustment is less conservative here than BH
  across independent pairwise tests. It still assumes normal,
  equal-variance residuals across groups (the ANOVA omnibus test), which is
  a stretch at n=3/group, so it is treated as a second opinion alongside the
  Welch's t-test version above rather than a replacement. The ssGSEA version
  of this figure is saved as
  `figures/ExtData5e_bh4_boxplot_by_emtstate_dunnett.svg`.

## Visualization

- **Cross-method score correlation**: Pearson correlation matrix comparing
  `EMT_KS_TumorSig` and the `Epi.tum`/`Mes.tum` GSVA scores; GSVA vs. ssGSEA
  scatter plots for all 9 gene sets, one panel per set.
- **BH4-focused comparisons**: raw log2 expression of each of BH4's 7 member
  genes vs. the EMT KS score.
- **Heatmaps** (`ComplexHeatmap`):
  - BH4 gene log2FC heatmap (DESeq2 + apeglm, described above), no sample
    clustering, genes/states in fixed order, significance-flagged.
  - Gene-set-score and BH4 gene-level heatmaps, samples ordered by
    `EMT_State` and z-scored per row, with `EMT_State` as a column color bar
    and no sample clustering (columns/rows both fixed-order):
    - Gene-set-score heatmap: all 7 gene-set scores for GSVA; the 5 pathway
      gene sets only (`Epi.tum`/`Mes.tum` excluded) for ssGSEA.
    - BH4 gene-level heatmap: raw log2 expression of BH4's 7 member genes
      (gene-level values, not a gene-set score), shown both z-scored per row
      and as raw log2 values.
- **Bar charts**: one panel per BH4 gene, bars grouped by `EMT_State` with
  one bar per replicate, shown both as per-gene z-scores and as raw log2
  expression.
- **Boxplots**: one panel per gene set and scoring method (GSVA, ssGSEA,
  UCell), box + jitter by `EMT_State`, y-axis limited to a tight
  whisker-based range (`boxplot.stats()`). The BH4-specific variants
  additionally overlay either the BH-adjusted pairwise Welch's t-test
  brackets or the ANOVA + Dunnett's test brackets described above.
- **Dot plot**: mean expression / percent-expressed summary of BH4's 7
  genes by `EMT_State`, colored by min-max normalized mean expression
  (`pctExp` was not used as the color channel, as it saturates near 100%
  for bulk TPM data with little dropout).

## Software

Analysis was performed in R (see `SUM149PT_EMT_scoring_plotting_sessionInfo.txt`
for the full session/package version listing), primarily using `tidyverse`,
`DESeq2` (with `apeglm` shrinkage), `GSVA`, `UCell`/`Seurat`,
`ComplexHeatmap`, `ggpubr`, and `multcomp` (Dunnett's test). Scoring
methodology mirrors the companion `TCGA_BRCA_EMT_GeneSet_Scoring.Rmd`
analysis, with grouping here based on each sample's own annotated
`EMT_State` label rather than a data-driven tertile split, since this
dataset (unlike the TCGA cohort) already has ground-truth EMT-spectrum
labels.
