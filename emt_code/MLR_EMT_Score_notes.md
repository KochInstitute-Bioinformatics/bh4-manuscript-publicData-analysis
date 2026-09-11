# MLR EMT Scoring: Porting Notes

## Background

The MLR (ordinal multinomial logistic regression) method scores samples along
an Epithelial (0) → Hybrid (1) → Mesenchymal (2) continuum, based on
`George et al., 2017` (method) and `Chakraborty et al., 2020` (comparative
benchmark against KS-score and 76GS). The original implementation is MATLAB
code at
[priyanka8993/EMT_score_calculation](https://github.com/priyanka8993/EMT_score_calculation/tree/master/MLR_Matlab_Code),
built around a fixed, pre-trained ordinal logistic model:

- Predictors: `CLDN7` expression, and the `VIM`/`CDH1` expression ratio
- Coefficients (`B1`): 2 intercepts (thresholds) + 2 slopes, fit on **NCI60
  microarray data** (Affymetrix HG-U133 Plus 2.0 / GPL570)
- Prediction: 3-class ordinal logistic regression (`mnrval(..., 'model',
  'ordinal')` in MATLAB), collapsed to a continuous 0–2 score using the
  hybrid-class probability as a distance from whichever pure class (E or M)
  is more probable

## The normalization problem

The original pipeline normalizes a new dataset onto the NCI60 microarray
scale via a **mean-shift**: it averages a 20-gene reference panel in the new
data, compares that average to the same panel's average in NCI60 (matched
via the GPL570 probe annotation table), and subtracts the difference from
the new dataset's values.

This doesn't transfer to RNA-seq/scRNA-seq input:
- The 20-gene mean-shift assumes a simple additive relationship between two
  *microarray* datasets — not defensible across fundamentally different
  assay types with different noise models and dynamic ranges. Specifically:
  - **Different measurement physics.** Array intensity saturates at the high
    end (probe binding capacity) and is noise-floor-limited at the low end
    (background cross-hybridization), giving a narrower dynamic range
    (~2-3 orders of magnitude) than RNA-seq counts (~4-5 orders of
    magnitude). A correction sized for one part of that range isn't
    constant across the whole range, so a single global shift can't
    represent it.
  - **Different noise models.** Array noise is roughly log-normal
    (multiplicative, becomes additive after log2 transform); RNA-seq/
    scRNA-seq counts are Poisson/negative-binomial, with variance scaling
    with the mean in a way arrays don't exhibit. A location-shift implicitly
    assumes the same kind of noise is just being recentered, not a
    structurally different error process.
  - **Zero-inflation has no array equivalent.** This is what actually broke
    downstream (see the CDH1 fix below): CDH1 was ~64% zero in raw scRNA-seq
    counts — dropout, a detection failure, not a low-but-real value. A
    mean-shift can only move a distribution's center; it can't recover
    signal that was never captured, and shifting a pile of hard zeros just
    distorts the non-zero values instead.
  - **Panel-average generalization is weaker across assay types.** Even
    array-to-array, using one 20-gene panel's average offset to correct the
    whole transcriptome is already an approximation. Across assay types it's
    a bigger leap: array-specific artifacts (probe cross-reactivity, GC
    bias) have no RNA-seq counterpart, so there's no reason to expect the
    panel's offset represents the transcriptome-wide offset.

  This is why quantile mapping was adopted instead (below): it only
  requires each gene's rank/percentile to correspond between the new data's
  distribution and NCI60's, not the same absolute scale or noise model — a
  weaker, better-justified assumption across assay types (at the cost of
  the whole-cohort-shift blindness discussed later in this document).
- The required GPL570 probe annotation file wasn't included in the repo and
  wasn't needed for the approach adopted here. **Correction (2026-08-05):**
  this was assumed to be an impractically large download at the time, but
  turned out not to be true — the summarized annotation file
  (`GPL570.annot.gz`) is only ~8MB and was later fetched directly from NCBI
  GEO's FTP to look up a handful of housekeeping-gene probe IDs for the
  housekeeping-anchor method below. The real obstacle was never file size;
  it was that a probe-to-symbol crosswalk is a bigger, separate problem
  from a location-shift normalization's own defensibility (previous bullet)
  — for the original quantile-mapping approach, no crosswalk was needed at
  all, so this was never actually blocking.

## What was verified directly (not assumed)

`RelevantData.mat` was inspected directly via Python/`scipy.io.loadmat`
(MATLAB itself unavailable for scripting file I/O conveniently, but
`module load matlab/2016b` **is** available on this cluster and was used for
validation — see below). Confirmed contents relevant to this port:

- `B1 = [-7.87139071, 0.04129755, 1.35705632, -1.95661141]`
- `DataNCI60` (54703 probes × 60 NCI60 cell lines) and `LabelsNCI60` (probe
  IDs) — **probe-indexed, not gene-symbol-indexed**
- `GeneList1` confirms the predictor probes: CLDN7 = `202790_at`, VIM =
  `201426_s_at`, CDH1 = `201131_s_at`
- All 3 probes matched to exactly one row each in `LabelsNCI60`; extracted
  their 60-cell-line values directly (1 NaN cell line per probe) — written
  to [`NCI60_MLR_reference.csv`](NCI60_MLR_reference.csv). This sidesteps
  the GPL570 dependency entirely: only these 2 predictor genes' NCI60
  reference values are needed, not a full probe-to-symbol crosswalk.
- All 25 row-position-critical genes from the original pipeline's
  `genes_for_EMT_score.txt` were confirmed present in Wu2021's `Genes.txt`
  (not a blocker for this dataset).

**Formula validation:** wrote a small MATLAB test script
(`module load matlab/2016b`) calling `mnrval(B1, X, 'model', 'ordinal')` on
a 10-point grid spanning the NCI60 predictor ranges. Reverse-engineered the
cumulative-logit formula below and confirmed it reproduces MATLAB's class
probabilities to within <1e-4 across the whole grid:

```
eta    = B1[3] * CLDN7 + B1[4] * (VIM/CDH1 ratio)
P(Y<=1) = plogis(B1[1] + eta)      # class 1 = Epithelial
P(Y<=2) = plogis(B1[2] + eta)
P_Epi    = P(Y<=1)
P_Hybrid = P(Y<=2) - P(Y<=1)
P_Mes    = 1 - P(Y<=2)
score    = ifelse(P_Epi > P_Mes, P_Hybrid, 2 - P_Hybrid)
```

## The adopted approach: quantile mapping instead of mean-shift

Each predictor (`CLDN7`, `VIM`, `CDH1`) is individually quantile-mapped onto
its NCI60 training distribution: a sample's percentile rank for that gene,
*within its own dataset*, is substituted with the NCI60 value at that same
percentile. This is technology-agnostic (applies the same way to bulk,
pseudobulk, and single-cell input) and needs nothing beyond the NCI60
reference CSV above.

Implementation: [`MLR_EMT_Score.R`](MLR_EMT_Score.R) — `MLRScore(cldn7, vim,
cdh1)`.

### Fix: map VIM and CDH1 *before* taking their ratio, not after

First version computed the `VIM/CDH1` ratio on raw values, then
quantile-mapped the ratio. On the single-cell data (`Wu_04`), CDH1 is
broadly undetected (~64% zero even in raw counts), so most cells produced
`Inf` or `NaN` ratios before any correction was applied — these all
quantile-mapped to the extreme top of the NCI60 ratio distribution (or to
`NA`), regardless of true biology. Result: `EMT_MLR_Score` saturated near
2.0 for most cells, with >14,000 `NA`s.

**Why it broke this way.** CDH1 dropout in scRNA-seq (a detection failure,
not necessarily true biological absence) splits into two cases at the raw
`vim_raw / cdh1_raw` division:

- `vim_raw > 0, cdh1_raw == 0` → ratio = `Inf`
- `vim_raw == 0, cdh1_raw == 0` → ratio = `0/0 = NaN`

`.quantileMap()` (`MLR_EMT_Score.R:44-51`) checks `is.finite()` only on the
*reference* vector, not on `x` — infinite values in `x` are kept and ranked
normally, so every dropout-driven `Inf` lands at the very top of the NCI60
ratio distribution ("maximally mesenchymal"), regardless of how much real
VIM signal that cell actually had. `NaN` satisfies R's `is.na()`, so it's
treated like a missing input and passed straight through as `NA`. With ~64%
CDH1 dropout, that's enough `Inf`s to saturate most scores near 2.0, and
enough simultaneous VIM/CDH1 zeros to produce the `NA` pileup.

**Fix:** quantile-map `VIM` and `CDH1` *separately* first (each onto its own
NCI60 reference distribution, which is continuous and never near zero),
then take the ratio of the two mapped values. This removed the raw-scale
division-by-near-zero problem and the NA pileup.

**Why mapping first works.** NCI60's reference values are continuous log2
microarray intensities — never zero, never near-zero, since array
hybridization intensity has no dropout floor the way scRNA-seq counts do.
Mapping raw CDH1 (including the pile of technical zeros) onto that
reference means every value, even all the tied zeros, lands on some
**finite, non-zero** point of the NCI60 CDH1 distribution: all
`cdh1_raw == 0` cells are tied for the lowest rank, `rank(...,
ties.method = "average")` gives them all the same averaged rank, and they
all map to the exact same (low but finite) NCI60 quantile value — same for
VIM. With both operands of the ratio guaranteed finite and bounded away
from zero, `vim_mapped / cdh1_mapped` can never produce `Inf` or `NaN`; the
division-by-near-zero problem is eliminated structurally rather than
patched.

**What this does and doesn't fix.** It removes the *artifact*: a cell no
longer gets snapped to an extreme score purely because of a detection
failure in one gene. It does **not** recover lost biological information —
a true dropout cell and a "barely detected" cell both collapse to the same
tied rank and the same mapped CDH1 value, since the method can't
distinguish "no CDH1 mRNA" from "CDH1 mRNA present but not captured by
sequencing." Those cells still get pushed toward the low-CDH1 / high-ratio
(more mesenchymal-leaning) end, just now as a bounded, defensible rank
position instead of an unbounded artifact. Any `NA`s remaining in the
output now come only from genuinely missing input values (`NA` in the raw
data itself, not zeros) — the expected behavior of `.quantileMap()`, not
the pileup that motivated this fix.

## Current findings

- **Single-cell (`Wu_04`)**: after the fix, no more NA pileup, but scores
  are still skewed toward Mesenchymal/Hybrid, in tension with the KS score
  (which calls most cells Epithelial).
- **Pseudobulk (`Wu_05`)**: most of the 20 patients score as Hybrid
  (`EMT_MLR_Score` near 1) rather than clearly Epithelial, despite these
  being solid breast tumor samples.
- **Gene-level sanity check (`mlr-vs-genes-scatter` in `Wu_05`)**: plotted
  `EMT_MLR_Score` against its own inputs. All three relationships go in the
  direction the model's coefficients predict:
  - `CLDN7` vs. score: negative (higher CLDN7 → more Epithelial, as expected)
  - `CDH1` vs. score: negative (higher CDH1 → smaller ratio → more Epithelial)
  - `VIM` vs. score: positive (higher VIM → larger ratio → more Mesenchymal)

  This is good evidence the R port itself is behaving correctly — the
  relationships aren't inverted or broken.

## Working interpretation

The directionally-correct-but-overall-shifted pattern points to a
**calibration effect of quantile mapping against NCI60**, not an
implementation bug or proof that MLR "doesn't work" on RNA-seq data.

**ECDF diagnostic (`nci60-ecdf-diagnostic` chunk in `Wu_05.Rmd`):** plots
NCI60's empirical CDF for `CLDN7`/`VIM`/`CDH1` against each of the 20
patients' quantile-mapped values. This refined an earlier, less precise
version of this explanation (originally framed as "this cohort's real range
is narrower than NCI60's and shifted toward the mesenchymal end"). What the
diagnostic actually demonstrates is more fundamental:

`.quantileMap()` only encodes each patient's *rank within this 20-patient
cohort* — by construction, the lowest-ranked patient always lands near
NCI60's low percentile, the highest-ranked near its high percentile, and
critically, **the median patient always lands at NCI60's median value**,
regardless of what this cohort's true absolute biological level actually
is. The diagnostic plot confirms this directly: the mapped points spread
roughly evenly across NCI60's full percentile range rather than clustering
in one shifted slice, which is exactly what within-cohort rank-based
quantile mapping guarantees regardless of the input data's real biology.

This means the method has **no way to represent "this whole cohort is more
epithelial than NCI60's median"** — it can only express relative ordering
*within* the cohort being scored, never an absolute shift of the whole
cohort relative to the reference population. If NCI60's own median
CLDN7/VIM/CDH1 profile sits close to the model's Hybrid decision region
(plausible, since NCI60 spans full-spectrum cancer cell lines including
frankly mesenchymal ones), then any new cohort — however genuinely
epithelial in absolute terms — will have its median patient pulled toward
Hybrid purely as an artifact of the mapping, not as a biological finding.

**Practical takeaway:** treat `EMT_MLR_Score`'s relative ordering within a
cohort as informative, but don't over-interpret its absolute position on the
published 0–2 (E/Hybrid/M) scale as a literal biological classification for
this data type — the quantile-mapping approach used here is structurally
unable to capture whole-cohort shifts relative to NCI60.

## Alternative calibration methods (pseudobulk, Atlas_04)

Once scoring moved to donor-level pseudobulk (`Atlas_04_Malignant_pseudobulk.Rmd`),
CDH1 dropout is far less of a concern than in single-cell, so it became
worth directly testing whether an alternative to quantile mapping improves
concordance with the KS score -- the whole-cohort-shift blind spot above
was the leading suspect for the poor concordance observed (donors reading
MLR≈1-2 despite a low, clearly-Epithelial KS score).

Both alternatives are implemented in `MLR_EMT_Score.R` alongside the
original `MLRScore()`, compared side by side in `Atlas_04`'s
`mlr-vs-ks-scatter-methods` chunk:

- **`MLRScoreZ()` — robust z-score rescaling.** Replaces `.quantileMap()`'s
  within-cohort rank with a median/MAD location-scale transform onto NCI60.
  Still centers this cohort's median predictor on NCI60's median (same
  whole-cohort-shift blindness as quantile mapping), but smoother and less
  prone to tie-driven saturation with a ~136-donor cohort. Low risk, direct
  swap-in.
- **`MLRScoreHK()` — housekeeping-gene anchor.** Estimates a single
  additive offset from the median difference between this cohort's and
  NCI60's housekeeping-gene expression (`ACTB`, `GAPDH`, `B2M`, `TBP`,
  `PPIA`, pooled), then subtracts that one constant directly from raw
  CLDN7/VIM/CDH1 -- no rank- or distribution-based mapping at all. This is
  the one alternative that can actually represent a genuine whole-cohort
  shift relative to NCI60, since the correction doesn't depend on this
  cohort's own predictor-gene distribution the way quantile/z-score mapping
  do. Tradeoffs: (1) still assumes a single additive constant, estimated
  from 5 stable genes, transfers meaningfully to CLDN7/VIM/CDH1's own
  expression range and to this assay's noise model -- i.e. it's a narrower,
  more principled version of the original pipeline's mean-shift (previous
  section), not a fix for the "different assay types" defensibility problem
  discussed there, just a smaller, better-justified panel to estimate the
  shift from; (2) unlike quantile/z-score mapping, nothing bounds the
  corrected CDH1 away from zero, so a donor with very low CDH1 can still
  produce an unstable ratio the way raw single-cell scoring originally did.

**Getting the housekeeping probe IDs:** `RelevantData.mat`'s `GeneList1/2/3`
only cover the model's own predictor/ratio genes (confirmed by loading the
file directly in MATLAB, `module load matlab/2016b`) — no general
probe-to-symbol crosswalk is bundled with it. Probe IDs for the 5
housekeeping genes were instead looked up in `GPL570.annot.gz` (fetched
fresh from NCBI GEO's FTP — see the correction above), matched on the exact
`Gene symbol` column, then their NCI60 values extracted from
`RelevantData.mat`'s `DataNCI60`/`LabelsNCI60` the same way as the original
3 predictor probes. Reference values are saved to
`NCI60_housekeeping_reference.csv`; the annotation file itself isn't needed
again and wasn't added to the repo.

### Result (2026-08-05): both alternatives broke, quantile mapping remains primary

Ran on the actual `Atlas_04` pseudobulk cohort. Both alternatives produced
clearly wrong output, for reasons traced directly to NCI60's own reference
statistics (`NCI60_MLR_reference.csv`/`NCI60_housekeeping_reference.csv`) —
not implementation bugs, but real design flaws in each calibration:

- **`MLRScoreZ()` — most donors saturated at score≈2, including the
  lowest-KS (most clearly Epithelial) ones.** Cause: NCI60's own CLDN7
  reference distribution is nearly degenerate — `median=2.336, mad=0.036`,
  with 59/60 lines packed into a razor-thin band and only a couple of real
  outliers pulling the max to 10.08. A MAD-based z-score transform breaks on
  a distribution shaped like that: `median(ref) + z*mad(ref)` maps almost
  any input z onto a value within ~0.04 of 2.34, crushing CLDN7 to a
  near-constant regardless of a donor's real CLDN7 level — effectively
  deleting CLDN7 from the model. With CLDN7's contribution to `eta` pinned,
  the score becomes almost entirely ratio-driven, and NCI60's own typical
  VIM/CDH1 ratio (median 13.68/3.71 ≈ 3.7) already sits well past this
  model's Mesenchymal threshold once CLDN7 is fixed — so nearly any donor
  near or above the *cohort's own* median ratio gets pushed to score≈2,
  independent of true biology. Quantile mapping doesn't have this failure
  mode because it never uses NCI60's spread as a divisor.
- **`MLRScoreHK()` — every donor scored near 0 (max ~0.015).** Cause: pooled
  NCI60 housekeeping median is ~13.98 (log2 microarray intensity for
  highly-expressed genes), while pseudobulk's housekeeping values are on
  Seurat's natural-log RNA-seq scale, which runs much lower for the same
  genes -- producing a large offset. Applying that one offset identically to
  CLDN7, VIM, *and* CDH1 breaks two ways at once: (1) adding the same large
  constant to both VIM and CDH1 crushes their ratio toward 1 (ratios aren't
  shift-invariant -- only differences are), destroying the VIM/CDH1 contrast
  the model depends on; (2) CLDN7 gets shifted far outside the ~2-10 range
  it was ever trained on. Working through the formula with representative
  values, both effects saturate `P_Epi≈1` for nearly every donor.

**Common thread:** both alternatives reintroduce sensitivity to *absolute*
expression scale (via NCI60's spread, or via raw log-unit magnitude), which
is exactly what quantile mapping was built to avoid by only ever using
within-cohort rank. That's why quantile mapping -- despite the documented
whole-cohort-shift blindness -- remains the better-behaved method here.
**Decision: keep `MLRScore()` (quantile mapping) as the primary/reported
method.** `MLRScoreZ()`/`MLRScoreHK()` are left in `MLR_EMT_Score.R` and
`Atlas_04` as a documented negative result, not for routine use.

## Open items

- **Not yet done:** applying `MLRScore()` to the bulk RNA-seq data
  (technology-agnostic, should work unchanged once CLDN7/VIM/CDH1 vectors
  are available from that dataset) — the same median-matching caveat above
  will apply there too.
- MATLAB was only needed for one-time formula/reference-value validation
  (originally the 3 predictor probes, now also the 5 housekeeping probes)
  — routine scoring going forward is pure R, no MATLAB dependency.

## Files

| File | Purpose |
|---|---|
| `MLR_EMT_Score.R` | `MLRScore()` (quantile mapping), `MLRScoreZ()` (robust z-score rescaling), `MLRScoreHK()` (housekeeping anchor) -- all wrapping shared ordinal-logistic scoring |
| `NCI60_MLR_reference.csv` | NCI60 CLDN7/VIM/CDH1 values (60 cell lines), extracted from `RelevantData.mat` |
| `NCI60_housekeeping_reference.csv` | NCI60 ACTB/GAPDH/B2M/TBP/PPIA values (60 cell lines), extracted from `RelevantData.mat` via GPL570.annot.gz probe lookup |
| `MLR_EMT_Score_notes.md` | This document |
