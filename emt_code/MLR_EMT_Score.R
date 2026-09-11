# MLR ordinal-logistic EMT scoring (George et al., 2017; Chakraborty et al.,
# 2020), ported from the MATLAB implementation at
# https://github.com/priyanka8993/EMT_score_calculation/tree/master/MLR_Matlab_Code
# (MLR3.m + EMT_automated.m).
#
# The published model's 4 coefficients (`MLR.B1`) and its 3-class ordinal-
# logistic prediction formula were fit on NCI60 microarray data and can't be
# applied to RNA-seq-scale expression values as-is -- the model's fixed
# intercepts are calibrated to the microarray platform's numeric range.
# Rather than reproduce the original pipeline's platform mean-shift
# normalization (which assumes a shared microarray-to-microarray
# relationship that doesn't hold for bulk/pseudobulk/single-cell RNA-seq),
# each predictor is instead quantile-mapped onto its NCI60 training
# distribution before scoring: a sample's percentile rank for a given
# predictor, within its own dataset, is substituted with the NCI60 training
# value at that same percentile. This corrects for scale AND shape
# differences between platforms/assays, applies the same way regardless of
# input technology, and needs nothing beyond the two predictor genes' NCI60
# reference values (no external probe-annotation file required, unlike the
# original pipeline's GPL570-based normalization).
#
# The ordinal-logistic formula and `MLR.B1` below were reverse-engineered
# from MATLAB's own `mnrval(B1, X, 'model', 'ordinal')` and validated
# directly against it on a 10-point test grid spanning the NCI60 predictor
# ranges (2026-08-03) -- max abs difference in class probabilities was <1e-4.
#
# MLRScore() (quantile mapping, above) is the default/baseline calibration.
# Two alternatives are also defined below, added to test whether either
# improves concordance with the KS score on donor-level pseudobulk data:
# MLRScoreZ() (robust z-score rescaling) and MLRScoreHK() (housekeeping-gene
# anchor). See MLR_EMT_Score_notes.md for the reasoning behind each and why
# quantile/z-score mapping structurally can't detect a whole-cohort shift
# relative to NCI60.

# B1 = [intercept_1, intercept_2, slope_CLDN7, slope_VIM/CDH1_ratio], read
# directly from RelevantData.mat in the MLR_Matlab_Code repo.
MLR.B1 <- c(-7.87139071184979, 0.0412975537779503, 1.35705632225475, -1.95661140686887)

# NCI60 reference values for the model's two predictor probes (CLDN7:
# 202790_at; VIM: 201426_s_at; CDH1: 201131_s_at), extracted directly from
# RelevantData.mat -- used only as the target distribution for quantile
# mapping below, never fed into the model directly.
MLR.NCI60Reference.file <- "emt_code/NCI60_MLR_reference.csv"

#' Map each finite element of `x` onto the corresponding quantile of `ref`,
#' using x's own within-vector percentile rank (ties averaged). `NA` in `x`
#' stays `NA`; non-finite values in `ref` are dropped before computing
#' quantiles. Infinite values in `x` are kept and ranked normally (e.g. a
#' VIM/CDH1 ratio of Inf from CDH1 == 0 is a real, maximally-mesenchymal
#' signal, not a missing value) -- ties among them just map to the same
#' top-of-distribution reference value.
.quantileMap <- function(x, ref) {
  ref <- ref[is.finite(ref)]
  out <- rep(NA_real_, length(x))
  ok  <- !is.na(x)
  p   <- (rank(x[ok], ties.method = "average") - 0.5) / sum(ok)
  out[ok] <- unname(quantile(ref, probs = p, na.rm = TRUE, type = 7))
  out
}

#' Shared ordinal-logistic scoring step for all MLRScore* variants below --
#' takes already-calibrated (mapped/corrected) predictor values and applies
#' the published model's fixed coefficients.
.MLROrdinalScore <- function(cldn7_mapped, ratio_mapped) {
  eta   <- MLR.B1[3] * cldn7_mapped + MLR.B1[4] * ratio_mapped
  P_le1 <- plogis(MLR.B1[1] + eta)
  P_le2 <- plogis(MLR.B1[2] + eta)

  P_Epi    <- P_le1
  P_Hybrid <- P_le2 - P_le1
  P_Mes    <- 1 - P_le2

  # matches EMT_automated.m's ScoreEMT3: distance from whichever pure
  # endpoint (Epi or Mes) has higher probability, using P_Hybrid as that
  # distance -- 0/2 when P_Hybrid is ~0 (confidently pure), 1 when P_Hybrid
  # dominates (confidently hybrid)
  score <- ifelse(P_Epi > P_Mes, P_Hybrid, 2 - P_Hybrid)

  data.frame(P_Epi = P_Epi, P_Hybrid = P_Hybrid, P_Mes = P_Mes, EMT_MLR_Score = score)
}

#' MLR EMT score (George et al., 2017 ordinal-logistic method), calibrated
#' via quantile mapping onto NCI60
#'
#' @param cldn7 numeric vector of CLDN7 expression
#' @param vim   numeric vector of VIM expression, same units/scale as cldn7 and cdh1
#' @param cdh1  numeric vector of CDH1 expression
#' @return data.frame (row-aligned to the inputs) with columns P_Epi,
#'   P_Hybrid, P_Mes, and EMT_MLR_Score (0 = pure epithelial, 1 = maximally
#'   hybrid, 2 = pure mesenchymal)
MLRScore <- function(cldn7, vim, cdh1) {
  stopifnot(length(cldn7) == length(vim), length(vim) == length(cdh1))

  nci60_ref <- read.csv(MLR.NCI60Reference.file)

  cldn7_mapped <- .quantileMap(cldn7, nci60_ref$CLDN7)

  # VIM and CDH1 are quantile-mapped individually, *before* taking their
  # ratio -- computing the ratio first on raw single-cell values divides by
  # a huge pile of near-zero/zero CDH1 (dropout is common for this gene in
  # this dataset), producing mostly Inf/NaN that then floods the mapped
  # ratio's top tail regardless of a cell's actual biology. NCI60's own
  # CDH1/VIM values are continuous log2 microarray intensities (never near
  # zero), so mapping each gene onto its own NCI60 distribution first means
  # the ratio is only ever taken between two well-behaved, non-zero-
  # inflated values.
  vim_mapped   <- .quantileMap(vim, nci60_ref$VIM)
  cdh1_mapped  <- .quantileMap(cdh1, nci60_ref$CDH1)
  ratio_mapped <- vim_mapped / cdh1_mapped

  .MLROrdinalScore(cldn7_mapped, ratio_mapped)
}

#' Map `x` onto `ref`'s scale via a robust z-score (median/MAD) transform,
#' rather than `.quantileMap()`'s within-vector rank. Still centers `x`'s
#' own median on `ref`'s median (same whole-cohort-shift blind spot as
#' quantile mapping -- see MLR_EMT_Score_notes.md), but the transform is
#' smooth rather than rank-based, so it isn't forced to spread tied/
#' near-tied values across `ref`'s full range the way quantile mapping can
#' with a small cohort. `NA` in `x` stays `NA`.
.zscoreMap <- function(x, ref) {
  ref <- ref[is.finite(ref)]
  z <- (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE)
  median(ref) + z * mad(ref)
}

#' MLR EMT score, calibrated via robust z-score rescaling onto NCI60 instead
#' of quantile mapping -- see MLRScore() for arguments/return value.
MLRScoreZ <- function(cldn7, vim, cdh1) {
  stopifnot(length(cldn7) == length(vim), length(vim) == length(cdh1))

  nci60_ref <- read.csv(MLR.NCI60Reference.file)

  cldn7_mapped <- .zscoreMap(cldn7, nci60_ref$CLDN7)
  vim_mapped   <- .zscoreMap(vim, nci60_ref$VIM)
  cdh1_mapped  <- .zscoreMap(cdh1, nci60_ref$CDH1)
  ratio_mapped <- vim_mapped / cdh1_mapped

  .MLROrdinalScore(cldn7_mapped, ratio_mapped)
}

# Housekeeping-gene NCI60 reference (ACTB, GAPDH, B2M, TBP, PPIA), extracted
# from RelevantData.mat the same way as MLR.NCI60Reference.file -- probe IDs
# looked up via the GPL570 annotation file (GPL570.annot.gz, ~8MB, much
# smaller than the original pipeline's normalization assumed -- see
# MLR_EMT_Score_notes.md).
MLR.NCI60Housekeeping.file <- "emt_code/NCI60_housekeeping_reference.csv"

#' MLR EMT score, calibrated via a housekeeping-gene anchor instead of
#' quantile/z-score mapping. Estimates a single additive offset (median
#' expression across ACTB/GAPDH/B2M/TBP/PPIA, pooled) between this cohort
#' and NCI60, then subtracts that one constant directly from CLDN7/VIM/CDH1's
#' raw values -- no per-predictor-gene mapping at all. Unlike quantile/
#' z-score mapping, this *can* represent a genuine whole-cohort shift
#' relative to NCI60, since the correction doesn't depend on this cohort's
#' own predictor-gene distribution. Tradeoffs: assumes the housekeeping
#' genes' offset transfers to CLDN7/VIM/CDH1's own expression range, and
#' (like the original pipeline's mean-shift, which this narrows to a small
#' definitionally-stable panel instead of a possibly EMT-associated 20-gene
#' panel) still assumes a single additive constant is meaningful across
#' assay types with different noise models/dynamic range -- see
#' MLR_EMT_Score_notes.md. Because it doesn't bound predictors away from
#' zero the way quantile/z-score mapping do, a donor with corrected CDH1
#' very close to zero can still produce an unstable ratio.
#'
#' @param housekeeping data.frame/matrix of raw housekeeping-gene expression
#'   (same units as cldn7/vim/cdh1), one row per sample, columns matching
#'   the gene columns of MLR.NCI60Housekeeping.file (ACTB, GAPDH, B2M, TBP, PPIA)
MLRScoreHK <- function(cldn7, vim, cdh1, housekeeping) {
  stopifnot(length(cldn7) == length(vim), length(vim) == length(cdh1))

  nci60_ref <- read.csv(MLR.NCI60Reference.file)
  nci60_hk  <- read.csv(MLR.NCI60Housekeeping.file)

  hk_genes <- setdiff(colnames(nci60_hk), "nci60_index")
  stopifnot(all(hk_genes %in% colnames(housekeeping)))

  # pool all samples x all housekeeping genes into one distribution on each
  # side -- a single global offset, not five separate per-gene ones, since
  # the goal is one correction constant representing "how much higher/lower
  # does this cohort run compared to NCI60"
  offset <- median(unlist(housekeeping[, hk_genes]), na.rm = TRUE) -
    median(unlist(nci60_hk[, hk_genes]), na.rm = TRUE)

  cldn7_mapped <- cldn7 - offset
  vim_mapped   <- vim - offset
  cdh1_mapped  <- cdh1 - offset
  ratio_mapped <- vim_mapped / cdh1_mapped

  result <- .MLROrdinalScore(cldn7_mapped, ratio_mapped)
  result$housekeeping_offset <- offset
  result
}
