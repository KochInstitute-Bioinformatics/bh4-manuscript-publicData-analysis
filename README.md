# BH4 paper — figure-generating code and data

This directory is a stand-alone, lightweight collection of the code and
input data needed to regenerate the figure panels for the BH4 paper. Each
dataset lives in its own directory with an R Notebook (`.Rmd`) that scores
and plots that dataset, and a `figures/` subdirectory holding the SVG
panels pulled into the paper.

## Layout

```
090826_BH4_paper/
├── singularity_Rstudio_451.sh   # Luria-specific launcher for an RStudio
│                                 # Server session (see "Running the
│                                 # analyses" below)
├── engaging_launch_Rstudio.sh   # Engaging-specific launcher for an
│                                 # RStudio Server session (see "Running
│                                 # the analyses" below)
├── R/                            # Local installs of packages missing from
│                                 # the container image (picked up
│                                 # automatically as R's default per-user
│                                 # library path -- see "Running the
│                                 # analyses"). Not tracked in git --
│                                 # reinstalled on first use.
├── Rcode/                        # Shared R utility code (e.g. gene-symbol
│                                 # resolution helpers used by Chen's prep
│                                 # notebook)
├── emt_code/                     # EMT-spectrum scoring code and marker/
│                                 # signature gene lists, shared across
│                                 # every dataset directory below
├── Baugh/                        # Dataset directory (see below)
├── Chen/                         # Dataset directory
├── TCGA/                         # Dataset directory
└── SUM149PT/                     # Dataset directory
```

Each dataset directory follows the same target pattern:

- A **prep** notebook that turns raw input data into the (typically large,
  gitignored) intermediate `.rds` checkpoint the dataset's scoring notebook
  starts from.
- A **scoring/plotting** `.Rmd` notebook that loads that checkpoint, scores
  it (EMT status and/or gene-set activity), and produces figure panels.
- A `figures/` subdirectory with the SVG panels pulled into the paper.
- A `*_sessionInfo.txt` file recording the R session/package versions used
  the last time each notebook was run.
- A `*_Methods.md` file documenting the processing/scoring methodology and
  the statistics behind every figure, in prose. All four directories now
  have one.

Three of the four directories (Baugh, Chen, TCGA) now split prep from
scoring/plotting into separate notebooks; SUM149PT still does both in one
notebook (see "Per-dataset status" below for why that's fine there). See
"Known gaps" for the fix history and what's intentionally left as-is.

## Running the analyses

Notebooks are run inside a container via an RStudio Server launcher script.
Two are provided, one per cluster this project has run on — pick whichever
matches where you're working; both do the same job (start an RStudio Server
session inside the project's container and print SSH-tunnel instructions to
connect to it) but are written for different clusters and so aren't
interchangeable:

| Script | Cluster | Container image | Container runtime |
| --- | --- | --- | --- |
| `singularity_Rstudio_451.sh` | MIT Luria | `docker://bumproo/r4_5_3_singlecell_bulk_rnaseq:latest` | Singularity |
| `engaging_launch_Rstudio.sh` | MIT Engaging (ORCD) | `docker://bumproo/bulk_r451:v1` | Apptainer |

Each is self-contained (SLURM `#SBATCH` directives, `module load` calls,
and the SSH-tunnel/login instructions it prints when launched are specific
to that cluster) — see the script itself for exact connection steps. Both
launch the same way: `sbatch <script>`, then follow the tunnel/login
instructions written to the job's `slurm.*.out` log. Each dataset directory
has its own `.Rproj` so that relative paths inside its `.Rmd` resolve
correctly when opened as its own RStudio project (working directory = that
dataset directory).

`R/` holds local installs of packages the container image lacks
(`multcomp`, `sandwich`, `TH.data` as of this writing), but **is not
tracked in git** (~150 binary package files, easily reinstalled — not
worth carrying in the repo). Both scripts launch their container with
`$HOME` bound to this directory (`-H $PWD:/home/rstudio`), and
`R/x86_64-pc-linux-gnu-library/4.5/` is R's own default per-user library
location (`~/R/<platform>-library/<R version>`), so once populated, R
finds it there automatically with nothing extra to configure. On a fresh
checkout, install these once from inside an RStudio session (they'll land
in the right place automatically):

```r
install.packages(c("multcomp", "sandwich", "TH.data"))
```

This only keeps working as long as the container's R version/platform
matches that path (`4.5`, `x86_64-pc-linux-gnu`) — true of both images
above (`r4_5_3_...` and `bulk_r451` are both R 4.5.x) — if either image is
ever upgraded to a different R version, packages will land in a new
`R/x86_64-pc-linux-gnu-library/<new version>/` directory instead.

Both scripts hardcode the same plaintext RStudio Server login password by
design (single-user session, reached only via an SSH tunnel) — this is a
deliberate choice, not an oversight, but worth a reminder to reconsider
before either script is shared outside the lab or published to Zenodo.

## Figure panels

| Dataset | Panel(s) | File(s) in `<dataset>/figures/` |
| --- | --- | --- |
| Baugh | Fig1a, Fig1b, Fig1c | `Fig1a_EMT_KS_Score_FeaturePlot_chakroborty_RdBu.svg`, `Fig1b_EMT_KS_Score_byCondition_boxplot_chakroborty.svg`, `Fig1c_EMT_Stem_Epi_Dotplot_byCluster_res0.4.svg` |
| Chen | Fig5b, ExtData5a–d | `Fig5b_ssgsea_BH4_Chen.svg`, `ExtData5a_ssgsea_Glutathione_Chen.svg`, `ExtData5b_ssgsea_CoQ10_Chen.svg`, `ExtData5c_ssgsea_Thioredoxin_Chen.svg`, `ExtData5d_ssgsea_Peroxiredoxin_Chen.svg` |
| TCGA | Fig5c | `Fig5c_ssgsea_BH4_TCGA_BRCA.svg` |
| SUM149PT | ExtData5e, ExtData5f | `ExtData5e_bh4_boxplot_by_emtstate_dunnett.svg`, `ExtData5f_heatmap_bh4genes_lfc_sig.svg` |

## Large files / data availability

Several `.rds` objects in this directory (pre-processed Seurat objects,
pseudobulk objects, scored objects) are large (tens of MB to several
hundred MB) and are **not tracked in git** (see `.gitignore`) to keep this
repository lightweight. Likewise untracked: `SUM149PT/from_geo/` (raw GEO
download), each dataset's own `data/` (Baugh's raw 10x directories,
Chen's raw h5ad, TCGA's Firehose matrix/metadata), `R/` (locally
installed packages — see "Running the analyses"), and every rendered
`*.nb.html` notebook preview (regenerated by re-knitting the source
`.Rmd`). To reproduce a dataset directory from scratch:

1. Obtain that directory's raw input data (see each notebook's early
   "Load ..." chunks for what's expected, and where it should come from).
2. Run the directory's prep notebook (if it has one), then its
   scoring/plotting notebook, top to bottom inside the provided container.

## Per-dataset status

- **Chen** is the model for how this should work.
  `ChenAtlas_malignant_pseudobulk_prep.Rmd` is not meant to be run
  routinely -- it's a one-shot notebook that downloads the raw source atlas
  (a CZ CELLxGENE h5ad, fetched via `download.file()` into `Chen/data/`
  on first run -- gitignored, since it's ~1-2GB and re-fetchable from the
  URL in the notebook), imports it, subsets to
  malignant cells, and aggregates to the donor-level pseudobulk checkpoint
  (`Atlas_malignant_pseudobulk_donor.rds`) that
  `ChenAtlas_malignant_pseudobulk.Rmd` then loads to do the actual EMT/
  gene-set scoring and figure generation. Nothing multi-GB is checkpointed
  along the way in the prep notebook -- intermediate Seurat objects are
  dropped once no longer needed. `*_Methods.md` is present and accurate.
  Not yet actually re-run end to end since being written, but per the
  user that verification isn't needed (see "Known gaps").

- **Baugh** now has `baugh_BH4paper_prep.Rmd`, which consolidates five
  upstream notebooks (per-sample import/doublet-calling → integration →
  annotation → tumor subsetting → single-cell KS scoring) into one file,
  with each stage guarded by `if (file.exists(...))` so re-running only
  redoes stages whose checkpoint is missing. Its output filenames line up
  correctly with what `baugh_BH4paper_figures.Rmd` expects, and
  `doublet_called_rds/` (Stage 1's output) is now populated.
  `baugh_BH4paper_figures_Methods.md` now documents both notebooks. Raw
  data now lives at `Baugh/data/` (GEO accession GSE261252, from Baugh et
  al. 2024 — see the prep notebook's "Raw data" section and `Known gaps`
  below), with `data/rename.sh` documenting the GSM-to-sample-directory
  file renaming needed before `Read10X()` can read it.

- **TCGA** now has `TCGA_BRCA_prep.Rmd`, which builds a Seurat object from
  the Broad Firehose batch-corrected pan-cancer matrix and clinical
  metadata and subsets it to BRCA; its save-filename bug and its
  undefined-`tcga`-object bug (see prior review) have both been fixed/
  removed by the user. `TCGA_BRCA_EMT_GeneSet_Scoring_Methods.md` now
  documents both notebooks. The original GDC/Broad Firehose hosting pages
  for its raw inputs are dead; provenance is now tracked instead via a
  citation to Liu et al. 2018 (*Cell* 173(2), 400–416.e11), added to
  prep's new "Raw data" section — see "Fixed since the last review."

- **SUM149PT** does prep and scoring/plotting in one notebook
  (`SUM149PT_EMT_scoring_plotting.Rmd`) rather than splitting them out. This
  is reasonable given how little prep work is involved: the notebook reads
  three already-processed GEO supplementary files directly (normalized
  TPM matrix, raw counts, gene annotation) with no heavy
  import/integration/subsetting step to isolate into its own notebook.
  `nf-core_fetchngs.sh`/`ids.csv`/`from_geo/results/` are kept
  intentionally, for anyone who wants to obtain and reprocess the raw
  reads instead of using the GEO-supplied processed matrices — not dead
  files to remove.

## Known gaps

This section tracks known issues to resolve before this directory is
considered publication/Zenodo-ready.

**Fixed since the last review:**

- **Baugh prep's Stage 5 file paths didn't resolve.**
  `baugh_BH4paper_prep.Rmd`'s Stage 5 chunk (single-cell KS scoring) read
  `source("emt_code/Chakraborty_et_al_KS_code_vCAW.R")`,
  `read.xlsx("emt_code/EM_gene_signature_tumor_KS.xlsx")`,
  `read.xlsx("Mouse_Human_Orthologs_OnetoOne.xlsx")`, and
  `write.xlsx(..., "emt_code/EM_gene_signature_tumor_mouse_KS.xlsx")` —
  all without a `../` prefix. Since this notebook's working directory is
  `Baugh/` (its own `.Rproj`) and `emt_code/`/
  `Mouse_Human_Orthologs_OnetoOne.xlsx` live at the repo root, none of
  those paths resolved. Fixed by adding `../` to all four, matching what
  `baugh_BH4paper_figures.Rmd` already does correctly. (The prep
  notebook's own "Summary" section previously claimed the reverse — that
  the figures notebook's `../emt_code/...` paths were the broken ones —
  which was no longer accurate given where these files actually sit;
  that note has been corrected too.)
- TCGA's `saveRDS()`/`readRDS()` filename mismatch (`TCGA_BRCA_test.rds`
  vs. `TCGA_BRCA_seurat.rds`) — fixed by the user.
- TCGA prep's undefined-`tcga`-object bug (the "Check whether additional
  per-sample normalization is needed" section) — the user removed that
  section entirely; it was prep work already done elsewhere and not
  necessary in this notebook.
- **Baugh's raw-data provenance.** The user copied the 8 raw 10x outs
  directories into `Baugh/data/` (scRNA-seq from Baugh et al. 2024,
  *Clin Exp Metastasis* 41, 733–746,
  <https://doi.org/10.1007/s10585-024-10289-z>; GEO accession
  [GSE261252](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261252)),
  with `data/rename.sh` mapping GEO's per-GSM filenames to the standard
  10x names `Read10X()` expects. `baugh_BH4paper_prep.Rmd`'s Stage 1 now
  documents this (citation, accession, and the `rename.sh` step) and reads
  from `data/` (previously `../data/`, which never resolved to anything —
  see the now-fixed Stage 5 path bug above for the same root cause).
  `Baugh/data/` and `TCGA/data/` (previously untracked by no `.gitignore`
  rule at all — a latent risk of committing ~2GB of raw TCGA input) are
  now both gitignored, matching `SUM149PT/from_geo/` and the root-level
  `data/`; `Baugh/data/rename.sh` is kept as tracked provenance code via a
  `!` negation.
- **TCGA's raw-data provenance.** The GDC/Broad Firehose pages that
  originally hosted `geneexpressionmatrix_uniqueIDs.txt`
  (`EBPlusPlusAdjustPANCAN` batch-corrected matrix) and
  `Updated_Metadata_fixed.txt` (matched clinical metadata) are no longer
  live, so no working download link exists. Per the user, citing the
  associated publication is sufficient provenance tracking here:
  Liu, J., Lichtenberg, T., Hoadley, K.A. et al. An Integrated TCGA
  Pan-Cancer Clinical Data Resource to Drive High-Quality Survival
  Outcome Analytics. *Cell* 173(2), 400–416.e11 (2018).
  <https://doi.org/10.1016/j.cell.2018.02.052>. Added to
  `TCGA_BRCA_prep.Rmd`'s new "Raw data" section.
- **SUM149PT's raw-data provenance.** `SUM149PT_EMT_scoring_plotting.Rmd`
  and its Methods.md now cite the source publication: Brown, M.S.,
  Abdollahi, B., Wilkins, O.M. et al. Phenotypic heterogeneity driven by
  plasticity of the intermediate EMT state governs disease progression
  and metastasis in breast cancer. *Sci. Adv.* 8, eabj8002 (2022).
  <https://doi.org/10.1126/sciadv.abj8002>. The raw data
  (`from_geo/GSE172609_*.tsv`, `Human.GRCh38.p13.annot.tsv`) comes from
  GEO accession GSE172609, the RNA-seq SubSeries of that study's
  SuperSeries, GSE172613.
- **Chen consolidated onto the same per-dataset `data/` convention.**
  The user added `Chen/data/` and updated
  `ChenAtlas_malignant_pseudobulk_prep.Rmd` to download its raw h5ad
  there instead of a shared repo-root `data/` (nothing had been
  downloaded yet, so this was a zero-cost change). A leftover `../data/`
  reference in the notebook's prose (the code paths were already
  correctly updated) has been fixed to match. All four dataset
  directories now use the same self-contained `<Dataset>/data/` layout;
  the root-level `/data/` `.gitignore` rule has been replaced with
  `Chen/data/` alongside `Baugh/data/` and `TCGA/data/`.
- **`.gitignore` re-reviewed for large/primary files.** Two gaps closed:
  `R/` (the local package-install library — ~150 binary files, previously
  tracked in full) is no longer tracked at all — install `multcomp`,
  `sandwich`, `TH.data` from inside R instead (see "Running the
  analyses"); and every rendered `*.nb.html` notebook preview (up to
  ~5MB each, all seven currently on disk) is now ignored too, since
  they're fully regenerated by re-knitting the source `.Rmd`. Everything
  else audited clean: all raw per-dataset `data/` directories and
  `SUM149PT/from_geo/` were already covered, and `.Rproj.user/`/
  `.Rhistory`/etc. already matched inside every dataset subdirectory.

**Methods documentation:**

All four dataset directories now have a `*_Methods.md` file (Baugh's
`baugh_BH4paper_figures_Methods.md` added 2026-09-11, alongside TCGA's),
each with a "Statistical comparisons" section checked against that
directory's figures notebook code directly (not just its prose):
correlation via `stat_cor()`, pairwise Wilcoxon/BH, Welch's t-tests/BH,
ANOVA+Dunnett's, or (Baugh) a single unadjusted Welch's t-test, as
applicable to each. `baugh_BH4paper_prep.Rmd`'s header no longer links to
the nonexistent `METHODS.md`/`baugh_notes.md` — it now points at
`baugh_BH4paper_figures_Methods.md`.

**Kept intentionally (not gaps):**

- `SUM149PT/nf-core_fetchngs.sh`, `ids.csv`, and `from_geo/results/` are
  kept on purpose, for anyone who wants to obtain and process the raw
  reads themselves rather than use the GEO-supplied processed matrices
  that `SUM149PT_EMT_scoring_plotting.Rmd` actually reads. The unused
  `library(tximport)` call in that notebook is being removed separately.
- Baugh's and Chen's prep notebooks have no `<name>_prep_sessionInfo.txt`
  in their directories, meaning they haven't been executed since being
  authored/consolidated (every stage is guarded by
  `if (file.exists(...))` against pre-existing checkpoints, so re-running
  them now would just skip every stage without exercising the new code).
  Per the user, a real from-scratch verification run isn't needed — not
  treated as a gap. (TCGA's prep notebook *has* been run — its
  sessionInfo file is dated today — which is how its filename mismatch
  was caught.)

## License

See `LICENSE` (MIT).
