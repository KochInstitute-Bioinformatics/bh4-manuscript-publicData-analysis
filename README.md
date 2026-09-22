# BH4 paper — figure-generating code and data

This directory is a stand-alone, lightweight collection of the code and
input data needed to regenerate the figure panels for the BH4 paper. Each
dataset lives in its own directory under `datasets/` with an R Notebook
(`.Rmd`) that scores and plots that dataset, and a `figures/` subdirectory
holding the PDF/SVG panels pulled into the paper.

## Layout

```
090826_BH4_paper/
├── singularity_Rstudio_451.sh   # Luria-specific launcher for an RStudio
│                                 # Server session (see "Running the
│                                 # analyses" below)
├── engaging_launch_Rstudio.sh   # Engaging-specific launcher for an
│                                 # RStudio Server session (see "Running
│                                 # the analyses" below)
├── pack_input_data.sh            # Collects every dataset's .rds/data/
│                                 # from_geo/ into input_data.tar.gz for
│                                 # distribution outside git (see "Large
│                                 # files / data availability")
├── restore_input_data.sh         # Un-tars that bundle back onto
│                                 # datasets/ after a fresh clone
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
│                                 # every dataset directory below (each
│                                 # reaches it via a same-directory
│                                 # `emt_code`/`Rcode` symlink -- see below)
├── manuscript/                   # Manuscript drafts (figure layouts,
│                                 # text) -- not analysis code or data
└── datasets/
    ├── Baugh2024_GSE261252/               # mouse mammary tumor model
    ├── BRCA_atlas_CellxGene_Chen/          # CZ CELLxGENE BRCA atlas
    ├── BRCA_TCGA/                          # TCGA BRCA bulk RNA-seq
    ├── TNBC_cellline_SUM149PT_GSE172609/   # SUM149PT EMT-spectrum panel
    ├── PDAC_mouse_35952360/                # Pitter et al., mouse PDAC scRNA-seq
    ├── PDAC_human_38702773/                # Park et al., human PDAC scRNA-seq
    ├── PDAC_human_35952360/                # Lin et al. 2020, human PDAC scRNA-seq
    ├── SKCM_TCGA/                          # TCGA SKCM bulk RNA-seq
    ├── Neuroblastoma_human_25150838_SEQC_GSE62564/         # SEQC NB cohort, bulk
    ├── Neuroblastoma_human_40603900_Treehouse_GSE294351/   # Treehouse NB cohort, bulk
    └── Neuroblastoma_cellline_30948783_GSE90803_GSE116890/ # NB ADRN/MES cell lines, microarray
```

Each dataset directory:

- Is its own self-contained RStudio project (own `.Rproj`), and has a
  local `Rcode`/`emt_code` symlink pointing at the repo-root shared code
  (`ln -s ../../Rcode Rcode`, `ln -s ../../emt_code emt_code`) rather than
  a `../`-relative path -- so notebooks always `source("emt_code/...")`
  regardless of how deeply this directory is nested.
- Has a **scoring/plotting** `.Rmd` notebook (some also have a separate
  **prep** notebook -- see below) that scores a dataset (EMT status and/or
  gene-set activity) and produces figure panels.
- Has a `figures/` subdirectory with the PDF/SVG panels pulled into the
  paper.
- Has a `*_sessionInfo.txt` file recording the R session/package versions
  used the last time each notebook was run.
- Has a `*_Methods.md` file documenting the processing/scoring
  methodology and the statistics behind every figure, in prose.

3 of the 11 directories (`Baugh2024_GSE261252`, `BRCA_atlas_CellxGene_Chen`,
`BRCA_TCGA`) split prep from scoring/plotting into separate notebooks; the
rest (`TNBC_cellline_SUM149PT_GSE172609` and the 7 public-dataset
directories added 2026-09-22) combine both in one notebook, since each
reads an already-prepared checkpoint with no heavy import/integration
step of its own to isolate (see "Per-dataset status" below). See "Known
gaps" for the fix history and what's intentionally left as-is.

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
instructions written to the job's `slurm.*.out` log. Each `datasets/*`
directory has its own `.Rproj` so that relative paths inside its `.Rmd`
resolve correctly when opened as its own RStudio project (working
directory = that dataset directory); shared code is reached via that
directory's own `emt_code`/`Rcode` symlink (see "Layout" above), not a
`../` path, so it doesn't matter how deeply a dataset directory is
nested under `datasets/`.

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
| Baugh2024_GSE261252 | Fig1a, Fig1b, Fig1c | `Fig1a_EMT_KS_Score_FeaturePlot_chakroborty_RdBu.svg`, `Fig1b_EMT_KS_Score_byCondition_boxplot_chakroborty.svg`, `Fig1c_EMT_Stem_Epi_Dotplot_byCluster_res0.4.svg` |
| BRCA_atlas_CellxGene_Chen | Fig5b, ExtData5a–d | `Fig5b_ssgsea_BH4_Chen.svg`, `ExtData5a_ssgsea_Glutathione_Chen.svg`, `ExtData5b_ssgsea_CoQ10_Chen.svg`, `ExtData5c_ssgsea_Thioredoxin_Chen.svg`, `ExtData5d_ssgsea_Peroxiredoxin_Chen.svg` |
| BRCA_TCGA | Fig5c | `Fig5c_ssgsea_BH4_TCGA_BRCA.svg` |
| TNBC_cellline_SUM149PT_GSE172609 | ExtData5e, ExtData5f | `ExtData5e_bh4_boxplot_by_emtstate_dunnett.svg`, `ExtData5f_heatmap_bh4genes_lfc_sig.svg` |
| PDAC_mouse_35952360 | Fig5d, Fig5e, ExtData6a, ExtData6c | `5d_UMAP_BH4_score.pdf/.svg`, `5e_Boxplot_byCellState_BH4_score.pdf/.svg`, `Extended_6a_UMAP_byCellState_GeneExpression.pdf/.svg` (blocked -- see Known Gaps), `Extended_6c_Boxplot_byCellState_GeneExpression.pdf/.svg` |
| PDAC_human_38702773 | Fig5f | `5f_Boxplot_byCellState_BH4_score.pdf/.svg` |
| PDAC_human_35952360 | ExtData6d | `Extended_6d_Boxplot_byCellState_BH4_score.pdf/.svg` |
| SKCM_TCGA | Fig5g | `5g_Boxplot_byCellState_BH4_score.pdf/.svg` |
| Neuroblastoma_human_25150838_SEQC_GSE62564 | Fig5h | `5h_scatterplot_BH4vsMESADRN_score.pdf/.svg` |
| Neuroblastoma_human_40603900_Treehouse_GSE294351 | Fig5i | `5i_scatterplot_BH4vsMESADRN_score.pdf/.svg` |
| Neuroblastoma_cellline_30948783_GSE90803_GSE116890 | ExtData4a, ExtData4b, ExtData4c | `Extended_4ab_barplot_GCH1_GeneExpression.pdf/.svg`, `Extended_4c_lineplot_GCH1_GeneExpression.pdf/.svg` |

The last 7 rows were split out of a single monolithic script on
2026-09-22 (see "Known gaps") and haven't yet been re-knit in the
container to confirm the `.pdf`/`.svg` pair actually regenerates
byte-for-byte identically to the PDFs already on disk from the original
run -- see "Per-dataset status".

## Large files / data availability

Several `.rds` objects in this directory (pre-processed Seurat objects,
pseudobulk objects, scored objects) are large (tens of MB to several
hundred MB) and are **not tracked in git** (see `.gitignore`, blanket
`*.rds` rule) to keep this repository lightweight. Likewise untracked:
`datasets/TNBC_cellline_SUM149PT_GSE172609/from_geo/` (raw GEO download),
each of the original 4 datasets' own `data/` (Baugh's raw 10x
directories, Chen's raw h5ad, TCGA's Firehose matrix/metadata), `R/`
(locally installed packages — see "Running the analyses"), and every
rendered `*.nb.html` notebook preview (regenerated by re-knitting the
source `.Rmd`). The 7 newer public-dataset directories have no `data/`
subdirectory of their own -- their input `.rds` checkpoints sit at the
dataset directory's top level (matching `BRCA_TCGA`'s
`TCGA_BRCA_seurat.rds` convention) and are covered by the same blanket
`*.rds` rule.

**This means a fresh `git clone` of this repo is not, by itself, enough
to run any notebook** -- every `.rds` checkpoint and each dataset's raw
`data/`/`from_geo/` directory has to come from somewhere else. Two ways
to get them:

- **Fastest: `input_data.tar.gz`.** `pack_input_data.sh` collects every
  `.rds` checkpoint and raw `data/`/`from_geo/` directory across all 11
  `datasets/*` into `input_data/` (mirroring each file's exact
  `datasets/<name>/...` path) and tars it as `input_data.tar.gz`, hosted
  on Google Drive:
  <https://drive.google.com/file/d/1x3JLTMNisNe02R4xyrH1we21wqEEObnL/view?usp=sharing>
  (~5GB; GitHub itself never gets more than the two scripts below). After
  cloning, download that tarball, extract it at the repo root
  (`tar -xzf input_data.tar.gz`), and run `./restore_input_data.sh` -- it
  copies `input_data/datasets/` onto `datasets/`, merging into the
  existing tree. Re-run `pack_input_data.sh` whenever a dataset's data
  changes, to refresh the tarball for the next person (and re-upload it
  to the same Drive location, or wherever replaces it).
- **From scratch, per dataset.** Obtain that directory's raw input data
  (see each notebook's early "Load ..." chunks for what's expected, and
  where it should come from -- also each `*_Methods.md`'s "Data source"
  section), then run the directory's prep notebook (if it has one),
  then its scoring/plotting notebook, top to bottom inside the provided
  container. This is the only option for a dataset whose original hosting
  page has since gone dark (e.g. `BRCA_TCGA`'s Firehose matrix -- see
  "Known gaps") and provenance is tracked via citation instead.

## Per-dataset status

- **BRCA_atlas_CellxGene_Chen** is the model for how this should work.
  `ChenAtlas_malignant_pseudobulk_prep.Rmd` is not meant to be run
  routinely -- it's a one-shot notebook that downloads the raw source atlas
  (a CZ CELLxGENE h5ad, fetched via `download.file()` into this
  directory's `data/` on first run -- gitignored, since it's ~1-2GB and
  re-fetchable from the URL in the notebook), imports it, subsets to
  malignant cells, and aggregates to the donor-level pseudobulk checkpoint
  (`Atlas_malignant_pseudobulk_donor.rds`) that
  `ChenAtlas_malignant_pseudobulk.Rmd` then loads to do the actual EMT/
  gene-set scoring and figure generation. Nothing multi-GB is checkpointed
  along the way in the prep notebook -- intermediate Seurat objects are
  dropped once no longer needed. `*_Methods.md` is present and accurate.
  Not yet actually re-run end to end since being written, but per the
  user that verification isn't needed (see "Known gaps").

- **Baugh2024_GSE261252** has `baugh_BH4paper_prep.Rmd`, which
  consolidates five upstream notebooks (per-sample import/doublet-calling
  → integration → annotation → tumor subsetting → single-cell KS scoring)
  into one file, with each stage guarded by `if (file.exists(...))` so
  re-running only redoes stages whose checkpoint is missing. Its output
  filenames line up correctly with what `baugh_BH4paper_figures.Rmd`
  expects, and `doublet_called_rds/` (Stage 1's output) is now populated.
  `baugh_BH4paper_figures_Methods.md` documents both notebooks. Raw data
  lives at this directory's `data/` (GEO accession GSE261252, from Baugh
  et al. 2024 — see the prep notebook's "Raw data" section and `Known
  gaps` below), with `data/rename.sh` documenting the GSM-to-sample-
  directory file renaming needed before `Read10X()` can read it.

- **BRCA_TCGA** has `TCGA_BRCA_prep.Rmd`, which builds a Seurat object
  from the Broad Firehose batch-corrected pan-cancer matrix and clinical
  metadata and subsets it to BRCA; its save-filename bug and its
  undefined-`tcga`-object bug (see prior review) have both been fixed/
  removed by the user. `TCGA_BRCA_EMT_GeneSet_Scoring_Methods.md`
  documents both notebooks. The original GDC/Broad Firehose hosting pages
  for its raw inputs are dead; provenance is now tracked instead via a
  citation to Liu et al. 2018 (*Cell* 173(2), 400–416.e11), added to
  prep's new "Raw data" section — see "Fixed since the last review."

- **TNBC_cellline_SUM149PT_GSE172609** does prep and scoring/plotting in
  one notebook (`SUM149PT_EMT_scoring_plotting.Rmd`) rather than
  splitting them out. This is reasonable given how little prep work is
  involved: the notebook reads three already-processed GEO supplementary
  files directly (normalized TPM matrix, raw counts, gene annotation)
  with no heavy import/integration/subsetting step to isolate into its
  own notebook. `nf-core_fetchngs.sh`/`ids.csv`/`from_geo/results/` are
  kept intentionally, for anyone who wants to obtain and reprocess the
  raw reads instead of using the GEO-supplied processed matrices — not
  dead files to remove.

- **The 7 public-dataset directories** (`PDAC_mouse_35952360`,
  `PDAC_human_38702773`, `PDAC_human_35952360`, `SKCM_TCGA`,
  `Neuroblastoma_human_25150838_SEQC_GSE62564`,
  `Neuroblastoma_human_40603900_Treehouse_GSE294351`,
  `Neuroblastoma_cellline_30948783_GSE90803_GSE116890`) were split out of
  a single monolithic script on 2026-09-22 (each into its own combined
  scoring/plotting notebook, same rationale as SUM149PT above -- every
  input `.rds` was already a prepared checkpoint, not raw data needing
  its own prep stage). All of them have PDF figures already on disk from
  the original combined-script run, but **none have been re-knit as
  their new, split-out notebook** -- their `*_sessionInfo.txt` files are
  placeholders, not real output, since no R/container access was
  available during the split. Before relying on these for the paper:
  open each `.Rproj` in the container and re-knit top-to-bottom to
  confirm the split didn't introduce a regression, and to generate the
  paired `.svg` alongside the existing `.pdf`. See "Known gaps" for two
  correctness issues found in the source script during the split.

## Known gaps

This section tracks known issues to resolve before this directory is
considered publication/Zenodo-ready.

**Note on paths below:** entries in this section predate the 2026-09-22
`datasets/` reorg (Baugh → `Baugh2024_GSE261252`, Chen →
`BRCA_atlas_CellxGene_Chen`, TCGA → `BRCA_TCGA`, SUM149PT →
`TNBC_cellline_SUM149PT_GSE172609`, all moved under `datasets/`) and
reference the old top-level paths, which were correct at the time each
fix was made.

**2026-09-22 reorg:**

- **Directory rename + `datasets/` consolidation.** The 4 original
  dataset directories were renamed to a more informative
  `<disease/tissue>_<species>_<identifier>` convention and moved under
  `datasets/`; each now reaches shared `Rcode`/`emt_code` via a
  same-directory symlink instead of a `../`-relative path (see "Layout").
  7 more public datasets, previously sitting in an untracked
  `tmp/WH_BH4pathway_final/` as one monolithic script, were split into
  their own dataset directories under the same convention (see
  "Per-dataset status").
- **`data/rename.sh` negation rule never actually worked.**
  `.gitignore`'s `Baugh/data/` + `!Baugh/data/rename.sh` pattern silently
  failed to track `rename.sh` -- git can't re-include a file whose parent
  directory is excluded via a trailing-slash pattern. Despite the
  README's prior claim that this file was "kept as tracked provenance
  code," it was never actually committed. Fixed by ignoring the
  directory's *contents* (`datasets/*/data/*`) instead of the directory
  itself, which lets the negation apply; `rename.sh` is now genuinely
  tracked.
- **PDAC_mouse_35952360's Extended 6a panel references an undefined
  `p1`.** The source script's Data1 section called
  `plot_grid(p1, p3, rel_widths = c(1.5, 6))`, but `p1` was never defined
  anywhere in that section -- the only `p1` in the whole source script is
  an unrelated local variable in the (later-running) cell-line Data7
  section. Per the repo owner's decision, this was carried over verbatim
  (not silently patched) -- see that dataset's Methods.md. Will error on
  a fresh knit until `p1` is defined.
- **PDAC_human_35952360 (Lin) read a nonexistent
  `state_signatures.xlsx`.** No such file exists anywhere in the source
  script's exported bundle. Since Data2 (Park, `PDAC_human_38702773`)
  runs the identical Raghavan-rule classification code against
  `34890551_Raghavan_state_signatures.xlsx`, and both datasets share the
  same `ic_threshold`/signature columns, this was repointed at that same
  file per the repo owner's confirmation -- see that dataset's
  Methods.md.
- **`sessionInfo.txt` placeholders for all 7 newly split datasets.** No
  R/Singularity/Apptainer was available in the environment that
  performed the split, so none of the 7 new notebooks have actually been
  knit yet -- their `sessionInfo.txt` files are stub text explaining
  this, not real output. Real ones (and confirmation the split didn't
  break anything, and the paired `.svg` exports) require opening each
  `.Rproj` in the container and re-knitting.

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
  `data/`; `Baugh/data/rename.sh` was *intended* to be kept as tracked
  provenance code via a `!` negation, but that negation silently never
  worked until the 2026-09-22 reorg fixed it — see "2026-09-22 reorg"
  above.
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
