# Shared setup for the 7 BH4-pathway public-dataset notebooks
# (datasets/PDAC_mouse_35952360, PDAC_human_38702773, PDAC_human_35952360,
# SKCM_TCGA, Neuroblastoma_human_25150838_SEQC_GSE62564,
# Neuroblastoma_human_40603900_Treehouse_GSE294351,
# Neuroblastoma_cellline_30948783_GSE90803_GSE116890).
#
# Extracted verbatim from the original monolithic
# tmp/WH_BH4pathway_final/src/Analysis_script.Rmd (lines 9-216) during the
# datasets/ reorg, so all 7 notebooks share one copy instead of inlining it
# 7 times. Each notebook does: source("emt_code/bh4_pathway_shared_setup.R")

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(rstatix)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(Matrix)
  library(UCell)
  library(GSVA)
  library(cowplot)
  library(readxl)
})

# figure-theme for ggplot
figure_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 6, margin = margin(b = 2)),
    plot.subtitle = element_text(hjust = 0.5, size = 5, margin = margin(b = 2)),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 5, margin = margin(l = 0, unit = "pt")),
    legend.title = element_text(size = 5, margin = margin(b = 1, unit = "pt"), face = "bold"),
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.margin = margin(0.2, 0, 0.2, 0),
    plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")
  )
panel_height <- 1.2
calc_fig_height <- function(n_panels, ncol, panel_h = panel_height) {
  ceiling(n_panels / ncol) * panel_h + 1
}

# ComplexHeatmap defaults (raster backend set for a Cairo-free environment)
ht_defaults <- list(
  show_column_dend = TRUE, column_dend_side = "bottom",
  column_names_side = c("top"), column_names_rot = 45,
  row_title_gp = gpar(fontsize = 5), row_title_rot = 0,
  column_title_gp = gpar(fontsize = 5), column_title_rot = 0,
  row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(5, "mm"), row_dend_gp = gpar(lwd = 0.5),
  column_dend_height = unit(3, "mm"), column_dend_gp = gpar(lwd = 0.5),
  row_gap = unit(0.2, "mm"), column_gap = unit(0.2, "mm"),
  border_gp = gpar(col = "black", lwd = 0.3),
  use_raster = TRUE, raster_quality = 4, raster_device = "png", raster_by_magick = FALSE,
  heatmap_legend_param = list(
    direction = "horizontal", legend_width = unit(1, "cm"),
    grid_height = unit(2, "mm"), tick_length = unit(0.5, "mm"),
    title_gp = gpar(fontsize = 5, fontface = "bold"), labels_gp = gpar(fontsize = 5))
)
make_legend_params <- function(df, cat_params, cont_params) {
  setNames(lapply(colnames(df), function(col) {
    if (is.numeric(df[[col]])) cont_params else cat_params
  }), colnames(df))
}
legend_params_cat  <- list(title_gp = gpar(fontsize = 5, fontface = "bold"), labels_gp = gpar(fontsize = 5),
                           grid_height = unit(2, "mm"), grid_width = unit(2, "mm"), ncol = 2)
legend_params_cont <- list(title_gp = gpar(fontsize = 5, fontface = "bold"), labels_gp = gpar(fontsize = 5),
                           grid_width = unit(2, "mm"), tick_length = unit(0.5, "mm"),
                           direction = "vertical", legend_height = unit(1.2, "cm"))

save_fig <- function(p, name, width_in, height_in) {
  ggsave(file.path(current.outputfolder, paste0(name, ".pdf")), p,
         width = width_in, height = height_in, units = "in", bg = "white")
  ggsave(file.path(current.outputfolder, paste0(name, ".svg")), p,
         width = width_in, height = height_in, units = "in", bg = "white")
  cat("  Saved:", name, ".pdf /", name, ".svg\n")
  print(p)
}

make_boxplot_state <- function(box_long, quantities, title_str, out_pdf, levels, fills = state_color, ncol = 3) {
  make_single <- function(q) {
    gdf <- box_long %>% filter(Quantity == q, !is.na(cancer_state)) %>%
      mutate(cancer_state = factor(cancer_state, levels = levels))
    y_lab <- unique(gdf$y_label)
    # Stats only on states with >= 2 cells (drop empty/singleton groups so
    # wilcox.test never receives a 0-observation group).
    grp_n     <- table(gdf$cancer_state)
    test_levs <- names(grp_n)[grp_n >= 2]
    has_test  <- length(test_levs) > 1
    gdf_stat  <- gdf %>% filter(cancer_state %in% test_levs) %>%
      mutate(cancer_state = droplevels(factor(cancer_state, levels = levels)))
    kw <- if (has_test) kruskal.test(value ~ cancer_state, data = gdf_stat) else list(p.value = NA_real_)
    kw_label <- if (!is.finite(kw$p.value)) "KW n/a" else if (kw$p.value < 0.001) sprintf("KW p = %.2e", kw$p.value) else sprintf("KW p = %.3f", kw$p.value)

    wstats <- tapply(gdf$value, gdf$cancer_state, function(v) boxplot.stats(v)$stats)
    wstats <- wstats[!vapply(wstats, is.null, logical(1))]
    wu <- max(vapply(wstats, `[`, numeric(1), 5), na.rm = TRUE)
    wl <- min(vapply(wstats, `[`, numeric(1), 1), na.rm = TRUE)
    yr <- wu - wl; if (!is.finite(yr) || yr == 0) yr <- 1

    # x-axis labels: state name + per-state cell count
    short_state <- function(v) { v[v == "Mesenchymal_I"] <- "Mes_I"; v[v == "Mesenchymal_II"] <- "Mes_II"; v }
    n_by   <- table(factor(gdf$cancer_state, levels = levels))
    x_labs <- setNames(sprintf("%s\n(n=%s)", short_state(levels), format(as.integer(n_by[levels]), big.mark = ",")), levels)
    p <- ggplot(gdf) +
      geom_boxplot(aes(x = cancer_state, y = value, fill = cancer_state),
                   width = 0.6, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.4) +
      scale_fill_manual(values = fills) +
      scale_x_discrete(labels = x_labs) +
      labs(title = paste0(q, "\n", kw_label), x = NULL, y = y_lab) +
      figure_theme + theme(legend.position = "none")

    y_top <- wu
    if (has_test && is.finite(kw$p.value) && kw$p.value < 0.05) {
      pw <- gdf_stat %>% rstatix::wilcox_test(value ~ cancer_state, p.adjust.method = "BH")
      if (!"p.adj" %in% names(pw)) pw$p.adj <- p.adjust(pw$p, method = "BH")
      pw <- pw[is.finite(pw$p.adj) & pw$p.adj < 0.05, , drop = FALSE]   # label only significant comparisons
      if (nrow(pw) > 0) {
        span <- abs(match(pw$group2, levels) - match(pw$group1, levels))
        pw <- pw[order(span), ]
        step <- 0.14 * yr
        pw$y.position <- wu + step * seq_len(nrow(pw))
        pw$p.label <- ifelse(pw$p.adj < 0.001, sprintf("%.1e", pw$p.adj), sprintf("%.3f", pw$p.adj))
        p <- p + ggpubr::stat_pvalue_manual(pw, label = "p.label", tip.length = 0.01,
               bracket.size = 0.3, label.size = 1.6)
        y_top <- max(pw$y.position) + step
      }
    }
    p + coord_cartesian(ylim = c(wl - 0.05 * yr, y_top + 0.05 * yr), clip = "off")
  }
  panels <- lapply(quantities, make_single)
  p <- wrap_plots(panels, ncol = ncol)
  legend_df <- data.frame(cancer_state = factor(levels, levels = levels), x = seq_along(levels), y = 0)
  legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = cancer_state)) +
    geom_col(width = 0.1) + scale_fill_manual(values = fills) +
    theme_void() + theme(legend.position = "bottom", legend.title = element_blank(),
                         legend.text = element_text(size = 5))
  final_p <- legend_plot / p + plot_layout(heights = c(0.05, 1)) +
    plot_annotation(title = title_str,
      theme = theme(plot.title = element_text(face = "bold", size = 6, hjust = 0.5)))
  ggsave(out_pdf, final_p, width = 6, height = calc_fig_height(length(quantities), ncol),
         bg = "white", limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".svg", out_pdf), final_p, width = 6, height = calc_fig_height(length(quantities), ncol),
         bg = "white", limitsize = FALSE)
  cat("  Saved:", basename(out_pdf), "/", basename(sub("\\.pdf$", ".svg", out_pdf)), "\n")
  print(final_p)
}

# read excel sheet
read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

# ── Helper: extract gene expression from Affymetrix raw matrix ───────────────
extract_affy_genes <- function(raw_matrix, probeset_map, gene_probeset_map, sample_meta) {
  rows <- list()
  for (gene in names(gene_probeset_map)) {
    ps <- gene_probeset_map[gene]
    if (ps %in% rownames(raw_matrix)) {
      for (i in 1:nrow(sample_meta)) {
        gsm <- sample_meta$GSM[i]
        if (gsm %in% colnames(raw_matrix)) {
          rows[[length(rows) + 1]] <- data.frame(
            Gene = gene,
            Probeset = ps,
            GSM = gsm,
            MAS5_log2 = as.numeric(raw_matrix[ps, gsm]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

# ── BH4 signatures (human) — 3 variants, names match the 3CA notebook ──
bh4_sigs <- list(
  BH4_denovo_synthesis                        = c("GCH1", "PTS", "SPR"),
  BH4_recycling_regeneration                  = c("QDPR", "DHFR", "PCBD1", "PCBD2"),
  BH4_denovo_synthesis_recycling_regeneration = c("GCH1", "PTS", "SPR","QDPR", "DHFR", "PCBD1", "PCBD2")
)
bh4_score_names <- names(bh4_sigs)
bh4_full        <- "BH4_denovo_synthesis_recycling_regeneration"    # composite (full pathway)
bh4_genes       <- c("GCH1", "PTS", "SPR", "QDPR", "PCBD1", "PCBD2", "DHFR")  # gene-panel order

bh4_sigs_mouse <- list(
  BH4_denovo_synthesis                        = c("Gch1", "Pts", "Spr"),
  BH4_recycling_regeneration                  = c("Qdpr", "Dhfr", "Pcbd1", "Pcbd2"),
  BH4_denovo_synthesis_recycling_regeneration = c("Gch1", "Pts", "Spr", "Qdpr", "Dhfr", "Pcbd1", "Pcbd2")
)
bh4_genes_mouse       <- c("Gch1", "Pts", "Spr", "Qdpr", "Pcbd1", "Pcbd2", "Dhfr")  # gene-panel order
