# =============================================================================
# 02_preprocessing_QC.R
#
# TFM: Transcriptomic response of human induced neurons to ionizing radiation
#
# Purpose:
#   Perform low-expression filtering, TMM normalization, voom transformation
#   and quality-control analyses for the 54-library RNA-seq dataset.
#
# Input:
#   data/processed/count_matrix_54.rds
#   data/processed/sample_metadata_54.csv
#
# Outputs:
#   data/processed/dge_54_filtered_TMM.rds
#   data/processed/voom_54_QC.rds
#   data/processed/logCPM_54.rds
#   results/QC_summary.csv
#   results/library_QC_metrics.csv
#   results/QC_flagged_libraries.csv
#   results/PCA_coordinates.csv
#   figures/QC/
# =============================================================================


# =============================================================================
# 1. REQUIRED PACKAGES
# =============================================================================

required <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "patchwork"
)

missing_packages <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(edgeR)
library(limma)
library(ggplot2)
library(patchwork)

# =============================================================================
# 2. OUTPUT DIRECTORIES
# =============================================================================

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/QC",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 3. LOAD INPUT DATA
# =============================================================================

counts <- readRDS(
  "data/processed/count_matrix_54.rds"
)

meta <- read.csv(
  "data/processed/sample_metadata_54.csv",
  stringsAsFactors = FALSE
)


# =============================================================================
# 4. DEFINE EXPERIMENTAL VARIABLES
# =============================================================================

meta$treatment <- factor(
  meta$treatment,
  levels = c("CTL", "IR")
)

meta$batch <- factor(meta$batch)

meta$donor <- factor(meta$donor)

meta$donor_batch <- interaction(
  meta$donor,
  meta$batch,
  drop = TRUE
)


# =============================================================================
# 5. INPUT CHECKS
# =============================================================================

stopifnot(
  nrow(counts) == 78428L,
  ncol(counts) == 54L,
  nrow(meta) == 54L,
  ncol(counts) == nrow(meta),
  identical(
    colnames(counts),
    meta$count_column
  ),
  sum(meta$treatment == "CTL") == 27L,
  sum(meta$treatment == "IR") == 27L,
  nlevels(meta$donor) == 17L,
  nlevels(meta$batch) == 10L,
  nlevels(meta$donor_batch) == 27L
)


# =============================================================================
# 6. LIBRARY-LEVEL QC BEFORE FILTERING
# =============================================================================

library_size <- colSums(counts)

detected_genes <- colSums(
  counts > 0
)

library_qc <- data.frame(
  GSM = meta$GSM,
  donor = meta$donor,
  batch = meta$batch,
  treatment = meta$treatment,
  library_size = library_size,
  library_size_millions = library_size / 1e6,
  detected_genes = detected_genes,
  stringsAsFactors = FALSE
)


# =============================================================================
# 7. EXPLORATORY 1.5 x IQR OUTLIER CRITERION
# =============================================================================

iqr_limits <- function(x) {

  q1 <- unname(
    quantile(
      x,
      0.25,
      na.rm = TRUE
    )
  )

  q3 <- unname(
    quantile(
      x,
      0.75,
      na.rm = TRUE
    )
  )

  iqr <- q3 - q1

  c(
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr
  )
}


lib_limits <- iqr_limits(
  library_qc$library_size
)

gene_limits <- iqr_limits(
  library_qc$detected_genes
)


library_qc$outlier_library_size <- (
  library_qc$library_size < lib_limits[["lower"]] |
    library_qc$library_size > lib_limits[["upper"]]
)

library_qc$outlier_detected_genes <- (
  library_qc$detected_genes < gene_limits[["lower"]] |
    library_qc$detected_genes > gene_limits[["upper"]]
)

library_qc$outlier_any <- (
  library_qc$outlier_library_size |
    library_qc$outlier_detected_genes
)


flagged <- library_qc[
  which(library_qc$outlier_any),
  c(
    "GSM",
    "donor",
    "batch",
    "treatment",
    "library_size_millions",
    "detected_genes",
    "outlier_library_size",
    "outlier_detected_genes"
  ),
  drop = FALSE
]


write.csv(
  library_qc,
  "results/library_QC_metrics.csv",
  row.names = FALSE
)

write.csv(
  flagged,
  "results/QC_flagged_libraries.csv",
  row.names = FALSE
)


# =============================================================================
# 8. LOW-EXPRESSION FILTERING
# =============================================================================
#
# The donor-batch pairing structure is used here to guide filterByExpr().
# The final differential-expression model is fitted separately in script 03.
# =============================================================================

filter_design <- model.matrix(
  ~ donor_batch + treatment,
  data = meta
)

dge <- DGEList(
  counts = counts,
  samples = meta
)

keep <- filterByExpr(
  dge,
  design = filter_design
)

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]


# =============================================================================
# 9. TMM NORMALIZATION
# =============================================================================

dge <- calcNormFactors(
  dge,
  method = "TMM"
)


# =============================================================================
# 10. FILTERING CHECKPOINT
# =============================================================================

if (nrow(dge) != 22439L) {

  warning(
    "Filtered gene count is ",
    nrow(dge),
    "; the final TFM reports 22,439. ",
    "Check package versions and sample mapping."
  )
}


# =============================================================================
# 11. DESIGN MATRIX FOR VOOM QC
# =============================================================================
#
# The PCA used in the final TFM was calculated from the voom-transformed
# expression matrix (v$E), rather than directly from log2-CPM.
#
# Batch and treatment are included in the design used to estimate the
# mean-variance relationship. No expression values are removed or residualized
# for PCA; PCA is performed directly on v$E.
#
# Donor correlation is estimated later in 03_differential_expression.R.
# =============================================================================

design_qc <- model.matrix(
  ~ batch + treatment,
  data = meta
)


# =============================================================================
# 12. VOOM TRANSFORMATION
# =============================================================================

v <- voom(
  dge,
  design = design_qc,
  plot = FALSE
)


# =============================================================================
# 13. LOG2-CPM OBJECT
# =============================================================================
#
# Retained as a useful QC object and for compatibility with downstream
# exploratory analyses. PCA below is intentionally based on v$E.
# =============================================================================

log_cpm <- cpm(
  dge,
  log = TRUE,
  prior.count = 0.5
)


# =============================================================================
# 14. PRINCIPAL COMPONENT ANALYSIS ON VOOM EXPRESSION
# =============================================================================

pca <- prcomp(
  t(v$E),
  scale. = FALSE
)

var_exp <- 100 * (
  pca$sdev^2 /
    sum(pca$sdev^2)
)


# Reproducibility checkpoint from the final TFM:
# PC1 ~ 11.9 %
# PC2 ~  9.8 %

if (
  round(var_exp[1], 1) != 11.9 ||
  round(var_exp[2], 1) != 9.8
) {

  warning(
    sprintf(
      paste0(
        "PCA variance differs from the final TFM: ",
        "PC1 = %.3f%%, PC2 = %.3f%%."
      ),
      var_exp[1],
      var_exp[2]
    )
  )
}


# =============================================================================
# 15. PCA DATA FRAME
# =============================================================================

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  treatment = meta$treatment,
  batch = meta$batch,
  donor = meta$donor,
  GSM = meta$GSM,
  stringsAsFactors = FALSE
)

pca_df$QC_flagged <- (
  pca_df$GSM %in% flagged$GSM
)


write.csv(
  pca_df,
  "results/PCA_coordinates.csv",
  row.names = FALSE
)


# =============================================================================
# 16. PCA BY TREATMENT
# =============================================================================

p_treatment <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    fill = treatment
  )
) +

  geom_point(
    shape = 21,
    size = 3,
    alpha = 0.75
  ) +

  labs(
    title = "Principal component analysis",
    subtitle = "Human induced neurons, 336 h",
    x = sprintf(
      "PC1 (%.1f%%)",
      var_exp[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_exp[2]
    ),
    fill = "Condition"
  ) +

  theme_classic(
    base_size = 12
  )


print(p_treatment)


ggsave(
  "figures/QC/PCA_by_treatment.png",
  plot = p_treatment,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/QC/PCA_by_treatment.pdf",
  plot = p_treatment,
  width = 7,
  height = 5
)


# =============================================================================
# 17. PCA BY BATCH
# =============================================================================

p_batch <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = batch
  )
) +

  geom_point(
    size = 3,
    alpha = 0.8
  ) +

  labs(
    title = "PCA by experimental batch",
    subtitle = "Human induced neurons, 336 h",
    x = sprintf(
      "PC1 (%.1f%%)",
      var_exp[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_exp[2]
    ),
    colour = "Batch"
  ) +

  theme_classic(
    base_size = 12
  )


print(p_batch)


ggsave(
  "figures/QC/PCA_by_batch.png",
  plot = p_batch,
  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# 18. PCA BY DONOR
# =============================================================================

p_donor <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = donor
  )
) +

  geom_point(
    size = 3,
    alpha = 0.8
  ) +

  labs(
    title = "PCA by donor",
    subtitle = "Human induced neurons, 336 h",
    x = sprintf(
      "PC1 (%.1f%%)",
      var_exp[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_exp[2]
    ),
    colour = "Donor"
  ) +

  theme_classic(
    base_size = 12
  )


print(p_donor)


ggsave(
  "figures/QC/PCA_by_donor.png",
  plot = p_donor,
  width = 9,
  height = 6,
  dpi = 300
)


# =============================================================================
# 19. PCA WITH QC-FLAGGED LIBRARIES
# =============================================================================

p_flagged <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2
  )
) +

  geom_point(
    aes(
      fill = treatment
    ),
    shape = 21,
    size = 3,
    alpha = 0.7
  ) +

  geom_point(
    data = pca_df[
      pca_df$QC_flagged,
      ,
      drop = FALSE
    ],
    shape = 21,
    size = 5,
    stroke = 1.3
  ) +

  geom_text(
    data = pca_df[
      pca_df$QC_flagged,
      ,
      drop = FALSE
    ],
    aes(
      label = GSM
    ),
    nudge_y = 5,
    size = 3
  ) +

  labs(
    title = "PCA and libraries flagged by exploratory QC",
    subtitle = "1.5 x IQR criterion",
    x = sprintf(
      "PC1 (%.1f%%)",
      var_exp[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      var_exp[2]
    ),
    fill = "Condition"
  ) +

  theme_classic(
    base_size = 12
  )


print(p_flagged)


ggsave(
  "figures/QC/PCA_QC_flagged_libraries.png",
  plot = p_flagged,
  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# 20. LIBRARY-SIZE QC FIGURE
# =============================================================================

p_library_size <- ggplot(
  library_qc,
  aes(
    x = reorder(
      GSM,
      library_size_millions
    ),
    y = library_size_millions,
    fill = treatment
  )
) +

  geom_col(
    width = 0.8
  ) +

  geom_hline(
    yintercept = lib_limits[["lower"]] / 1e6,
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept = lib_limits[["upper"]] / 1e6,
    linetype = "dashed"
  ) +

  labs(
    title = "RNA-seq library sizes",
    subtitle = "Dashed lines indicate 1.5 x IQR limits",
    x = "Library",
    y = "Library size (million counts)",
    fill = "Condition"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    )
  )


print(p_library_size)


ggsave(
  "figures/QC/library_size_QC.png",
  plot = p_library_size,
  width = 11,
  height = 6,
  dpi = 300
)


# =============================================================================
# 21. DETECTED-GENE QC FIGURE
# =============================================================================

p_detected <- ggplot(
  library_qc,
  aes(
    x = reorder(
      GSM,
      detected_genes
    ),
    y = detected_genes,
    fill = treatment
  )
) +

  geom_col(
    width = 0.8
  ) +

  geom_hline(
    yintercept = gene_limits[["lower"]],
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept = gene_limits[["upper"]],
    linetype = "dashed"
  ) +

  labs(
    title = "Detected genes per RNA-seq library",
    subtitle = "Dashed lines indicate 1.5 x IQR limits",
    x = "Library",
    y = "Detected genes",
    fill = "Condition"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    )
  )


print(p_detected)


ggsave(
  "figures/QC/detected_genes_QC.png",
  plot = p_detected,
  width = 11,
  height = 6,
  dpi = 300
)


# =============================================================================
# 22. VOOM-NORMALIZED EXPRESSION DISTRIBUTIONS
# =============================================================================

expression_df <- data.frame(
  expression = as.vector(v$E),
  library = rep(
    colnames(v$E),
    each = nrow(v$E)
  )
)


p_expression <- ggplot(
  expression_df,
  aes(
    x = library,
    y = expression
  )
) +

  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.25
  ) +

  labs(
    title = "Expression distributions after TMM normalization and voom",
    x = "RNA-seq library",
    y = "Voom log2 expression"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


print(p_expression)


ggsave(
  "figures/QC/expression_distribution_TMM_voom.png",
  plot = p_expression,
  width = 10,
  height = 6,
  dpi = 300
)
# =============================================================================
# 23. COMPOSITE QC FIGURE
# =============================================================================
#
# Composite representation corresponding to the QC figure used in the TFM:
#   A) Library size
#   B) Number of detected genes
#   C) PCA highlighting libraries flagged by the exploratory 1.5 x IQR rule
#
# The PCA is based on the voom-transformed expression matrix (v$E).
# =============================================================================


# -----------------------------------------------------------------------------
# Prepare ordered library-level QC data
# -----------------------------------------------------------------------------

qc_plot_df <- library_qc[
  order(library_qc$library_size_millions),
  ,
  drop = FALSE
]

qc_plot_df$library_order_size <- seq_len(
  nrow(qc_plot_df)
)


genes_plot_df <- library_qc[
  order(library_qc$detected_genes),
  ,
  drop = FALSE
]

genes_plot_df$library_order_genes <- seq_len(
  nrow(genes_plot_df)
)


# -----------------------------------------------------------------------------
# Panel A: library size
# -----------------------------------------------------------------------------

p_qc_A <- ggplot(
  qc_plot_df,
  aes(
    x = library_order_size,
    y = library_size_millions
  )
) +

  geom_point(
    aes(
      colour = ifelse(
        outlier_any,
        "Flagged",
        "Other"
      )
    ),
    size = 2.2
  ) +

  geom_text(
    data = qc_plot_df[
      qc_plot_df$outlier_any,
      ,
      drop = FALSE
    ],
    aes(
      label = GSM
    ),
    nudge_y = 0.45,
    size = 2.5,
    check_overlap = TRUE
  ) +

  scale_colour_manual(
    values = c(
      "Other" = "#F8766D",
      "Flagged" = "#00BFC4"
    ),
    labels = c(
      "Other" = "Resto",
      "Flagged" = "Señalada"
    )
  ) +

  labs(
    title = "Tamaño de biblioteca",
    x = "Bibliotecas ordenadas",
    y = "Conteos totales (millones)",
    colour = NULL
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      face = "bold"
    )
  )


# -----------------------------------------------------------------------------
# Panel B: detected genes
# -----------------------------------------------------------------------------

p_qc_B <- ggplot(
  genes_plot_df,
  aes(
    x = library_order_genes,
    y = detected_genes
  )
) +

  geom_point(
    aes(
      colour = ifelse(
        outlier_any,
        "Flagged",
        "Other"
      )
    ),
    size = 2.2
  ) +

  geom_text(
    data = genes_plot_df[
      genes_plot_df$outlier_any,
      ,
      drop = FALSE
    ],
    aes(
      label = GSM
    ),
    nudge_y = 120,
    size = 2.5,
    check_overlap = TRUE
  ) +

  scale_colour_manual(
    values = c(
      "Other" = "#F8766D",
      "Flagged" = "#00BFC4"
    ),
    labels = c(
      "Other" = "Resto",
      "Flagged" = "Señalada"
    )
  ) +

  labs(
    title = "Genes detectados por biblioteca",
    x = "Bibliotecas ordenadas",
    y = "Número de genes detectados",
    colour = NULL
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(
      face = "bold"
    )
  )


# -----------------------------------------------------------------------------
# Panel C: PCA highlighting QC-flagged libraries
# -----------------------------------------------------------------------------

p_qc_C <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2
  )
) +

  geom_point(
    colour = "grey60",
    size = 2.1,
    alpha = 0.9
  ) +

  geom_point(
    data = pca_df[
      pca_df$QC_flagged,
      ,
      drop = FALSE
    ],
    shape = 21,
    fill = "white",
    colour = "black",
    size = 4,
    stroke = 1
  ) +

  geom_text(
    data = pca_df[
      pca_df$QC_flagged,
      ,
      drop = FALSE
    ],
    aes(
      label = GSM
    ),
    nudge_x = -4,
    nudge_y = 5,
    size = 2.6,
    check_overlap = FALSE
  ) +

  labs(
    title = "Evaluación de posibles valores atípicos",
    x = sprintf(
      "PC1 (%.1f %%)",
      var_exp[1]
    ),
    y = sprintf(
      "PC2 (%.1f %%)",
      var_exp[2]
    )
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


# -----------------------------------------------------------------------------
# Combine panels
# -----------------------------------------------------------------------------

p_qc_composite <- (
  p_qc_A | p_qc_B
) / p_qc_C +

  plot_layout(
    heights = c(
      1,
      1.15
    )
  ) +

  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 12
      )
    )
  )


# Display in RStudio
print(p_qc_composite)


# -----------------------------------------------------------------------------
# Save composite QC figure
# -----------------------------------------------------------------------------

ggsave(
  "figures/QC/QC_composite_TFM.png",
  plot = p_qc_composite,
  width = 13,
  height = 7,
  dpi = 300
)

ggsave(
  "figures/QC/QC_composite_TFM.pdf",
  plot = p_qc_composite,
  width = 13,
  height = 7
)

# =============================================================================
# 24. QC SUMMARY TABLE
# =============================================================================

qc_summary <- data.frame(

  metric = c(
    "initial_genes",
    "filtered_genes",
    "libraries",
    "CTL",
    "IR",
    "donors",
    "batches",
    "donor_batch_pairs",
    "PC1_variance_percent",
    "PC2_variance_percent",
    "median_library_size_millions",
    "median_detected_genes",
    "flagged_libraries"
  ),

  value = c(
    nrow(counts),
    nrow(dge),
    ncol(dge),
    sum(meta$treatment == "CTL"),
    sum(meta$treatment == "IR"),
    nlevels(meta$donor),
    nlevels(meta$batch),
    nlevels(meta$donor_batch),
    var_exp[1],
    var_exp[2],
    median(library_qc$library_size_millions),
    median(library_qc$detected_genes),
    nrow(flagged)
  )
)


write.csv(
  qc_summary,
  "results/QC_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 25. SAVE PROCESSED OBJECTS
# =============================================================================

saveRDS(
  dge,
  "data/processed/dge_54_filtered_TMM.rds"
)

saveRDS(
  v,
  "data/processed/voom_54_QC.rds"
)

saveRDS(
  log_cpm,
  "data/processed/logCPM_54.rds"
)


# =============================================================================
# 26. SESSION INFORMATION
# =============================================================================

capture.output(
  sessionInfo(),
  file = "results/sessionInfo_QC.txt"
)


# =============================================================================
# 27. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("TFM PREPROCESSING AND QC SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Initial genes:",
  nrow(counts),
  "\n"
)

cat(
  "Genes after filterByExpr:",
  nrow(dge),
  "\n"
)

cat(
  "Libraries:",
  ncol(dge),
  "\n"
)

cat(
  "CTL:",
  sum(meta$treatment == "CTL"),
  "| IR:",
  sum(meta$treatment == "IR"),
  "\n"
)

cat(
  "Donors:",
  nlevels(meta$donor),
  "| Batches:",
  nlevels(meta$batch),
  "| Donor-batch pairs:",
  nlevels(meta$donor_batch),
  "\n\n"
)


cat(
  sprintf(
    "PCA variance (voom): PC1 %.3f%% | PC2 %.3f%%\n",
    var_exp[1],
    var_exp[2]
  )
)

cat(
  sprintf(
    "PCA variance (rounded): PC1 %.1f%% | PC2 %.1f%%\n\n",
    var_exp[1],
    var_exp[2]
  )
)


cat("Library-size QC\n")

cat(
  sprintf(
    "Range: %.2f - %.2f million counts\n",
    min(library_qc$library_size_millions),
    max(library_qc$library_size_millions)
  )
)

cat(
  sprintf(
    "Median: %.2f million counts\n\n",
    median(library_qc$library_size_millions)
  )
)


cat("Detected genes per library\n")

cat(
  "Range:",
  min(library_qc$detected_genes),
  "-",
  max(library_qc$detected_genes),
  "\n"
)

cat(
  "Median:",
  median(library_qc$detected_genes),
  "\n\n"
)


cat("1.5 x IQR limits\n")

cat(
  sprintf(
    "Library size: %.2f - %.2f million counts\n",
    lib_limits[["lower"]] / 1e6,
    lib_limits[["upper"]] / 1e6
  )
)

cat(
  sprintf(
    "Detected genes: %.0f - %.0f\n\n",
    gene_limits[["lower"]],
    gene_limits[["upper"]]
  )
)


cat(
  "Libraries flagged by 1.5 x IQR criterion:",
  nrow(flagged),
  "\n"
)

print(flagged)


cat("\nExpected TFM PCA checkpoint:\n")
cat("  PC1: 11.9%\n")
cat("  PC2:  9.8%\n")


cat("\nFigures saved in:\n")
cat("  figures/QC/\n")


cat("\n============================================================\n")
cat("QC COMPLETE\n")
cat("============================================================\n")
