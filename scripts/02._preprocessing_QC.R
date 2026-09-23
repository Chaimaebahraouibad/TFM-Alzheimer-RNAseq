# =============================================================================
# 02_preprocessing_QC.R
# TFM: Transcriptomic response of human induced neurons to ionizing radiation
#
# Purpose:
#   - Low-expression filtering using filterByExpr
#   - TMM normalization
#   - PCA quality control
#   - Library-size and detected-gene QC
#   - Exploratory 1.5 x IQR outlier screening
#   - PCA visualization by treatment, batch and donor
#   - PCA visualization of QC-flagged libraries
#
# Input:
#   data/processed/count_matrix_54.rds
#   data/processed/sample_metadata_54.csv
#
# Output:
#   results/
#   figures/
#   data/processed/dge_54_filtered_TMM.rds
#   data/processed/logCPM_54.rds
# =============================================================================


# =============================================================================
# 0. PACKAGES AND DIRECTORIES
# =============================================================================

required <- c("edgeR", "ggplot2")

missing_packages <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(edgeR)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 1. LOAD THE 54-LIBRARY ANALYSIS SET
# =============================================================================

counts <- readRDS(
  "data/processed/count_matrix_54.rds"
)

meta <- read.csv(
  "data/processed/sample_metadata_54.csv",
  stringsAsFactors = FALSE
)

meta$treatment <- factor(
  meta$treatment,
  levels = c("CTL", "IR")
)

meta$batch <- factor(meta$batch)
meta$donor <- factor(meta$donor)

stopifnot(
  nrow(counts) == 78428L,
  ncol(counts) == 54L,
  nrow(meta) == 54L,
  ncol(counts) == nrow(meta),
  all(colnames(counts) == meta$count_column),
  sum(meta$treatment == "CTL") == 27L,
  sum(meta$treatment == "IR") == 27L,
  nlevels(meta$donor) == 17L,
  nlevels(meta$batch) == 10L
)


# =============================================================================
# 2. DEFINE DONOR-BATCH PAIRS
# =============================================================================

# The 54 libraries correspond to 27 donor-batch combinations.
# Each donor-batch combination contains one CTL and one IR library.

meta$donor_batch <- interaction(
  meta$donor,
  meta$batch,
  drop = TRUE
)

pair_table <- table(
  meta$donor_batch,
  meta$treatment
)

stopifnot(
  nlevels(meta$donor_batch) == 27L,
  all(pair_table == 1L)
)


# =============================================================================
# 3. LOW-EXPRESSION FILTERING
# =============================================================================

# The experimental design is considered during filterByExpr.
#
# donor_batch is used here only to guide the low-expression filtering.
# The final differential-expression model is implemented separately in
# 03_differential_expression.R.

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

# TFM checkpoint:
# 78,428 initial gene IDs -> 22,439 genes after filterByExpr.

if (nrow(dge) != 22439L) {
  warning(
    "Filtered gene count is ",
    nrow(dge),
    "; the final TFM reports 22,439. ",
    "Check package versions and sample mapping before proceeding."
  )
}


# =============================================================================
# 4. TMM NORMALIZATION
# =============================================================================

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

saveRDS(
  dge,
  "data/processed/dge_54_filtered_TMM.rds"
)


# =============================================================================
# 5. LOG2-CPM FOR EXPLORATORY QC
# =============================================================================

# log2-CPM values are used here for exploratory PCA visualization.

log_cpm <- cpm(
  dge,
  log = TRUE,
  prior.count = 0.5
)

saveRDS(
  log_cpm,
  "data/processed/logCPM_54.rds"
)


# =============================================================================
# 6. PRINCIPAL COMPONENT ANALYSIS
# =============================================================================

pca <- prcomp(
  t(log_cpm),
  scale. = FALSE
)

var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  GSM = meta$GSM,
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  treatment = meta$treatment,
  batch = meta$batch,
  donor = meta$donor,
  stringsAsFactors = FALSE
)

write.csv(
  pca_df,
  "results/PCA_coordinates.csv",
  row.names = FALSE
)


# =============================================================================
# 7. PCA BY EXPERIMENTAL CONDITION
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

ggsave(
  "figures/PCA_by_treatment.png",
  p_treatment,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/PCA_by_treatment.pdf",
  p_treatment,
  width = 7,
  height = 5
)


# =============================================================================
# 8. PCA BY EXPERIMENTAL BATCH
# =============================================================================

p_batch <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    fill = batch
  )
) +
  geom_point(
    shape = 21,
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
    fill = "Batch"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  "figures/PCA_by_batch.png",
  p_batch,
  width = 7,
  height = 5,
  dpi = 300
)


# =============================================================================
# 9. PCA BY DONOR
# =============================================================================

p_donor <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    fill = donor
  )
) +
  geom_point(
    shape = 21,
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
    fill = "Donor"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  "figures/PCA_by_donor.png",
  p_donor,
  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# 10. LIBRARY-SIZE AND DETECTED-GENE QC
# =============================================================================

# These metrics are calculated from the original 54-library count matrix
# before low-expression filtering.

library_size <- colSums(counts)

detected_genes <- colSums(
  counts > 0
)

library_qc <- data.frame(
  GSM = meta$GSM,
  donor = meta$donor,
  batch = meta$batch,
  treatment = meta$treatment,
  library_size = as.numeric(library_size),
  library_size_millions = as.numeric(library_size) / 1e6,
  detected_genes = as.numeric(detected_genes),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(library_qc) == 54L,
  !anyNA(library_qc$library_size),
  !anyNA(library_qc$detected_genes)
)


# =============================================================================
# 11. EXPLORATORY 1.5 x IQR OUTLIER SCREENING
# =============================================================================

# Return explicitly named numeric limits.
# names = FALSE prevents quantile names such as "25%" and "75%" from being
# propagated into the resulting vector.

iqr_limits <- function(x) {

  q <- quantile(
    x,
    probs = c(0.25, 0.75),
    na.rm = TRUE,
    names = FALSE
  )

  iqr_value <- q[2] - q[1]

  c(
    lower = q[1] - 1.5 * iqr_value,
    upper = q[2] + 1.5 * iqr_value
  )
}


lib_limits <- iqr_limits(
  library_qc$library_size
)

gene_limits <- iqr_limits(
  library_qc$detected_genes
)

stopifnot(
  !anyNA(lib_limits),
  !anyNA(gene_limits)
)


library_qc$outlier_library_size <-
  library_qc$library_size < lib_limits[["lower"]] |
  library_qc$library_size > lib_limits[["upper"]]


library_qc$outlier_detected_genes <-
  library_qc$detected_genes < gene_limits[["lower"]] |
  library_qc$detected_genes > gene_limits[["upper"]]


library_qc$outlier_any <-
  library_qc$outlier_library_size |
  library_qc$outlier_detected_genes


# There should never be NA values in the QC flags.
stopifnot(
  !anyNA(library_qc$outlier_library_size),
  !anyNA(library_qc$outlier_detected_genes),
  !anyNA(library_qc$outlier_any)
)


write.csv(
  library_qc,
  "results/library_QC.csv",
  row.names = FALSE
)


# =============================================================================
# 12. EXTRACT FLAGGED LIBRARIES
# =============================================================================

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
  flagged,
  "results/QC_flagged_libraries.csv",
  row.names = FALSE
)


# =============================================================================
# 13. PCA WITH QC-FLAGGED LIBRARIES
# =============================================================================

# Match QC information to PCA coordinates using the sample/count identifier.

qc_match <- match(
  pca_df$sample,
  meta$count_column
)

stopifnot(
  !anyNA(qc_match)
)

pca_df$outlier <- library_qc$outlier_any[
  qc_match
]

stopifnot(
  !anyNA(pca_df$outlier)
)


p_outliers <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(
    aes(fill = treatment),
    shape = 21,
    size = 3,
    alpha = 0.65
  ) +
  geom_point(
    data = subset(
      pca_df,
      outlier
    ),
    shape = 21,
    size = 5,
    stroke = 1.2
  ) +
  geom_text(
    data = subset(
      pca_df,
      outlier
    ),
    aes(label = GSM),
    nudge_y = 5,
    size = 3,
    check_overlap = TRUE
  ) +
  labs(
    title = "PCA with exploratory QC flags",
    subtitle = "Libraries flagged using the 1.5 x IQR criterion",
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

ggsave(
  "figures/PCA_QC_flagged_libraries.png",
  p_outliers,
  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# 14. QC SUMMARY TABLE
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
    "library_size_min_million",
    "library_size_max_million",
    "library_size_median_million",
    "detected_genes_min",
    "detected_genes_max",
    "detected_genes_median",
    "libraries_flagged_IQR"
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
    min(library_qc$library_size_millions),
    max(library_qc$library_size_millions),
    median(library_qc$library_size_millions),
    min(library_qc$detected_genes),
    max(library_qc$detected_genes),
    median(library_qc$detected_genes),
    sum(library_qc$outlier_any)
  )
)

write.csv(
  qc_summary,
  "results/QC_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 15. PRINT VALIDATION RESULTS
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
    "PCA variance: PC1 %.1f%% | PC2 %.1f%%\n\n",
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


cat("\n============================================================\n")
cat("QC COMPLETE\n")
cat("============================================================\n")
