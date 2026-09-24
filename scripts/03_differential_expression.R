# =============================================================================
# 03_differential_expression.R
# TFM: Transcriptomic response of human induced neurons to ionizing radiation
#
# Purpose:
#   Differential-expression analysis between irradiated (IR) and control (CTL)
#   human induced neurons using limma-voom.
#
# Main model:
#   - Experimental batch included as a fixed effect
#   - Within-donor dependence estimated with duplicateCorrelation
#   - Contrast of interest: IR vs CTL
#
# Input:
#   data/processed/dge_54_filtered_TMM.rds
#   data/processed/sample_metadata_54.csv
#
# Output:
#   results/DE_results_all.csv
#   results/DE_significant_FDR.csv
#   results/DEGs_FDR_logFC.csv
#   results/DE_summary.csv
#   results/selected_genes_DE.csv
#   data/processed/voom_main_model.rds
#   data/processed/fit_main_model.rds
# =============================================================================


# =============================================================================
# 0. PACKAGES AND DIRECTORIES
# =============================================================================

required <- c("edgeR", "limma")

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
library(limma)

dir.create("results", showWarnings = FALSE)
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 1. LOAD FILTERED AND TMM-NORMALIZED DATA
# =============================================================================

dge <- readRDS(
  "data/processed/dge_54_filtered_TMM.rds"
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


# =============================================================================
# 2. VALIDATE INPUT
# =============================================================================

stopifnot(
  nrow(dge) == 22439L,
  ncol(dge) == 54L,
  nrow(meta) == 54L,
  all(colnames(dge$counts) == meta$count_column),
  sum(meta$treatment == "CTL") == 27L,
  sum(meta$treatment == "IR") == 27L,
  nlevels(meta$donor) == 17L,
  nlevels(meta$batch) == 10L
)


# =============================================================================
# 3. DESIGN MATRIX
# =============================================================================

# Main TFM model:
#
#   expression ~ batch + treatment
#
# Batch is included as a fixed effect.
# Dependence between observations from the same donor is handled below
# using duplicateCorrelation.

design <- model.matrix(
  ~ batch + treatment,
  data = meta
)

colnames(design) <- make.names(
  colnames(design)
)

cat("\nDesign matrix dimensions:",
    nrow(design), "x", ncol(design), "\n")

cat("\nDesign matrix coefficients:\n")
print(colnames(design))


# =============================================================================
# 4. INITIAL VOOM TRANSFORMATION
# =============================================================================

# First voom pass before estimating within-donor correlation.

v0 <- voom(
  dge,
  design = design,
  plot = FALSE
)


# =============================================================================
# 5. ESTIMATE WITHIN-DONOR CORRELATION
# =============================================================================

corfit <- duplicateCorrelation(
  v0,
  design = design,
  block = meta$donor
)

consensus_correlation <- corfit$consensus.correlation

cat(
  sprintf(
    "\nEstimated within-donor correlation: %.6f\n",
    consensus_correlation
  )
)


# =============================================================================
# 6. FINAL VOOM TRANSFORMATION WITH DONOR CORRELATION
# =============================================================================

v <- voom(
  dge,
  design = design,
  block = meta$donor,
  correlation = consensus_correlation,
  plot = FALSE
)

saveRDS(
  v,
  "data/processed/voom_main_model.rds"
)


# =============================================================================
# 7. FIT LINEAR MODEL
# =============================================================================

fit <- lmFit(
  v,
  design = design,
  block = meta$donor,
  correlation = consensus_correlation
)

fit <- eBayes(fit)

saveRDS(
  fit,
  "data/processed/fit_main_model.rds"
)


# =============================================================================
# 8. IDENTIFY TREATMENT COEFFICIENT
# =============================================================================

cat("\nAvailable coefficients:\n")
print(colnames(fit$coefficients))

treatment_coef <- grep(
  "^treatmentIR$",
  colnames(fit$coefficients)
)

if (length(treatment_coef) != 1L) {
  stop(
    "Could not uniquely identify the IR vs CTL treatment coefficient."
  )
}

cat(
  "\nTreatment coefficient:",
  colnames(fit$coefficients)[treatment_coef],
  "\n"
)


# =============================================================================
# 9. EXTRACT ALL DIFFERENTIAL-EXPRESSION RESULTS
# =============================================================================

de_all <- topTable(
  fit,
  coef = treatment_coef,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

de_all$gene_id <- rownames(de_all)

de_all <- de_all[
  ,
  c(
    "gene_id",
    setdiff(
      colnames(de_all),
      "gene_id"
    )
  )
]

write.csv(
  de_all,
  "results/DE_results_all.csv",
  row.names = FALSE
)


# =============================================================================
# 10. FDR < 0.05
# =============================================================================

de_fdr <- de_all[
  !is.na(de_all$adj.P.Val) &
    de_all$adj.P.Val < 0.05,
  ,
  drop = FALSE
]

write.csv(
  de_fdr,
  "results/DE_significant_FDR.csv",
  row.names = FALSE
)


# =============================================================================
# 11. FDR < 0.05 AND |log2FC| >= 1
# =============================================================================

degs <- de_all[
  !is.na(de_all$adj.P.Val) &
    de_all$adj.P.Val < 0.05 &
    abs(de_all$logFC) >= 1,
  ,
  drop = FALSE
]

degs$direction <- ifelse(
  degs$logFC > 0,
  "Up_in_IR",
  "Down_in_IR"
)

write.csv(
  degs,
  "results/DEGs_FDR_logFC.csv",
  row.names = FALSE
)


# =============================================================================
# 12. COUNT DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

n_fdr <- nrow(de_fdr)

n_degs <- nrow(degs)

n_up <- sum(
  degs$logFC >= 1
)

n_down <- sum(
  degs$logFC <= -1
)


# =============================================================================
# 13. SAVE ANALYSIS SUMMARY
# =============================================================================

de_summary <- data.frame(
  metric = c(
    "genes_tested",
    "within_donor_correlation",
    "FDR_lt_0.05",
    "FDR_lt_0.05_and_abs_log2FC_ge_1",
    "upregulated_in_IR",
    "downregulated_in_IR"
  ),

  value = c(
    nrow(de_all),
    consensus_correlation,
    n_fdr,
    n_degs,
    n_up,
    n_down
  )
)

write.csv(
  de_summary,
  "results/DE_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 14. CHECK GENES HIGHLIGHTED IN THE TFM
# =============================================================================

# The count matrix uses Ensembl identifiers. If gene symbols are already
# present in the limma annotation, use them here. Otherwise this section
# reports that symbol annotation must be added in a later step.

selected_symbols <- c(
  "CDKN1A",
  "DDB2",
  "GDF15",
  "CXCL5",
  "GRIA4"
)

symbol_columns <- intersect(
  c(
    "gene_symbol",
    "symbol",
    "SYMBOL",
    "Gene.symbol"
  ),
  colnames(de_all)
)

if (length(symbol_columns) > 0) {

  symbol_col <- symbol_columns[1]

  selected_results <- de_all[
    de_all[[symbol_col]] %in% selected_symbols,
    ,
    drop = FALSE
  ]

  write.csv(
    selected_results,
    "results/selected_genes_DE.csv",
    row.names = FALSE
  )

} else {

  cat(
    "\nNOTE: Gene symbols are not currently attached to the DGE object.\n",
    "The five TFM genes will be validated after Ensembl-to-symbol annotation.\n",
    sep = ""
  )
}


# =============================================================================
# 15. VALIDATION AGAINST TFM CHECKPOINTS
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("TFM DIFFERENTIAL-EXPRESSION SUMMARY\n")
cat("============================================================\n\n")

cat(
  sprintf(
    "Within-donor correlation: %.3f\n",
    consensus_correlation
  )
)

cat(
  "Genes tested:",
  nrow(de_all),
  "\n"
)

cat(
  "Significant genes (FDR < 0.05):",
  n_fdr,
  "\n"
)

cat(
  "DEGs (FDR < 0.05 & |log2FC| >= 1):",
  n_degs,
  "\n"
)

cat(
  "Upregulated in IR:",
  n_up,
  "\n"
)

cat(
  "Downregulated in IR:",
  n_down,
  "\n\n"
)


# =============================================================================
# 16. CHECK EXPECTED TFM RESULTS
# =============================================================================

expected_correlation <- 0.537
expected_fdr <- 7706L
expected_degs <- 466L
expected_up <- 220L
expected_down <- 246L


if (abs(consensus_correlation - expected_correlation) > 0.01) {

  warning(
    "Within-donor correlation differs from the TFM value (~0.537): ",
    round(consensus_correlation, 4)
  )
}


if (n_fdr != expected_fdr) {

  warning(
    "Number of genes with FDR < 0.05 differs from TFM: ",
    n_fdr,
    " vs 7,706."
  )
}


if (n_degs != expected_degs) {

  warning(
    "Number of DEGs differs from TFM: ",
    n_degs,
    " vs 466."
  )
}


if (n_up != expected_up) {

  warning(
    "Number of upregulated DEGs differs from TFM: ",
    n_up,
    " vs 220."
  )
}


if (n_down != expected_down) {

  warning(
    "Number of downregulated DEGs differs from TFM: ",
    n_down,
    " vs 246."
  )
}


cat("Expected TFM checkpoints:\n")
cat("  Within-donor correlation: ~0.537\n")
cat("  FDR < 0.05:              7,706 genes\n")
cat("  FDR + |log2FC| >= 1:       466 genes\n")
cat("  Upregulated in IR:          220 genes\n")
cat("  Downregulated in IR:        246 genes\n")

cat("\n============================================================\n")
cat("DIFFERENTIAL-EXPRESSION ANALYSIS COMPLETE\n")
cat("============================================================\n")
# =============================================================================
# 17. GENE ANNOTATION: ENSEMBL -> GENE SYMBOL
# =============================================================================

annotation_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_annotation_packages <- annotation_packages[
  !vapply(
    annotation_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_annotation_packages) > 0) {
  stop(
    "Install missing Bioconductor packages before continuing: ",
    paste(missing_annotation_packages, collapse = ", "),
    "\n\nUse:\n",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) ",
    "install.packages('BiocManager')\n",
    "BiocManager::install(c('AnnotationDbi', 'org.Hs.eg.db'))"
  )
}

library(AnnotationDbi)
library(org.Hs.eg.db)


# Remove Ensembl version suffix if present.
# Example:
# ENSG00000141510.18 -> ENSG00000141510

de_all$ensembl_id <- sub(
  "\\..*$",
  "",
  de_all$gene_id
)


# Map Ensembl identifiers to HGNC gene symbols.

de_all$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = de_all$ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)


cat("\nGene annotation completed.\n")
cat(
  "Genes with SYMBOL:",
  sum(!is.na(de_all$gene_symbol)),
  "/",
  nrow(de_all),
  "\n"
)


# =============================================================================
# 18. REBUILD SIGNIFICANT TABLES WITH GENE SYMBOLS
# =============================================================================

de_fdr <- de_all[
  !is.na(de_all$adj.P.Val) &
    de_all$adj.P.Val < 0.05,
  ,
  drop = FALSE
]


degs <- de_all[
  !is.na(de_all$adj.P.Val) &
    de_all$adj.P.Val < 0.05 &
    abs(de_all$logFC) >= 1,
  ,
  drop = FALSE
]


degs$direction <- ifelse(
  degs$logFC > 0,
  "Up_in_IR",
  "Down_in_IR"
)


# Save annotated tables

write.csv(
  de_all,
  "results/DE_results_all_annotated.csv",
  row.names = FALSE
)

write.csv(
  de_fdr,
  "results/DE_significant_FDR_annotated.csv",
  row.names = FALSE
)

write.csv(
  degs,
  "results/DEGs_FDR_logFC_annotated.csv",
  row.names = FALSE
)


# =============================================================================
# 19. CHECK GSEA ANNOTATION CHECKPOINT
# =============================================================================

# The TFM reports a ranked list of 17,919 unique mapped gene symbols
# after Ensembl -> SYMBOL conversion and removal of unmapped/duplicate IDs.

gsea_annotation <- de_all[
  !is.na(de_all$gene_symbol) &
    de_all$gene_symbol != "",
  ,
  drop = FALSE
]


# Sort by absolute moderated t statistic so that, when duplicate symbols
# occur, the gene entry carrying the strongest statistical signal is retained.

gsea_annotation <- gsea_annotation[
  order(abs(gsea_annotation$t), decreasing = TRUE),
  ,
  drop = FALSE
]


gsea_annotation <- gsea_annotation[
  !duplicated(gsea_annotation$gene_symbol),
  ,
  drop = FALSE
]


cat(
  "Unique mapped gene symbols:",
  nrow(gsea_annotation),
  "\n"
)


if (nrow(gsea_annotation) != 17919L) {

  warning(
    "Unique mapped gene-symbol count is ",
    nrow(gsea_annotation),
    "; the TFM GSEA ranking reports 17,919. ",
    "This may reflect annotation-database version differences."
  )
}


# Save ranking input for the next GSEA script.

gsea_ranking <- data.frame(
  gene_symbol = gsea_annotation$gene_symbol,
  t = gsea_annotation$t,
  logFC = gsea_annotation$logFC,
  adj.P.Val = gsea_annotation$adj.P.Val
)


gsea_ranking <- gsea_ranking[
  order(gsea_ranking$t, decreasing = TRUE),
]


write.csv(
  gsea_ranking,
  "results/GSEA_ranked_genes.csv",
  row.names = FALSE
)


# =============================================================================
# 20. SELECTED GENES
# =============================================================================

selected_symbols <- c(
  "CDKN1A",
  "DDB2",
  "GDF15",
  "CXCL5",
  "GRIA4"
)


selected_results <- de_all[
  de_all$gene_symbol %in% selected_symbols,
  ,
  drop = FALSE
]


# =============================================================================
# 21. 95% CONFIDENCE INTERVALS
# =============================================================================

# Confidence intervals are calculated from:
#
#   beta +/- t_(df, 0.975) * SE
#
# limma stores the unscaled standard error and posterior residual variance.

coef_name <- colnames(fit$coefficients)[treatment_coef]

coef_values <- fit$coefficients[
  ,
  treatment_coef
]

standard_errors <- fit$stdev.unscaled[
  ,
  treatment_coef
] * sqrt(fit$s2.post)


critical_t <- qt(
  0.975,
  df = fit$df.total
)


ci_lower <- coef_values -
  critical_t * standard_errors

ci_upper <- coef_values +
  critical_t * standard_errors


ci_table <- data.frame(
  gene_id = rownames(fit$coefficients),
  CI95_lower = ci_lower,
  CI95_upper = ci_upper,
  stringsAsFactors = FALSE
)


selected_results <- merge(
  selected_results,
  ci_table,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)


# Restore requested gene order

selected_results$gene_symbol <- factor(
  selected_results$gene_symbol,
  levels = selected_symbols
)

selected_results <- selected_results[
  order(selected_results$gene_symbol),
  ,
  drop = FALSE
]

selected_results$gene_symbol <- as.character(
  selected_results$gene_symbol
)


write.csv(
  selected_results,
  "results/selected_genes_DE.csv",
  row.names = FALSE
)


cat("\nSelected genes:\n")

print(
  selected_results[
    ,
    c(
      "gene_symbol",
      "logFC",
      "CI95_lower",
      "CI95_upper",
      "adj.P.Val"
    )
  ],
  row.names = FALSE
)


# =============================================================================
# 22. CHECK SELECTED GENES AGAINST TFM
# =============================================================================

expected_selected <- data.frame(
  gene_symbol = c(
    "CDKN1A",
    "DDB2",
    "GDF15",
    "CXCL5",
    "GRIA4"
  ),

  expected_logFC = c(
    1.701,
    1.282,
    2.307,
    2.218,
    -1.529
  )
)


selected_check <- merge(
  expected_selected,
  selected_results[
    ,
    c(
      "gene_symbol",
      "logFC",
      "CI95_lower",
      "CI95_upper"
    )
  ],
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)


selected_check$difference <- (
  selected_check$logFC -
    selected_check$expected_logFC
)


write.csv(
  selected_check,
  "results/selected_genes_validation.csv",
  row.names = FALSE
)


cat("\nSelected-gene validation:\n")
print(
  selected_check,
  row.names = FALSE
)


# =============================================================================
# 23. VOLCANO PLOT
# =============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 before generating figures.")
}

library(ggplot2)

dir.create(
  "figures",
  showWarnings = FALSE
)


volcano_df <- de_all


volcano_df$category <- "Not significant"

volcano_df$category[
  volcano_df$adj.P.Val < 0.05 &
    volcano_df$logFC >= 1
] <- "Up in IR"

volcano_df$category[
  volcano_df$adj.P.Val < 0.05 &
    volcano_df$logFC <= -1
] <- "Down in IR"


# Avoid Inf when FDR is numerically zero.

volcano_df$minus_log10_FDR <- -log10(
  pmax(
    volcano_df$adj.P.Val,
    .Machine$double.xmin
  )
)


# Label a small set of highly changed significant genes.

label_df <- degs[
  !is.na(degs$gene_symbol) &
    degs$gene_symbol != "",
  ,
  drop = FALSE
]

label_df <- label_df[
  order(abs(label_df$logFC), decreasing = TRUE),
  ,
  drop = FALSE
]

label_df <- head(
  label_df,
  12
)

label_df$minus_log10_FDR <- -log10(
  pmax(
    label_df$adj.P.Val,
    .Machine$double.xmin
  )
)


p_volcano <- ggplot(
  volcano_df,
  aes(
    x = logFC,
    y = minus_log10_FDR
  )
) +

  geom_point(
    aes(color = category),
    alpha = 0.55,
    size = 1.4
  ) +

  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +

  labs(
    title = "Differential expression: IR vs CTL",
    subtitle = "Human induced neurons at 336 h",
    x = expression(log[2] * " fold change (IR vs CTL)"),
    y = expression(-log[10] * "(FDR)"),
    color = NULL
  ) +

  theme_classic(
    base_size = 12
  )


# Use ggrepel when available to avoid overlapping labels.

if (requireNamespace("ggrepel", quietly = TRUE)) {

  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        label = gene_symbol
      ),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    )

} else {

  p_volcano <- p_volcano +
    geom_text(
      data = label_df,
      aes(
        label = gene_symbol
      ),
      size = 3,
      vjust = -0.5,
      show.legend = FALSE
    )
}


ggsave(
  "figures/volcano_IR_vs_CTL.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/volcano_IR_vs_CTL.pdf",
  p_volcano,
  width = 8,
  height = 6
)


# =============================================================================
# 24. HEATMAP: TOP 30 DEGs
# =============================================================================
#
# Reproduce the TFM heatmap using the 15 restrictive DEGs with the highest
# log2FC and the 15 restrictive DEGs with the lowest log2FC.
#
# Expression values correspond to voom-normalized expression and are
# standardized per gene using row-wise Z-scores.
# =============================================================================


# -----------------------------------------------------------------------------
# 24.1 Select 15 most upregulated and 15 most downregulated DEGs
# -----------------------------------------------------------------------------

top_up_heat <- degs[
  order(
    degs$logFC,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

top_up_heat <- head(
  top_up_heat,
  15
)


top_down_heat <- degs[
  order(
    degs$logFC,
    decreasing = FALSE
  ),
  ,
  drop = FALSE
]

top_down_heat <- head(
  top_down_heat,
  15
)


heatmap_genes <- rbind(
  top_up_heat,
  top_down_heat
)


stopifnot(
  nrow(heatmap_genes) == 30L
)


# -----------------------------------------------------------------------------
# 24.2 Prepare gene labels
# -----------------------------------------------------------------------------
#
# Use gene symbol when available. If no symbol is available, use the Ensembl
# identifier without its version suffix.
# -----------------------------------------------------------------------------

heatmap_labels <- ifelse(
  !is.na(heatmap_genes$gene_symbol) &
    heatmap_genes$gene_symbol != "",
  heatmap_genes$gene_symbol,
  heatmap_genes$ensembl_id
)

heatmap_labels <- make.unique(
  heatmap_labels
)


# -----------------------------------------------------------------------------
# 24.3 Extract voom-normalized expression
# -----------------------------------------------------------------------------

heatmap_expression <- v$E[
  heatmap_genes$gene_id,
  ,
  drop = FALSE
]

rownames(heatmap_expression) <-
  heatmap_labels


# -----------------------------------------------------------------------------
# 24.4 Standardize expression per gene
# -----------------------------------------------------------------------------

heatmap_z <- t(
  scale(
    t(heatmap_expression)
  )
)


# -----------------------------------------------------------------------------
# 24.5 Treatment annotation
# -----------------------------------------------------------------------------

heatmap_annotation <- data.frame(
  Tratamiento = meta$treatment
)

rownames(heatmap_annotation) <-
  colnames(heatmap_z)

heatmap_annotation <- heatmap_annotation[
  colnames(heatmap_z),
  ,
  drop = FALSE
]


# -----------------------------------------------------------------------------
# 24.6 Treatment colours
# -----------------------------------------------------------------------------

annotation_colors <- list(
  Tratamiento = c(
    CTL = "#00BFC4",
    IR = "#F8766D"
  )
)


# -----------------------------------------------------------------------------
# 24.7 Generate heatmap
# -----------------------------------------------------------------------------

if (requireNamespace("pheatmap", quietly = TRUE)) {

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    main = "Patrón de expresión de los 30 DEGs más destacados",
    filename = "figures/heatmap_top30_DEGs.png",
    width = 9,
    height = 7
  )


  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    main = "Patrón de expresión de los 30 DEGs más destacados",
    filename = "figures/heatmap_top30_DEGs.pdf",
    width = 9,
    height = 7
  )


  cat(
    "\nHeatmap generated using 15 upregulated + 15 downregulated DEGs.\n"
  )
  cat("PNG: figures/heatmap_top30_DEGs.png\n")
  cat("PDF: figures/heatmap_top30_DEGs.pdf\n")


} else {

  warning(
    "Package 'pheatmap' is not installed. ",
    "Heatmap was skipped. Install it with install.packages('pheatmap')."
  )
}


# =============================================================================
# 25. SELECTED-GENE EXPRESSION DATA
# =============================================================================

# Save normalized expression values for the five representative genes.
# This can later be used for boxplots or supplementary figures.

selected_expression_ids <- selected_results$gene_id

selected_expression <- v$E[
  selected_expression_ids,
  ,
  drop = FALSE
]


rownames(selected_expression) <-
  selected_results$gene_symbol


selected_expression_long <- data.frame(
  gene_symbol = rep(
    rownames(selected_expression),
    each = ncol(selected_expression)
  ),

  sample = rep(
    colnames(selected_expression),
    times = nrow(selected_expression)
  ),

  expression = as.vector(
    t(selected_expression)
  ),

  stringsAsFactors = FALSE
)


sample_meta_small <- data.frame(
  sample = meta$count_column,
  GSM = meta$GSM,
  treatment = meta$treatment,
  donor = meta$donor,
  batch = meta$batch,
  stringsAsFactors = FALSE
)


selected_expression_long <- merge(
  selected_expression_long,
  sample_meta_small,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)


write.csv(
  selected_expression_long,
  "results/selected_genes_expression.csv",
  row.names = FALSE
)


# =============================================================================
# 26. FINAL CHECKPOINT
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("ANNOTATION AND VISUALIZATION SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Genes tested:",
  nrow(de_all),
  "\n"
)

cat(
  "Genes with SYMBOL:",
  sum(!is.na(de_all$gene_symbol)),
  "\n"
)

cat(
  "Unique mapped SYMBOLs for GSEA:",
  nrow(gsea_annotation),
  "\n\n"
)

cat(
  "FDR < 0.05:",
  nrow(de_fdr),
  "\n"
)

cat(
  "DEGs:",
  nrow(degs),
  "\n"
)

cat(
  "Up in IR:",
  sum(degs$logFC >= 1),
  "\n"
)

cat(
  "Down in IR:",
  sum(degs$logFC <= -1),
  "\n\n"
)

cat("Selected genes:\n")

print(
  selected_results[
    ,
    c(
      "gene_symbol",
      "logFC",
      "CI95_lower",
      "CI95_upper"
    )
  ],
  row.names = FALSE
)

cat("\nFiles generated:\n")
cat("  results/DE_results_all_annotated.csv\n")
cat("  results/DE_significant_FDR_annotated.csv\n")
cat("  results/DEGs_FDR_logFC_annotated.csv\n")
cat("  results/selected_genes_DE.csv\n")
cat("  results/selected_genes_validation.csv\n")
cat("  results/selected_genes_expression.csv\n")
cat("  results/GSEA_ranked_genes.csv\n")
cat("  figures/volcano_IR_vs_CTL.png\n")
cat("  figures/volcano_IR_vs_CTL.pdf\n")
cat("  figures/heatmap_top30_DEGs.png\n")
cat("  figures/heatmap_top30_DEGs.pdf\n")

cat("\n============================================================\n")
cat("03 DIFFERENTIAL EXPRESSION COMPLETE\n")
cat("============================================================\n")
# =============================================================================
# 27. DISPLAY SAVED FIGURES IN RSTUDIO
# =============================================================================
# Force the principal figures produced by this script to the active graphics
# device. Saved PNG/PDF files are unaffected.

if (interactive()) {

  if (exists("p_volcano")) {
    print(p_volcano)
  }

  if (exists("heatmap_z") && requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(
      heatmap_z,
      annotation_col = heatmap_annotation,
      annotation_colors = annotation_colors,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      fontsize_row = 8,
      main = "Patrón de expresión de los 30 DEGs más destacados"
    )
  }

  cat("\nFigures sent to the active RStudio Plots device.\n")
  cat("Use the Plots history arrows to move between them.\n")
}
