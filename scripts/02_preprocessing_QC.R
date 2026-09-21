# 02_preprocessing_QC.R
# Purpose: low-expression filtering, TMM normalization and PCA QC for the
# 54-library TFM analysis set.

required <- c("edgeR", "ggplot2")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Install missing packages before running this script: ",
       paste(required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)], collapse = ", "))
}

library(edgeR)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

counts <- readRDS("data/processed/count_matrix_54.rds")
meta <- read.csv("data/processed/sample_metadata_54.csv", stringsAsFactors = FALSE)
meta$treatment <- factor(meta$treatment, levels = c("CTL", "IR"))
meta$batch <- factor(meta$batch)
meta$donor <- factor(meta$donor)

stopifnot(ncol(counts) == nrow(meta), colnames(counts) == meta$count_column)

# The filtering step in the TFM considered the experimental design.
# This donor-batch pairing design is used ONLY to guide filterByExpr; the final
# differential-expression model is reconstructed in 03_differential_expression.R.
meta$donor_batch <- interaction(meta$donor, meta$batch, drop = TRUE)
filter_design <- model.matrix(~ donor_batch + treatment, data = meta)

dge <- DGEList(counts = counts, samples = meta)
keep <- filterByExpr(dge, design = filter_design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# TFM checkpoint: 78,428 initial gene IDs -> 22,439 after filterByExpr.
qc_summary <- data.frame(
  metric = c("initial_genes", "filtered_genes", "libraries", "CTL", "IR", "donors", "batches"),
  value = c(nrow(counts), nrow(dge), ncol(dge), sum(meta$treatment == "CTL"),
            sum(meta$treatment == "IR"), nlevels(meta$donor), nlevels(meta$batch))
)
write.csv(qc_summary, "results/QC_summary.csv", row.names = FALSE)

if (nrow(dge) != 22439L) {
  warning("Filtered gene count is ", nrow(dge),
          "; the final TFM reports 22,439. Check package versions and sample mapping before proceeding.")
}

# PCA described in TFM methodology: log2-CPM with prior.count = 0.5.
log_cpm <- cpm(dge, log = TRUE, prior.count = 0.5)
pca <- prcomp(t(log_cpm), scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  treatment = meta$treatment,
  batch = meta$batch,
  donor = meta$donor
)
write.csv(pca_df, "results/PCA_coordinates.csv", row.names = FALSE)

p <- ggplot(pca_df, aes(PC1, PC2, fill = treatment)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  labs(
    title = "Principal component analysis",
    subtitle = "Human induced neurons, 336 h",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    fill = "Condition"
  ) +
  theme_classic(base_size = 12)

ggsave("figures/PCA_54_libraries.pdf", p, width = 7, height = 5)
ggsave("figures/PCA_54_libraries.png", p, width = 7, height = 5, dpi = 300)

saveRDS(dge, "data/processed/dge_54_filtered_TMM.rds")
saveRDS(log_cpm, "data/processed/logCPM_54.rds")

cat("Initial genes:", nrow(counts), "\n")
cat("Genes after filterByExpr:", nrow(dge), "\n")
cat(sprintf("PCA variance: PC1 %.1f%% | PC2 %.1f%%\n", var_exp[1], var_exp[2]))
