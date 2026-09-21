# 01_data_preparation.R
# TFM: Transcriptomic response of human induced neurons to ionizing radiation
# Purpose: download GSE329677 processed counts and reconstruct the 54-library
# analysis set used in the final TFM (iNs, CTL/IR, 336 h).

required <- c("GEOquery", "data.table")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Install missing packages before running this script: ",
       paste(required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)], collapse = ", "))
}

library(GEOquery)
library(data.table)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# 1. Download processed gene-count matrix deposited with GSE329677
counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329677/suppl/",
  "GSE329677_salmon.merged.gene_counts_ALL_samples.tsv.gz"
)
counts_file <- "data/raw/GSE329677_counts.tsv.gz"
if (!file.exists(counts_file)) {
  download.file(counts_url, counts_file, mode = "wb")
}
counts <- fread(counts_file)

# The downloaded table contains 78,428 genes and 232 columns:
# one gene_id column plus 231 RNA-seq libraries.
stopifnot(
  nrow(counts) == 78428L,
  ncol(counts) == 232L
)

# 2. Obtain GEO sample metadata
# The TFM analysis uses induced neurons at 336 h in CTL and IR conditions.
gse <- getGEO("GSE329677", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

is_in <- grepl("iN", pheno$title)
is_336 <- grepl("time 336", pheno$title)
is_ctl_ir <- grepl(", (CTL|IR),", pheno$title)
pheno_54 <- pheno[is_in & is_336 & is_ctl_ir, , drop = FALSE]

# 3. Parse experimental variables from GEO titles
pheno_54$donor <- sub(".*donor ([^,]+).*", "\\1", pheno_54$title)
pheno_54$treatment <- ifelse(grepl(", CTL,", pheno_54$title), "CTL", "IR")
pheno_54$batch <- sub(".*batch ([^,]+).*", "\\1", pheno_54$title)
pheno_54$time_h <- 336L

# 4. Match GEO samples to columns in the processed count matrix.
# In GEO metadata, the processed count-column identifier is embedded in description.
count_names <- names(counts)[-1]
stopifnot(length(count_names) == 231L)
pheno_54$count_column <- vapply(pheno_54$description, function(desc) {
  hits <- count_names[vapply(count_names, function(nm) grepl(nm, desc, fixed = TRUE), logical(1))]
  if (length(hits) == 1L) hits else NA_character_
}, character(1))

if (anyNA(pheno_54$count_column)) {
  stop("Could not uniquely map every selected GEO sample to a count-matrix column.")
}
if (!all(pheno_54$count_column %in% names(counts))) {
  stop("At least one selected count column is absent from the downloaded matrix.")
}

# 5. Reorder and define factors used downstream
pheno_54$treatment <- factor(pheno_54$treatment, levels = c("CTL", "IR"))
pheno_54$donor <- factor(pheno_54$donor)
pheno_54$batch <- factor(pheno_54$batch)

# The final TFM contains 54 libraries: 27 CTL + 27 IR, 17 donor IDs, 10 batches.
stopifnot(
  nrow(pheno_54) == 54L,
  all(table(pheno_54$treatment) == c(CTL = 27L, IR = 27L)),
  nlevels(pheno_54$donor) == 17L,
  nlevels(pheno_54$batch) == 10L
)

# 27 donor-batch combinations, each containing one CTL and one IR library.
pheno_54$donor_batch <- interaction(pheno_54$donor, pheno_54$batch, drop = TRUE)
pair_table <- table(pheno_54$donor_batch, pheno_54$treatment)
stopifnot(nrow(pair_table) == 27L, all(pair_table == 1L))

# 6. Build the analysis count matrix
counts_54 <- counts[, c("gene_id", pheno_54$count_column), with = FALSE]
count_matrix_54 <- as.matrix(counts_54[, -1, with = FALSE])
rownames(count_matrix_54) <- counts_54$gene_id
storage.mode(count_matrix_54) <- "numeric"

# Save compact reproducible inputs for subsequent scripts.
metadata_out <- data.frame(
  GSM = pheno_54$geo_accession,
  donor = as.character(pheno_54$donor),
  batch = as.character(pheno_54$batch),
  treatment = as.character(pheno_54$treatment),
  time_h = pheno_54$time_h,
  count_column = pheno_54$count_column,
  stringsAsFactors = FALSE
)
write.csv(metadata_out, "data/processed/sample_metadata_54.csv", row.names = FALSE)
saveRDS(count_matrix_54, "data/processed/count_matrix_54.rds")

cat("Prepared", ncol(count_matrix_54), "libraries and", nrow(count_matrix_54), "genes.\n")
cat("CTL/IR:", paste(names(table(pheno_54$treatment)), table(pheno_54$treatment), collapse = "; "), "\n")
cat("Donors:", nlevels(pheno_54$donor), "| Batches:", nlevels(pheno_54$batch), "| Donor-batch pairs:", nlevels(pheno_54$donor_batch), "\n")
