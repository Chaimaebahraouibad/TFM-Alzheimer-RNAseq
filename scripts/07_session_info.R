# =============================================================================
# 07_session_info.R
#
# TFM: Transcriptomic response of human induced neurons to ionizing radiation
#
# Purpose:
#   Record the R environment and the versions of the packages used in the
#   transcriptomic reanalysis, so that the computational environment can be
#   documented in the GitHub repository.
#
# Run this script at the end of the analysis.
#
# Outputs:
#   results/reproducibility/sessionInfo.txt
#   results/reproducibility/package_versions.csv
# =============================================================================


# =============================================================================
# 1. OUTPUT DIRECTORY
# =============================================================================

output_dir <- "results/reproducibility"

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# =============================================================================
# 2. PACKAGES USED ACROSS SCRIPTS 01-06
# =============================================================================

packages_used <- c(
    "GEOquery",
    "data.table",
    "edgeR",
    "limma",
    "ggplot2",
    "patchwork",
    "AnnotationDbi",
    "org.Hs.eg.db",
    "fgsea",
    "msigdbr"
)


# =============================================================================
# 3. CHECK PACKAGE AVAILABILITY
# =============================================================================

installed <- vapply(
    packages_used,
    requireNamespace,
    logical(1),
    quietly = TRUE
)

if (any(!installed)) {
    warning(
        "The following packages used in the analysis are not installed in ",
        "the current R environment: ",
        paste(packages_used[!installed], collapse = ", "),
        ". Their versions will be recorded as NA."
    )
}


# =============================================================================
# 4. SAVE PACKAGE VERSIONS
# =============================================================================

package_versions <- data.frame(
    package = packages_used,
    version = vapply(
        packages_used,
        function(pkg) {
            if (requireNamespace(pkg, quietly = TRUE)) {
                as.character(packageVersion(pkg))
            } else {
                NA_character_
            }
        },
        character(1)
    ),
    stringsAsFactors = FALSE
)

write.csv(
    package_versions,
    file.path(
        output_dir,
        "package_versions.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# 5. SAVE COMPLETE R SESSION INFORMATION
# =============================================================================

session_file <- file.path(
    output_dir,
    "sessionInfo.txt"
)

session_lines <- capture.output(
    sessionInfo()
)

writeLines(
    session_lines,
    con = session_file
)


# =============================================================================
# 6. OPTIONAL REPRODUCIBILITY SUMMARY
# =============================================================================

summary_file <- file.path(
    output_dir,
    "reproducibility_summary.txt"
)

summary_lines <- c(
    "TFM computational environment",
    "================================",
    "",
    paste("R version:", R.version.string),
    paste("Platform:", R.version$platform),
    "",
    "Packages used in scripts 01-06:",
    paste0(
        "  ",
        package_versions$package,
        ": ",
        ifelse(
            is.na(package_versions$version),
            "not installed in current environment",
            package_versions$version
        )
    ),
    "",
    "Complete session information is available in sessionInfo.txt."
)

writeLines(
    summary_lines,
    con = summary_file
)


# =============================================================================
# 7. CONSOLE SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("07 REPRODUCIBILITY INFORMATION\n")
cat("============================================================\n\n")

cat("R version:\n")
cat(" ", R.version.string, "\n\n")

cat("Package versions:\n")
print(
    package_versions,
    row.names = FALSE
)

cat("\nFiles generated:\n")
cat(" ", session_file, "\n")
cat(
    " ",
    file.path(output_dir, "package_versions.csv"),
    "\n"
)
cat(" ", summary_file, "\n")

cat("\n============================================================\n")
cat("07 REPRODUCIBILITY INFORMATION COMPLETE\n")
cat("============================================================\n")
