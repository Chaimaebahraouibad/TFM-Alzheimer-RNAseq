# =============================================================================
# 06_original_study_comparison.R
#
# TFM: Transcriptomic response of human induced neurons to ionizing radiation
#
# Purpose:
#   Prepare the results used for the qualitative comparison between the
#   present transcriptomic reanalysis and the original study.
#
# IMPORTANT:
#   This script does NOT attempt to reproduce the original DESeq2 analysis.
#
#   The comparison is qualitative because the original study and the present
#   reanalysis differ in statistical method, model specification, DEG
#   threshold and functional-enrichment strategy.
#
# Inputs:
#   results/DE_results_all_annotated.csv
#   results/GSEA/Hallmark_GSEA_all.csv
#
# Outputs:
#   results/original_study_comparison/
#   figures/original_study_comparison/
# =============================================================================


# =============================================================================
# 1. REQUIRED PACKAGES
# =============================================================================

required <- c(
    "ggplot2"
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
        paste(
            missing_packages,
            collapse = ", "
        )
    )
}

library(ggplot2)


# =============================================================================
# 2. OUTPUT DIRECTORIES
# =============================================================================

dir.create(
    "results/original_study_comparison",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "figures/original_study_comparison",
    recursive = TRUE,
    showWarnings = FALSE
)


# =============================================================================
# 3. INPUT FILES
# =============================================================================

de_file <- "results/DE_results_all_annotated.csv"

hallmark_file <- "results/GSEA/Hallmark_GSEA_all.csv"


if (!file.exists(de_file)) {

    stop(
        "Cannot find ",
        de_file,
        ". Run 03_differential_expression.R first."
    )
}


if (!file.exists(hallmark_file)) {

    stop(
        "Cannot find ",
        hallmark_file,
        ". Run 04_GSEA.R first."
    )
}


# =============================================================================
# 4. LOAD RESULTS
# =============================================================================

de <- read.csv(
    de_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


hallmark <- read.csv(
    hallmark_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("06 ORIGINAL-STUDY COMPARISON\n")
cat("============================================================\n\n")


cat(
    "DE table:",
    nrow(de),
    "genes\n"
)


cat(
    "Hallmark table:",
    nrow(hallmark),
    "pathways\n\n"
)


# =============================================================================
# 5. VALIDATE DE TABLE
# =============================================================================

required_de_columns <- c(
    "gene_id",
    "gene_symbol",
    "logFC",
    "P.Value",
    "adj.P.Val"
)


missing_de_columns <- setdiff(
    required_de_columns,
    names(de)
)


if (length(missing_de_columns) > 0L) {

    stop(
        "DE table does not contain the expected columns: ",
        paste(
            missing_de_columns,
            collapse = ", "
        )
    )
}


if (nrow(de) != 22439L) {

    warning(
        "Expected 22,439 genes, but found ",
        nrow(de),
        "."
    )
}


# =============================================================================
# 6. DEFINE DIFFERENTIAL-EXPRESSION RESULTS
# =============================================================================

de$significant_FDR <- (
    !is.na(de$adj.P.Val) &
    de$adj.P.Val < 0.05
)


de$DEG_restrictive <- (
    !is.na(de$adj.P.Val) &
    de$adj.P.Val < 0.05 &
    abs(de$logFC) >= 1
)


de$direction <- "Not_restrictive_DEG"


de$direction[
    de$DEG_restrictive &
    de$logFC > 0
] <- "Up_in_IR"


de$direction[
    de$DEG_restrictive &
    de$logFC < 0
] <- "Down_in_IR"


n_tested <- nrow(de)


n_fdr <- sum(
    de$significant_FDR,
    na.rm = TRUE
)


n_degs <- sum(
    de$DEG_restrictive,
    na.rm = TRUE
)


n_up <- sum(
    de$direction == "Up_in_IR",
    na.rm = TRUE
)


n_down <- sum(
    de$direction == "Down_in_IR",
    na.rm = TRUE
)


cat("Differential-expression summary:\n")

cat(
    "  Genes tested:",
    n_tested,
    "\n"
)

cat(
    "  FDR < 0.05:",
    n_fdr,
    "\n"
)

cat(
    "  Restrictive DEGs:",
    n_degs,
    "\n"
)

cat(
    "  Upregulated in IR:",
    n_up,
    "\n"
)

cat(
    "  Downregulated in IR:",
    n_down,
    "\n\n"
)


# =============================================================================
# 7. REPRESENTATIVE GENES
# =============================================================================
#
# Genes used in the interpretation of the transcriptomic response.
# MMP3 is also included because it links the transcriptomic and structural
# sections of the TFM.
# =============================================================================

representative_symbols <- c(
    "CDKN1A",
    "DDB2",
    "GDF15",
    "CXCL5",
    "GRIA4",
    "MMP3"
)


representative_genes <- de[
    !is.na(de$gene_symbol) &
    de$gene_symbol %in%
        representative_symbols,
    ,
    drop = FALSE
]


# In case several Ensembl entries map to the same symbol,
# keep the entry with the smallest adjusted p-value.

representative_genes <- representative_genes[
    order(
        representative_genes$gene_symbol,
        representative_genes$adj.P.Val
    ),
    ,
    drop = FALSE
]


representative_genes <- representative_genes[
    !duplicated(
        representative_genes$gene_symbol
    ),
    ,
    drop = FALSE
]


representative_genes$gene_symbol <- factor(
    representative_genes$gene_symbol,
    levels = representative_symbols
)


representative_genes <- representative_genes[
    order(
        representative_genes$gene_symbol
    ),
    ,
    drop = FALSE
]


representative_genes$gene_symbol <- as.character(
    representative_genes$gene_symbol
)


representative_output_columns <- intersect(
    c(
        "gene_id",
        "ensembl_id",
        "gene_symbol",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B",
        "direction"
    ),
    names(representative_genes)
)


write.csv(
    representative_genes[
        ,
        representative_output_columns,
        drop = FALSE
    ],
    "results/original_study_comparison/representative_genes_reanalysis.csv",
    row.names = FALSE
)


cat("Representative genes found:\n")

print(
    representative_genes[
        ,
        intersect(
            c(
                "gene_symbol",
                "logFC",
                "adj.P.Val",
                "direction"
            ),
            names(representative_genes)
        ),
        drop = FALSE
    ],
    row.names = FALSE
)


# =============================================================================
# 8. SAFE EXTRACTION FUNCTION
# =============================================================================
#
# Returns one numeric value from a one-row data frame.
# If the value cannot be extracted unambiguously, NA is returned.
# =============================================================================

get_single_value <- function(
    df,
    column
) {

    if (
        is.data.frame(df) &&
        nrow(df) == 1L &&
        column %in% names(df)
    ) {

        value <- suppressWarnings(
            as.numeric(
                df[[column]][1]
            )
        )

        if (
            length(value) == 1L &&
            is.finite(value)
        ) {

            return(value)
        }
    }

    return(NA_real_)
}


# =============================================================================
# 9. CDKN1A RESULT
# =============================================================================

cdkn1a <- representative_genes[
    representative_genes$gene_symbol ==
        "CDKN1A",
    ,
    drop = FALSE
]


cdkn1a_logfc <- get_single_value(
    cdkn1a,
    "logFC"
)


cdkn1a_fdr <- get_single_value(
    cdkn1a,
    "adj.P.Val"
)


if (nrow(cdkn1a) != 1L) {

    warning(
        "CDKN1A could not be uniquely identified."
    )
}


# =============================================================================
# 10. VALIDATE HALLMARK TABLE
# =============================================================================

required_hallmark_columns <- c(
    "pathway",
    "NES",
    "padj"
)


missing_hallmark_columns <- setdiff(
    required_hallmark_columns,
    names(hallmark)
)


if (length(missing_hallmark_columns) > 0L) {

    stop(
        "Hallmark GSEA table does not contain the expected columns: ",
        paste(
            missing_hallmark_columns,
            collapse = ", "
        )
    )
}


# =============================================================================
# 11. SELECT HALLMARK PROCESSES
# =============================================================================

hallmark_interest_names <- c(
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_G2M_CHECKPOINT"
)


hallmark_interest <- hallmark[
    hallmark$pathway %in%
        hallmark_interest_names,
    ,
    drop = FALSE
]


hallmark_interest$pathway <- factor(
    hallmark_interest$pathway,
    levels = hallmark_interest_names
)


hallmark_interest <- hallmark_interest[
    order(
        hallmark_interest$pathway
    ),
    ,
    drop = FALSE
]


hallmark_interest$pathway <- as.character(
    hallmark_interest$pathway
)


write.csv(
    hallmark_interest,
    paste0(
        "results/original_study_comparison/",
        "Hallmark_processes_for_original_study_comparison.csv"
    ),
    row.names = FALSE
)


cat("\nSelected Hallmark pathways:\n")

print(
    hallmark_interest[
        ,
        intersect(
            c(
                "pathway",
                "NES",
                "padj"
            ),
            names(hallmark_interest)
        ),
        drop = FALSE
    ],
    row.names = FALSE
)


# =============================================================================
# 12. EXTRACT SPECIFIC HALLMARK RESULTS
# =============================================================================

p53_result <- hallmark[
    hallmark$pathway ==
        "HALLMARK_P53_PATHWAY",
    ,
    drop = FALSE
]


tnfa_result <- hallmark[
    hallmark$pathway ==
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    ,
    drop = FALSE
]


p53_nes <- get_single_value(
    p53_result,
    "NES"
)


p53_fdr <- get_single_value(
    p53_result,
    "padj"
)


tnfa_nes <- get_single_value(
    tnfa_result,
    "NES"
)


tnfa_fdr <- get_single_value(
    tnfa_result,
    "padj"
)


# =============================================================================
# 13. METHODOLOGICAL COMPARISON
# =============================================================================
#
# Descriptive comparison of the analytical strategies.
#
# No statistical equivalence between the analyses is implied.
# =============================================================================

method_comparison <- data.frame(

    characteristic = c(
        "Differential-expression method",
        "Statistical design",
        "Treatment of donor dependence",
        "Adjusted-p-value threshold",
        "Absolute log2FC threshold",
        "Functional analysis"
    ),

    original_study = c(
        "DESeq2",
        "~ batch + condition",
        "Not explicitly incorporated in the design described for comparison",
        "< 0.05",
        "> 0.5",
        "Functional analysis of increased and decreased genes"
    ),

    present_reanalysis = c(
        "limma-voom",
        "~ batch + treatment",
        "duplicateCorrelation with donor as blocking factor",
        "< 0.05",
        ">= 1",
        "GSEA using the complete ranked gene list"
    ),

    stringsAsFactors = FALSE
)


write.csv(
    method_comparison,
    paste0(
        "results/original_study_comparison/",
        "methodological_comparison.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# 14. BIOLOGICAL COMPARISON FRAMEWORK
# =============================================================================
#
# This table organizes the qualitative biological comparison used in the TFM.
#
# It does NOT imply that the gene sets or enrichment statistics from the
# original study and the present reanalysis are directly equivalent.
# =============================================================================

biological_comparison <- data.frame(

    biological_theme = c(
        "p53-associated response",
        "Inflammatory / NF-kB signalling",
        "DNA-damage response",
        "Apoptosis",
        "Oxidative / metabolic response",
        "Synaptic / neuronal processes"
    ),

    present_reanalysis_evidence = c(
        "Positive HALLMARK_P53_PATHWAY enrichment and increased CDKN1A expression",
        "Positive HALLMARK_TNFA_SIGNALING_VIA_NFKB enrichment",
        "Positive HALLMARK_DNA_REPAIR enrichment and increased DDB2 expression",
        "Positive HALLMARK_APOPTOSIS enrichment",
        "Positive HALLMARK_OXIDATIVE_PHOSPHORYLATION enrichment",
        "Several neuronal and synaptic processes showed relative enrichment toward CTL"
    ),

    comparison_scope = rep(
        "Qualitative comparison of biological trends",
        6
    ),

    stringsAsFactors = FALSE
)


write.csv(
    biological_comparison,
    paste0(
        "results/original_study_comparison/",
        "biological_comparison_framework.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# 15. FIGURE: REPRESENTATIVE GENE EFFECT SIZES
# =============================================================================
#
# IMPORTANT:
#
# This figure contains effect sizes from the PRESENT REANALYSIS ONLY.
# It is not a quantitative effect-size comparison with the original study.
# =============================================================================

plot_symbols <- c(
    "CDKN1A",
    "DDB2",
    "GDF15",
    "CXCL5",
    "GRIA4"
)


plot_genes <- representative_genes[
    representative_genes$gene_symbol %in%
        plot_symbols,
    ,
    drop = FALSE
]


plot_genes <- plot_genes[
    is.finite(plot_genes$logFC),
    ,
    drop = FALSE
]


if (nrow(plot_genes) > 0L) {

    plot_genes$gene_symbol <- factor(
        plot_genes$gene_symbol,
        levels = c(
            "GRIA4",
            "DDB2",
            "CDKN1A",
            "CXCL5",
            "GDF15"
        )
    )


    p_genes <- ggplot(
        plot_genes,
        aes(
            x = gene_symbol,
            y = logFC
        )
    ) +

        geom_col(
            width = 0.7
        ) +

        geom_hline(
            yintercept = 0,
            linetype = "dashed"
        ) +

        coord_flip() +

        labs(
            title =
                "Genes representativos de la respuesta a irradiación",

            subtitle =
                "Estimaciones del presente reanálisis",

            x = NULL,

            y = expression(
                log[2] * "FC (IR vs CTL)"
            )
        ) +

        theme_classic(
            base_size = 12
        )


    print(
        p_genes
    )


    ggsave(
        paste0(
            "figures/original_study_comparison/",
            "representative_gene_logFC.png"
        ),
        plot = p_genes,
        width = 7,
        height = 5,
        dpi = 300
    )


    ggsave(
        paste0(
            "figures/original_study_comparison/",
            "representative_gene_logFC.pdf"
        ),
        plot = p_genes,
        width = 7,
        height = 5
    )
}


# =============================================================================
# 16. FIGURE: SELECTED HALLMARK PROCESSES
# =============================================================================

hallmark_plot <- hallmark_interest[
    is.finite(hallmark_interest$NES) &
    !is.na(hallmark_interest$padj),
    ,
    drop = FALSE
]


if (nrow(hallmark_plot) > 0L) {

    hallmark_plot$short_name <- hallmark_plot$pathway


    hallmark_plot$short_name <- sub(
        "^HALLMARK_",
        "",
        hallmark_plot$short_name
    )


    hallmark_plot$short_name <- gsub(
        "_",
        " ",
        hallmark_plot$short_name
    )


    hallmark_plot <- hallmark_plot[
        order(
            hallmark_plot$NES
        ),
        ,
        drop = FALSE
    ]


    hallmark_plot$short_name <- factor(
        hallmark_plot$short_name,
        levels = hallmark_plot$short_name
    )


    p_hallmark <- ggplot(
        hallmark_plot,
        aes(
            x = short_name,
            y = NES
        )
    ) +

        geom_col(
            width = 0.7
        ) +

        geom_hline(
            yintercept = 0,
            linetype = "dashed"
        ) +

        coord_flip() +

        labs(
            title =
                "Procesos Hallmark relevantes para la comparación",

            subtitle =
                "GSEA del presente reanálisis",

            x = NULL,

            y =
                "Normalized Enrichment Score (NES)"
        ) +

        theme_classic(
            base_size = 12
        )


    print(
        p_hallmark
    )


    ggsave(
        paste0(
            "figures/original_study_comparison/",
            "Hallmark_comparison_processes.png"
        ),
        plot = p_hallmark,
        width = 8,
        height = 6,
        dpi = 300
    )


    ggsave(
        paste0(
            "figures/original_study_comparison/",
            "Hallmark_comparison_processes.pdf"
        ),
        plot = p_hallmark,
        width = 8,
        height = 6
    )
}


# =============================================================================
# 17. REANALYSIS SUMMARY
# =============================================================================
#
# All values are converted explicitly to numeric scalars before constructing
# the data frame. This avoids failures caused by zero-length objects or
# data-frame columns.
# =============================================================================

reanalysis_summary <- data.frame(

    metric = c(
        "Genes tested",
        "Genes with FDR < 0.05",
        "Restrictive DEGs",
        "Upregulated DEGs in IR",
        "Downregulated DEGs in IR",
        "CDKN1A log2FC",
        "CDKN1A FDR",
        "P53 Hallmark NES",
        "P53 Hallmark FDR",
        "TNFA/NFKB Hallmark NES",
        "TNFA/NFKB Hallmark FDR"
    ),

    value = as.numeric(
        c(
            n_tested,
            n_fdr,
            n_degs,
            n_up,
            n_down,
            cdkn1a_logfc,
            cdkn1a_fdr,
            p53_nes,
            p53_fdr,
            tnfa_nes,
            tnfa_fdr
        )
    ),

    stringsAsFactors = FALSE
)


stopifnot(
    nrow(reanalysis_summary) == 11L,
    all(
        c(
            "metric",
            "value"
        ) %in%
            names(reanalysis_summary)
    )
)


write.csv(
    reanalysis_summary,
    paste0(
        "results/original_study_comparison/",
        "reanalysis_summary.csv"
    ),
    row.names = FALSE
)


cat("\nReanalysis summary:\n")

print(
    reanalysis_summary,
    row.names = FALSE
)


# =============================================================================
# 18. TEXT REPORT
# =============================================================================

report_file <- paste0(
    "results/original_study_comparison/",
    "original_study_comparison_summary.txt"
)


report_connection <- file(
    report_file,
    open = "wt"
)


writeLines(
    c(
        "============================================================",
        "COMPARISON WITH ORIGINAL STUDY",
        "============================================================",
        "",
        "IMPORTANT",
        paste(
            "This comparison is qualitative and does not attempt",
            "to reproduce the original DESeq2 analysis gene by gene."
        ),
        "",
        "PRESENT REANALYSIS",
        paste(
            "Genes tested:",
            n_tested
        ),
        paste(
            "FDR < 0.05:",
            n_fdr
        ),
        paste(
            "Restrictive DEGs:",
            n_degs
        ),
        paste(
            "Upregulated in IR:",
            n_up
        ),
        paste(
            "Downregulated in IR:",
            n_down
        ),
        "",
        paste(
            "CDKN1A log2FC:",
            format(
                cdkn1a_logfc,
                digits = 6
            )
        ),
        paste(
            "CDKN1A FDR:",
            format(
                cdkn1a_fdr,
                scientific = TRUE,
                digits = 6
            )
        ),
        "",
        paste(
            "HALLMARK P53 NES:",
            format(
                p53_nes,
                digits = 6
            )
        ),
        paste(
            "HALLMARK P53 FDR:",
            format(
                p53_fdr,
                scientific = TRUE,
                digits = 6
            )
        ),
        "",
        paste(
            "HALLMARK TNFA/NFKB NES:",
            format(
                tnfa_nes,
                digits = 6
            )
        ),
        paste(
            "HALLMARK TNFA/NFKB FDR:",
            format(
                tnfa_fdr,
                scientific = TRUE,
                digits = 6
            )
        ),
        "",
        "============================================================",
        "END OF REPORT",
        "============================================================"
    ),
    con = report_connection
)


close(
    report_connection
)


# =============================================================================
# 19. CHECK EXPECTED REANALYSIS RESULTS
# =============================================================================
#
# These checks refer to results already established in scripts 03 and 04.
# They are used to detect accidental changes in upstream analyses.
# =============================================================================

expected_tested <- 22439L
expected_fdr <- 7706L
expected_degs <- 466L
expected_up <- 220L
expected_down <- 246L


if (n_tested != expected_tested) {

    warning(
        "Genes tested differ from expected value: ",
        n_tested,
        " vs 22,439."
    )
}


if (n_fdr != expected_fdr) {

    warning(
        "FDR-significant gene count differs from expected value: ",
        n_fdr,
        " vs 7,706."
    )
}


if (n_degs != expected_degs) {

    warning(
        "Restrictive DEG count differs from expected value: ",
        n_degs,
        " vs 466."
    )
}


if (n_up != expected_up) {

    warning(
        "Upregulated DEG count differs from expected value: ",
        n_up,
        " vs 220."
    )
}


if (n_down != expected_down) {

    warning(
        "Downregulated DEG count differs from expected value: ",
        n_down,
        " vs 246."
    )
}


# =============================================================================
# 20. SESSION INFORMATION
# =============================================================================

capture.output(
    sessionInfo(),
    file = paste0(
        "results/original_study_comparison/",
        "sessionInfo_comparison.txt"
    )
)


# =============================================================================
# 21. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("TFM ORIGINAL-STUDY COMPARISON SUMMARY\n")
cat("============================================================\n\n")


cat("Present reanalysis:\n")


cat(
    "  Genes tested:",
    n_tested,
    "\n"
)


cat(
    "  FDR < 0.05:",
    n_fdr,
    "\n"
)


cat(
    "  Restrictive DEGs:",
    n_degs,
    "\n"
)


cat(
    "  Upregulated in IR:",
    n_up,
    "\n"
)


cat(
    "  Downregulated in IR:",
    n_down,
    "\n\n"
)


cat("Representative result:\n")


if (is.finite(cdkn1a_logfc)) {

    cat(
        sprintf(
            "  CDKN1A: log2FC = %.3f | FDR = %.3e\n\n",
            cdkn1a_logfc,
            cdkn1a_fdr
        )
    )

} else {

    cat(
        "  CDKN1A: result not available\n\n"
    )
}


cat("Selected functional results:\n")


if (is.finite(p53_nes)) {

    cat(
        sprintf(
            "  HALLMARK P53: NES = %.3f | FDR = %.3e\n",
            p53_nes,
            p53_fdr
        )
    )

} else {

    cat(
        "  HALLMARK P53: result not available\n"
    )
}


if (is.finite(tnfa_nes)) {

    cat(
        sprintf(
            "  HALLMARK TNFA/NFKB: NES = %.3f | FDR = %.3e\n",
            tnfa_nes,
            tnfa_fdr
        )
    )

} else {

    cat(
        "  HALLMARK TNFA/NFKB: result not available\n"
    )
}


cat("\nMethodological context:\n")


cat(
    paste0(
        "  Original study: DESeq2, ~ batch + condition, ",
        "adjusted p < 0.05, |log2FC| > 0.5\n"
    )
)


cat(
    paste0(
        "  Present reanalysis: limma-voom + duplicateCorrelation, ",
        "FDR < 0.05, |log2FC| >= 1\n"
    )
)


cat("\nComparison scope:\n")


cat(
    paste0(
        "  Qualitative comparison of biological trends. ",
        "No exact DEG-list overlap is calculated because the analytical ",
        "strategies and DEG thresholds differ.\n"
    )
)


cat("\nFiles generated under:\n")

cat(
    "  results/original_study_comparison/\n"
)

cat(
    "  figures/original_study_comparison/\n"
)


cat("\n============================================================\n")
cat("06 ORIGINAL-STUDY COMPARISON COMPLETE\n")
cat("============================================================\n")
