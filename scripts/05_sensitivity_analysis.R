# =============================================================================
# 05_sensitivity_analysis.R
#
# TFM: Transcriptomic response of human induced neurons to ionizing radiation
#
# Purpose:
#   Evaluate the robustness of the differential-expression results using an
#   alternative statistical model in which donor is included as a fixed effect.
#
# Main model (03):
#   expression ~ batch + treatment
#   within-donor dependence handled with duplicateCorrelation
#
# Sensitivity model (05):
#   expression ~ donor + treatment
#
# Comparisons:
#   - genome-wide log2FC correlation
#   - restrictive DEGs (FDR < 0.05 and |log2FC| >= 1)
#   - DEG overlap
#   - direction concordance
#   - log2FC correlation among shared DEGs
#   - selected genes discussed in the TFM
#
# IMPORTANT:
#   The Venn diagram reports absolute counts only.
#   Percentages are calculated separately because percentages relative to the
#   union of both DEG sets can be misleading for interpretation.
# =============================================================================


# =============================================================================
# 1. REQUIRED PACKAGES
# =============================================================================

required <- c(
    "edgeR",
    "limma",
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

library(edgeR)
library(limma)
library(ggplot2)


# =============================================================================
# 2. OUTPUT DIRECTORIES
# =============================================================================

dir.create(
    "results/sensitivity",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "figures/sensitivity",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "data/processed",
    recursive = TRUE,
    showWarnings = FALSE
)


# =============================================================================
# 3. INPUT FILES
# =============================================================================

dge_file <- "data/processed/dge_54_filtered_TMM.rds"

meta_file <- "data/processed/sample_metadata_54.csv"

main_file <- "results/DE_results_all_annotated.csv"


if (!file.exists(dge_file)) {

    stop(
        "Cannot find ",
        dge_file,
        ". Run 02_preprocessing_QC.R first."
    )
}


if (!file.exists(meta_file)) {

    stop(
        "Cannot find ",
        meta_file,
        ". Run 01_data_preparation.R first."
    )
}


if (!file.exists(main_file)) {

    stop(
        "Cannot find ",
        main_file,
        ". Run 03_differential_expression.R first."
    )
}


# =============================================================================
# 4. LOAD DATA
# =============================================================================

dge <- readRDS(
    dge_file
)


meta <- read.csv(
    meta_file,
    stringsAsFactors = FALSE
)


main <- read.csv(
    main_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


# =============================================================================
# 5. PREPARE METADATA
# =============================================================================

meta$treatment <- factor(
    meta$treatment,
    levels = c(
        "CTL",
        "IR"
    )
)


meta$donor <- factor(
    meta$donor
)


meta$batch <- factor(
    meta$batch
)


# =============================================================================
# 6. VALIDATE INPUT DATA
# =============================================================================

stopifnot(
    nrow(dge) == 22439L,
    ncol(dge) == 54L,
    nrow(meta) == 54L,
    all(
        colnames(dge$counts) ==
            meta$count_column
    ),
    sum(meta$treatment == "CTL") == 27L,
    sum(meta$treatment == "IR") == 27L,
    nlevels(meta$donor) == 17L,
    nlevels(meta$batch) == 10L
)


required_main_columns <- c(
    "gene_id",
    "logFC",
    "adj.P.Val"
)


if (!all(
    required_main_columns %in%
        names(main)
)) {

    stop(
        "Main DE table does not contain the expected columns: ",
        paste(
            required_main_columns,
            collapse = ", "
        )
    )
}


cat("\n")
cat("============================================================\n")
cat("SENSITIVITY ANALYSIS INPUT\n")
cat("============================================================\n\n")


cat(
    "Genes:",
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
    "\n\n"
)


# =============================================================================
# 7. SENSITIVITY DESIGN MATRIX
# =============================================================================
#
# Alternative model:
#
#     expression ~ donor + treatment
#
# Donor is included directly as a fixed effect.
#
# duplicateCorrelation is NOT used in this model.
# =============================================================================

design_sensitivity <- model.matrix(
    ~ donor + treatment,
    data = meta
)


colnames(design_sensitivity) <- make.names(
    colnames(design_sensitivity)
)


cat(
    "Sensitivity design dimensions:",
    nrow(design_sensitivity),
    "x",
    ncol(design_sensitivity),
    "\n"
)


cat(
    "\nSensitivity model coefficients:\n"
)


print(
    colnames(design_sensitivity)
)


# =============================================================================
# 8. CHECK DESIGN MATRIX
# =============================================================================

design_rank <- qr(
    design_sensitivity
)$rank


cat(
    "\nDesign rank:",
    design_rank,
    "/",
    ncol(design_sensitivity),
    "\n"
)


if (
    design_rank !=
    ncol(design_sensitivity)
) {

    stop(
        "Sensitivity design matrix is not full rank."
    )
}


# =============================================================================
# 9. VOOM TRANSFORMATION
# =============================================================================

v_sensitivity <- voom(
    dge,
    design = design_sensitivity,
    plot = FALSE
)


saveRDS(
    v_sensitivity,
    "data/processed/voom_sensitivity_model.rds"
)


# =============================================================================
# 10. FIT SENSITIVITY MODEL
# =============================================================================

fit_sensitivity <- lmFit(
    v_sensitivity,
    design_sensitivity
)


fit_sensitivity <- eBayes(
    fit_sensitivity
)


saveRDS(
    fit_sensitivity,
    "data/processed/fit_sensitivity_model.rds"
)


# =============================================================================
# 11. IDENTIFY TREATMENT COEFFICIENT
# =============================================================================

treatment_coef <- grep(
    "^treatmentIR$",
    colnames(
        fit_sensitivity$coefficients
    )
)


if (length(treatment_coef) != 1L) {

    stop(
        "Could not uniquely identify treatmentIR in sensitivity model."
    )
}


cat(
    "\nTreatment coefficient:",
    colnames(
        fit_sensitivity$coefficients
    )[treatment_coef],
    "\n"
)


# =============================================================================
# 12. EXTRACT SENSITIVITY RESULTS
# =============================================================================

sensitivity <- topTable(
    fit_sensitivity,
    coef = treatment_coef,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
)


sensitivity$gene_id <- rownames(
    sensitivity
)


sensitivity <- sensitivity[
    ,
    c(
        "gene_id",
        setdiff(
            names(sensitivity),
            "gene_id"
        )
    ),
    drop = FALSE
]


# =============================================================================
# 13. ADD GENE ANNOTATION
# =============================================================================

annotation_columns <- intersect(
    c(
        "gene_id",
        "ensembl_id",
        "gene_symbol"
    ),
    names(main)
)


annotation_table <- unique(
    main[
        ,
        annotation_columns,
        drop = FALSE
    ]
)


sensitivity_order <- sensitivity$gene_id


sensitivity <- merge(
    sensitivity,
    annotation_table,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


sensitivity <- sensitivity[
    match(
        sensitivity_order,
        sensitivity$gene_id
    ),
    ,
    drop = FALSE
]


rownames(sensitivity) <- NULL


write.csv(
    sensitivity,
    "results/sensitivity/DE_sensitivity_all.csv",
    row.names = FALSE
)


# =============================================================================
# 14. DEFINE SIGNIFICANT GENES IN SENSITIVITY MODEL
# =============================================================================

sensitivity_fdr <- sensitivity[
    !is.na(sensitivity$adj.P.Val) &
        sensitivity$adj.P.Val < 0.05,
    ,
    drop = FALSE
]


sensitivity_degs <- sensitivity[
    !is.na(sensitivity$adj.P.Val) &
        sensitivity$adj.P.Val < 0.05 &
        abs(sensitivity$logFC) >= 1,
    ,
    drop = FALSE
]


sensitivity_degs$direction <- ifelse(
    sensitivity_degs$logFC > 0,
    "Up_in_IR",
    "Down_in_IR"
)


write.csv(
    sensitivity_fdr,
    "results/sensitivity/DE_sensitivity_FDR.csv",
    row.names = FALSE
)


write.csv(
    sensitivity_degs,
    "results/sensitivity/DEGs_sensitivity.csv",
    row.names = FALSE
)


# =============================================================================
# 15. DEFINE MAIN-MODEL DEGs
# =============================================================================

main_degs <- main[
    !is.na(main$adj.P.Val) &
        main$adj.P.Val < 0.05 &
        abs(main$logFC) >= 1,
    ,
    drop = FALSE
]


main_degs$direction <- ifelse(
    main_degs$logFC > 0,
    "Up_in_IR",
    "Down_in_IR"
)


# =============================================================================
# 16. MATCH ALL GENES BETWEEN MODELS
# =============================================================================

main_comparison_columns <- c(
    "gene_id",
    "logFC",
    "adj.P.Val"
)


if ("gene_symbol" %in% names(main)) {

    main_comparison_columns <- c(
        "gene_id",
        "gene_symbol",
        "logFC",
        "adj.P.Val"
    )
}


comparison <- merge(
    main[
        ,
        main_comparison_columns,
        drop = FALSE
    ],
    sensitivity[
        ,
        c(
            "gene_id",
            "logFC",
            "adj.P.Val"
        ),
        drop = FALSE
    ],
    by = "gene_id",
    suffixes = c(
        "_main",
        "_sensitivity"
    ),
    all = FALSE
)


comparison <- comparison[
    is.finite(comparison$logFC_main) &
        is.finite(
            comparison$logFC_sensitivity
        ),
    ,
    drop = FALSE
]


cat(
    "\nGenes comparable between models:",
    nrow(comparison),
    "\n"
)


# =============================================================================
# 17. GENOME-WIDE LOG2FC CORRELATION
# =============================================================================

global_logfc_correlation <- cor(
    comparison$logFC_main,
    comparison$logFC_sensitivity,
    method = "pearson",
    use = "complete.obs"
)


cat(
    sprintf(
        "Genome-wide log2FC correlation: r = %.3f\n",
        global_logfc_correlation
    )
)


# =============================================================================
# 18. DEG OVERLAP
# =============================================================================

main_deg_ids <- unique(
    main_degs$gene_id
)


sensitivity_deg_ids <- unique(
    sensitivity_degs$gene_id
)


shared_ids <- intersect(
    main_deg_ids,
    sensitivity_deg_ids
)


main_only_ids <- setdiff(
    main_deg_ids,
    sensitivity_deg_ids
)


sensitivity_only_ids <- setdiff(
    sensitivity_deg_ids,
    main_deg_ids
)


n_main <- length(
    main_deg_ids
)


n_sensitivity <- length(
    sensitivity_deg_ids
)


n_shared <- length(
    shared_ids
)


n_main_only <- length(
    main_only_ids
)


n_sensitivity_only <- length(
    sensitivity_only_ids
)


# =============================================================================
# 19. SHARED DEG TABLE
# =============================================================================

shared_degs <- comparison[
    comparison$gene_id %in%
        shared_ids,
    ,
    drop = FALSE
]


shared_degs$direction_main <- ifelse(
    shared_degs$logFC_main > 0,
    "Up_in_IR",
    "Down_in_IR"
)


shared_degs$direction_sensitivity <- ifelse(
    shared_degs$logFC_sensitivity > 0,
    "Up_in_IR",
    "Down_in_IR"
)


shared_degs$same_direction <- (
    shared_degs$direction_main ==
        shared_degs$direction_sensitivity
)


n_same_direction <- sum(
    shared_degs$same_direction
)


direction_concordance_pct <- (
    100 *
        n_same_direction /
        nrow(shared_degs)
)


shared_logfc_correlation <- cor(
    shared_degs$logFC_main,
    shared_degs$logFC_sensitivity,
    method = "pearson",
    use = "complete.obs"
)


write.csv(
    shared_degs,
    "results/sensitivity/shared_DEGs.csv",
    row.names = FALSE
)


# =============================================================================
# 20. EXCLUSIVE DEG TABLES
# =============================================================================

main_only <- main_degs[
    main_degs$gene_id %in%
        main_only_ids,
    ,
    drop = FALSE
]


sensitivity_only <- sensitivity_degs[
    sensitivity_degs$gene_id %in%
        sensitivity_only_ids,
    ,
    drop = FALSE
]


write.csv(
    main_only,
    "results/sensitivity/DEGs_main_only.csv",
    row.names = FALSE
)


write.csv(
    sensitivity_only,
    "results/sensitivity/DEGs_sensitivity_only.csv",
    row.names = FALSE
)


# =============================================================================
# 21. OVERLAP PERCENTAGES
# =============================================================================
#
# IMPORTANT:
#
# These percentages use each model as its own denominator.
#
# Main-model overlap:
#     shared / all main-model DEGs
#
# Sensitivity-model overlap:
#     shared / all sensitivity-model DEGs
#
# They should NOT be confused with percentages relative to the union of the
# two sets.
# =============================================================================

pct_main_shared <- (
    100 *
        n_shared /
        n_main
)


pct_sensitivity_shared <- (
    100 *
        n_shared /
        n_sensitivity
)


# =============================================================================
# 22. SELECTED TFM GENES
# =============================================================================

selected_symbols <- c(
    "CDKN1A",
    "DDB2",
    "GDF15",
    "CXCL5",
    "GRIA4"
)


if ("gene_symbol" %in% names(comparison)) {

    selected_comparison <- comparison[
        comparison$gene_symbol %in%
            selected_symbols,
        ,
        drop = FALSE
    ]


    selected_comparison$gene_symbol <- factor(
        selected_comparison$gene_symbol,
        levels = selected_symbols
    )


    selected_comparison <- selected_comparison[
        order(
            selected_comparison$gene_symbol
        ),
        ,
        drop = FALSE
    ]


    selected_comparison$gene_symbol <- as.character(
        selected_comparison$gene_symbol
    )


    write.csv(
        selected_comparison,
        "results/sensitivity/selected_genes_model_comparison.csv",
        row.names = FALSE
    )

} else {

    selected_comparison <- data.frame()

    warning(
        "gene_symbol is not available in the main DE table; ",
        "selected-gene comparison was skipped."
    )
}


# =============================================================================
# 23. SAVE COMPLETE MODEL COMPARISON
# =============================================================================

comparison$DEG_main <- (
    comparison$adj.P.Val_main < 0.05 &
        abs(comparison$logFC_main) >= 1
)


comparison$DEG_sensitivity <- (
    comparison$adj.P.Val_sensitivity < 0.05 &
        abs(
            comparison$logFC_sensitivity
        ) >= 1
)


write.csv(
    comparison,
    "results/sensitivity/model_comparison_all_genes.csv",
    row.names = FALSE
)


# =============================================================================
# 24. SUMMARY TABLE
# =============================================================================

summary_table <- data.frame(

    metric = c(
        "genes_tested_sensitivity",
        "genes_comparable",
        "sensitivity_FDR_lt_0.05",
        "main_DEGs",
        "sensitivity_DEGs",
        "shared_DEGs",
        "main_only_DEGs",
        "sensitivity_only_DEGs",
        "shared_pct_of_main",
        "shared_pct_of_sensitivity",
        "genome_wide_logFC_correlation",
        "shared_DEG_logFC_correlation",
        "shared_DEGs_same_direction",
        "shared_DEG_direction_concordance_pct"
    ),

    value = c(
        nrow(sensitivity),
        nrow(comparison),
        nrow(sensitivity_fdr),
        n_main,
        n_sensitivity,
        n_shared,
        n_main_only,
        n_sensitivity_only,
        pct_main_shared,
        pct_sensitivity_shared,
        global_logfc_correlation,
        shared_logfc_correlation,
        n_same_direction,
        direction_concordance_pct
    ),

    stringsAsFactors = FALSE
)


write.csv(
    summary_table,
    "results/sensitivity/sensitivity_summary.csv",
    row.names = FALSE
)


# =============================================================================
# 25. FIGURE: GENOME-WIDE LOG2FC CONCORDANCE
# =============================================================================

p_concordance <- ggplot(
    comparison,
    aes(
        x = logFC_main,
        y = logFC_sensitivity
    )
) +

    geom_point(
        alpha = 0.35,
        size = 1.2
    ) +

    geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed"
    ) +

    labs(
        title = "Concordancia entre modelos de expresión diferencial",

        subtitle = sprintf(
            "Correlación de Pearson: r = %.3f",
            global_logfc_correlation
        ),

        x = expression(
            log[2] * "FC - modelo principal"
        ),

        y = expression(
            log[2] * "FC - modelo alternativo"
        )
    ) +

    theme_classic(
        base_size = 12
    )


print(
    p_concordance
)


ggsave(
    "figures/sensitivity/logFC_model_concordance.png",
    plot = p_concordance,
    width = 7,
    height = 6,
    dpi = 300
)


ggsave(
    "figures/sensitivity/logFC_model_concordance.pdf",
    plot = p_concordance,
    width = 7,
    height = 6
)


# =============================================================================
# 26. FIGURE: SHARED DEG CONCORDANCE
# =============================================================================

p_shared <- ggplot(
    shared_degs,
    aes(
        x = logFC_main,
        y = logFC_sensitivity
    )
) +

    geom_point(
        alpha = 0.65,
        size = 1.8
    ) +

    geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed"
    ) +

    labs(
        title = "Concordancia de los DEGs compartidos",

        subtitle = sprintf(
            "Pearson r = %.3f | Misma dirección = %.1f%%",
            shared_logfc_correlation,
            direction_concordance_pct
        ),

        x = expression(
            log[2] * "FC - modelo principal"
        ),

        y = expression(
            log[2] * "FC - modelo alternativo"
        )
    ) +

    theme_classic(
        base_size = 12
    )


print(
    p_shared
)


ggsave(
    "figures/sensitivity/shared_DEG_logFC_concordance.png",
    plot = p_shared,
    width = 7,
    height = 6,
    dpi = 300
)


ggsave(
    "figures/sensitivity/shared_DEG_logFC_concordance.pdf",
    plot = p_shared,
    width = 7,
    height = 6
)


# =============================================================================
# 27. FIGURE: DEG OVERLAP
# =============================================================================
#
# IMPORTANT:
#
# The Venn diagram contains ABSOLUTE COUNTS ONLY.
#
# Percentages are deliberately omitted because ggVennDiagram normally
# calculates percentages relative to the union of both sets.
#
# The percentages relevant to the TFM are calculated separately:
#
#     shared / main DEGs
#
# and
#
#     shared / sensitivity DEGs
#
# This avoids confusion between different denominators.
# =============================================================================

if (
    requireNamespace(
        "ggVennDiagram",
        quietly = TRUE
    )
) {

    venn_list <- list(

        "Modelo principal" =
            main_deg_ids,

        "Modelo alternativo" =
            sensitivity_deg_ids
    )


    p_venn <- ggVennDiagram::ggVennDiagram(
        venn_list,
        label = "count",
        label_alpha = 0,
        edge_size = 0.8,
        set_size = 4
    ) +

        labs(
            title = "Solapamiento de DEGs entre modelos"
        ) +

        coord_cartesian(
            clip = "off"
        ) +

        theme_void() +

        theme(

            plot.title = element_text(
                hjust = 0.5,
                size = 14,
                face = "bold",
                margin = margin(
                    b = 15
                )
            ),

            plot.margin = margin(
                t = 40,
                r = 140,
                b = 40,
                l = 140
            )
        )


    print(
        p_venn
    )


    ggsave(
        "figures/sensitivity/DEG_overlap_venn.png",
        plot = p_venn,
        width = 11,
        height = 7,
        dpi = 300,
        limitsize = FALSE
    )


    ggsave(
        "figures/sensitivity/DEG_overlap_venn.pdf",
        plot = p_venn,
        width = 11,
        height = 7,
        limitsize = FALSE
    )

} else {

    # -------------------------------------------------------------------------
    # Fallback figure if ggVennDiagram is not installed
    # -------------------------------------------------------------------------

    overlap_df <- data.frame(

        category = factor(

            c(
                "Solo principal",
                "Compartidos",
                "Solo alternativo"
            ),

            levels = c(
                "Solo principal",
                "Compartidos",
                "Solo alternativo"
            )
        ),

        genes = c(
            n_main_only,
            n_shared,
            n_sensitivity_only
        )
    )


    p_overlap <- ggplot(
        overlap_df,
        aes(
            x = category,
            y = genes
        )
    ) +

        geom_col(
            width = 0.65
        ) +

        geom_text(
            aes(
                label = genes
            ),
            vjust = -0.4,
            size = 4.5
        ) +

        labs(
            title = "Solapamiento de DEGs entre modelos",
            x = NULL,
            y = "Número de DEGs"
        ) +

        theme_classic(
            base_size = 12
        )


    print(
        p_overlap
    )


    ggsave(
        "figures/sensitivity/DEG_overlap.png",
        plot = p_overlap,
        width = 8,
        height = 5,
        dpi = 300
    )


    ggsave(
        "figures/sensitivity/DEG_overlap.pdf",
        plot = p_overlap,
        width = 8,
        height = 5
    )
}


# =============================================================================
# 28. SESSION INFORMATION
# =============================================================================

capture.output(
    sessionInfo(),
    file =
        "results/sensitivity/sessionInfo_sensitivity.txt"
)


# =============================================================================
# 29. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("TFM SENSITIVITY ANALYSIS SUMMARY\n")
cat("============================================================\n\n")


cat(
    "Sensitivity model: expression ~ donor + treatment\n\n"
)


cat(
    "Genes tested:",
    nrow(sensitivity),
    "\n"
)


cat(
    "Genes comparable between models:",
    nrow(comparison),
    "\n\n"
)


cat(
    "Genome-wide log2FC correlation:",
    sprintf(
        "r = %.3f",
        global_logfc_correlation
    ),
    "\n\n"
)


cat(
    "Restrictive DEGs:\n"
)


cat(
    "  Main model:",
    n_main,
    "\n"
)


cat(
    "  Sensitivity model:",
    n_sensitivity,
    "\n"
)


cat(
    "  Shared:",
    n_shared,
    "\n"
)


cat(
    "  Main only:",
    n_main_only,
    "\n"
)


cat(
    "  Sensitivity only:",
    n_sensitivity_only,
    "\n\n"
)


cat(
    sprintf(
        "Shared DEGs represent %.1f%% of main-model DEGs\n",
        pct_main_shared
    )
)


cat(
    sprintf(
        "Shared DEGs represent %.1f%% of sensitivity-model DEGs\n\n",
        pct_sensitivity_shared
    )
)


cat(
    sprintf(
        "Shared-DEG log2FC correlation: r = %.3f\n",
        shared_logfc_correlation
    )
)


cat(
    sprintf(
        "Shared DEGs with same direction: %d/%d (%.1f%%)\n\n",
        n_same_direction,
        n_shared,
        direction_concordance_pct
    )
)


if (nrow(selected_comparison) > 0L) {

    cat(
        "Selected TFM genes:\n"
    )


    print(
        selected_comparison[
            ,
            c(
                "gene_symbol",
                "logFC_main",
                "adj.P.Val_main",
                "logFC_sensitivity",
                "adj.P.Val_sensitivity"
            ),
            drop = FALSE
        ],
        row.names = FALSE
    )
}


cat("\n")
cat("Interpretation of overlap percentages:\n")


cat(
    sprintf(
        "  %d / %d = %.1f%% of main-model DEGs reproduced\n",
        n_shared,
        n_main,
        pct_main_shared
    )
)


cat(
    sprintf(
        "  %d / %d = %.1f%% of sensitivity-model DEGs shared\n",
        n_shared,
        n_sensitivity,
        pct_sensitivity_shared
    )
)


cat(
    "\nNote: percentages relative to the union of both DEG sets ",
    "are intentionally not reported.\n",
    sep = ""
)


cat("\nFiles generated under:\n")
cat("  results/sensitivity/\n")
cat("  figures/sensitivity/\n")


cat("\n============================================================\n")
cat("05 SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================================\n")
