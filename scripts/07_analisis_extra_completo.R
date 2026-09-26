# =============================================================================
# TFM - SCRIPT DE ANALISIS
# Respuesta transcriptomica de neuronas humanas inducidas a radiacion ionizante
# Reconstruido del registro de comandos de RStudio: sin prompts ni salidas de consola.
# =============================================================================


    data.frame(
        Proceso = unname(pathways[out$pathway]),
        NES = out$NES,
        FDR = out$padj,
        Coleccion = collection,
        stringsAsFactors = FALSE
    )
}

integrated_hallmark_df <- make_integrated_df(
    hallmark_gsea,
    integrated_hallmark,
    "Hallmark"
)

integrated_go_df <- make_integrated_df(
    gobp_gsea,
    integrated_go,
    "GO"
)

integrated_kegg_df <- make_integrated_df(
    kegg_gsea,
    integrated_kegg,
    "KEGG"
)

funcional <- rbind(
    integrated_hallmark_df,
    integrated_go_df,
    integrated_kegg_df
)

funcional$Direccion <- ifelse(
    funcional$NES > 0,
    "Enriquecida hacia IR",
    "Enriquecida hacia CTL"
)

funcional$Significancia <- -log10(funcional$FDR)

funcional$Proceso <- factor(
    funcional$Proceso,
    levels = funcional$Proceso[order(funcional$NES)]
)

p_integrated <- ggplot(
    funcional,
    aes(
        x = NES,
        y = Proceso,
        color = Direccion,
        shape = Coleccion,
        size = Significancia
    )
) +
    geom_point(alpha = 0.9) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
    ) +
    scale_color_manual(
        values = c(
            "Enriquecida hacia CTL" = "#F8766D",
            "Enriquecida hacia IR" = "#00BFC4"
        )
    ) +
    scale_shape_manual(
        values = c(
            Hallmark = 17,
            GO = 16,
            KEGG = 15
        )
    ) +
    labs(
        title = "Perfil funcional integrado de la respuesta a IR",
        subtitle = "Hallmark, GO y KEGG en neuronas inducidas a 336 h",
        x = "NES",
        y = NULL,
        color = "Dirección",
        shape = "Colección",
        size = expression(-log[10](FDR))
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
    )

ggsave(
    "figures/GSEA/integrated_functional_profile.pdf",
    p_integrated,
    width = 12,
    height = 5.5
)
ggsave(
    "figures/GSEA/integrated_functional_profile.png",
    p_integrated,
    width = 12,
    height = 5.5,
    dpi = 300
)


# =============================================================================
# 21. SAVE R OBJECTS
# =============================================================================

saveRDS(
    ranks,
    "results/GSEA/GSEA_ranked_genes.rds"
)

saveRDS(
    hallmark_gsea,
    "results/GSEA/Hallmark_GSEA.rds"
)

saveRDS(
    gobp_gsea,
    "results/GSEA/GO_BP_GSEA.rds"
)

saveRDS(
    kegg_gsea,
    "results/GSEA/KEGG_Legacy_GSEA.rds"
)


# =============================================================================
# 22. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("TFM GSEA SUMMARY\n")
cat("============================================================\n\n")


cat(
    "Genes in ranked list:",
    length(ranks),
    "\n\n"
)


cat("HALLMARK\n")

cat(
    "  Tested:",
    nrow(hallmark_gsea),
    "\n"
)

cat(
    "  Significant:",
    nrow(hallmark_sig),
    "\n"
)

cat(
    "  Positive / IR:",
    sum(hallmark_sig$NES > 0),
    "\n"
)

cat(
    "  Negative / CTL:",
    sum(hallmark_sig$NES < 0),
    "\n\n"
)


cat("GO:BP\n")

cat(
    "  Tested:",
    nrow(gobp_gsea),
    "\n"
)

cat(
    "  Significant:",
    nrow(gobp_sig),
    "\n"
)

cat(
    "  Positive / IR:",
    nrow(gobp_positive),
    "\n"
)

cat(
    "  Negative / CTL:",
    nrow(gobp_negative),
    "\n\n"
)


cat("KEGG LEGACY\n")

cat(
    "  Tested:",
    nrow(kegg_gsea),
    "\n"
)

cat(
    "  Significant:",
    nrow(kegg_sig),
    "\n"
)

cat(
    "  Positive / IR:",
    nrow(kegg_positive),
    "\n"
)

cat(
    "  Negative / CTL:",
    nrow(kegg_negative),
    "\n\n"
)


cat("Selected Hallmark pathways:\n")

print(
    hallmark_interest[
        ,
        c(
            "pathway",
            "size",
            "ES",
            "NES",
            "padj"
        )
    ],
    row.names = FALSE
)


cat("\nLeading-edge sizes:\n")

print(
    hallmark_leading_edge_table[
        ,
        c(
            "pathway",
            "leading_edge_n"
        )
    ],
    row.names = FALSE
)


cat("\nLeading-edge overlap:\n")

print(
    leading_edge_overlap[
        ,
        c(
            "comparison",
            "shared_genes_n"
        )
    ],
    row.names = FALSE
)


# =============================================================================
# 23. REPRODUCIBILITY CHECKPOINTS
# =============================================================================

expected_ranked_genes <- 17919L
expected_hallmark <- 50L
expected_gobp <- 3825L
expected_kegg <- 165L


cat("\nExpected reproducibility checkpoints:\n")

cat(
    "  Ranked genes:              17,919\n"
)

cat(
    "  Hallmark tested:               50\n"
)

cat(
    "  GO:BP terms evaluated:      3,825\n"
)

cat(
    "  KEGG pathways evaluated:      165\n"
)

cat(
    "  Hallmark p53 NES:           ~3.25\n"
)

cat(
    "  Hallmark DNA repair NES:    ~2.34\n"
)

cat(
    "  Hallmark G2/M NES:         ~-2.09\n"
)


if (length(ranks) != expected_ranked_genes) {

    warning(
        "Ranked gene count differs from expected value: ",
        length(ranks),
        " vs 17,919."
    )
}


if (nrow(hallmark_gsea) != expected_hallmark) {

    warning(
        "Hallmark pathway count differs from expected value: ",
        nrow(hallmark_gsea),
        " vs 50."
    )
}


if (nrow(gobp_gsea) != expected_gobp) {

    warning(
        "GO:BP term count differs from expected value: ",
        nrow(gobp_gsea),
        " vs 3,825. Check MSigDB/msigdbr version."
    )
}


if (nrow(kegg_gsea) != expected_kegg) {

    warning(
        "KEGG pathway count differs from expected value: ",
        nrow(kegg_gsea),
        " vs 165. Check MSigDB/msigdbr version."
    )
}


# =============================================================================
# 24. CHECK KEY HALLMARK NES VALUES
# =============================================================================

get_nes <- function(pathway_name) {

    value <- hallmark_gsea$NES[
        hallmark_gsea$pathway ==
            pathway_name
    ]

    if (length(value) != 1L) {
        return(NA_real_)
    }

    value
}


p53_nes <- get_nes(
    "HALLMARK_P53_PATHWAY"
)

dna_repair_nes <- get_nes(
    "HALLMARK_DNA_REPAIR"
)

g2m_nes <- get_nes(
    "HALLMARK_G2M_CHECKPOINT"
)


if (
    !is.na(p53_nes) &&
    abs(p53_nes - 3.25) > 0.05
) {

    warning(
        "P53 pathway NES differs from expected value: ",
        round(p53_nes, 3)
    )
}


if (
    !is.na(dna_repair_nes) &&
    abs(dna_repair_nes - 2.34) > 0.05
) {

    warning(
        "DNA repair NES differs from expected value: ",
        round(dna_repair_nes, 3)
    )
}


if (
    !is.na(g2m_nes) &&
    abs(g2m_nes - (-2.09)) > 0.05
) {

    warning(
        "G2/M checkpoint NES differs from expected value: ",
        round(g2m_nes, 3)
    )
}


# =============================================================================
# 25. SESSION INFORMATION
# =============================================================================

capture.output(
    sessionInfo(),
    file = "results/GSEA/sessionInfo_GSEA.txt"
)


# =============================================================================
# 26. COMPLETE
# =============================================================================

cat("\nFiles generated under:\n")
cat("  results/GSEA/\n")
cat("  figures/GSEA/\n")


cat("\n============================================================\n")
cat("04 GSEA COMPLETE\n")
cat("============================================================\n")

print(p_gobp)
print(p_kegg)
print(p_integrated)
print (p_p53)
print(p_tnfa)
# Curvas GSEA
print(p_p53)
print(p_dna_repair)
print(p_tnfa)

# Figuras resumen
print(p_gobp)
print(p_kegg)
print(p_integrated)
c(
    P53 = exists("p_p53"),
    DNA_repair = exists("p_dna_repair"),
    TNFa_NFkB = exists("p_tnfa"),
    GO_BP = exists("p_gobp"),
    KEGG = exists("p_kegg"),
    Integrada = exists("p_integrated")
)
list.files(
    "figures/GSEA",
    pattern = "\\.(png|pdf)$",
    full.names = FALSE
)
file.remove(
    "figures/GSEA/Hallmark_selected_NES.png",
    "figures/GSEA/Hallmark_selected_NES.pdf",
    "figures/GSEA/Hallmark_top15_NES.png",
    "figures/GSEA/Hallmark_top15_NES.pdf"
)

list.files("figures/GSEA")
# Figura 8 — P53
print(p_p53)
print(p_p53)
print(p_dna_repair)
print(p_tnfa)
print(p_gobp)
print(p_kegg)
print(p_integrated)
shell.exec(
    normalizePath("figures/GSEA")
)
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
source("scripts/05_sensitivity_analysis.R")
source("scripts/05_sensitivity_analysis.R")
list.files("figures/sensitivity")
# Concordancia global de log2FC
print(p_concordance)

# Concordancia de los 362 DEGs compartidos
print(p_shared)

# Solapamiento de DEGs
print(p_overlap)

# Venn
print(p_venn)
print (p_concordance)
print(p_shared)
print(p_overlap)
print(p_venn)
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
# 27. ADDITIONAL TFM FIGURES
# =============================================================================
#
# Figures:
#   1. DEG direction counts
#   2. log2FC of representative genes
#   3. Expression of representative genes
#   4. Direction consistency across matched CTL-IR comparisons
#   5. Neuronal / non-neuronal marker heatmap
#
# =============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 with install.packages('ggplot2').")
}

library(ggplot2)

dir.create(
    "figures",
    recursive = TRUE,
    showWarnings = FALSE
)


# =============================================================================
# 27.1 NUMBER OF DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

deg_counts <- data.frame(
    Direction = factor(
        c(
            "Mayor expresión en IR",
            "Menor expresión en IR"
        ),
        levels = c(
            "Mayor expresión en IR",
            "Menor expresión en IR"
        )
    ),
    n = c(n_up, n_down)
)


p_deg_counts <- ggplot(
    deg_counts,
    aes(
        x = Direction,
        y = n,
        fill = Direction
    )
) +
    geom_col(
        width = 0.62,
        show.legend = FALSE
    ) +
    geom_text(
        aes(label = n),
        vjust = -0.45,
        size = 5,
        fontface = "bold"
    ) +
    scale_fill_manual(
        values = c(
            "Mayor expresión en IR" = "#00BFC4",
            "Menor expresión en IR" = "#F8766D"
        )
    ) +
    labs(
        title = "Genes diferencialmente expresados tras la irradiación",
        subtitle = sprintf(
            "%d DEGs (FDR < 0,05 y |log2FC| ≥ 1)",
            n_degs
        ),
        x = NULL,
        y = "Número de DEGs"
    ) +
    expand_limits(
        y = max(deg_counts$n) * 1.10
    ) +
    theme_classic(
        base_size = 13
    ) +
    theme(
        plot.title = element_text(
            face = "bold",
            hjust = 0.5
        ),
        plot.subtitle = element_text(
            hjust = 0.5
        ),
        axis.text.x = element_text(
            size = 11
        )
    )


ggsave(
    "figures/DEG_direction_counts.png",
    p_deg_counts,
    width = 7,
    height = 5,
    dpi = 300
)

ggsave(
    "figures/DEG_direction_counts.pdf",
    p_deg_counts,
    width = 7,
    height = 5
)


# =============================================================================
# 27.2 REPRESENTATIVE GENES: log2 FOLD CHANGE
# =============================================================================

representative_symbols <- c(
    "GDF15",
    "CXCL5",
    "CDKN1A",
    "DDB2",
    "GRIA4"
)


representative_genes <- selected_results[
    match(
        representative_symbols,
        selected_results$gene_symbol
    ),
    ,
    drop = FALSE
]


if (any(is.na(representative_genes$gene_symbol))) {
    stop(
        "One or more representative genes could not be recovered."
    )
}


representative_genes$gene_symbol <- factor(
    representative_genes$gene_symbol,
    levels = rev(representative_symbols)
)


representative_genes$functional_category <- c(
    "Daño/estrés",
    "Inflamación",
    "Daño/estrés",
    "Daño/estrés",
    "Función neuronal"
)


representative_genes$functional_category <- factor(
    representative_genes$functional_category,
    levels = c(
        "Daño/estrés",
        "Función neuronal",
        "Inflamación"
    )
)


p_selected_logfc <- ggplot(
    representative_genes,
    aes(
        x = logFC,
        y = gene_symbol,
        fill = functional_category
    )
) +
    geom_col(
        width = 0.72
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    scale_fill_manual(
        values = c(
            "Daño/estrés" = "#F8766D",
            "Función neuronal" = "#00BA38",
            "Inflamación" = "#619CFF"
        )
    ) +
    labs(
        title = "Cambio de expresión en genes representativos",
        subtitle = "IR frente a CTL en neuronas inducidas a 336 h",
        x = "log2 Fold Change (IR vs CTL)",
        y = NULL,
        fill = "Categoría funcional"
    ) +
    theme_classic(
        base_size = 13
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        legend.position = "right"
    )


ggsave(
    "figures/selected_genes_logFC.png",
    p_selected_logfc,
    width = 7.5,
    height = 5,
    dpi = 300
)

ggsave(
    "figures/selected_genes_logFC.pdf",
    p_selected_logfc,
    width = 7.5,
    height = 5
)


# =============================================================================
# 27.3 EXPRESSION OF REPRESENTATIVE GENES
# =============================================================================

expression_plot_df <- selected_expression_long

expression_plot_df$gene_symbol <- factor(
    expression_plot_df$gene_symbol,
    levels = c(
        "CDKN1A",
        "DDB2",
        "GDF15",
        "CXCL5",
        "GRIA4"
    )
)


expression_plot_df$treatment <- factor(
    expression_plot_df$treatment,
    levels = c(
        "CTL",
        "IR"
    )
)


p_selected_expression <- ggplot(
    expression_plot_df,
    aes(
        x = treatment,
        y = expression,
        fill = treatment
    )
) +
    geom_boxplot(
        width = 0.55,
        outlier.shape = NA,
        alpha = 0.75
    ) +
    geom_jitter(
        width = 0.12,
        size = 1.5,
        alpha = 0.75
    ) +
    facet_wrap(
        ~ gene_symbol,
        scales = "free_y",
        nrow = 1
    ) +
    scale_fill_manual(
        values = c(
            CTL = "#F8766D",
            IR = "#00BFC4"
        )
    ) +
    labs(
        title = "Expresión de genes representativos de la respuesta a IR",
        subtitle = "Neuronas inducidas a 336 h",
        x = NULL,
        y = "Expresión normalizada (voom)",
        fill = "Tratamiento"
    ) +
    theme_classic(
        base_size = 12
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold"
        ),
        legend.position = "right"
    )


ggsave(
    "figures/selected_genes_expression.png",
    p_selected_expression,
    width = 10,
    height = 5.5,
    dpi = 300
)

ggsave(
    "figures/selected_genes_expression.pdf",
    p_selected_expression,
    width = 10,
    height = 5.5
)


# =============================================================================
# 27.4 DIRECTION CONSISTENCY ACROSS MATCHED CTL-IR COMPARISONS
# =============================================================================
#
# The metadata are used to identify donor-batch combinations containing
# both CTL and IR libraries.
#
# For each gene:
#
#   difference = IR expression - CTL expression
#
# The percentage indicates the proportion of matched comparisons with the
# same direction as the global differential-expression result.
#
# =============================================================================

consistency_symbols <- c(
    "CDKN1A",
    "DDB2",
    "GDF15",
    "CXCL5",
    "GRIA4",
    "RELN"
)


consistency_ids <- de_all$gene_id[
    match(
        consistency_symbols,
        de_all$gene_symbol
    )
]


names(consistency_ids) <- consistency_symbols


if (any(is.na(consistency_ids))) {

    missing_consistency <- names(
        consistency_ids
    )[is.na(consistency_ids)]

    stop(
        "Genes missing for direction-consistency analysis: ",
        paste(
            missing_consistency,
            collapse = ", "
        )
    )
}


consistency_expression <- v$E[
    consistency_ids,
    ,
    drop = FALSE
]


rownames(consistency_expression) <-
    consistency_symbols


consistency_meta <- data.frame(
    sample = meta$count_column,
    donor = as.character(meta$donor),
    batch = as.character(meta$batch),
    treatment = as.character(meta$treatment),
    stringsAsFactors = FALSE
)


consistency_meta$pair_id <- paste(
    consistency_meta$donor,
    consistency_meta$batch,
    sep = "__"
)


pair_table <- table(
    consistency_meta$pair_id,
    consistency_meta$treatment
)


valid_pairs <- rownames(pair_table)[
    pair_table[, "CTL"] >= 1 &
        pair_table[, "IR"] >= 1
]


pair_differences <- lapply(
    valid_pairs,
    function(pair_id) {

        pair_meta <- consistency_meta[
            consistency_meta$pair_id == pair_id,
            ,
            drop = FALSE
        ]

        ctl_samples <- pair_meta$sample[
            pair_meta$treatment == "CTL"
        ]

        ir_samples <- pair_meta$sample[
            pair_meta$treatment == "IR"
        ]

        ctl_mean <- rowMeans(
            consistency_expression[
                ,
                ctl_samples,
                drop = FALSE
            ]
        )

        ir_mean <- rowMeans(
            consistency_expression[
                ,
                ir_samples,
                drop = FALSE
            ]
        )

        ir_mean - ctl_mean
    }
)


pair_differences <- do.call(
    cbind,
    pair_differences
)


colnames(pair_differences) <- valid_pairs


global_logfc <- de_all$logFC[
    match(
        consistency_symbols,
        de_all$gene_symbol
    )
]


same_direction <- matrix(
    FALSE,
    nrow = length(consistency_symbols),
    ncol = ncol(pair_differences)
)


for (i in seq_along(consistency_symbols)) {

    if (global_logfc[i] > 0) {

        same_direction[i, ] <-
            pair_differences[i, ] > 0

    } else {

        same_direction[i, ] <-
            pair_differences[i, ] < 0
    }
}


n_same <- rowSums(
    same_direction,
    na.rm = TRUE
)


n_pairs <- ncol(pair_differences)


consistency_df <- data.frame(
    gene_symbol = consistency_symbols,
    n_same = n_same,
    n_pairs = n_pairs,
    percent = 100 * n_same / n_pairs,
    global_direction = ifelse(
        global_logfc > 0,
        "Mayor expresión en IR",
        "Menor expresión en IR"
    ),
    stringsAsFactors = FALSE
)


consistency_df$gene_symbol <- factor(
    consistency_df$gene_symbol,
    levels = rev(consistency_symbols)
)


consistency_df$label <- paste0(
    consistency_df$n_same,
    "/",
    consistency_df$n_pairs
)


p_direction_consistency <- ggplot(
    consistency_df,
    aes(
        x = percent,
        y = gene_symbol,
        fill = global_direction
    )
) +
    geom_col(
        width = 0.65
    ) +
    geom_text(
        aes(label = label),
        hjust = 1.15,
        fontface = "bold",
        size = 3.8
    ) +
    scale_fill_manual(
        values = c(
            "Mayor expresión en IR" = "#00BFC4",
            "Menor expresión en IR" = "#F8766D"
        )
    ) +
    scale_x_continuous(
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100)
    ) +
    labs(
        title = "Consistencia de la dirección del cambio de expresión",
        subtitle = paste0(
            n_pairs,
            " comparaciones CTL–IR correspondientes a combinaciones donante–lote"
        ),
        x = "Comparaciones con la misma dirección del cambio (%)",
        y = NULL,
        fill = "Dirección global"
    ) +
    theme_classic(
        base_size = 12
    ) +
    theme(
        plot.title = element_text(
            face = "bold",
            hjust = 0.5
        ),
        plot.subtitle = element_text(
            hjust = 0.5
        ),
        legend.position = "top"
    )


ggsave(
    "figures/gene_direction_consistency.png",
    p_direction_consistency,
    width = 8,
    height = 6,
    dpi = 300
)

ggsave(
    "figures/gene_direction_consistency.pdf",
    p_direction_consistency,
    width = 8,
    height = 6
)


write.csv(
    consistency_df,
    "results/gene_direction_consistency.csv",
    row.names = FALSE
)


# =============================================================================
# 27.5 NEURONAL / NON-NEURONAL MARKER HEATMAP
# =============================================================================

marker_symbols <- c(
    # Neuronal markers
    "MAP2",
    "RBFOX3",
    "TUBB3",
    "NEFL",
    "NEFM",
    "SNAP25",
    "SYN1",
    "SYT1",
    "DLG4",

    # Extracellular matrix / non-neuronal-associated markers
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "VIM",
    "FN1",

    # Representative response genes
    "MMP3",
    "CXCL5",
    "PSG1"
)


marker_ids <- de_all$gene_id[
    match(
        marker_symbols,
        de_all$gene_symbol
    )
]


names(marker_ids) <- marker_symbols


available_markers <- marker_symbols[
    !is.na(marker_ids)
]


missing_markers <- marker_symbols[
    is.na(marker_ids)
]


if (length(missing_markers) > 0) {

    warning(
        "Markers not found and omitted from heatmap: ",
        paste(
            missing_markers,
            collapse = ", "
        )
    )
}


marker_expression <- v$E[
    marker_ids[
        available_markers
    ],
    ,
    drop = FALSE
]


rownames(marker_expression) <-
    available_markers


marker_z <- t(
    scale(
        t(marker_expression)
    )
)


marker_annotation <- data.frame(
    Condición = meta$treatment
)


rownames(marker_annotation) <-
    meta$count_column


marker_annotation <- marker_annotation[
    colnames(marker_z),
    ,
    drop = FALSE
]


marker_annotation_colors <- list(
    Condición = c(
        CTL = "#F8766D",
        IR = "#00BFC4"
    )
)


if (requireNamespace("pheatmap", quietly = TRUE)) {

    pheatmap::pheatmap(
        marker_z,
        annotation_col = marker_annotation,
        annotation_colors = marker_annotation_colors,
        cluster_rows = FALSE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        show_rownames = TRUE,
        fontsize_row = 9,
        border_color = NA,
        main = paste(
            "Marcadores neuronales y genes asociados",
            "a componentes no neuronales"
        ),
        filename =
            "figures/neuronal_non_neuronal_markers.png",
        width = 11,
        height = 6.5
    )


    pheatmap::pheatmap(
        marker_z,
        annotation_col = marker_annotation,
        annotation_colors = marker_annotation_colors,
        cluster_rows = FALSE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        show_rownames = TRUE,
        fontsize_row = 9,
        border_color = NA,
        main = paste(
            "Marcadores neuronales y genes asociados",
            "a componentes no neuronales"
        ),
        filename =
            "figures/neuronal_non_neuronal_markers.pdf",
        width = 11,
        height = 6.5
    )

} else {

    warning(
        "Package 'pheatmap' is not installed. ",
        "Marker heatmap was skipped."
    )
}


# =============================================================================
# 27.6 FIGURE CHECK
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("ADDITIONAL TFM FIGURES\n")
cat("============================================================\n\n")

cat(
    sprintf(
        "DEG counts: %d up in IR | %d down in IR\n",
        n_up,
        n_down
    )
)

cat(
    "Matched donor-batch comparisons:",
    n_pairs,
    "\n\n"
)

cat("Figures generated:\n")

additional_figures <- c(
    "DEG_direction_counts",
    "selected_genes_logFC",
    "selected_genes_expression",
    "gene_direction_consistency",
    "neuronal_non_neuronal_markers"
)

for (fig in additional_figures) {

    png_file <- file.path(
        "figures",
        paste0(fig, ".png")
    )

    pdf_file <- file.path(
        "figures",
        paste0(fig, ".pdf")
    )

    cat(
        sprintf(
            "  %-40s PNG: %-5s | PDF: %-5s\n",
            fig,
            ifelse(file.exists(png_file), "OK", "MISSING"),
            ifelse(file.exists(pdf_file), "OK", "MISSING")
        )
    )
}

cat("\n============================================================\n")
cat("ADDITIONAL FIGURES COMPLETE\n")
cat("============================================================\n")
source("scripts/03_differential_expression.R")
list.files(
    "figures",
    pattern = "\\.(png|pdf)$"
)
print(p_deg_counts)
print(p_selected_logfc)
print(p_selected_expression)
print(p_direction_consistency)
# =============================================================================
# 27.7 FORCE FIGURES TO DISPLAY IN RSTUDIO
# =============================================================================

cat("\nShowing additional TFM figures in RStudio Plots...\n")

# Figure 1 — DEG counts
if (exists("p_deg_counts")) {
    print(p_deg_counts)
}

# Figure 2 — Representative genes log2FC
if (exists("p_selected_logfc")) {
    print(p_selected_logfc)
}

# Figure 3 — Representative-gene expression
if (exists("p_selected_expression")) {
    print(p_selected_expression)
}

# Figure 4 — Direction consistency
if (exists("p_direction_consistency")) {
    print(p_direction_consistency)
}

# Figure 5 — Neuronal / non-neuronal marker heatmap
if (
    exists("marker_z") &&
    exists("marker_annotation") &&
    requireNamespace("pheatmap", quietly = TRUE)
) {

    p_marker_heatmap <- pheatmap::pheatmap(
        marker_z,
        annotation_col = marker_annotation,
        annotation_colors = marker_annotation_colors,
        cluster_rows = FALSE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        show_rownames = TRUE,
        fontsize_row = 9,
        border_color = NA,
        main = paste(
            "Marcadores neuronales y genes asociados",
            "a componentes no neuronales"
        ),
        silent = TRUE
    )

    grid::grid.newpage()
    grid::grid.draw(p_marker_heatmap$gtable)
}

cat("\nAll additional figures sent to the Plots device.\n")
# =============================================================================
# INTERACTIVE FIGURE VIEWER
# =============================================================================

if (interactive()) {

    cat("\n")
    cat("============================================================\n")
    cat("FIGURE VIEWER\n")
    cat("Press ENTER in the Console to show the next figure.\n")
    cat("============================================================\n")

    readline("\nENTER -> DEG direction counts")
    print(p_deg_counts)

    readline("\nENTER -> Representative genes log2FC")
    print(p_selected_logfc)

    readline("\nENTER -> Representative-gene expression")
    print(p_selected_expression)

    readline("\nENTER -> Direction consistency")
    print(p_direction_consistency)

    readline("\nENTER -> Neuronal/non-neuronal marker heatmap")

    p_marker_heatmap <- pheatmap::pheatmap(
        marker_z,
        annotation_col = marker_annotation,
        annotation_colors = marker_annotation_colors,
        cluster_rows = FALSE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        show_rownames = TRUE,
        fontsize_row = 9,
        border_color = NA,
        main = paste(
            "Marcadores neuronales y genes asociados",
            "a componentes no neuronales"
        ),
        silent = TRUE
    )

    grid::grid.newpage()
    grid::grid.draw(p_marker_heatmap$gtable)

    cat("\nAll figures displayed.\n")
}
print(p_deg_counts)
print(p_deg_counts)
print(p_selected_logfc)
print(p_selected_expression)
print(p_direction_consistency)
# =============================================================================
# 28. FORCE ALL TFM FIGURES TO RSTUDIO PLOTS
# =============================================================================

# Close file graphics devices that may still be open
while (!is.null(dev.list())) {
    dev.off()
}

# Open a normal RStudio graphics device
if (Sys.getenv("RSTUDIO") == "1") {
    dev.new()
}

# -------------------------------------------------------------------------
# 1. DEG direction
# -------------------------------------------------------------------------

print(p_deg_counts)

# -------------------------------------------------------------------------
# 2. Representative log2FC
# -------------------------------------------------------------------------

print(p_selected_logfc)

# -------------------------------------------------------------------------
# 3. Representative expression
# -------------------------------------------------------------------------

print(p_selected_expression)

# -------------------------------------------------------------------------
# 4. Direction consistency
# -------------------------------------------------------------------------

print(p_direction_consistency)

# -------------------------------------------------------------------------
# 5. Marker heatmap
# -------------------------------------------------------------------------

p_marker_heatmap <- pheatmap::pheatmap(
    marker_z,
    annotation_col = marker_annotation,
    annotation_colors = marker_annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 9,
    border_color = NA,
    main = paste(
        "Marcadores neuronales y genes asociados",
        "a componentes no neuronales"
    )
)

cat("\n")
cat("============================================================\n")
cat("ALL TFM FIGURES SENT TO THE RSTUDIO PLOTS PANEL\n")
cat("Use the <- and -> arrows in Plots to move between figures.\n")
cat("============================================================\n")

# =============================================================================
# 28. FORCE ALL TFM FIGURES TO RSTUDIO PLOTS
# =============================================================================

# Close file graphics devices that may still be open
while (!is.null(dev.list())) {
    dev.off()
}

# Open a normal RStudio graphics device
if (Sys.getenv("RSTUDIO") == "1") {
    dev.new()
}

# -------------------------------------------------------------------------
# 1. DEG direction
# -------------------------------------------------------------------------

print(p_deg_counts)

# -------------------------------------------------------------------------
# 2. Representative log2FC
# -------------------------------------------------------------------------

print(p_selected_logfc)

# -------------------------------------------------------------------------
# 3. Representative expression
# -------------------------------------------------------------------------

print(p_selected_expression)

# -------------------------------------------------------------------------
# 4. Direction consistency
# -------------------------------------------------------------------------

print(p_direction_consistency)

# -------------------------------------------------------------------------
# 5. Marker heatmap
# -------------------------------------------------------------------------

p_marker_heatmap <- pheatmap::pheatmap(
    marker_z,
    annotation_col = marker_annotation,
    annotation_colors = marker_annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 9,
    border_color = NA,
    main = paste(
        "Marcadores neuronales y genes asociados",
        "a componentes no neuronales"
    )
)

cat("\n")
cat("============================================================\n")
cat("ALL TFM FIGURES SENT TO THE RSTUDIO PLOTS PANEL\n")
cat("Use the <- and -> arrows in Plots to move between figures.\n")
cat("============================================================\n")
list.files(
    "figures",
    pattern = "\\.(png|pdf)$"
)
print(p_concordance)
print(p_shared)
print(p_overlap)
print(p_venn)
list.files("figures/sensitivity")
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
print(p_overlap)
print(p_venn)
print(p_concordance)
print(p_shared)
normalizePath("figures/sensitivity")
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
list.files("figures/original_study_comparison")
print(p_genes)
print(p_hallmark)
shell.exec(normalizePath("figures/sensitivity"))
cat("sensitivity:", nrow(sensitivity), "\n")
cat("comparison:", nrow(comparison), "\n")

cat("\nNA en modelo principal:\n")
print(colSums(is.na(comparison)))

cat("\nDimensiones:\n")
print(dim(sensitivity))
print(dim(comparison))
files <- list.files(
    "scripts",
    pattern = "\\.(R|r|Rmd)$",
    full.names = TRUE,
    recursive = TRUE
)

hits <- lapply(files, function(f) {
    x <- readLines(f, warn = FALSE)
    i <- grep("21713|21\\.713", x)
    if (length(i)) {
        data.frame(
            file = f,
            line = i,
            text = x[i]
        )
    }
})

do.call(rbind, hits)
list.files("results", recursive = TRUE)
# Buscar 21713 / 21.713 en TODOS los archivos del proyecto

files <- list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
)

files <- files[file.info(files)$isdir == FALSE]

hits <- lapply(files, function(f) {

    txt <- tryCatch(
        readLines(f, warn = FALSE),
        error = function(e) NULL
    )

    if (is.null(txt)) return(NULL)

    idx <- grep(
        "21713|21\\.713",
        txt
    )

    if (length(idx) > 0) {
        data.frame(
            archivo = f,
            linea = idx,
            contenido = txt[idx],
            stringsAsFactors = FALSE
        )
    }
})
hits <- do.call(rbind, hits)

print(hits)
# ============================================================
# BUSCAR EL ORIGEN DOCUMENTAL DE 21.713
# ============================================================

files_text <- list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(R|Rmd|md|txt|log)$",
    ignore.case = TRUE
)

hits_21713 <- lapply(files_text, function(f) {

    x <- tryCatch(
        readLines(f, warn = FALSE),
        error = function(e) NULL
    )

    if (is.null(x)) return(NULL)

    i <- grep(
        "21713|21\\.713",
        x
    )

    if (length(i) == 0) return(NULL)

    data.frame(
        archivo = f,
        linea = i,
        texto = x[i],
        stringsAsFactors = FALSE
    )
})

hits_21713 <- Filter(Negate(is.null), hits_21713)

if (length(hits_21713) == 0) {
    cat("No se encontró 21.713 en scripts/documentación.\n")
} else {
    print(
        do.call(rbind, hits_21713),
        row.names = FALSE
    )
}
cat("\n===== SCRIPT 06 =====\n")

x06 <- readLines(
    "scripts/06_original_study_comparison.R",
    warn = FALSE
)
history(max.show = Inf)
files <- list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(R|Rmd)$"
)

for (f in files) {
    x <- readLines(f, warn = FALSE)

    i <- grep(
        "Hughes|intersect|merge|common",
        x,
        ignore.case = TRUE
    )

    if (length(i) > 0) {
        cat("\n---", f, "---\n")
        cat(x[i], sep = "\n")
    }
}
cat("\n===== ARCHIVOS RELACIONADOS CON HUGHES / GSEA / ORIGINAL =====\n\n")

all_files <- list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE
)

grep(
    "hughes|original|gsea|gene.set|deg|differential",
    all_files,
    value = TRUE,
    ignore.case = TRUE
)
cat("\n===== OBJETOS DEL SCRIPT 06 =====\n\n")

source("scripts/06_original_study_comparison.R")
intersect(genes_hughes, genes_chaimae)
