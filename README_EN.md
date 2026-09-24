[🇪🇸 Versión en español](README.md)

# Master's Thesis - Transcriptomic Analysis of the Response of Human Neurons to Ionizing Radiation-Induced DNA Damage

[🇪🇸 Versión en español](README.md)

## Description

This repository contains the scripts, graphical results, and documentation associated with the Master's Thesis developed as part of the Master's Degree in Bioinformatics at Universidad Internacional de Valencia (VIU).

The project consists of a transcriptomic reanalysis of RNA-seq data from human induced neurons exposed to ionizing radiation, with the aim of characterizing gene expression changes and biological processes associated with the DNA damage response.

The analysis includes quality control, differential gene expression analysis, functional enrichment using Gene Set Enrichment Analysis (GSEA), statistical sensitivity analysis, qualitative comparison with the original study, and an exploratory structural characterization of MMP3.

## Dataset

The data were obtained from the publicly available GSE329677 dataset deposited in the NCBI Gene Expression Omnibus (GEO).

The main analysis focuses on human induced neurons (iNs) corresponding to control (CTL) and irradiated (IR) conditions at 336 hours.

## Analysis Workflow

The analysis workflow comprises:

1. Data and metadata preparation and verification.
2. Quality control of the RNA-seq libraries.
3. Selection of samples corresponding to the main analysis.
4. Gene expression normalization and statistical modelling using limma-voom.
5. Differential expression analysis between IR and CTL.
6. Functional enrichment analysis using GSEA with Hallmark, Gene Ontology Biological Process, and KEGG gene-set collections.
7. Sensitivity analysis using an alternative statistical model specification.
8. Qualitative comparison with the results of the original study.
9. Exploratory structural characterization of MMP3 and the SPI ligand.

## Differential Expression Analysis

The main differential expression analysis was performed using limma-voom.

The statistical model included batch as a fixed effect and accounted for the dependence between observations from the same donor using `duplicateCorrelation`.

Differentially expressed genes (DEGs) were defined using the following criteria:

- FDR < 0.05
- |log2FC| ≥ 1

Using these criteria, 466 differentially expressed genes were identified:

- 220 genes showed higher expression in IR.
- 246 genes showed lower expression in IR.

## Functional Enrichment Analysis

Gene Set Enrichment Analysis (GSEA) was performed using a ranked list of 17,919 genes.

Three functional gene-set collections were analyzed:

- Hallmark
- Gene Ontology Biological Process (GO:BP)
- KEGG

The results showed enrichment toward the irradiated condition of processes associated with the DNA damage response, p53 signaling, DNA repair, inflammatory signaling, and oxidative and metabolic processes.

In contrast, relative enrichment toward the control condition was observed for processes associated with neuronal and synaptic functions.

## Sensitivity Analysis

The robustness of the differential expression results was evaluated using an alternative statistical model specification.

The main model identified 466 DEGs, whereas the sensitivity model identified 423 DEGs.

A total of 362 DEGs were shared between both models, representing:

- 77.7% of the DEGs identified by the main model.
- 85.6% of the DEGs identified by the sensitivity model.

All 362 shared DEGs retained the same direction of differential expression.

The genome-wide log2FC correlation between the two models was:

- r = 0.966

Among the 362 shared DEGs, the log2FC correlation was:

- r = 0.996

## Comparison with the Original Study

The results of the reanalysis were qualitatively compared with those reported in the original study.

Because of differences in statistical methodology, DEG selection criteria, and functional enrichment strategies, the comparison was primarily performed at the level of transcriptomic trends and biological processes rather than as an exact reproduction of the original DEG lists.

## Exploratory Structural Characterization of MMP3

As an exploratory extension of the transcriptomic analysis, MMP3 was selected for structural characterization using the experimental human MMP3 structure PDB 1D8F and the SPI ligand.

Molecular docking was performed using AutoDock Vina via SwissDock, and structural inspection was carried out using UCSF Chimera.

Nine docking poses were generated. The best-scoring pose had an AutoDock Vina score of -9.032 kcal/mol and was selected for subsequent structural inspection.

The selected pose was examined in relation to residues located within the MMP3 binding site and the catalytic Zn801 ion.

This structural characterization was performed independently of the R-based RNA-seq pipeline and should be interpreted as an exploratory in silico analysis. Docking scores alone do not constitute experimental evidence of binding affinity or inhibitory activity.

Structural documentation is available in:

`figures/MMP3_docking/`


## Repository Structure

```text
TFM/
├── scripts/
│   ├── 01_...
│   ├── 02_...
│   ├── 03_differential_expression.R
│   ├── 04_...
│   ├── 05_sensitivity_analysis.R
│   └── 06_original_study_comparison.R
│
├── figures/
│   ├── QC/
│   ├── GSEA/
│   ├── sensitivity/
│   ├── original_study_comparison/
│   └── MMP3_docking/
│
├── README.md
└── README_EN.md
```

## Reproducibility

The scripts included in the `scripts/` directory document the transcriptomic analysis workflow used in this Master's Thesis.

Figures generated during the different stages of the analysis are organized in the `figures/` directory.

Package dependencies required for each stage of the analysis are specified in the corresponding scripts.

The exploratory structural characterization of MMP3 was performed using external structural bioinformatics tools and is therefore documented separately from the R-based transcriptomic pipeline.

## Main Tools

The transcriptomic analysis was primarily performed in R. Specific package dependencies are documented in the corresponding analysis scripts.

Main tools and packages include:

- R
- limma
- edgeR
- fgsea
- msigdbr
- ggplot2
- pheatmap
- UCSF Chimera
- AutoDock Vina / SwissDock

## Author

**Chaimae Bahraoui Badouch**  
Master's Degree in Bioinformatics  
Universidad Internacional de Valencia (VIU)





