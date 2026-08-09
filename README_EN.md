# TFM - Bioinformatics Analysis of Alzheimer's Disease

## Description

This repository contains the scripts and documentation associated with the Master's Thesis at the Universidad Internacional de Valencia, focused on the bioinformatic analysis of RNA-seq data related to Alzheimer's disease.

## Dataset

The data were obtained from the **GSE329677** dataset, available through NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329677

## Workflow

1. Sample selection, including a total of 16 samples at 0h (control) and 336h (irradiated and control).
2. Quality control using FastQC.
3. Overall quality assessment using MultiQC.
4. Read preprocessing.
5. Alignment to the reference genome.
6. Gene expression quantification.
7. Differential gene expression analysis.
8. Biological interpretation of the results.

## Analyses included:

* Sample selection
* Quality control of RNA-seq reads using FastQC

## Dataset:

The analyses were performed using RNA-seq data from the **GSE329677** dataset, available through NCBI GEO.

## Scripts:

01_fastqc.sh — Quality control of RNA-seq reads using FastQC

## Tools

* Linux
* Conda
* FastQC 0.12.1
* MultiQC 1.19

## Author

Chaimae Bahraoui Badouch
Master's Degree in Bioinformatics
Universidad Internacional de Valencia (VIU)
