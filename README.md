[🇬🇧 English version](README_EN.md)

# TFM - Análisis bioinformático sobre el Alzheimer

## Descripción

Este repositorio contiene los scripts y documentación asociados al Trabajo Fin de Máster asociado a la UNIVERSIDAD INTERNACIONAL DE VALENCIA, centrado en el análisis bioinformático de datos de RNA-seq, sobre Alzheimer. 

## Dataset

Los datos proceden del dataset GSE329677 disponible en NCBI GEO. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329677

## Flujo de trabajo

1. Selección de las muestras, un total de 16 muestras 0h (control) y 336h (irradiadas y control)
2. Control de calidad mediante FastQC.
3. Evaluación conjunta mediante MultiQC.
4. Preprocesamiento de las lecturas.
6. Alineamiento contra el genoma de referencia.
7. Cuantificación de la expresión génica.
8. Análisis de expresión diferencial.
9. Interpretación biológica de los resultados.

## Análisis incluidos:
Selección de muestras
Control de calidad de las lecturas de RNA-seq mediante FastQC

## Dataset:
Los análisis se realizaron utilizando datos de RNA-seq del conjunto de datos GSE329677, disponible a través de NCBI GEO.

## Scripts:
01_fastqc.sh — Control de calidad de las lecturas de RNA-seq mediante FastQC

## Herramientas

- Linux
- Conda
- FastQC 0.12.1
- MultiQC 1.19

## Autora

Chaimae Bahraoui Badouch
Máster en Bioinformática
Universidad Internacional de Valencia (VIU)
