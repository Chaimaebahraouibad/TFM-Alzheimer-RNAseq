# TFM - Análisis bioinformático de RNA-seq

## Descripción

Este repositorio contiene los scripts y documentación asociados al Trabajo Fin de Máster, centrado en el análisis bioinformático de datos de RNA-seq.

## Dataset

Los datos proceden del dataset GSE329677 disponible en NCBI GEO.

## Flujo de trabajo

1. Selección de las muestras.
2. Control de calidad mediante FastQC.
3. Evaluación conjunta mediante MultiQC.
4. Preprocesamiento de las lecturas, si fuese necesario.
5. Alineamiento contra el genoma de referencia.
6. Cuantificación de la expresión génica.
7. Análisis de expresión diferencial.
8. Interpretación biológica de los resultados.

## Herramientas

- Linux
- Conda
- FastQC 0.12.1
- MultiQC 1.19
