[🇪🇸 Versión en español](README.md)

# TFM - Análisis transcriptómico de la respuesta de neuronas humanas al daño en el ADN inducido por radiación ionizante

[🇬🇧 English version](README_EN.md)

## Descripción

Este repositorio contiene los scripts, resultados gráficos y documentación asociados al Trabajo Fin de Máster del Máster Universitario en Bioinformática de la Universidad Internacional de Valencia (VIU).

El trabajo realiza un reanálisis transcriptómico de datos de RNA-seq de neuronas humanas inducidas expuestas a radiación ionizante, con el objetivo de caracterizar los cambios de expresión génica y los procesos biológicos asociados a la respuesta al daño en el ADN.

El análisis incluye control de calidad, expresión diferencial, enriquecimiento funcional mediante GSEA, análisis de sensibilidad del modelo estadístico, comparación cualitativa con el estudio original y una caracterización estructural exploratoria de MMP3.

## Dataset

Los datos proceden del conjunto GSE329677, disponible públicamente en NCBI Gene Expression Omnibus (GEO).

El análisis principal se centra en neuronas humanas inducidas (iNs) correspondientes a las condiciones control (CTL) e irradiada (IR) a las 336 horas.

## Flujo de trabajo

El flujo de análisis comprende:

1. Preparación y comprobación de los datos y metadatos.
2. Control de calidad de las bibliotecas de RNA-seq.
3. Selección de las muestras correspondientes al análisis principal.
4. Normalización y modelado de la expresión génica mediante limma-voom.
5. Análisis de expresión diferencial entre IR y CTL.
6. Análisis de enriquecimiento funcional mediante GSEA utilizando las colecciones Hallmark, Gene Ontology Biological Process y KEGG.
7. Análisis de sensibilidad mediante una especificación estadística alternativa.
8. Comparación cualitativa con los resultados del estudio original.
9. Caracterización estructural exploratoria de MMP3 y el ligando SPI.

## Análisis de expresión diferencial

El análisis principal se realizó mediante limma-voom.

El modelo estadístico incorporó el lote como efecto fijo y tuvo en cuenta la dependencia entre observaciones procedentes del mismo donante mediante `duplicateCorrelation`.

Los genes diferencialmente expresados se definieron utilizando:

- FDR < 0,05
- |log2FC| ≥ 1

Con estos criterios se identificaron 466 genes diferencialmente expresados:

- 220 con mayor expresión en IR.
- 246 con menor expresión en IR.

## Enriquecimiento funcional

Se realizó Gene Set Enrichment Analysis (GSEA) sobre una lista ordenada de 17.919 genes.

Se analizaron tres colecciones funcionales:

- Hallmark
- Gene Ontology Biological Process (GO:BP)
- KEGG

Los resultados mostraron enriquecimiento hacia la condición irradiada de procesos relacionados con la respuesta al daño del ADN, señalización de p53, reparación del ADN, respuesta inflamatoria y procesos metabólicos y oxidativos.

Por el contrario, hacia la condición control se observó enriquecimiento relativo de procesos relacionados con funciones neuronales y sinápticas.

## Análisis de sensibilidad

Se evaluó la robustez de los resultados mediante una especificación estadística alternativa.

El modelo principal identificó 466 DEGs y el modelo de sensibilidad 423 DEGs, de los cuales 362 fueron compartidos.

Los 362 genes compartidos conservaron la misma dirección del cambio de expresión.

## Comparación con el estudio original

Los resultados del reanálisis se compararon cualitativamente con los del estudio original.

Debido a las diferencias en los métodos estadísticos, los criterios de selección de DEGs y las estrategias de enriquecimiento funcional, la comparación se planteó principalmente a nivel de tendencias transcriptómicas y procesos biológicos y no como una reproducción exacta de las listas de genes diferencialmente expresados.

## Caracterización estructural exploratoria de MMP3

Como extensión exploratoria del análisis transcriptómico, se realizó una caracterización estructural de MMP3 utilizando la estructura experimental humana PDB 1D8F y el ligando SPI.

El docking molecular se realizó mediante AutoDock Vina a través de SwissDock y la inspección estructural mediante UCSF Chimera.

Esta caracterización se realizó de forma independiente del pipeline de RNA-seq en R y debe interpretarse como un análisis exploratorio *in silico*, no como evidencia experimental de afinidad o actividad inhibitoria.

La documentación y la figura estructural se encuentran en:

`figures/MMP3_docking/`




## Estructura del repositorio

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


## Herramientas principales

El análisis transcriptómico se realizó principalmente en R. Las dependencias específicas utilizadas en cada etapa se indican en los scripts correspondientes.

Entre las principales herramientas y paquetes utilizados se encuentran:

- R
- limma
- edgeR
- fgsea
- msigdbr
- ggplot2
- pheatmap
- UCSF Chimera
- AutoDock Vina / SwissDock

## Autora

Chaimae Bahraoui Badouch
Máster Universitario en Bioinformática
Universidad Internacional de Valencia (VIU)
