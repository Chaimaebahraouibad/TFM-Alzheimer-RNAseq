# Exploratory structural characterization of MMP3

This directory contains the structural outputs from the exploratory molecular docking analysis of human MMP3 with the SPI ligand.

The analysis used the experimental human MMP3 structure PDB 1D8F, which contains SPI as a co-crystallized ligand. Chain B was selected as the receptor and SPI was extracted from the same structure for docking.

Molecular docking was performed using AutoDock Vina via SwissDock. The search space was defined around the known SPI binding site.

Nine docking poses were generated, with scores ranging from -9.032 to -6.944 kcal/mol. The best-scoring pose (-9.032 kcal/mol) was selected for subsequent structural inspection.

The selected pose was examined using UCSF Chimera. Interactions involving LEU664 and ALA665 were inspected, together with the position of the ligand relative to the catalytic Zn801. A distance of 2.297 Å was measured between Zn801 and the OB atom of SPI.

This structural analysis was performed independently of the R-based RNA-seq pipeline and should be interpreted as an exploratory in silico characterization. Docking scores alone do not constitute experimental evidence of binding affinity or inhibitory activity.
