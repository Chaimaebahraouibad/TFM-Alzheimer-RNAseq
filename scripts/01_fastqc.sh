#!/bin/bash

# Control de calidad de las lecturas RNA-seq mediante FastQC

mkdir -p fastqc

fastqc -t 8 -o fastqc *.fastq.gz
