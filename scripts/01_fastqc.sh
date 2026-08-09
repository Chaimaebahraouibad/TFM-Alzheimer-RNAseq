#!/bin/bash

# Quality control of RNA-seq reads using FastQC

mkdir -p fastqc

fastqc -t 8 -o fastqc *.fastq.gz
