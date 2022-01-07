# ChIP-seq-standard

This repository contains a customised data analysis pipeline that facilitates simultaneous analysis of ChIP-seq data at ribosomal DNA (rDNA) repeats and genome-wide (For paired-end short read - Illumina data). 

### Table of contents 
1. Outline
2. Software installation
3. How to run - example dataset
  
     3.1 Quality check of raw fastq files - FastQC/MultiQC
  
     3.2 Read alignment, processing and post-alignment quality check - Bowtie2, Samtools and BamQC/MultiQC

     3.3 Coverage track generation and visualisation - deepTools and IGV
  
     3.4 Peak calling (rDNA) - Bamtools and MACS2
  
     3.5 Peak calling (genome-wide) - MACS2
  
     3.6 Differential binding analysis - R Bioconductor package 'diffbind'
 
     3.7 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'
  

### 1. Outline

This pipeline was run on an HPC (high-performance computing) system based on CentOS (Linux). R based analysis was carried out in RStudio/MacOS Catalina.
The 'Scripts' folder contains template bash scripts and example Snakemake workflows.

![alt text](https://github.com/tudumanne/ChIP-seq-standard/files/7828134/Picture.1.pdf)

2. Software installation 

The required software/command-line tools were installed via conda on Linux. 
(Miniconda https://docs.conda.io/en/latest/miniconda.html)


3. How to run - example dataset

The folder 'example dataset' contains 9 ChIPed and 9 input control samples, 3 biological replicates per each stage (WT, PreM and Mal).
These files contain a subset of reads sequenced on MiSeq platform.  
  
3.1 Quality check of raw fastq files - FastQC/MultiQC
  
3.2 Read alignment, processing and post-alignment quality check - Bowtie2, Samtools and BamQC/MultiQC

3.3 Coverage track generation and visualisation - deepTools and IGV
  
3.4 Peak calling (rDNA) - Bamtools and MACS2
  
3.5 Peak calling (genome-wide) - MACS2
  
3.6 Differential binding analysis - R Bioconductor package 'diffbind'
 
3.7 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'

