# ChIP-seq-standard

This repository contains a customised data analysis pipeline that facilitates simultaneous analysis of ChIP-seq data at ribosomal DNA (rDNA) repeats and genome-wide (For paired-end short read - Illumina data). 

Table of contents 
1. Overview of the pipeline
2. Software installation
3. How to run - example dataset
  
     3.1 Quality check of raw fastq files - FastQC/MultiQC
  
     3.2 Read alignment - Bowtie2
  
     3.3 Processing aligned files - Samtools
  
     3.4 Quality check of BAM files - BamQC/MultiQC and deepTools
  
     3.5 Peak calling (rDNA) - Bamtools and MACS2
  
     3.6 Peak calling (genome-wide) - MACS2
  
     3.7 Differential binding analysis - R Bioconductor package 'diffbind'
 
     3.8 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'
  
     3.9 Coverage track generation and visualisation - deepTools and IGV

1. Overview of the pipeline 

[Picture 1.pdf](https://github.com/tudumanne/ChIP-seq-standard/files/7828134/Picture.1.pdf)
