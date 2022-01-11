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

First part of the pipeline (3.1-3.5) was run on an HPC (high-performance computing) system based on CentOS (Linux). R based analyses (3.6-3.7) were carried out in RStudio installed on MacOS Catalina.

The 'Scripts' folder contains template bash scripts and example Snakemake workflows.

![alt text](https://github.com/tudumanne/ChIP-seq-standard/files/7828134/Picture.1.pdf)

2. Software installation 

The required software/command-line tools were installed via conda on Linux. 
(Miniconda https://docs.conda.io/en/latest/miniconda.html)

conda env create -f environment.yml

3. How to run - example dataset

The folder 'example dataset' contains 9 ChIPed (H3K4me3) and 9 input control samples, 3 biological replicates per each stage (WT, PreM and Mal).
These files contain a subset of reads sequenced on MiSeq platform.  
  
3.1 Quality check of raw fastq files - FastQC/MultiQC

fastqc -o fastqc --extract --dir fastqc_chip --format fastq h3k4me3_*.fastq.gz
fastqc -o fastqc --extract --dir fastqc_input --format fastq input_*.fastq.gz

multiqc fastqc_chip/

3.2 Read alignment, processing and post-alignment quality check - Bowtie2, Samtools and BamQC/MultiQC

bowtie2 -x index -1 input_{i}_rep{i}_R1.fastq.gz -2 input_{i}_rep{i}_R2.fastq.gz -S input_{i}_rep{i}.sam
samtools view -S -b H3K9ac_WT_rep01.sam > H3K9ac_WT_rep01.bam
samtools sort H3K9ac_WT_rep01.bam > H3K9ac_WT_rep01_sorted.bam
samtools index H3K9ac_WT_rep01_sorted.bam

3.3 Coverage track generation and visualisation - deepTools and IGV

multiBamSummary bins --bamfiles H3K9ac_*_sorted.bam -o H3K9ac_bam.npz --binSize 1000 --extendReads --centerReads

plotCorrelation --corData H3K9ac_rdna_bam.npz --corMethod spearman --whatToPlot heatmap --plotFile correlation_heatmap_H3K9ac_rdna.pdf --skipZeros --plotTitle "Spearman Correlation of read counts - H3K9ac (rDNA only)" --plotFileFormat pdf --plotNumbers
plotPCA --corData H3K9ac_rdna_bam.npz --plotFile pca_H3K9ac_rdna.pdf --plotTitle "PCA plot - H3K9ac (rDNA only)" --plotFileFormat pdf

plotFingerprint -b H2A.Z_WT_rep01_sorted.bam -plot fingerprint_WT_H2A.Z.pdf --extendReads --centerReads Input --plotFileFormat pdf --plotTitle "WT - H2A.Z"


computeMatrix scale-regions -S H3K9ac_WT_rep03_ratio.bw -R Mus_musculus.GRCm38.102.chr.gtf --regionBodyLength 2000 -b 5000 -a 5000 -o H3K9ac_WT_metagene.mat.gz --skipZeros --missingDataAsZero --metagene
plotHeatmap -m H3K9ac_WT_metagene.mat.gz --dpi 300 -o H3K9ac_WT_metagene.pdf

bamCompare -b1 H3K9ac_WT_rep01_sorted.REF_BK000964.3.bam -b2 input_WT_rdna_merged.bam -o H3K9ac_WT_rep01_rdna.bw --scaleFactorsMethod None --operation ratio --binSize 10 --normalizeUsing RPKM --smoothLength 25 --extendReads --effectiveGenomeSize 45000 

3.4 Peak calling (rDNA) - Bamtools and MACS2


  
3.5 Peak calling (genome-wide) - MACS2
  
3.6 Differential binding analysis - R Bioconductor package 'diffbind'
 
3.7 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'

