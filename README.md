# ChIP-seq-standard

This repository contains a customised data analysis pipeline that facilitates simultaneous analysis of ChIP-seq data at ribosomal DNA (rDNA) repeats and genome-wide (For paired-end short read - Illumina data). 

### Table of contents 
1. [Pipeline overview](https://github.com/tudumanne/ChIP-seq-standard/blob/main/README.md#pipeline-overview)
2. [Software installation](#
3. How to run - example dataset
  
     3.1 Quality check of raw fastq files - FastQC/MultiQC
  
     3.2 Read alignment, processing and post-alignment quality check - Bowtie2, Samtools and BamQC/MultiQC

     3.3 Coverage track generation and visualisation - deepTools and IGV
  
     3.4 Peak calling (rDNA) - Bamtools and MACS2
  
     3.5 Peak calling (genome-wide) - MACS2
  
     3.6 Differential binding analysis - R Bioconductor package 'diffbind'
 
     3.7 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'
  

### 1. Pipeline overview

First part of the pipeline (3.1-3.5) was run on an HPC (high-performance computing) system based on CentOS (Linux). R based analyses (3.6-3.7) were carried out in RStudio installed on MacOS Catalina.

The 'Scripts' folder contains template bash scripts and example Snakemake workflows.

![alt text](https://github.com/tudumanne/ChIP-seq-standard/files/7828134/Picture.1.pdf)

### Software installation 

The required software/command-line tools were installed via conda on Linux. 

(Miniconda https://docs.conda.io/en/latest/miniconda.html)

```console
conda env create -n chip-seq -f environment.yml
```

### How to run - example dataset

The folder 'example dataset' contains 9 ChIPed (H3K4me3) and 9 input control samples, 3 biological replicates per each stage (WT, PreM and Mal).

These files contain a subset of reads.  
  
3.1 Quality check of raw fastq files - FastQC/MultiQC

```console
fastqc -o fastqc --extract --dir fastqc_chip --format fastq h3k4me3_*.fastq.gz
fastqc -o fastqc --extract --dir fastqc_input --format fastq input_*.fastq.gz

multiqc fastqc_chip/
```
3.2 Read alignment, processing and post-alignment quality check - Bowtie2, Samtools and BamQC/MultiQC

Read alignment using Bowtie2, SAM to BAM conversion, sorting and indexing the BAM files using Samtools

```console
bowtie2 -x index -1 {sample}_R1.fastq.gz -2 {sample}_R2.fastq.gz | samtools view -bS > {sample}.bam
samtools sort {sample}.bam > {sample}_sorted.bam
samtools index {sample}_sorted.bam
```

Quality check of BAM files

```console
qualimap bamqc -bam file.bam -outdir qualimap_results -outformat pdf
multiqc fastqc/
```

3.3 Coverage track generation and visualisation - deepTools and IGV

Fingerprint plot (determines how well the signal in the ChIP sample can be differentiated from the background distribution)

```console
plotFingerprint -b {sample}_sorted.bam -plot {sample}.pdf --extendReads --centerReads --plotFileFormat pdf --plotTitle "sample"
```

Correlation plot and PCA analysis (demonstrate the overall similarity between two or more aligned sequencing files based on read coverage within genomic regions)

```console
multiBamSummary bins --bamfiles *_sorted.bam -o {ChIP target}.npz --binSize 1000 --extendReads --centerReads
plotCorrelation --corData {ChIP target}.npz --corMethod spearman --whatToPlot heatmap --plotFile correlation_heatmap.pdf --skipZeros --plotFileFormat pdf --plotNumbers --plotTitle "sample"
plotPCA --corData {ChIP target}.npz --plotFile pca_plot.pdf --plotTitle "sample" --plotFileFormat pdf
```

Coverage of 

bamCoverage -b {sample}_sorted.bam -o /g/data/ff84/TU_ChIP_seq02_august2020/raw_fastq_files/bam_files/input_WT_rep01_coverage.bw --outFileFormat bigwig --binSize 10 --smoothLength 25 --extendReads --ignoreForNormalization chrX chrY BK000964.3 chrMT -p 12 


bamCompare -b1 H3K9ac_WT_rep01_sorted.REF_BK000964.3.bam -b2 input_WT_rdna_merged.bam -o H3K9ac_WT_rep01_rdna.bw --scaleFactorsMethod None --operation ratio --binSize 10 --normalizeUsing RPKM --smoothLength 25 --extendReads --effectiveGenomeSize 45000 

computeMatrix scale-regions -S H3K9ac_WT_rep03_ratio.bw -R Mus_musculus.GRCm38.102.chr.gtf --regionBodyLength 2000 -b 5000 -a 5000 -o H3K9ac_WT_metagene.mat.gz --skipZeros --missingDataAsZero --metagene
plotHeatmap -m H3K9ac_WT_metagene.mat.gz --dpi 300 -o H3K9ac_WT_metagene.pdf



3.4 Peak calling (rDNA) - Bamtools and MACS2

bamtools split -in /g/data/ff84/TU_ChIP_seq01_april2020/bam_files/sample_11_sorted.bam -reference

macs2 callpeak -t /g/data/ff84/TU_ChIP_seq03_march2021/bam_files/H3_WT_rep01_sorted.REF_BK000964.3.bam -c /g/data/ff84/TU_ChIP_seq03_march2021/bam_files/input_WT_merged.bam -n /g/data/ff84/TU_ChIP_seq03_march2021/bam_files/H3_WT_rep01_rdna_keepdup -g 4.5e4 -f BAMPE --outdir /g/data/ff84/TU_ChIP_seq03_march2021/bam_files/H3_WT_rep01_rDNA -q 0.01 --keep-dup all

  
3.5 Peak calling (genome-wide) - MACS2

macs2 callpeak -t H3K27me3_WT_rep01_sorted.bam -c WT_input_merged.bam -g mm -n WT_H3K27me3_rep01_keepdup -q 0.01 --cutoff -f BAMPE --keep-dup all
  
3.6 Differential binding analysis - R Bioconductor package 'diffbind'


 
3.7 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'

