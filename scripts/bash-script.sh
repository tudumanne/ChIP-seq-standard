#!/bin/bash

# usage - chip-seq standard pipeline

#change the working directory
cd ~/chip-seq/experiment_01

#quality check of raw fastq files
fastqc -o fastqc --extract --dir fastqc_chip --format fastq h3k4me3_*.fastq.gz
fastqc -o fastqc --extract --dir fastqc_input --format fastq input_*.fastq.gz
multiqc fastqc_chip/

#read alignment and file processing
bowtie2 -x index -1 {sample}_R1.fastq.gz -2 {sample}_R2.fastq.gz | samtools view -bS > {sample}.bam
samtools sort {sample}.bam > {sample}_sorted.bam
samtools index {sample}_sorted.bam

#quality check of bam files
qualimap bamqc -bam {sample}.bam -outdir qualimap_results -outformat pdf
multiqc qualimap_results/

#deeptools fingerprint plot
plotFingerprint -b {sample}_sorted.bam -plot {sample}.pdf --extendReads --centerReads --plotFileFormat pdf --plotTitle "sample_name"

#deeptools correlation and pca
multiBamSummary bins --bamfiles *_sorted.bam -o {ChIP_target}.npz --binSize 1000 --extendReads --centerReads
plotCorrelation --corData {ChIP_target}.npz --corMethod spearman --whatToPlot heatmap --plotFile correlation_heatmap.pdf --skipZeros --plotFileFormat pdf --plotNumbers --plotTitle "ChIP_target_name"
plotPCA --corData {ChIP_target}.npz --plotFile pca_plot.pdf --plotTitle "ChIP_target_name" --plotFileFormat pdf

#deeptools coverage track
bamCoverage -b {sample}_sorted.bam -o {sample}_coverage.bw --outFileFormat bigwig --binSize 10 --smoothLength 25 --extendReads --ignoreForNormalization chrX chrY BK000964.3 chrMT

#deeptools normalised _coverage
bamCompare -b1 {sample}_sorted.bam -b2 {input}_merged.bam -o {sample}_ratio.bw --scaleFactorsMethod None --operation ratio --binSize 10 --normalizeUsing RPKM --smoothLength 25 --extendReads --ignoreForNormalization chrX chrY BK000964.3 chrMT

#deeptools score matrix, heatmap and profile
computeMatrix scale-regions -S {sample}*.bw -R Mus_musculus.GRCm38.102.chr.gtf --regionBodyLength 2000 -b 5000 -a 5000 -o {sample}.mat.gz --skipZeros --missingDataAsZero --metagene
plotHeatmap -m {sample}.mat.gz --dpi 300 -o {sample}.pdf

#macs2 peak calling
macs2 {sample}_sorted.bam -c {input}_merged.bam -g mm -n sample -q 0.01 --cutoff -f BAMPE --keep-dup all

#Extract rDNA as a separate chromosome
bamtools split -in {sample}_sorted.bam -reference

#macs2 peak calling rDNA only 
macs2 callpeak -t {sample}.REF_BK000964.3.bam -c {input}_merged.bam -n sample -g 4.5e4 -f BAMPE --outdir sample -q 0.01 --keep-dup all

