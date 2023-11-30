# Tools used in RNA-seq analysis

## parallel-fastq-dump 
efficiently downloads DNA sequencing data from a high-throughput sequencer in parallel, converting it into FASTQ format commonly used for downstream analysis, speeding up data retrieval for large-scale genomics projects.

##  trim-galore
automates quality control for both DNA and RNA sequence data by trimming low-quality bases and adapter sequences. This enhances data quality, benefiting subsequent analyses like genome assembly or variant calling.

## fastqc 
https://github.com/utsumidaisuke/invmed/tree/main/1_basic_rnaseq#:~:text=Quality%20chceck-,fastqc,-mkdir%20qc%0Afastqc  
assesses the quality of sequencing data from DNA or RNA experiments. It generates reports on various metrics like read quality, GC content, and sequence duplication, helping researchers identify issues and make informed decisions during data processing and analysis.

## hisat2 
is a bioinformatics tool for aligning DNA sequencing reads to a reference genome. It accurately maps reads, identifying their origin within the genome, a crucial step in tasks like RNA-seq analysis to study gene expression.

## samtools
processes and analyzes DNA and RNA sequencing data, particularly in the context of next-generation sequencing. It performs tasks like file format conversion, indexing, and variant calling, aiding researchers in genomics research and variant identification.

## stringtie
is a swift and effective RNA-Seq assembler that employs a unique network flow algorithm and offers optional de novo assembly. It quantifies full-length transcripts, including splice variants, using short read alignments and longer sequences assembled from them.

## deseq2
is a bioinformatics tool for differential gene expression analysis. It identifies genes that are significantly differentially expressed between two or more experimental conditions in RNA-seq data, helping researchers understand how gene expression varies across conditions, such as disease vs. control samples.