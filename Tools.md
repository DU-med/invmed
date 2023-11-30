# Tools used in RNA-seq analysis

## [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump) 
efficiently downloads DNA sequencing data from a high-throughput sequencer in parallel, converting it into FASTQ format commonly used for downstream analysis, speeding up data retrieval for large-scale genomics projects.

## [trim-galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)
automates quality control for both DNA and RNA sequence data by trimming low-quality bases and adapter sequences. This enhances data quality, benefiting subsequent analyses like genome assembly or variant calling.

## [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
assesses the quality of sequencing data from DNA or RNA experiments. It generates reports on various metrics like read quality, GC content, and sequence duplication, helping researchers identify issues and make informed decisions during data processing and analysis.

## [HISAT2](https://daehwankimlab.github.io/hisat2/manual/)  
is a bioinformatics tool for aligning DNA sequencing reads to a reference genome. It accurately maps reads, identifying their origin within the genome, a crucial step in tasks like RNA-seq analysis to study gene expression.

## [samtools](https://www.htslib.org/doc/samtools.html)
processes and analyzes DNA and RNA sequencing data, particularly in the context of next-generation sequencing. It performs tasks like file format conversion, indexing, and variant calling, aiding researchers in genomics research and variant identification.

## [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
is a swift and effective RNA-Seq assembler that employs a unique network flow algorithm and offers optional de novo assembly. It quantifies full-length transcripts, including splice variants, using short read alignments and longer sequences assembled from them.

## [deseq2](https://bioconductor.org/packages/3.17/bioc/html/DESeq2.html)
is a bioinformatics tool for differential gene expression analysis. It identifies genes that are significantly differentially expressed between two or more experimental conditions in RNA-seq data, helping researchers understand how gene expression varies across conditions, such as disease vs. control samples.