# Basic of RNA-seq analysis

## Environment setup
```
mamba create -n rnaseq python=3.10
```

## Installtion of necessary tools
```
mamba install -c bioconda parallel-fastq-dump -y
mamba install -c bioconda trim-galore -y
mamba install -c bioconda fastqc -y
mamba install -c bioconda hisat2 -y
mamba install -c bioconda stringtie -y
mamba install -c bioconda bioconductor-deseq2 -y
```

## Reference file preparation
Download all the necessary reference files from NCBI database  
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
```
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000001405.40.zip" -H "Accept: application/zip"
unzip GCF_000001405.40.zip "ncbi_dataset/data/GCF_000001405.40/*" -d ref
```


## Output directories  
Creation of directories for output data
```
for i in SRR11309003 SRR11309004 SRR11309005 SRR11309006
do
mkdir -p analysis/$i
done
```

## Example fastq files
Download public RNA-seq data  
https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP252863  
```
parallel-fastq-dump --sra-id SRR11309003 --threads 4 --outdir fastq --split-files --gzip
parallel-fastq-dump --sra-id SRR11309004 --threads 4 --outdir fastq --split-files --gzip
parallel-fastq-dump --sra-id SRR11309005 --threads 4 --outdir fastq --split-files --gzip
parallel-fastq-dump --sra-id SRR11309006 --threads 4 --outdir fastq --split-files --gzip
```


## RNA-seq workflow (SRR11309003 as an example)
<img src="fig/RNAseqWorkflow.png" width='300'>

**Move to working derectory**  
```
cd analysis/SRR11309003
```

**Quality chceck**  
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
```
mkdir qc
fastqc -t 10 -o qc  SRR11309003_1.fastq.gz SRR11309003_2.fastq.gz
```

**Adaptor trimming**  
trim_galore  
```
```

**Alignment**  
HISAT2  
```
```

**Read count**  
StringTie  
```
```

**DEG(differentially expressed genes)**  
DESeq2(R)
```
```

## Overview of the result