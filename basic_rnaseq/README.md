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
Download necessary reference files from NCBI database  
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
```
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000001405.40.zip" -H "Accept: application/zip"
```

## Example fastq files
Download public RNA-seq data  
https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP252863  
```
parallel-fastq-dump --sra-id SRR11309003 --threads 4 --outdir fastq/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11309004 --threads 4 --outdir fastq/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11309005 --threads 4 --outdir fastq/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11309006 --threads 4 --outdir fastq/ --split-files --gzip
```

## RNA-seq workflow
<img src="fig/RNAseqWorkflow.png" width='500'>
