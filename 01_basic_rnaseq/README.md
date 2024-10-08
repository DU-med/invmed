# Basic of RNA-seq analysis

## Environment setup
```
mamba create -n rnaseq python=3.10
mamba activate rnaseq
```

## Installtion of necessary tools
```
mamba install bioconda::parallel-fastq-dump -y
mamba install bioconda::trim-galore -y
mamba install bioconda::fastqc -y
mamba install bioconda::hisat2 -y
mamba install bioconda::samtools -y
mamba install bioconda::stringtie -y
mamba install bioconda::bioconductor-deseq2 -y
mamba install bioconda::tpmcalculator -y
mamba install anaconda::jupyter -y
mamba install anaconda::pandas -y
mamba install conda-forge::matplotlib -y
mamba install conda-forge::nbclassic -y
mamba install conda-forge::r-irkernel -y
```

## Reference file preparation
Download all the necessary reference files from NCBI database (takes approximately 20min)
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
```
wget -r -l1 -np -nH --cut-dirs=7 -R "index.html*" -P ref ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
gunzip ref/GCF_000001405.40_GRCh38.p14_genomic.fna.gz ref/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz ref/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
```

## Genome index file
Indexing genome file for HISAT2
```
mkdir index
hisat2-build -p 20 ref/GCF_000001405.40_GRCh38.p14_genomic.fna index/human_genome
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
HEK 293 cell: SRR11309003 SRR11309004  
HBEC 5i cell: SRR11309005 SRR11309006  

## Output directories  
Creation of directories for output data
```
for i in SRR11309003 SRR11309004 SRR11309005 SRR11309006
do
mkdir -p analysis/$i
done
```

## RNA-seq workflow  
<img src="fig/RNAseqWorkflow.png" width='300'>

## 1️⃣ Step 1 : Processing fastq in each sample directory  (SRR11309003 as an example)
### Move to the working derectory
```
cd analysis/SRR11309003
```

### Quality chceck
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
```
mkdir qc
fastqc -t 4 -o qc ../../fastq/SRR11309003_1.fastq.gz ../../fastq/SRR11309003_2.fastq.gz
```

### Adaptor trimming
[trim-galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)
```
mkdir trimmed_fastq
trim_galore -j 20 --paired ../../fastq/SRR11309003_1.fastq.gz ../../fastq/SRR11309003_2.fastq.gz -o trimmed_fastq
```

### Alignment
[HISAT2](https://daehwankimlab.github.io/hisat2/manual/)  
```
hisat2 -p 16 -x ../../index/human_genome -1 trimmed_fastq/SRR11309003_1_val_1.fq.gz -2 trimmed_fastq/SRR11309003_2_val_2.fq.gz -S SRR11309003.sam 
```

### SAM and BAM file processing
[samtools](https://www.htslib.org/doc/samtools.html)
```
samtools sort SRR11309003.sam -@ 10 -O bam -o SRR11309003.sort.bam 
samtools index SRR11309003.sort.bam
rm SRR11309003.sam
```

### Read count
[StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
```
stringtie -p 10 -e -G ../../ref/GCF_000001405.40_GRCh38.p14_genomic.gff -o SRR11309003.gtf -A  SRR11309003.table SRR11309003.sort.bam
```

### TPMCalculator
[TPMCalculator](https://github.com/ncbi/TPMCalculator)
```
TPMCalculator -g ../../ref/GCF_000001405.40_GRCh38.p14_genomic.gtf -b SRR11309003.sort.bam > output
```

## 2️⃣ Step 2: Integration and analysis of the results in Step 1

### Move to analysis/ directory
Go to analysis/ directory where subdirectories of each sample are located  
<img src="fig/Tree.png" width='300'>

### Create read count table (gene_count_matrix.csv) for DESeq2
[prepDE.py3](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#deseq)
```
wget -c https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
chmod +x prepDE.py3
python ./prepDE.py
```

### DEG(differentially expressed genes) analysis
[DESeq2(R)](https://github.com/thelovelab/DESeq2)  
**ipynb**: analysis/DESeq2.ipynb
```
# Variable settings
in_f <- "gene_count_matrix.csv"        # input readcount data
out_f1 <- "results/result.txt"        # output file for the result
out_f2 <- "results/result.png"        # output image file for M-A plot
param_G1 <- 2        # sample number of group 1
param_G2 <- 2        # sample number of group2 
param_FDR <- 0.05        # false discovery rate (FDR) threshold
param_fig <- c(400, 380)        # figure size of M-A plot

# Loading necessary R package
library(DESeq2)        # loading DESeq2 package

# Loading  readc ount data
data <- read.table(in_f, header=TRUE, row.names=1, sep=",", quote="")

# Preprocessing
data.cl <- c(rep(1, param_G1), rep(2, param_G2))
colData <- data.frame(condition=as.factor(data.cl))
d <- DESeqDataSetFromMatrix(countData=data, colData=colData, design=~condition)

# Defferentially expressed genes analysis (DESeq2)
d <- DESeq(d)                          # DESeq2 excution
tmp <- results(d)                      # assign the result 
p.value <- tmp$pvalue                  # assign the p-value
p.value[is.na(p.value)] <- 1           # replace NA by 1
q.value <- tmp$padj                    # assgin the adjusted p-value in q.value
q.value[is.na(q.value)] <- 1           # replace NA by 1
ranking <- rank(p.value)               # assign the ranking based on p-value
log2.FC <- tmp$log2FoldChange        # assign the log2FC
sum(q.value < param_FDR)               # the number of genes (q.value < param_FDR)
sum(p.adjust(p.value, method="BH") < param_FDR)        # the number of genes (q.value < param_FDR, BH-method)

# saving the result into the output file
tmp <- cbind(rownames(data), data, p.value, q.value, ranking, log2.FC)
write.table(tmp, out_f1, sep="\t", append=F, quote=F, row.names=F) 

# saving MA-plot
png(out_f2, pointsize=13, width=param_fig[1], height=param_fig[2])
plotMA(d)
```

## Overview of the result
**ipynb**: analysis/result_overview.ipynb
```
# importing pandas library
import pandas as pd

# Loading the DESeq2 result file 
df = pd.read_csv('results/result.txt', sep='\t')

# Sorting by "ranking"
df = df.sort_values('ranking')

# Showing the dataframe
display(df)
```
