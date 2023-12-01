# Example command
# excute this script in each sample directory
# bash rnaseq.sh <sample id>
# bash rnaseq.sh SRR1554536

# directory setting
mkdir qc trimmed_fastq

# quality control
fastqc -t 4 -o qc ../../fastq/$1\_1.fastq.gz ../../fastq/$1\_2.fastq.gz

# trimming
trim_galore -j 20 --paired ../../fastq/$1\_1.fastq.gz ../../fastq/$1\_2.fastq.gz -o trimmed_fastq

# mapping
hisat2 -p 16 -x ../../index/human_genome -1 trimmed_fastq/$1\_1_val_1.fq.gz -2 trimmed_fastq/$1\_2_val_2.fq.gz -S $1\.sam 

# convert sam to bam
samtools sort $1\.sam -@ 10 -O bam -o $1\.sort.bam 
samtools index $1\.sort.bam
rm $1\.sam

# read counting
stringtie -p 10 -G ../../ref/ncbi_dataset/data/GCF_000001405.40/genomic.gff -o $1\.gtf -A  $1\.table $1\.sort.bam
