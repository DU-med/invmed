# Slamdunk Setup Guide

### Reference  Repository
[GitHub Repository for Slamdunk](https://github.com/t-neumann/slamdunk)

## 1. Environment setup
### Creating Slandunk Analysis Environment 
```
mamba create -n slamdunk python=3.8
```

### Activating the Environment
```
mamba activate slamdunk
```

### Installing slamdunk
```
mamba install -c bioconda slamdunk -y
```

### Cloning the Slamdunk GitHub Repository
```
git clone https://github.com/t-neumann/slamdunk.git
```

### Test Run
```
cd slamdunk
slamdunk all -r slamdunk/test/data/ref.fa -b slamdunk/test/data/actb.bed -o slamdunk/test/data/output -rl 100 -mbq 27 -5 0 slamdunk/test/data/reads.fq
```

## 2. Preparing Reference files
### Downloading the Mouse Reference Genome 
**NCBI**: GRCm39 (GCF_000001635.27_GRCm39_genomic.fna.gz)  
[NCBI Download Link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/)
```
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
```  
  
**UCSC**: mm39 (mm39.fa.gz)  
[UCSC Download Link](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/)
```
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
```  
  
**Ensembl**: GRCm39 (Mus_musculus.GRCm39.dna.primary_assembly.fa.gz)     
- DO NOT USE FILES BELOW  
Mus_musculus.GRCm39.dna_rm.primary_assembly.fa.gz  
Mus_musculus.GRCm39.dna_rm.toplevel.fa.gz  
Mus_musculus.GRCm39.dna.toplevel.fa.gz  
Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz  
Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz  
[Ensembl Download Link](https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/)  
```
wget -c https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

```

**Gencode**: GRCm39 (GRCm39.primary_assembly.genome.fa.gz)  
[Gencode Download Link](https://www.gencodegenes.org/mouse/)  
```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz
```

