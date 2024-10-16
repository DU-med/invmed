# Slamdunk setup

### Reference site
[GitHub](https://github.com/DU-med/invmed)

## 1. Environment setup
### Creation of Slandunk analysis environment 
```
mamba create -n slamdunk python=3.8
```

### Environment activation
```
mamba activate slamdunk
```

### Installtion of slamdunk
```
mamba install -c bioconda slamdunk -y
```

### Cloning github site of Slamdunk
```
git clone https://github.com/t-neumann/slamdunk.git
```

### Test run
```
cd slamdunk
slamdunk all -r slamdunk/test/data/ref.fa -b slamdunk/test/data/actb.bed -o slamdunk/test/data/output -rl 100 -mbq 27 -5 0 slamdunk/test/data/reads.fq
```

## 2. Reference files prepartion
### mouse reference genome download
NCBI: GRCm39(GCF_000001635.27_GRCm39_genomic.fna.gz)  
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/  
  
UCSC: mm39(mm39.fa.gz)  
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/  
  
Ensembl: GRCm39(Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz)    
https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/  
