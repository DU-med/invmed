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
NCBI: GRCm39
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/