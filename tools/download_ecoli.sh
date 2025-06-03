#!/bin/bash
# filepath: tools/download_ecoli.sh

mkdir -p data
echo "Downloading E. coli reference genome..."
wget -O data/GCF_000005845.2_ASM584v2_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip -f data/GCF_000005845.2_ASM584v2_genomic.fna.gz
echo "Download complete: data/GCF_000005845.2_ASM584v2_genomic.fna"