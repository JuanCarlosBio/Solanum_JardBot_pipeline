#!/usr/bin/env bash

## Genoma de referencia del Solanum virginianum -> https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000787875.1/
url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/787/875/GCA_000787875.1_SME_r2.5.1/GCA_000787875.1_SME_r2.5.1_genomic.fna.gz"

genome=data/reference/genome.fna.gz 

wget -O $genome $url -q 
gzip -d $genome