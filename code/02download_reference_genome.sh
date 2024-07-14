#!/usr/bin/env bash

## Genoma de referencia del Solanum virginianum -> https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000787875.1/
url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL3.1/GCF_000188115.5_SL3.1_genomic.fna.gz"

genome=data/reference/genome.fna.gz 

wget -O $genome $url -q 
gzip -d $genome