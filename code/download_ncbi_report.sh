#!/usr/bin/env bash

# Descargar el reporte en ENBI
url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA556343&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp,bam_ftp&format=tsv&download=true&limit=0"

wget -p metadata/ -O metadata/report.tsv $url -q