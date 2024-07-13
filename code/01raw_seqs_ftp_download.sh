#!/usr/bin/env bash

## En este script de bash descargaremos los archivos FASTQ del archivo TSV de MEATDATA

## Para ello, seleccionamos la columna específica del archivo
ftp=$(cut -f 7 metadata/report.tsv | grep $1)

## En caso que el directorio de datos crudos no exista, lo creamos
if [[ ! -d data/raw ]];then
  mkdir data/raw
fi

## Descarga del archivo ftp
wget $ftp 
## Mover el archivo a raw data
mv $1 data/raw/