## Importamos la librería subprocess para usar comandos de la SHELL
import subprocess

## Creamos una lista de las carreras usando el report.tsv del ENBI
SAMPLES = subprocess.run(
  ["cut", "-f", "4", "metadata/report.tsv"], 
  capture_output = True,
  text = True
  ).stdout.strip().split('\n')[1:]

## Regla maestra que controlará que todos los outputs de las reglas están
## correctos
rule all:
    input:
      "metadata/report.tsv",
      ["data/raw/" + sample + ".fastq.gz" for sample in SAMPLES],
      ["data/processed/" + sample + "_fastp.fastq.gz" for sample in SAMPLES],
      "data/reference/genome.fna",
      ["data/mapped_reads/" + sample + ".sam" for sample in SAMPLES],
  
## Descargar el reporte del ENA/EBI
rule report:
  input:
    bash_script = "code/download_ncbi_report.sh"
  output:
    "metadata/report.tsv"
  conda:
    "code/environments/env.yml" 
  shell:
    """
    bash {input.bash_script} {output}
    mv report.tsv metadata/
    """

## Primera regla para descargar los datos de forma automática,
## this will take some time :)
rule download_files:
  input:
    metadatos = "metadata/report.tsv",
    bash_script = "code/01raw_seqs_ftp_download.sh"
  output:
      "data/raw/{filename}.fastq.gz" ## Esto es interesante, es una WildCard que selecciona todas las muestras de SAMPLES
  params:
    "{filename}.fastq.gz"
  conda:
    "code/environments/env.yml" 
  shell:
      """
      bash {input.bash_script} {params}
      """

## Procesado de los archivos con fastp
rule fastp_processing:
  input:
    "data/raw/{filename}.fastq.gz"
  output:
    "data/processed/{filename}_fastp.fastq.gz"
  log:
    "logs/fastp/{filename}_fastpinfo.out"
  conda:
    "code/environments/env.yml" 
  shell:
    """
    fastp -i {input} -o {output} 2> {log}
    rm fastp.html fastp.json 
    truncate -s 0 {input}
    """

### Descarga del genoma de referencia
rule reference_genome:
  input:
    bash_script = "code/02download_reference_genome.sh"
  output:
    "data/reference/genome.fna"
  conda:
    "code/environments/env.yml" 
  shell:
    """
    bash {input.bash_script}
    """

## Mapear con el genoma de referencia
rule ref_genome_mapping:
  input:
    ref_genome = "data/reference/genome.fna",
    reads = "data/processed/{filename}_fastp.fastq.gz"
  output:
    "data/mapped_reads/{filename}.sam"
  log:
    "logs/SAM/{filename}_infosam.out"
  conda:
    "code/environments/env.yml" 
  shell:
    """
    bwa index {input.ref_genome}
    bwa mem -a {input.ref_genome} {input.reads} > {output} 2> {log}
    """

