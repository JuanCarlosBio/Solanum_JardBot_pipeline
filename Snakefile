## Importamos la librería subprocess para usar comandos de la SHELL
import subprocess

## Creamos una lista de las carreras usando el report.tsv del ENBI
SAMPLES = subprocess.run(
  ["cut", "-f", "4", "metadata/report.tsv"], 
  capture_output = True,
  text = True
  ).stdout.strip().split('\n')[1:]

print(
f"""
Estas son las muestras que has seleccionado para el WorkFlow:
>>> {SAMPLES}
"""
)

## Regla maestra que controlará que todos los outputs de las reglas están
## correctos
rule all:
    input:
      "metadata/report.tsv",
      ["data/raw/" + sample + ".fastq.gz" for sample in SAMPLES],
      ["data/processed/" + sample + "_fastp.fastq.gz" for sample in SAMPLES],
      "data/reference/genome.fna",
      ["data/mapped_reads/" + sample + ".sam" for sample in SAMPLES],
      ["data/mapped_reads/BAM/" + sample + ".bam" for sample in SAMPLES],
      ["data/sortedBAM/" + sample + "_sorted.bam" for sample in SAMPLES],
      ["data/dedupBAM/" + sample + "_sorted_dedup.bam" for sample in SAMPLES],
      ["data/rgBAM/" + sample + "_sorted_dedup_rg.bam" for sample in SAMPLES],
      "data/reference/genome.dict",
      "data/reference/genome.fna.fai",
      ["data/GVCF/" + sample + ".gvcf" for sample in SAMPLES]
  
### Descargar el reporte del ENA/EBI
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
    sed -n '1,3p' report.tsv -i
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
    
## Transformar los archivos SAM a BAM
rule sam_to_bam:
  input:
    sam_file = "data/mapped_reads/{filename}.sam"
  output:
    bam_file = "data/mapped_reads/BAM/{filename}.bam",
  conda:
    "code/environments/env.yml"
  shell:
    """
    samtools view -Sb {input.sam_file} > {output.bam_file}
    """

## Ordenar los archivos BAM y obtener las estadísticas del alineamiento con el genoma de referencia
rule sorting_bam:
  input:
    bam_file = "data/mapped_reads/BAM/{filename}.bam"
  output:
    bam_file_sorted = "data/sortedBAM/{filename}_sorted.bam"
  log:
    "logs/BAMflagstats/{filename}.flagstats"
  conda:
    "code/environments/env.yml"
  shell:
    """
    samtools sort {input.bam_file} > {output.bam_file_sorted}

    samtools index {output.bam_file_sorted}

    samtools flagstats {output.bam_file_sorted} > {log} 
    """

## Eliminar duplicados 
rule delete_duplicates_bam_sorted:
  input:
    bam_file_sorted = "data/sortedBAM/{filename}_sorted.bam"
  output:
    bam_file_sorted_dedup = "data/dedupBAM/{filename}_sorted_dedup.bam"
  log:
    "logs/BAMdedup/{filename}.dedup"
  conda:
    "code/environments/env.yml"
  shell:
    """
    picard MarkDuplicates \
      --INPUT {input.bam_file_sorted} \
      --OUTPUT {output.bam_file_sorted_dedup} \
      --METRICS_FILE {log} \
      --ASSUME_SORTED True
    
    samtools index {output.bam_file_sorted_dedup} 
    """

## Añadir RGID RGLB RGPL RGPU RGSM no se muy bien porque los BAM no tienen esta
## info pero sin ella no puedo continuar con GATK :/
rule add_rg_bam:
  input:
    bam_file_sorted_dedup = "data/dedupBAM/{filename}_sorted_dedup.bam"
  output:
    bam_file_rg = "data/rgBAM/{filename}_sorted_dedup_rg.bam" 
  conda:
    "code/environments/env.yml"
  shell:
    """
    picard AddOrReplaceReadGroups \
      I={input.bam_file_sorted_dedup} \
      O={output.bam_file_rg} \
      RGID=1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=sample
    """

## Preparing GATK
rule ref_genome_gatk:
  input:
    reference_genome = "data/reference/genome.fna"
  output:
    dict_file = "data/reference/genome.dict",
    fai_file = "data/reference/genome.fna.fai"
  conda:
    "code/environments/env.yml"
  shell:
    """
    gatk CreateSequenceDictionary \
      -R {input.reference_genome} \
      -O {output.dict_file}
    
    samtools faidx {input.reference_genome}
    """

## Variant Calling usando GATK
rule variant_call_haplotypecaller:
  input:
    reference_genome = "data/reference/genome.fna",
    dict_file = "data/reference/genome.dict",
    fai_file = "data/reference/genome.fna.fai",
    bam_file_sorted = "data/rgBAM/{filename}_sorted_dedup_rg.bam"
  output:
    gvcf_file = "data/GVCF/{filename}.gvcf" 
  conda:
    "code/environments/env.yml"
  shell:
    """
    samtools index {input.bam_file_sorted}

    gatk HaplotypeCaller \
      -R  {input.reference_genome} \
      -I {input.bam_file_sorted} \
      -O {output.gvcf_file}   \
      -ERC GVCF   \
      --sample-name sample
    """
