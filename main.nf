#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = "data/*.fastq.gz"
params.outdir = "results"

workflow {

    Channel
        .fromPath(params.input)
        .ifEmpty { error "❌ No FASTQ files found in data/ folder" }
        .set { reads }

    fastqc(reads)
    cutadapt(reads)
}

process fastqc {

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path "*_fastqc.html"

    script:
    """
    fastqc $reads
    """
}

process cutadapt {

    publishDir "${params.outdir}/cutadapt", mode: 'copy'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fastq.gz"

    script:
    """
    cutadapt -a AGATCGGAAGAGC \
             -o ${reads.simpleName}_trimmed.fastq.gz \
             $reads
    """
}
