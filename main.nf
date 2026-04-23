#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Help
// ---------------------------------------------------------------------------
def helpMessage() {
    log.info """
    =========================================
    plsdb-search  ~  plasmid detection pipeline
    =========================================
    Usage:
        nextflow run main.nf --input samplesheet.csv --plsdb_sketch plsdb.msh

    Required:
        --input          Path to samplesheet CSV (columns: sample, fastq_1, fastq_2)
                         fastq_2 is optional for single-end data
        --plsdb_sketch   Path to PLSDB Mash sketch (.msh) from the META archive
                         Download: https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download

    Optional:
        --outdir         Output directory [default: results]
        --run_mapping    Also map reads to plasmid hits with BWA [default: false]
        --plsdb_fasta    PLSDB FASTA file, required when --run_mapping is set

    Mash screen thresholds:
        --mash_min_identity   Minimum identity to report [default: 0.90]
        --mash_max_pvalue     Maximum p-value to report  [default: 1e-5]
        --mash_threads        Threads for mash screen    [default: 4]

    fastp options:
        --fastp_min_length    Minimum read length after trimming [default: 50]
        --fastp_threads       Threads for fastp                  [default: 4]

    BWA/samtools options (only with --run_mapping):
        --bwa_threads         Threads for bwa mem [default: 4]
        --min_coverage_pct    Min % of plasmid bases covered to report [default: 50]
        --min_mean_depth      Min mean read depth to report            [default: 1]
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }

// ---------------------------------------------------------------------------
// Parameter validation
// ---------------------------------------------------------------------------
if (!params.input)        error "Please provide --input samplesheet"
if (!params.plsdb_sketch) error "Please provide --plsdb_sketch (Mash .msh file)"
if (params.run_mapping && !params.plsdb_fasta) {
    error "--run_mapping requires --plsdb_fasta"
}

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process FASTP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/fastp/${meta.id}", mode: 'copy'

    conda      'bioconda::fastp=0.23.4'
    container  'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    path "*.json",                                emit: json
    path "*.html",                                emit: html

    script:
    def prefix = meta.id
    def single_end = meta.single_end
    if (single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}.trimmed.fastq.gz \\
            --length_required ${params.fastp_min_length} \\
            --thread ${params.fastp_threads} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_R1.trimmed.fastq.gz \\
            --out2 ${prefix}_R2.trimmed.fastq.gz \\
            --length_required ${params.fastp_min_length} \\
            --thread ${params.fastp_threads} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html
        """
    }
}

process MASH_SCREEN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/mash/${meta.id}", mode: 'copy'

    conda     'bioconda::mash=2.3'
    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    input:
    tuple val(meta), path(reads)
    path sketch

    output:
    tuple val(meta), path("${meta.id}.mash_screen.txt"), emit: screen

    script:
    // Pass all reads as positional args — mash screen handles multiple files
    def reads_args = reads.collect { it.toString() }.join(' ')
    """
    mash screen \\
        -w \\
        -p ${params.mash_threads} \\
        ${sketch} \\
        ${reads_args} \\
        > ${meta.id}.mash_screen.txt
    """
}

process FILTER_MASH {
    tag "$meta.id"
    label 'process_single'
    publishDir "${params.outdir}/hits/${meta.id}", mode: 'copy'

    container 'docker.io/library/python:3.11-slim'

    input:
    tuple val(meta), path(screen)

    output:
    tuple val(meta), path("${meta.id}.hits.tsv"),         emit: hits
    tuple val(meta), path("${meta.id}.hit_accessions.txt"), emit: accessions

    script:
    """
    filter_mash.py \\
        --input ${screen} \\
        --min_identity ${params.mash_min_identity} \\
        --max_pvalue ${params.mash_max_pvalue} \\
        --output ${meta.id}.hits.tsv \\
        --accessions ${meta.id}.hit_accessions.txt
    """
}

// Only runs when --run_mapping is set
process EXTRACT_PLASMID_SEQS {
    tag "$meta.id"
    label 'process_single'
    publishDir "${params.outdir}/mapping/${meta.id}", mode: 'copy'

    conda     'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    input:
    tuple val(meta), path(accessions)
    path plsdb_fasta

    output:
    tuple val(meta), path("${meta.id}.plasmids.fasta"), emit: fasta

    script:
    """
    # Extract only the hit sequences to keep the index small
    samtools faidx ${plsdb_fasta}
    xargs samtools faidx ${plsdb_fasta} < ${accessions} > ${meta.id}.plasmids.fasta
    """
}

process BWA_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda     'bioconda::bwa=0.7.17'
    container 'quay.io/biocontainers/bwa:0.7.17--h7132678_9'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${fasta}.*"), emit: index

    script:
    """
    bwa index ${fasta}
    """
}

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda     'bioconda::bwa=0.7.17'
    container 'quay.io/biocontainers/bwa:0.7.17--h7132678_9'

    input:
    tuple val(meta), path(fasta), path(index_files), path(reads)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam

    script:
    def reads_args = reads instanceof List
        ? reads.collect { it.toString() }.join(' ')
        : reads.toString()
    """
    bwa mem \\
        -t ${params.bwa_threads} \\
        ${fasta} \\
        ${reads_args} \\
        > ${meta.id}.sam
    """
}

process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/mapping/${meta.id}", mode: 'copy'

    conda     'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam ${sam}
    samtools index ${meta.id}.sorted.bam
    """
}

process SAMTOOLS_COVERAGE {
    tag "$meta.id"
    label 'process_single'
    publishDir "${params.outdir}/coverage/${meta.id}", mode: 'copy'

    conda     'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.coverage.tsv"), emit: coverage

    script:
    """
    samtools coverage ${bam} > ${meta.id}.coverage.tsv
    """
}

process FILTER_COVERAGE {
    tag "$meta.id"
    label 'process_single'
    publishDir "${params.outdir}/coverage/${meta.id}", mode: 'copy'

    container 'docker.io/library/python:3.11-slim'

    input:
    tuple val(meta), path(coverage)

    output:
    tuple val(meta), path("${meta.id}.coverage.filtered.tsv"), emit: filtered

    script:
    """
    filter_coverage.py \\
        --input ${coverage} \\
        --min_coverage_pct ${params.min_coverage_pct} \\
        --min_mean_depth ${params.min_mean_depth} \\
        --output ${meta.id}.coverage.filtered.tsv
    """
}

process MULTIQC {
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    conda     'bioconda::multiqc=1.21'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    input:
    path fastp_jsons

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
    // Parse samplesheet
    ch_reads = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = []
            if (row.fastq_2 && row.fastq_2.trim()) {
                meta.single_end = false
                reads = [file(row.fastq_1), file(row.fastq_2)]
            } else {
                meta.single_end = true
                reads = [file(row.fastq_1)]
            }
            [meta, reads]
        }

    ch_sketch = Channel.fromPath(params.plsdb_sketch, checkIfExists: true)

    // QC + trim
    FASTP(ch_reads)

    // Mash screen trimmed reads against PLSDB
    MASH_SCREEN(FASTP.out.reads, ch_sketch.collect())

    // Filter hits by identity / p-value
    FILTER_MASH(MASH_SCREEN.out.screen)

    // MultiQC report across all samples
    MULTIQC(FASTP.out.json.collect())

    // Optional: map reads to plasmid hits for coverage stats
    if (params.run_mapping) {
        ch_plsdb_fasta = Channel.fromPath(params.plsdb_fasta, checkIfExists: true)

        EXTRACT_PLASMID_SEQS(
            FILTER_MASH.out.accessions,
            ch_plsdb_fasta.collect()
        )

        BWA_INDEX(EXTRACT_PLASMID_SEQS.out.fasta)

        // Join per-sample index with per-sample trimmed reads
        ch_bwa_input = BWA_INDEX.out.index
            .join(FASTP.out.reads, by: 0)

        BWA_MEM(ch_bwa_input)
        SAMTOOLS_SORT(BWA_MEM.out.sam)

        SAMTOOLS_COVERAGE(SAMTOOLS_SORT.out.bam)
        FILTER_COVERAGE(SAMTOOLS_COVERAGE.out.coverage)
    }
}
