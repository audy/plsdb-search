# plsdb-search

Nextflow pipeline for detecting plasmids in Illumina short-read data using [PLSDB2025](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/).

## How it works

1. **FASTP** — adapter trimming and QC
2. **Mash screen** — rapidly screens raw reads against the PLSDB Mash sketch to identify candidate plasmids (no assembly required)
3. **Filter** — removes hits below identity/p-value thresholds
4. **MultiQC** — aggregated QC report across all samples
5. *(optional)* **BWA + samtools** — maps trimmed reads back to hit plasmid sequences and computes per-plasmid coverage breadth and depth

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- Singularity (all tools run in containers; no other software needed on the host)
- PLSDB Mash sketch — download the META archive from the [PLSDB2025 download page](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download) and extract the `.msh` file

## Quick start

```bash
# Set a shared cache dir so images are only pulled/converted once
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache

nextflow run main.nf \
    --input samplesheet.csv \
    --plsdb_sketch /path/to/plsdb.msh
```

## Samplesheet format

CSV with columns `sample`, `fastq_1`, `fastq_2`. Leave `fastq_2` empty for single-end data.

```csv
sample,fastq_1,fastq_2
sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
sample2,/data/sample2_R1.fastq.gz,
```

## Parameters

| Parameter | Description | Default |
|---|---|---|
| `--input` | Samplesheet CSV | required |
| `--plsdb_sketch` | PLSDB Mash sketch (`.msh`) | required |
| `--outdir` | Output directory | `results` |
| `--run_mapping` | Also map reads to hits with BWA for coverage stats | `false` |
| `--plsdb_fasta` | PLSDB FASTA file (required when `--run_mapping`) | — |
| `--mash_min_identity` | Minimum Mash identity to report | `0.90` |
| `--mash_max_pvalue` | Maximum Mash p-value to report | `1e-5` |
| `--mash_threads` | Threads for `mash screen` | `4` |
| `--fastp_min_length` | Minimum read length after trimming | `50` |
| `--fastp_threads` | Threads for `fastp` | `4` |
| `--bwa_threads` | Threads for `bwa mem` | `4` |
| `--min_coverage_pct` | Min % of plasmid bases covered (mapping mode) | `50.0` |
| `--min_mean_depth` | Min mean read depth (mapping mode) | `1.0` |

## Outputs

```
results/
├── fastp/<sample>/          # trimmed reads, QC JSON + HTML
├── mash/<sample>/           # raw mash screen output
├── hits/<sample>/
│   ├── <sample>.hits.tsv           # filtered plasmid hits (identity, p-value, accession)
│   └── <sample>.hit_accessions.txt # accession list
├── mapping/<sample>/        # (--run_mapping) BAM + extracted plasmid FASTA
├── coverage/<sample>/       # (--run_mapping) per-plasmid coverage TSV
├── multiqc/                 # aggregated QC report
└── pipeline_info/           # Nextflow execution timeline, report, trace
```

## Running on HPC (Slurm)

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/singularity_cache

nextflow run main.nf \
    -profile slurm \
    --input samplesheet.csv \
    --plsdb_sketch /shared/db/plsdb.msh
```

The `slurm` profile sets `--cpus-per-task` automatically from each process's resource label.

## Mapping mode

When `--run_mapping` is set, the pipeline extracts the sequences of Mash hits from the PLSDB FASTA, maps trimmed reads against them with BWA, and reports per-plasmid coverage breadth and depth. Plasmids with `< 50%` breadth or `< 1×` depth are filtered out (adjust with `--min_coverage_pct` and `--min_mean_depth`).

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --plsdb_sketch /path/to/plsdb.msh \
    --run_mapping \
    --plsdb_fasta /path/to/plsdb.fasta
```
