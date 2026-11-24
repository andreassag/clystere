# Usage

## Samplesheet

Prepare a comma-separated samplesheet describing your input genomes:

| Column       | Required | Description                                                                          |
| ------------ | -------- | ------------------------------------------------------------------------------------ |
| `sample`     | Yes      | Unique sample identifier (used as the output directory name)                         |
| `genome`     | Yes      | Path to a GenBank or FASTA file (`.gbff`, `.gbk`, `.fa`, `.fna`, `.fasta`, or `.gz`) |
| `annotation` | No       | Path to an optional GFF3 annotation file (`.gff`, `.gff3`, or `.gz`)                 |
| `fastq_1`    | No       | Optional R1 or single-end FASTQ for BGC abundance/expression profiling               |
| `fastq_2`    | No       | Optional R2 paired-end FASTQ                                                         |
| `assay_type` | No       | Optional assay type: `metagenome` (default) or `metatranscriptome`                   |

```csv title="samplesheet.csv"
sample,genome,annotation
strain_A,data/strain_A.gbff.gz,
strain_B,data/strain_B.fna,
```

When a row provides a raw `.fa`/`.fna`/`.fasta` genome without an `annotation` GFF3, clystere automatically runs
`--annotation_tool` (`pyrodigal` by default, or `bakta`) before antiSMASH.

For multi-sample metagenomic or metatranscriptomic read cohorts, pass a separate `--reads_input` CSV:

```csv title="reads.csv"
sample,fastq_1,fastq_2,assay_type
site1_dna,reads/site1_dna_R1.fastq.gz,reads/site1_dna_R2.fastq.gz,metagenome
site1_rna,reads/site1_rna_R1.fastq.gz,reads/site1_rna_R2.fastq.gz,metatranscriptome
```

---

## Basic run

```bash
nextflow run andreassag/clystere \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile docker
```

By default, clystere runs antiSMASH 8, GECCO, and deepBGC for every sample, aggregates via `comBGC.py`, and generates
both antiSMASH tabulation tables and cross-tool consensus summaries (`consensus_loci.tsv` and
`genome_qc_bgc_summary.tsv`). To disable an individual predictor:

```bash
--gecco_run false
--deepbgc_run false
```

---

## Profile combinations

Profiles are combined with commas. The first group selects the software environment; the second selects the execution
backend.

```bash
# Apptainer on a local machine
-profile apptainer

# Docker on SLURM (elevated resource ceiling)
-profile docker,slurm

# Conda on an HPC without a scheduler integration
-profile conda,hpc

# Compact fast test dataset with Apptainer
-profile test,apptainer
```

| Profile                 | Type        | Description                                                       |
| ----------------------- | ----------- | ----------------------------------------------------------------- |
| `docker`                | Container   | Run all processes in Docker containers                            |
| `apptainer/singularity` | Container   | Run all processes in Apptainer/Singularity images                 |
| `podman`                | Container   | Run all processes in Podman containers                            |
| `conda`                 | Environment | Use per-process Conda environments                                |
| `hpc`                   | Resource    | Raise default CPU/memory/time ceilings for cluster nodes          |
| `slurm`                 | Executor    | Submit jobs to a SLURM scheduler (combine with `hpc`)             |
| `test`                  | Test        | Fast compact BGC-bearing genome slice dataset                     |
| `test_full`             | Test        | Full bundled _E. coli_ + _B. subtilis_ dataset with BiG-SCAPE run |

---

## Enabling optional analyses

### Tier 1: Preprocessing, MAG QC & Taxonomy

```bash
nextflow run andreassag/clystere \
    --input samplesheet.csv \
    --outdir results \
    --annotation_tool pyrodigal \
    --pyrodigal_mode meta \
    --quast_run \
    --checkm2_run \
    --gtdbtk_run \
    --gtdbtk_db /path/to/gtdbtk_release \
    -profile apptainer
```

To use Bakta instead of Pyrodigal on unannotated FASTA inputs, pass `--annotation_tool bakta` (automatically downloads
the `light` database unless `--bakta_db` or `--bakta_db_type full` is provided).

### Tier 2: SanntiS, Macrel + AMPcombi2 & RiPP Preset

```bash
nextflow run andreassag/clystere \
    --input samplesheet.csv \
    --outdir results \
    --sanntis_run \
    --amp_run \
    --ripp_preset \
    -profile apptainer
```

### Tier 3: Clustering, BiG-FAM Novelty, clinker & cblaster

BiG-SCAPE and BiG-SLiCE can be run individually **or together** in the same pipeline execution:

```bash
nextflow run andreassag/clystere \
    --input samplesheet.csv \
    --outdir results \
    --bigscape_run \
    --bigslice_run \
    --bigfam_db /path/to/bigfam_db \
    --clinker_run \
    --cblaster_run \
    -profile apptainer
```

- `--clinker_run` generates interactive HTML synteny plots for every multi-member GCF produced by BiG-SCAPE and/or
  BiG-SLiCE.
- `--cblaster_run` defaults to `--cblaster_mode local` (builds an offline cohort DIAMOND + SQLite DB from all input
  genomes unless `--cblaster_db` is given, or set `--cblaster_mode remote` for NCBI).

### Tier 4: Genome Dereplication & BGC/GCF Abundance Profiling

```bash
nextflow run andreassag/clystere \
    --input samplesheet.csv \
    --reads_input reads.csv \
    --outdir results \
    --dereplicate_genomes \
    --dereplicate_tool skani \
    --dereplicate_ani 99.0 \
    --abundance_mapper bowtie2 \
    -profile apptainer
```

---

## Resuming runs

Nextflow caches intermediate results in the `work/` directory. Use `-resume` to skip already completed steps after a
pipeline modification or failure:

```bash
nextflow run andreassag/clystere --input samplesheet.csv --outdir results -profile docker -resume
```
