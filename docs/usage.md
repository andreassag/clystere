# Usage

## Samplesheet

Prepare a comma-separated samplesheet describing your input genomes:

| Column       | Required | Description                                                                                                                       |
| ------------ | -------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `sample`     | Yes      | Unique sample identifier (used as the output directory name)                                                                      |
| `genome`     | Yes      | Path to a genome file (GenBank `.gbff`/`.gbk`, EMBL `.embl`, or FASTA `.fna`/`.fa`/`.fasta`). Gzip-compressed files are accepted. |
| `annotation` | No       | Path to a GFF3 annotation file. When provided, antiSMASH skips its built-in gene-finding step.                                    |

```csv title="samplesheet.csv"
sample,genome,annotation
strain_A,data/strain_A.gbff.gz,
strain_B,data/strain_B.fna,
strain_C,data/strain_C.fna.gz,data/strain_C.gff3
```

---

## Basic run

```bash
nextflow run exterex/clystere \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

---

## Profile combinations

Profiles are combined with commas. The first group selects the software environment; the second selects the execution
backend.

```bash
# Singularity on a local machine
-profile singularity

# Docker on SLURM (elevated resource ceiling)
-profile docker,slurm

# Conda on an HPC without a scheduler integration
-profile conda,hpc
```

| Profile                 | Type        | Description                                              |
| ----------------------- | ----------- | -------------------------------------------------------- |
| `docker`                | Container   | Run all processes in Docker containers                   |
| `apptainer/singularity` | Container   | Run all processes in Apptainer images                    |
| `podman`                | Container   | Run all processes in Podman containers                   |
| `conda`                 | Environment | Use per-process Conda environments                       |
| `hpc`                   | Resource    | Raise default CPU/memory/time ceilings for cluster nodes |
| `slurm`                 | Executor    | Submit jobs to a SLURM scheduler (combine with `hpc`)    |
| `test`                  | Test        | Preset params for the bundled example data               |

---

## Enabling optional analyses

### BiG-SCAPE

```bash
nextflow run exterex/clystere \
    --input samplesheet.csv \
    --outdir results \
    --run_bigscape \
    -profile docker
```

### KnownClusterBlast annotation

Adds `knownclusterblast_hit`, `knownclusterblast_accession`, and `knownclusterblast_similarity` columns to
`all_regions.tsv`.

```bash
--antismash_cb_knownclusters
```

### Per-contig BGC counts

```bash
--count_per_contig        # one row per contig instead of per assembly
--split_hybrids           # count each BGC type in a hybrid region individually
```

### Disabling minimal mode

antiSMASH's minimal mode is enabled by default. To run the full detection suite, disable it and selectively enable
modules:

```bash
--antismash_minimal false \
--antismash_enable_nrps_pks \
--antismash_enable_terpene \
--antismash_clusterhmmer
```

---

## Resuming runs

Nextflow caches intermediate results in the `work/` directory. Use `-resume` to skip already completed steps after a
pipeline modification or failure:

```bash
nextflow run exterex/clystere --input samplesheet.csv --outdir results -profile docker -resume
```

---

## Reusing existing antiSMASH results

If antiSMASH has already been run, set `--antismash_reuse_results` to skip re-annotation and go directly to tabulation
and BiG-SCAPE:

```bash
--antismash_reuse_results   # skip antiSMASH; use results already in --outdir/antismash/
```
