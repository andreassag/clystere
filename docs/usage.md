# Usage

## Samplesheet

Prepare a comma-separated samplesheet describing your input genomes:

| Column       | Required | Description                                                                          |
| ------------ | -------- | ------------------------------------------------------------------------------------ |
| `sample`     | Yes      | Unique sample identifier (used as the output directory name)                         |
| `genome`     | Yes      | Path to a GenBank or FASTA file (`.gbff`, `.gbk`, `.fa`, `.fna`, `.fasta`, or `.gz`) |
| `annotation` | No       | Path to an optional GFF3 annotation file (`.gff`, `.gff3`, or `.gz`)                 |

```csv title="samplesheet.csv"
sample,genome,annotation
strain_A,data/strain_A.gbff.gz,
strain_B,data/strain_B.fna,data/strain_B.gff3
```

---

## Basic run

```bash
nextflow run andreassag/clystere \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile docker
```

By default, clystere runs antiSMASH, GECCO, and deepBGC for every sample. To disable an individual predictor:

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

<!-- markdownlint-disable MD046 -->

=== "BiG-SCAPE"

    ```bash
    nextflow run andreassag/clystere \
        --input samplesheet.csv \
        --outdir results \
        --bigscape_run \
        -profile docker
    ```

=== "BiG-SLiCE"

    ```bash
    nextflow run andreassag/clystere \
        --input samplesheet.csv \
        --outdir results \
        --bigslice_run \
        -profile docker
    ```

<!-- markdownlint-enable MD046 -->

`--bigscape_run` and `--bigslice_run` are mutually exclusive. Both work with any enabled subset of predictors (antiSMASH
always runs; GECCO and deepBGC are optional).

When `--bigscape_run` is enabled, clystere runs `bigscape dereplicate` by default before clustering to collapse
redundant regions found across tools. You can disable or tune this behavior with:

```bash
--bigscape_dereplicate false
--bigscape_dereplicate_cutoff 0.8
```

When `--bigslice_run` is enabled and `--bigslice_db` is not supplied, `BIGSLICE_DOWNLOADDB` automatically downloads the
BiG-SLiCE HMM database into `--clystere_db_dir`.

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
nextflow run andreassag/clystere --input samplesheet.csv --outdir results -profile docker -resume
```
