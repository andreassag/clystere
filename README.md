# clystere

**clystere** is a Nextflow pipeline for automated biosynthetic gene cluster (BGC) discovery and comparative analysis. It
runs [antiSMASH](https://antismash.secondarymetabolites.org/), [GECCO](https://github.com/zellerlab/GECCO), and
[deepBGC](https://github.com/Merck/deepbgc) by default across a collection of genomes, unifies predictions with
[comBGC](https://github.com/tomrichtermeier/comBGC-Filter), and optionally groups resulting BGCs into gene cluster
families (GCFs) with [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) or
[BiG-SLiCE](https://github.com/medema-group/bigslice).

---

## Features

- Parallel antiSMASH + GECCO + deepBGC annotation across any number of genome assemblies or GenBank files
- comBGC-based unification of overlapping predictions from all three tools before clustering
- Per-region tabulation and per-genome BGC count summary
- Optional BiG-SCAPE or BiG-SLiCE clustering (mutually exclusive)
- Optional automatic `bigscape dereplicate` step before BiG-SCAPE clustering

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- One of: Docker, Singularity, Podman, or Conda

---

## Quick start

Run with Docker on the bundled example data

```bash
nextflow run andreassag/clystere \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile docker
```

---

## Inputs

### Samplesheet

A comma-separated file with the following columns:

| Column       | Required | Description                                              |
| ------------ | -------- | -------------------------------------------------------- |
| `sample`     | Yes      | Unique sample identifier                                 |
| `genome`     | Yes      | Path to a genome file (GenBank, EMBL, or FASTA)          |
| `annotation` | No       | Path to a GFF3 annotation file (suppresses gene-finding) |

```csv
sample,genome,annotation
strain_A,data/strain_A.gbff.gz,
strain_B,data/strain_B.fna,,
strain_C,data/strain_C.fna.gz,data/strain_C.gff3
```

### antiSMASH database

A pre-built antiSMASH database directory. If `--antismash_db` points to a missing or empty directory the pipeline will
download the database there automatically. The database is also resolvable via `antismash-download-databases`.

---

## Parameters

A full parameter reference is available in [`nextflow_schema.json`](nextflow_schema.json). Key parameters are summarised
below.

### Input / output

| Parameter        | Default   | Description                          |
| ---------------- | --------- | ------------------------------------ |
| `--input`        | —         | Path to samplesheet CSV (required)   |
| `--outdir`       | `results` | Directory for all pipeline outputs   |
| `--antismash_db` | —         | Path to antiSMASH database directory |

### antiSMASH

| Parameter                      | Default    | Description                                                         |
| ------------------------------ | ---------- | ------------------------------------------------------------------- |
| `--antismash_taxon`            | `bacteria` | Taxonomic scope (`bacteria` or `fungi`)                             |
| `--antismash_minimal`          | `true`     | Run in minimal mode; enable modules individually                    |
| `--antismash_cb_knownclusters` | `false`    | Run KnownClusterBlast; adds similarity columns to `all_regions.tsv` |
| `--antismash_genefinding_tool` | `prodigal` | Gene caller when no annotation is supplied                          |
| `--antismash_minlength`        | `1000`     | Minimum sequence length (bp)                                        |
| `--antismash_accept_failure`   | `false`    | Continue if antiSMASH fails for a sample                            |
| `--antismash_extra_args`       | `""`       | Arbitrary additional flags passed to antiSMASH                      |

### GECCO

| Parameter            | Default | Description                          |
| -------------------- | ------- | ------------------------------------ |
| `--gecco_run`        | `true`  | Run GECCO BGC prediction             |
| `--gecco_extra_args` | `""`    | Additional arguments passed to GECCO |

### deepBGC

| Parameter              | Default | Description                                                      |
| ---------------------- | ------- | ---------------------------------------------------------------- |
| `--deepbgc_run`        | `true`  | Run deepBGC prediction                                           |
| `--deepbgc_data_dir`   | —       | Path to deepBGC model/Pfam downloads (auto-downloaded if absent) |
| `--deepbgc_extra_args` | `""`    | Additional arguments passed to deepBGC                           |

### comBGC unification

| Parameter              | Default | Description                                  |
| ---------------------- | ------- | -------------------------------------------- |
| `--combgc_min_length`  | `3000`  | Minimum BGC length retained by comBGC        |
| `--combgc_contig_edge` | `2`     | Exclude BGCs close to contig edges in comBGC |

### BiG-SCAPE

BiG-SCAPE and BiG-SLiCE in clystere run on unified comBGC-filtered regions and require `--gecco_run true` and
`--deepbgc_run true`.

| Parameter                       | Default       | Description                                  |
| ------------------------------- | ------------- | -------------------------------------------- |
| `--bigscape_run`                | `false`       | Enable BiG-SCAPE GCF clustering              |
| `--bigscape_dereplicate`        | `true`        | Run `bigscape dereplicate` before clustering |
| `--bigscape_dereplicate_cutoff` | `0.8`         | Similarity cutoff for dereplication          |
| `--bigscape_gcf_cutoffs`        | `0.3 0.5 0.7` | Space-separated list of distance cutoffs     |
| `--bigscape_mix`                | `true`        | Combine all BGC classes into one network     |
| `--bigscape_include_singletons` | `true`        | Include singletons in the output             |

### BiG-SLiCE

| Parameter               | Default | Description                                                            |
| ----------------------- | ------- | ---------------------------------------------------------------------- |
| `--bigslice_run`        | `false` | Enable BiG-SLiCE clustering (mutually exclusive with `--bigscape_run`) |
| `--bigslice_extra_args` | `""`    | Additional arguments passed to BiG-SLiCE                               |
| `--bigslice_zip_output` | `false` | Compress BiG-SLiCE output directory                                    |

### Tabulation

| Parameter            | Default | Description                                        |
| -------------------- | ------- | -------------------------------------------------- |
| `--run_tabulation`   | `true`  | Generate `all_regions.tsv` and `region_counts.tsv` |
| `--count_per_contig` | `false` | Report counts per contig rather than per assembly  |
| `--split_hybrids`    | `false` | Count each product type in hybrid BGCs separately  |

---

## Outputs

```text
results/
├── antismash/
│   └── <sample>/          # Full antiSMASH output per genome
├── gecco/
│   └── <sample>/          # GECCO outputs per genome (+ BiG-SLiCE-compatible regions)
├── deepbgc/
│   └── <sample>/          # deepBGC outputs per genome (+ converted region GBKs)
├── combgc/
│   └── <sample>/
│       ├── combgc_summary.tsv
│       └── combined_regions/   # Unified representative region GBKs used for clustering
├── bigscape/              # BiG-SCAPE output (when --bigscape_run)
├── bigslice/              # BiG-SLiCE output (when --bigslice_run)
├── summary/
│   ├── all_regions.tsv    # One row per BGC region across all samples
│   └── region_counts.tsv  # BGC type counts per genome (or per contig)
└── pipeline_info/         # Execution timeline, report, trace, and DAG
```

### `all_regions.tsv`

| Column                         | Description                                                                                  |
| ------------------------------ | -------------------------------------------------------------------------------------------- |
| `file`                         | Source antiSMASH run (genome stem)                                                           |
| `record_id`                    | Sequence/contig identifier                                                                   |
| `region`                       | Region number within the record                                                              |
| `start` / `end`                | Genomic coordinates (bp)                                                                     |
| `contig_edge`                  | Whether the region extends to a contig boundary                                              |
| `product`                      | BGC product class(es)                                                                        |
| `knownclusterblast_hit`        | Top MIBiG hit description _(only when `--antismash_cb_knownclusters`)_                       |
| `knownclusterblast_accession`  | MIBiG accession _(only when `--antismash_cb_knownclusters`)_                                 |
| `knownclusterblast_similarity` | Similarity category: `low`, `medium`, or `high` _(only when `--antismash_cb_knownclusters`)_ |
| `record_desc`                  | Sequence description from the source file                                                    |

### `region_counts.tsv`

One row per genome (or per contig with `--count_per_contig`) with integer counts for each BGC product class detected,
plus `total_count` and `description` columns.

---

## Profiles

```bash
# Docker (default for local runs)
-profile docker

# Singularity (recommended for HPC)
-profile singularity

# Apptainer (recommended for HPC)
-profile apptainer

# Conda
-profile conda

# SLURM cluster — sets executor + raises resource ceilings
-profile singularity,slurm

# Generic HPC — raises resource ceilings without binding to a scheduler
-profile singularity,hpc
```

---

## Citations

Please cite the pipeline and its dependencies. See [CITATIONS.md](CITATIONS.md) for full references.

---

## Licence

[MIT](LICENSE)
