# clystere

**clystere** is a Nextflow pipeline for automated biosynthetic gene cluster (BGC) discovery and comparative analysis. It
orchestrates [antiSMASH](https://antismash.secondarymetabolites.org/) across a collection of genomes, optionally groups
the resulting BGCs into gene cluster families (GCFs) with [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE), and
produces ready-to-use summary tables for downstream statistical analysis.

---

## Features

- Parallel antiSMASH annotation across any number of genome assemblies or GenBank files
- Per-region tabulation and per-genome BGC count summary
- Optional BiG-SCAPE clustering

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- One of: Docker, Singularity, Podman, or Conda

---

## Quick start

Run with Docker on the bundled example data

```bash
nextflow run exterex/clystere \
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

### BiG-SCAPE

| Parameter                       | Default       | Description                              |
| ------------------------------- | ------------- | ---------------------------------------- |
| `--run_bigscape`                | `false`       | Enable BiG-SCAPE GCF clustering          |
| `--bigscape_gcf_cutoffs`        | `0.3 0.5 0.7` | Space-separated list of distance cutoffs |
| `--bigscape_mix`                | `true`        | Combine all BGC classes into one network |
| `--bigscape_include_singletons` | `true`        | Include singletons in the output         |

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
├── bigscape/              # BiG-SCAPE output (when --run_bigscape)
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
