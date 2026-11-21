# clystere

**clystere** is a Nextflow pipeline for automated biosynthetic gene cluster (BGC) discovery and comparative analysis. It
runs [antiSMASH](https://antismash.secondarymetabolites.org/), [GECCO](https://github.com/zellerlab/GECCO), and
[deepBGC](https://github.com/Merck/deepbgc) across a collection of genomes, aggregates predictions with
[comBGC](https://nf-co.re/funcscan), and optionally groups resulting BGCs into gene cluster families (GCFs) with
[BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) or [BiG-SLiCE](https://github.com/medema-group/bigslice).

---

## Features

- Parallel antiSMASH + GECCO + deepBGC annotation across pre-annotated GenBank genomes (with optional GFF3 annotation)
- Optional per-predictor toggles (`--gecco_run`, `--deepbgc_run`)
- `comBGC.py`-based summary aggregation and sample-prefixed region GBK staging before clustering
- Per-region tabulation (`all_regions.tsv`) and per-genome/per-contig BGC count summary (`region_counts.tsv`)
- MultiQC HTML report aggregating software versions and BGC metrics
- Optional BiG-SCAPE (with optional `bigscape dereplicate`) or BiG-SLiCE (with automatic `BIGSLICE_DOWNLOADDB`)
  clustering

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 25.10.0
- Java ≥ 17 (required by Nextflow)
- One of: Docker, Apptainer, Singularity, Podman, or Conda

---

## Quick start

Run with Docker or Apptainer on the bundled example data:

```bash
nextflow run andreassag/clystere \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile docker
```

Or run the fast compact test profile:

```bash
nextflow run andreassag/clystere -profile test,apptainer
```

---

## Inputs

### Samplesheet

A comma-separated file with the following columns:

| Column       | Required | Description                                                                |
| ------------ | -------- | -------------------------------------------------------------------------- |
| `sample`     | Yes      | Unique sample identifier                                                   |
| `genome`     | Yes      | Path to a GenBank/FASTA genome file (`.gbk`, `.gbff`, `.fa`, `.fna`, etc.) |
| `annotation` | No       | Path to an optional GFF3 annotation file (`.gff`, `.gff3`, `.gz`)          |

```csv
sample,genome,annotation
strain_A,data/strain_A.gbff.gz,
strain_B,data/strain_B.gbk,data/strain_B.gff3
```

### Reference databases

Reference databases for antiSMASH (`--antismash_db`), deepBGC (`--deepbgc_data_dir`), and BiG-SLiCE (`--bigslice_db`)
are downloaded automatically into `--clystere_db_dir` (`${launchDir}/.clystere/databases`) if not already present.

---

## Key Parameters

A full parameter reference is available in [`nextflow_schema.json`](nextflow_schema.json) and
[docs/parameters.md](docs/parameters.md).

### Input / output

| Parameter           | Default                            | Description                              |
| ------------------- | ---------------------------------- | ---------------------------------------- |
| `--input`           | —                                  | Path to samplesheet CSV (required)       |
| `--outdir`          | `results`                          | Directory for all pipeline outputs       |
| `--clystere_db_dir` | `${launchDir}/.clystere/databases` | Directory for cached reference databases |

### BGC Predictors

| Parameter                      | Default    | Description                                                         |
| ------------------------------ | ---------- | ------------------------------------------------------------------- |
| `--antismash_db`               | —          | Path to antiSMASH database directory                                |
| `--antismash_taxon`            | `bacteria` | Taxonomic scope (`bacteria` or `fungi`)                             |
| `--antismash_minimal`          | `true`     | Run in minimal mode; enable modules individually                    |
| `--antismash_cb_knownclusters` | `false`    | Run KnownClusterBlast; adds similarity columns to `all_regions.tsv` |
| `--gecco_run`                  | `true`     | Run GECCO BGC prediction                                            |
| `--deepbgc_run`                | `true`     | Run deepBGC BGC prediction                                          |
| `--deepbgc_data_dir`           | —          | Path to deepBGC downloads directory                                 |

### Clustering & Tabulation

| Parameter                | Default       | Description                                                            |
| ------------------------ | ------------- | ---------------------------------------------------------------------- |
| `--bigscape_run`         | `false`       | Enable BiG-SCAPE GCF clustering                                        |
| `--bigscape_dereplicate` | `true`        | Run `bigscape dereplicate` before clustering                           |
| `--bigscape_gcf_cutoffs` | `0.3 0.5 0.7` | Space-separated list of distance cutoffs                               |
| `--bigslice_run`         | `false`       | Enable BiG-SLiCE clustering (mutually exclusive with `--bigscape_run`) |
| `--bigslice_db`          | —             | Path to BiG-SLiCE HMM database (`bigslice-models`)                     |
| `--run_tabulation`       | `true`        | Generate `all_regions.tsv` and `region_counts.tsv` via `SUMMARY`       |
| `--count_per_contig`     | `false`       | Report counts per contig rather than per assembly                      |
| `--split_hybrids`        | `false`       | Count each product type in hybrid BGCs separately                      |

---

- Usage, parameter reference, and output formats are documented in [`docs/`](docs/index.md).

---

## Citations

Please cite the pipeline and its dependencies. See [CITATIONS.md](CITATIONS.md) for full references.

---

## Licence

[MIT](LICENSE)
