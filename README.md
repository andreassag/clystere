# clystere

**clystere** is a Nextflow pipeline for end-to-end biosynthetic gene cluster (BGC) and antimicrobial peptide (AMP)
discovery, assembly QC, taxonomic annotation, gene cluster family (GCF) clustering, novelty scoring, comparative
visualization, and metagenomic/metatranscriptomic abundance profiling.

It combines [antiSMASH 8](https://antismash.secondarymetabolites.org/), [GECCO](https://github.com/zellerlab/GECCO),
[deepBGC](https://github.com/Merck/deepbgc), and optional [SanntiS](https://github.com/Finn-Lab/SanntiS) +
[Macrel](https://github.com/BigDataBiology/macrel)/[AMPcombi](https://github.com/paleobiotechnology/AMPcombi), unifies
predictions with [comBGC](https://nf-co.re/funcscan) and cross-tool consensus scoring, and supports downstream
clustering and comparative genomics via [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE),
[BiG-SLiCE](https://github.com/medema-group/bigslice) (including BiG-FAM `--query`),
[clinker](https://github.com/gamcil/clinker), [cblaster](https://github.com/gamcil/cblaster), and
[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2)/[Minimap2](https://github.com/lh3/minimap2) +
[CoverM](https://github.com/wwood/CoverM).

---

## Features

- **Tier 1 — Input flexibility, MAG QC & taxonomy**:
  - Automatic gene-calling for raw unannotated FASTA genomes via `--annotation_tool` (`pyrodigal` default, or `bakta`)
    plus automatic GenBank-to-FASTA extraction when FASTA-requiring tools are enabled
  - Optional assembly contiguity QC via **QUAST** (`--quast_run`) and completeness/contamination via **CheckM2**
    (`--checkm2_run`, auto-downloads CheckM2 DB to `--clystere_db_dir`)
  - Optional **GTDB-Tk** (`--gtdbtk_run`) classification feeding 7-rank lineages directly into BiG-SLiCE's
    `taxonomy.tsv` and the summary/MultiQC tables
- **Tier 2 — Expanded BGC & peptide discovery**:
  - Parallel **antiSMASH 8**, **GECCO**, **deepBGC**, and optional **SanntiS** (`--sanntis_run`) integrated into
    `comBGC.py` and region GBK staging
  - Optional **Macrel + AMPcombi2** (`--amp_run`) antimicrobial peptide discovery
  - **RiPP & RRE-Finder preset** (`--ripp_preset`) enabling antiSMASH RiPP/RRE modules and extracting
    `ripp_candidates.tsv`
- **Tier 3 — Comparative genomics, novelty scoring & visualization**:
  - Co-execution support for **BiG-SCAPE v2** (`--bigscape_run`) and **BiG-SLiCE v2** (`--bigslice_run`), plus **BiG-FAM
    query mode** (`--bigfam_db`) emitting `bigfam_novelty_scores.tsv`
  - Optional **clinker** (`--clinker_run`) interactive HTML synteny plots per multi-member GCF and **cblaster**
    (`--cblaster_run`) offline cohort-wide (default) or external/remote gene-neighborhood searches
  - Default **Cross-Tool Consensus & Novelty Summary** emitting `consensus_loci.tsv`, `genome_qc_bgc_summary.tsv`, and
    MultiQC custom-content tables
- **Tier 4 — Metagenomic & metatranscriptomic ecological context**:
  - Optional pre-annotation **genome dereplication** (`--dereplicate_genomes`, `skani` default or `drep`)
  - Optional **BGC/GCF read mapping & abundance/expression quantification** (`--reads_input` or inline
    `fastq_1`/`fastq_2` columns) using `bowtie2` (default) or `minimap2` + `samtools` + `coverm`

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 26.04.0
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

| Column       | Required | Description                                                                           |
| ------------ | -------- | ------------------------------------------------------------------------------------- |
| `sample`     | Yes      | Unique sample identifier                                                              |
| `genome`     | Yes      | Path to a GenBank/FASTA genome file (`.gbk`, `.gbff`, `.fa`, `.fna`, `.fasta`, `.gz`) |
| `annotation` | No       | Path to an optional GFF3 annotation file (`.gff`, `.gff3`, `.gz`)                     |
| `fastq_1`    | No       | Optional R1 or single-end FASTQ for BGC abundance/expression profiling                |
| `fastq_2`    | No       | Optional R2 paired-end FASTQ                                                          |
| `assay_type` | No       | Optional assay type: `metagenome` (default) or `metatranscriptome`                    |

```csv
sample,genome,annotation
strain_A,data/strain_A.gbff.gz,
strain_B,data/strain_B.fna,
```

Additional cohort read samples can also be passed via `--reads_input` (`sample,fastq_1[,fastq_2,assay_type]`).

### Reference databases

Reference databases for antiSMASH (`--antismash_db`), deepBGC (`--deepbgc_data_dir`), BiG-SLiCE (`--bigslice_db`), Bakta
(`--bakta_db`), and CheckM2 (`--checkm2_db`) are downloaded automatically into `--clystere_db_dir`
(`${launchDir}/.clystere/databases`) if absent. GTDB-Tk requires `--gtdbtk_db` (or `--gtdbtk_download_db true`).

---

## Key Parameters

A full parameter reference is available in [`nextflow_schema.json`](nextflow_schema.json) and
[docs/parameters.md](docs/parameters.md).

### Input, Preprocessing & QC

| Parameter               | Default                            | Description                                                             |
| ----------------------- | ---------------------------------- | ----------------------------------------------------------------------- |
| `--input`               | —                                  | Path to samplesheet CSV (required)                                      |
| `--reads_input`         | `null`                             | Optional reads samplesheet CSV for BGC abundance/expression             |
| `--outdir`              | `results`                          | Directory for all pipeline outputs                                      |
| `--clystere_db_dir`     | `${launchDir}/.clystere/databases` | Directory for cached reference databases                                |
| `--annotation_tool`     | `pyrodigal`                        | Auto gene-caller for raw FASTA (`pyrodigal`, `bakta`, `none`)           |
| `--dereplicate_genomes` | `false`                            | Dereplicate redundant genomes via `--dereplicate_tool` (`skani`/`drep`) |
| `--quast_run`           | `false`                            | Run QUAST assembly QC                                                   |
| `--checkm2_run`         | `false`                            | Run CheckM2 completeness/contamination QC                               |
| `--gtdbtk_run`          | `false`                            | Run GTDB-Tk taxonomy classification (requires `--gtdbtk_db`)            |

### BGC & Peptide Discovery

| Parameter                      | Default | Description                                                         |
| ------------------------------ | ------- | ------------------------------------------------------------------- |
| `--antismash_minimal`          | `true`  | Run antiSMASH in minimal mode; enable modules individually          |
| `--antismash_cb_knownclusters` | `false` | Run KnownClusterBlast; adds similarity columns to `all_regions.tsv` |
| `--gecco_run`                  | `true`  | Run GECCO BGC prediction                                            |
| `--deepbgc_run`                | `true`  | Run deepBGC BGC prediction                                          |
| `--sanntis_run`                | `false` | Run SanntiS as a 4th neural-network BGC predictor                   |
| `--amp_run`                    | `false` | Run Macrel + AMPcombi2 antimicrobial peptide discovery              |
| `--ripp_preset`                | `false` | Enable antiSMASH RiPP/RRE modules and emit `ripp_candidates.tsv`    |

### Clustering, Comparative Genomics & Abundance

| Parameter            | Default   | Description                                                              |
| -------------------- | --------- | ------------------------------------------------------------------------ |
| `--bigscape_run`     | `false`   | Enable BiG-SCAPE v2 GCF clustering                                       |
| `--bigslice_run`     | `false`   | Enable BiG-SLiCE v2 GCF clustering (co-executable with `--bigscape_run`) |
| `--bigfam_db`        | `null`    | Path to BiG-FAM/BiG-SLiCE DB for `bigslice --query` novelty scoring      |
| `--clinker_run`      | `false`   | Generate interactive clinker HTML synteny plots for multi-member GCFs    |
| `--cblaster_run`     | `false`   | Run cblaster gene-neighborhood search (`local` cohort or `remote`)       |
| `--abundance_mapper` | `bowtie2` | Read aligner for BGC/GCF abundance (`bowtie2` or `minimap2`)             |
| `--run_tabulation`   | `true`    | Generate antiSMASH summaries + cross-tool `consensus_loci.tsv`           |

---

- Usage, parameter reference, and output formats are documented in [`docs/`](docs/index.md).

---

## Citations

Please cite the pipeline and its dependencies. See [CITATIONS.md](CITATIONS.md) for full references.

---

## Licence

[MIT](LICENSE)
