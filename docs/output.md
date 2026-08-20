# Output

All results are written to `--outdir` (default: `results`). The directory structure is:

```text
results/
├── antismash/
│   ├── <sample>/           # one directory per input genome
│   │   ├── <sample>.json   # antiSMASH JSON (parsed by tabulation scripts)
│   │   ├── <sample>.gbk    # annotated GenBank output
│   │   └── ...
│   └── ...
├── gecco/
│   ├── <sample>/           # GECCO output directory per genome
│   │   ├── *.clusters.tsv
│   │   ├── *.features.tsv
│   │   └── *.region*.gbk   # generated when BiG-SCAPE or BiG-SLiCE is enabled
│   └── ...
├── deepbgc/
│   ├── <sample>/           # deepBGC output directory per genome
│   │   ├── *.bgc.tsv
│   │   ├── *.full.gbk
│   │   └── *.region*.gbk   # converted to antiSMASH-like region files
│   └── ...
├── combgc/
│   ├── <sample>/
│   │   ├── combgc_summary.tsv
│   │   └── combined_regions/  # representative non-redundant region GBKs used for clustering
│   └── ...
├── bigscape/               # only when --bigscape_run
│   └── output_files/
│       ├── *.network       # GCF network files (one per class + mix)
│       ├── *.tsv           # cluster annotation tables
│       └── ...
├── bigslice/               # only when --bigslice_run
│   ├── result/             # SQLite database and analysis outputs
│   └── ...
├── summary/
│   ├── all_regions.tsv     # per-BGC-region table
│   └── region_counts.tsv   # per-genome BGC count table
├── multiqc/
│   ├── clystere-MultiQC-Report_multiqc_report.html
│   └── clystere-MultiQC-Report_multiqc_report_data/
└── pipeline_info/
    ├── execution_timeline_<timestamp>.html
    ├── execution_report_<timestamp>.html
    ├── execution_trace_<timestamp>.txt
    └── pipeline_dag_<timestamp>.html
```

---

## `summary/all_regions.tsv`

One row per biosynthetic region across all input genomes. Produced by `TABULATE_REGIONS`.

| Column                         | Description                                                                       |
| ------------------------------ | --------------------------------------------------------------------------------- |
| `file`                         | Input genome identifier (stem of the antiSMASH JSON filename)                     |
| `record_id`                    | Sequence record name (e.g. NCBI accession)                                        |
| `region`                       | Region number within the record                                                   |
| `start`                        | Start coordinate (bp, 0-based)                                                    |
| `end`                          | End coordinate (bp)                                                               |
| `contig_edge`                  | `True` if the region extends to the edge of the contig (potentially truncated)    |
| `product`                      | BGC product class(es), separated by `/` for hybrids                               |
| `knownclusterblast_hit`        | Top KnownClusterBlast hit name (only present when `--antismash_cb_knownclusters`) |
| `knownclusterblast_accession`  | MIBiG accession of the top hit                                                    |
| `knownclusterblast_similarity` | Similarity category: `high` (>75 %), `medium` (>50 %), or `low` (>15 %)           |
| `record_desc`                  | Full sequence record description                                                  |

!!! note The three `knownclusterblast_*` columns are only included when `--antismash_cb_knownclusters` is set.

---

## `summary/region_counts.tsv`

One row per genome assembly (or contig when `--count_per_contig`). Produced by `COUNT_REGIONS`.

| Column        | Description                                                                                      |
| ------------- | ------------------------------------------------------------------------------------------------ |
| `record`      | Genome assembly name (or `genome\|contig` when `--count_per_contig`)                             |
| `total_count` | Total number of BGC regions regardless of class                                                  |
| `<bgc_type>`  | Count of regions of that product type (one column per distinct type observed across all genomes) |
| `hybrid`      | Count of multi-product regions (absent when `--split_hybrids`)                                   |
| `description` | Sequence description; includes `[N total records]` suffix for multi-contig assemblies            |

---

## `antismash/<sample>/`

Raw antiSMASH output for each genome. The key file is `<sample>.json`, which is read by the tabulation scripts. Other
files (HTML, GenBank, SVG plots) are present depending on the antiSMASH flags used.

---

## `gecco/<sample>/`

Raw GECCO output for each genome. When `--bigscape_run` or `--bigslice_run` is enabled, clystere also runs
`gecco convert gbk --format bigslice` and publishes GECCO region GenBank files compatible with clustering tools.

---

## `deepbgc/<sample>/`

Raw deepBGC output for each genome (`*.bgc.tsv` and `*.full.gbk`). clystere also generates antiSMASH-like
`*.regionNNN.gbk` files using BiG-SLiCE's conversion script to ensure compatibility with BiG-SCAPE/BiG-SLiCE.

---

## `combgc/<sample>/`

Unified per-sample BGC selection generated from antiSMASH, GECCO, and deepBGC predictions. The `combined_regions/`
directory contains representative region GenBank files used as clustering input.

---

## `bigscape/output_files/`

Standard BiG-SCAPE output. The `.network` files are tab-separated edge lists suitable for import into Cytoscape or
Python `networkx`. One network is generated per BGC class plus a `mix` network when `--bigscape_mix` (default).

---

## `bigslice/`

Standard BiG-SLiCE output generated from the unified `combgc/*/combined_regions/` folders. The directory includes the
processed analysis results and SQLite-backed data used for downstream inspection.

---

## `multiqc/`

Aggregated MultiQC summary report (`clystere-MultiQC-Report_multiqc_report.html`) and associated data tables/plots.
Contains:

- Pipeline execution parameters and software tool versions.
- Unified BGC counts and region distributions across antiSMASH, GECCO, and deepBGC.

---

## `pipeline_info/`

Nextflow execution metadata generated per run:

| File                        | Description                                          |
| --------------------------- | ---------------------------------------------------- |
| `execution_timeline_*.html` | Gantt chart of process execution times               |
| `execution_report_*.html`   | Resource usage report (CPU, memory, I/O per process) |
| `execution_trace_*.txt`     | Raw per-task resource trace (parseable TSV)          |
| `pipeline_dag_*.html`       | Directed acyclic graph of the pipeline               |
