# Output

All results are written to `--outdir` (default: `results`). The directory structure is:

```text
results/
├── antismash/
│   └── <sample>/               # Full antiSMASH output per genome (<sample>.json, <sample>.gbk, *region*.gbk)
├── gecco/
│   └── <sample>/               # GECCO outputs per genome (*.clusters.tsv, *.features.tsv, *.region*.gbk)
├── deepbgc/
│   └── <sample>/               # deepBGC outputs per genome (*.bgc.tsv, *.bgc.gbk, *.full.gbk)
├── combgc/
│   └── <sample>/
│       ├── combgc_summary.tsv
│       └── <sample>_regions/   # Sample-prefixed region GBKs staged for clustering
├── bigscape/                   # only when --bigscape_run
│   ├── bigscape.db
│   └── output_files/           # GCF networks and annotation tables
├── bigscape_dereplicate/       # only when --bigscape_run and --bigscape_dereplicate
├── bigslice/                   # only when --bigslice_run
│   └── all_samples/result/     # SQLite database and tsv_export/
├── summary/
│   ├── all_regions.tsv         # per-BGC-region table
│   └── region_counts.tsv       # per-genome (or per-contig) BGC count table
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    ├── pipeline_dag.html
    └── params_*.json
```

---

## `summary/all_regions.tsv`

One row per biosynthetic region across all input genomes. Produced by `SUMMARY` (`tabulate_regions.py`).

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

!!! note "KnownClusterBlast columns" The three `knownclusterblast_*` columns are only included when
`--antismash_cb_knownclusters` is set.

---

## `summary/region_counts.tsv`

One row per genome assembly (or contig when `--count_per_contig`). Produced by `SUMMARY` (`count_regions.py`).

| Column        | Description                                                                                      |
| ------------- | ------------------------------------------------------------------------------------------------ |
| `record`      | Genome assembly name (or `genome\|contig` when `--count_per_contig`)                             |
| `total_count` | Total number of BGC regions regardless of class                                                  |
| `<bgc_type>`  | Count of regions of that product type (one column per distinct type observed across all genomes) |
| `hybrid`      | Count of multi-product regions (absent when `--split_hybrids`)                                   |
| `description` | Sequence description; includes `[N total records]` suffix for multi-contig assemblies            |

---

## `combgc/<sample>/`

Aggregated per-sample BGC summary (`combgc_summary.tsv`) generated across active predictors (`antiSMASH`, `GECCO`, and
`deepBGC`), along with `<sample>_regions/` containing sample-prefixed region GenBank files (`*.gbk`) staged for
clustering.

---

## `bigscape/` & `bigscape_dereplicate/`

Standard BiG-SCAPE v2 output generated from the staged `<sample>_regions/` directories, including `bigscape/bigscape.db`
and network/annotation files. When `--bigscape_dereplicate` is enabled (default), representative clusters are first
written to `bigscape_dereplicate/`.

---

## `bigslice/`

Standard BiG-SLiCE output generated from the staged region GBK files. Includes `result/data.db` and
`result/tsv_export/`.

---

## `multiqc/`

Aggregated MultiQC HTML report and associated data tables/plots summarizing software versions and `region_counts.tsv`.
