# Output

All results are written to `--outdir` (default: `results`). The directory structure is:

```text
results/
├── preprocessing/              # only when unannotated FASTA inputs are gene-called
│   ├── pyrodigal/              # Pyrodigal GFF3/FAA/FNA outputs
│   └── bakta/<sample>/         # Bakta GBFF/GFF3/FAA/TSV outputs
├── dereplication/              # only when --dereplicate_genomes
│   ├── dereplication_clusters.tsv
│   └── dereplicated_genome_ids.txt
├── qc/
│   ├── quast/                  # only when --quast_run
│   └── checkm2/                # only when --checkm2_run
├── taxonomy/
│   └── gtdbtk/                 # only when --gtdbtk_run
├── antismash/<sample>/         # Full antiSMASH 8 output per genome
├── gecco/<sample>/             # GECCO outputs per genome
├── deepbgc/<sample>/           # deepBGC outputs per genome
├── sanntis/<sample>/           # only when --sanntis_run (*.sanntis.full.gff)
├── amp/                        # only when --amp_run
│   ├── macrel/                 # Macrel smORF & AMP predictions
│   └── ampcombi/               # Ampcombi_summary.tsv & clustered AMPs
├── combgc/<sample>/
│   ├── combgc_summary.tsv      # Unified summary across up to 4 BGC predictors
│   └── <sample>_regions/       # Sample-prefixed region GBKs staged for clustering
├── bigscape/                   # only when --bigscape_run (bigscape.db, networks)
├── bigscape_dereplicate/       # only when --bigscape_run and --bigscape_dereplicate
├── bigslice/all_samples/result/# only when --bigslice_run (data.db, tsv_export/)
├── bigfam_query/               # only when --bigfam_db (bigfam_novelty_scores.tsv)
├── comparative/                # gcf_membership_summary.tsv
├── clinker/                    # only when --clinker_run (<gcf_id>.html, *_alignments.tsv)
├── cblaster/                   # only when --cblaster_run (cohort.dmnd/sqlite3, *_cblaster.tsv)
├── abundance/                  # only when --reads_input or fastq_1/2 provided
│   ├── bgc_catalog.fna
│   ├── bgc_catalog_metadata.tsv
│   ├── bgc_abundance_tpm.tsv
│   ├── bgc_abundance_rpkm.tsv
│   ├── bgc_covered_fraction.tsv
│   ├── gcf_abundance_tpm.tsv
│   ├── gcf_abundance_rpkm.tsv
│   ├── bgc_activity_summary.tsv
│   └── coverm/
├── summary/
│   ├── all_regions.tsv         # Per-region antiSMASH table
│   ├── region_counts.tsv       # Per-genome/contig antiSMASH counts
│   ├── ripp_candidates.tsv     # only when --ripp_preset or RiPP/RRE modules enabled
│   ├── consensus_loci.tsv      # Cross-tool unified BGC loci with support tiers & novelty
│   └── genome_qc_bgc_summary.tsv # Per-genome QC, taxonomy, predictor counts & BGC density
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

## `summary/consensus_loci.tsv`

One row per merged biosynthetic locus across all active BGC predictors (`antiSMASH`, `GECCO`, `deepBGC`, `SanntiS`),
clustered by coordinate interval Jaccard similarity (`>= --consensus_min_overlap`, default `0.25`):

- `support_tier`: `High (3-4 tools)`, `Medium (2 tools)`, or `Single-tool (1 tool)`
- `mean_jaccard`: Mean pairwise interval Jaccard overlap across supporting predictors
- Joins KnownClusterBlast hits, MIBiG IDs, BiG-SCAPE & BiG-SLiCE GCF IDs, and BiG-FAM `distance_d` / `is_putative_novel`
  flags

---

## `summary/genome_qc_bgc_summary.tsv`

One row per genome joining Tier 1 QC and taxonomy (`gtdb_classification`, `completeness`, `contamination`,
`assembly_length_bp`, `num_contigs`, `n50`), per-predictor BGC counts (`antismash_bgcs`, `gecco_bgcs`, `deepbgc_bgcs`,
`sanntis_bgcs`), consensus tier breakdowns, putative novel counts, assembly-normalized BGC density (`bgcs_per_mbp`,
`complete_bgc_fraction`, `contig_edge_rate`), and AMP candidate counts (`amp_count`).

---

## `summary/ripp_candidates.tsv`

Emitted when `--ripp_preset` or any antiSMASH RiPP/RRE module is active. Extracts `ripp_subclass`, `locus_tag`,
`rre_bitscore`, `leader_sequence`, `core_sequence`, `predicted_molecular_weight`, and `tailoring_modifications` from
antiSMASH 8 JSON outputs.

---

## `abundance/`

Emitted when metagenomic or metatranscriptomic reads are supplied via `--reads_input` or inline `fastq_1`/`fastq_2`
columns in `--input`:

- `bgc_abundance_tpm.tsv`, `bgc_abundance_rpkm.tsv`, `bgc_covered_fraction.tsv`: BGC-by-sample abundance matrices
- `gcf_abundance_tpm.tsv`, `gcf_abundance_rpkm.tsv`: GCF-level summed abundance matrices
- `bgc_activity_summary.tsv`: Cohort prevalence, mean DNA TPM, mean RNA TPM, and `rna_dna_tpm_ratio`

---

## `multiqc/`

Aggregated MultiQC HTML report summarizing Genome QC/Taxonomy & BGC Overview, Cross-Tool Consensus Support, antiSMASH
BGC counts, QUAST/CheckM2 metrics, AMP candidate counts, BGC abundance summary, and software versions.
