# Parameters

Full parameter reference for clystere. All parameters can be passed on the command line (`--parameter_name value`) or
set in a params YAML/JSON file (`-params-file params.yml`).

---

## Input / output

| Parameter            | Default   | Description                                                                |
| -------------------- | --------- | -------------------------------------------------------------------------- |
| `--input`            | `null`    | **Required.** Path to a samplesheet CSV.                                   |
| `--outdir`           | `results` | Output directory for published files.                                      |
| `--publish_dir_mode` | `copy`    | How Nextflow stages outputs: `copy`, `symlink`, `move`, `rellink`, `link`. |

---

## antiSMASH — basic options

| Parameter                      | Default    | Description                                                                  |
| ------------------------------ | ---------- | ---------------------------------------------------------------------------- |
| `--antismash_db`               | `null`     | Path to an antiSMASH database directory. Downloaded automatically if absent. |
| `--antismash_taxon`            | `bacteria` | Taxon preset: `bacteria` or `fungi`.                                         |
| `--antismash_genefinding_tool` | `prodigal` | Gene-finding tool: `prodigal`, `prodigal-m`, `none`, or `error`.             |
| `--antismash_minlength`        | `1000`     | Minimum sequence length (bp) to annotate.                                    |
| `--antismash_limit`            | `-1`       | Maximum number of sequences to annotate (`-1` = no limit).                   |
| `--antismash_minimal`          | `true`     | Run in minimal mode (core detection only).                                   |
| `--antismash_accept_failure`   | `false`    | Continue the pipeline even if antiSMASH exits non-zero.                      |
| `--antismash_reuse_results`    | `false`    | Skip antiSMASH and reuse existing results from `--outdir/antismash/`.        |
| `--antismash_extra_args`       | `''`       | Arbitrary extra arguments appended verbatim to the antiSMASH command.        |

---

## antiSMASH — detection modules

These flags enable additional detection modules. All require `--antismash_minimal false` (or set the flag alongside
`--antismash_minimal` to override it):

| Parameter                           | antiSMASH flag            | Description                        |
| ----------------------------------- | ------------------------- | ---------------------------------- |
| `--antismash_enable_nrps_pks`       | `--enable-nrps-pks`       | Full NRPS/PKS analysis             |
| `--antismash_enable_terpene`        | `--enable-terpene`        | Terpene pathway detection          |
| `--antismash_enable_lanthipeptides` | `--enable-lanthipeptides` | Lanthipeptide detection            |
| `--antismash_enable_lassopeptides`  | `--enable-lassopeptides`  | Lasso peptide detection            |
| `--antismash_enable_sactipeptides`  | `--enable-sactipeptides`  | Sactipeptide detection             |
| `--antismash_enable_thiopeptides`   | `--enable-thiopeptides`   | Thiopeptide detection              |
| `--antismash_enable_t2pks`          | `--enable-t2pks`          | Type II PKS detection              |
| `--antismash_enable_tta`            | `--enable-tta`            | TTA codon detection                |
| `--antismash_enable_genefunctions`  | `--enable-genefunctions`  | Gene function annotation           |
| `--antismash_enable_html`           | `--enable-html`           | Generate HTML reports              |
| `--antismash_fullhmmer`             | `--fullhmmer`             | Full-genome Pfam HMM scan          |
| `--antismash_clusterhmmer`          | `--clusterhmmer`          | Cluster-focused Pfam scan          |
| `--antismash_tigrfam`               | `--tigrfam`               | TIGRFam annotation                 |
| `--antismash_cassis`                | `--cassis`                | CASSIS promoter motif detection    |
| `--antismash_asf`                   | `--asf`                   | Active site finder                 |
| `--antismash_pfam2go`               | `--pfam2go`               | Pfam-to-GO mapping                 |
| `--antismash_rre`                   | `--rre`                   | RiPP recognition element detection |
| `--antismash_smcog_trees`           | `--smcog-trees`           | SMCOG phylogenetic trees           |
| `--antismash_tfbs`                  | `--tfbs`                  | Transcription factor binding sites |

---

## antiSMASH — ClusterBlast

| Parameter                           | Default | antiSMASH flag            | Description                                                            |
| ----------------------------------- | ------- | ------------------------- | ---------------------------------------------------------------------- |
| `--antismash_cc_mibig`              | `false` | `--cc-mibig`              | Compare clusters to MIBiG reference clusters                           |
| `--antismash_cb_general`            | `false` | `--cb-general`            | GeneralClusterBlast (NCBI database)                                    |
| `--antismash_cb_knownclusters`      | `false` | `--cb-knownclusters`      | KnownClusterBlast (MIBiG). Also adds hit columns to `all_regions.tsv`. |
| `--antismash_cb_subclusters`        | `false` | `--cb-subclusters`        | SubClusterBlast                                                        |
| `--antismash_cc_custom_dbs`         | `null`  | `--cc-custom-dbs`         | Comma-separated paths to custom ClusterCompare DB config files         |
| `--antismash_cb_nclusters`          | `null`  | `--cb-nclusters`          | Maximum clusters to report (default: 10)                               |
| `--antismash_cb_min_homology_scale` | `null`  | `--cb-min-homology-scale` | Minimum homology scale threshold (0.0–1.0)                             |

---

## antiSMASH — HMM detection

| Parameter                                                  | Default   | Description                                            |
| ---------------------------------------------------------- | --------- | ------------------------------------------------------ |
| `--antismash_hmmdetection_strictness`                      | `relaxed` | Detection strictness: `strict`, `relaxed`, or `loose`. |
| `--antismash_hmmdetection_fungal_cutoff_multiplier`        | `null`    | Override fungal HMM score cutoff multiplier.           |
| `--antismash_hmmdetection_fungal_neighbourhood_multiplier` | `null`    | Override fungal neighbourhood multiplier.              |
| `--antismash_hmmdetection_limit_to_rule_names`             | `null`    | Comma-separated rule names to restrict detection.      |
| `--antismash_hmmdetection_limit_to_rule_categories`        | `null`    | Comma-separated rule categories to restrict detection. |

---

## antiSMASH — output

| Parameter                        | Default | Description                                                             |
| -------------------------------- | ------- | ----------------------------------------------------------------------- |
| `--antismash_zip_output`         | `false` | `true` = force zip; `false` = force no-zip; `null` = antiSMASH default. |
| `--antismash_summary_gbk`        | `null`  | Toggle summary GenBank output (`true`/`false`/`null`).                  |
| `--antismash_region_gbks`        | `null`  | Toggle per-region GenBank files (`true`/`false`/`null`).                |
| `--antismash_allow_long_headers` | `true`  | Allow sequence headers longer than antiSMASH's default limit.           |

---

## GECCO

| Parameter            | Default | Description                                          |
| -------------------- | ------- | ---------------------------------------------------- |
| `--gecco_run`        | `true`  | Run GECCO for each input genome.                     |
| `--gecco_extra_args` | `''`    | Arbitrary extra arguments passed to GECCO.           |

---

## deepBGC

| Parameter              | Default | Description                                                           |
| ---------------------- | ------- | --------------------------------------------------------------------- |
| `--deepbgc_run`        | `true`  | Run deepBGC for each input genome.                                    |
| `--deepbgc_data_dir`   | `null`  | Path to deepBGC downloads/models. Auto-downloaded when absent/empty.  |
| `--deepbgc_extra_args` | `''`    | Arbitrary extra arguments passed to deepBGC.                          |

---

## comBGC unification

These parameters control filtering while combining antiSMASH, GECCO, and deepBGC predictions prior to BiG-SCAPE or
BiG-SLiCE.

| Parameter              | Default | Description                                       |
| ---------------------- | ------- | ------------------------------------------------- |
| `--combgc_min_length`  | `3000`  | Minimum BGC length retained by comBGC.            |
| `--combgc_contig_edge` | `2`     | Exclude BGCs near contig edges in comBGC.         |

---

## BiG-SCAPE

`--bigscape_run` and `--bigslice_run` are mutually exclusive.

BiG-SCAPE in clystere runs on unified comBGC-filtered regions and requires `--gecco_run true` and
`--deepbgc_run true`.

| Parameter                       | Default       | Description                                                                   |
| ------------------------------- | ------------- | ----------------------------------------------------------------------------- |
| `--bigscape_run`                | `false`       | Enable BiG-SCAPE GCF clustering.                                              |
| `--bigscape_dereplicate`        | `true`        | Run `bigscape dereplicate` before clustering.                                 |
| `--bigscape_dereplicate_cutoff` | `0.8`         | Similarity cutoff used for BiG-SCAPE dereplication.                           |
| `--bigscape_pfam_path`          | `null`        | Path to `Pfam-A.hmm`. Resolved automatically from `--antismash_db` if absent. |
| `--bigscape_gcf_cutoffs`        | `0.3 0.5 0.7` | Space-separated list of GCF similarity cutoffs.                               |
| `--bigscape_mix`                | `true`        | Include a mixed-class network alongside class-specific networks.              |
| `--bigscape_classify`           | `none`        | Classification mode: `none` or `classify`.                                    |
| `--bigscape_include_singletons` | `true`        | Include singleton BGCs (no GCF membership) in output.                         |
| `--bigscape_zip_output`         | `false`       | Compress BiG-SCAPE output directory.                                          |
| `--bigscape_extra_args`         | `''`          | Arbitrary extra arguments passed to BiG-SCAPE.                                |

## BiG-SLiCE

BiG-SLiCE requires its HMM model bundle. The pipeline bootstraps this bundle automatically during the first BiG-SLiCE
task execution.

BiG-SLiCE in clystere runs on unified comBGC-filtered regions and requires `--gecco_run true` and
`--deepbgc_run true`.

| Parameter               | Default | Description                                                            |
| ----------------------- | ------- | ---------------------------------------------------------------------- |
| `--bigslice_run`        | `false` | Enable BiG-SLiCE clustering. Mutually exclusive with `--bigscape_run`. |
| `--bigslice_zip_output` | `false` | Compress BiG-SLiCE output directory.                                   |
| `--bigslice_extra_args` | `''`    | Arbitrary extra arguments passed to BiG-SLiCE.                         |

---

## Tabulation

| Parameter            | Default | Description                                                                                           |
| -------------------- | ------- | ----------------------------------------------------------------------------------------------------- |
| `--run_tabulation`   | `true`  | Run `TABULATE_REGIONS` and `COUNT_REGIONS` after antiSMASH.                                           |
| `--count_per_contig` | `false` | Emit one row per contig in `region_counts.tsv` instead of per assembly.                               |
| `--split_hybrids`    | `false` | Count each product class in a hybrid BGC individually. `total_count` is always the number of regions. |

---

## Pipeline internals

| Parameter           | Default | Description                                                                |
| ------------------- | ------- | -------------------------------------------------------------------------- |
| `--max_memory`      | `14.GB` | Hard ceiling on memory per process (overridden by `hpc`/`slurm` profiles). |
| `--max_cpus`        | `8`     | Hard ceiling on CPUs per process.                                          |
| `--max_time`        | `240.h` | Hard ceiling on wall-time per process.                                     |
| `--monochrome_logs` | `false` | Disable colour in Nextflow log output.                                     |
| `--help`            | `false` | Print parameter help and exit.                                             |
| `--version`         | `false` | Print pipeline version and exit.                                           |
