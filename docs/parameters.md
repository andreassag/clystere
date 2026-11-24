# Parameters

Full parameter reference for clystere. All parameters can be passed on the command line (`--parameter_name value`) or
set in a params YAML/JSON file (`-params-file params.yml`).

---

## Input / output

| Parameter            | Default                            | Description                                                                |
| -------------------- | ---------------------------------- | -------------------------------------------------------------------------- |
| `--input`            | `null`                             | **Required.** Path to a genome samplesheet CSV.                            |
| `--reads_input`      | `null`                             | Optional reads samplesheet CSV (`sample,fastq_1[,fastq_2,assay_type]`).    |
| `--outdir`           | `results`                          | Output directory for published files.                                      |
| `--clystere_db_dir`  | `${launchDir}/.clystere/databases` | Workspace directory for downloading and caching reference databases.       |
| `--publish_dir_mode` | `copy`                             | How Nextflow stages outputs: `copy`, `symlink`, `move`, `rellink`, `link`. |

---

## Genome preprocessing & dereplication

| Parameter               | Default     | Description                                                                             |
| ----------------------- | ----------- | --------------------------------------------------------------------------------------- |
| `--annotation_tool`     | `pyrodigal` | Automatic gene caller for unannotated FASTA inputs: `pyrodigal`, `bakta`, or `none`.    |
| `--pyrodigal_mode`      | `single`    | Pyrodigal mode: `single` (isolates) or `meta` (metagenomic contigs/MAGs).               |
| `--bakta_db`            | `null`      | Path to existing Bakta DB. Downloaded automatically when `--annotation_tool bakta`.     |
| `--bakta_db_type`       | `light`     | Bakta DB type to download automatically: `light` or `full`.                             |
| `--dereplicate_genomes` | `false`     | Dereplicate redundant genomes/MAGs prior to BGC annotation.                             |
| `--dereplicate_tool`    | `skani`     | Dereplication engine: `skani` (default) or `drep`.                                      |
| `--dereplicate_ani`     | `99.0`      | ANI percentage threshold for dereplication (`99.0` strain-level, `95.0` species-level). |

---

## Assembly QC & taxonomy (Tier 1)

| Parameter                | Default     | Description                                                                                  |
| ------------------------ | ----------- | -------------------------------------------------------------------------------------------- |
| `--quast_run`            | `false`     | Run QUAST assembly contiguity QC (N50, L50, contig counts).                                  |
| `--checkm2_run`          | `false`     | Run CheckM2 completeness and contamination QC.                                               |
| `--checkm2_db`           | `null`      | Path to CheckM2 `.dmnd` file or directory. Auto-downloaded to `--clystere_db_dir` if absent. |
| `--checkm2_db_zenodo_id` | `14897628`  | Zenodo record ID for automatic CheckM2 database download.                                    |
| `--gtdbtk_run`           | `false`     | Run GTDB-Tk taxonomic classification and wire lineages into BiG-SLiCE and summary tables.    |
| `--gtdbtk_db`            | `null`      | Path to extracted GTDB-Tk reference package (required unless `--gtdbtk_download_db true`).   |
| `--gtdbtk_download_db`   | `false`     | Auto-download GTDB-Tk reference archive (~110 GB) if `--gtdbtk_db` is omitted.               |
| `--gtdbtk_db_url`        | GTDB latest | Archive URL used when `--gtdbtk_download_db` is enabled.                                     |
| `--gtdbtk_extra_args`    | `''`        | Additional CLI arguments appended to `gtdbtk classify_wf`.                                   |

---

## antiSMASH — basic options

| Parameter                | Default    | Description                                                                  |
| ------------------------ | ---------- | ---------------------------------------------------------------------------- |
| `--antismash_db`         | `null`     | Path to an antiSMASH database directory. Downloaded automatically if absent. |
| `--antismash_taxon`      | `bacteria` | Taxon preset: `bacteria` or `fungi`.                                         |
| `--antismash_minlength`  | `1000`     | Minimum sequence length (bp) to annotate.                                    |
| `--antismash_limit`      | `-1`       | Maximum number of sequences to annotate (`-1` = no limit).                   |
| `--antismash_minimal`    | `true`     | Run in minimal mode (core detection only).                                   |
| `--antismash_extra_args` | `''`       | Arbitrary extra arguments appended verbatim to the antiSMASH command.        |

---

## antiSMASH — detection modules

These flags enable additional detection modules alongside or instead of minimal mode:

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

## antiSMASH — ClusterBlast, HMM & output

| Parameter                             | Default   | Description                                                          |
| ------------------------------------- | --------- | -------------------------------------------------------------------- |
| `--antismash_cc_mibig`                | `false`   | Compare clusters to MIBiG reference clusters (`--cc-mibig`).         |
| `--antismash_cb_general`              | `false`   | GeneralClusterBlast (`--cb-general`).                                |
| `--antismash_cb_knownclusters`        | `false`   | KnownClusterBlast against MIBiG (adds columns to `all_regions.tsv`). |
| `--antismash_cb_subclusters`          | `false`   | SubClusterBlast (`--cb-subclusters`).                                |
| `--antismash_hmmdetection_strictness` | `relaxed` | Detection strictness: `strict`, `relaxed`, or `loose`.               |
| `--antismash_zip_output`              | `false`   | Archive antiSMASH directory as a zip file.                           |
| `--antismash_allow_long_headers`      | `true`    | Allow sequence headers longer than 80 characters.                    |

---

## GECCO, deepBGC & Tier 2 Expanded Discovery

| Parameter                | Default | Description                                                                                                                                |
| ------------------------ | ------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `--gecco_run`            | `true`  | Run GECCO for each input genome.                                                                                                           |
| `--gecco_extra_args`     | `''`    | Extra CLI arguments passed to GECCO.                                                                                                       |
| `--deepbgc_run`          | `true`  | Run deepBGC for each input genome.                                                                                                         |
| `--deepbgc_data_dir`     | `null`  | Path to deepBGC downloads/models. Auto-downloaded when absent.                                                                             |
| `--deepbgc_extra_args`   | `''`    | Extra CLI arguments passed to deepBGC.                                                                                                     |
| `--sanntis_run`          | `false` | Run SanntiS as a 4th neural-network BGC predictor integrated into `comBGC.py` and clustering.                                              |
| `--sanntis_interpro_dir` | `null`  | Optional directory of pre-computed per-sample InterProScan `.gff3`/`.tsv` files (`<sample_id>*`) passed to SanntiS via `--ip-file`.        |
| `--sanntis_extra_args`   | `''`    | Extra CLI arguments passed to SanntiS.                                                                                                     |
| `--amp_run`              | `false` | Run Macrel + AMPcombi2 antimicrobial peptide discovery.                                                                                    |
| `--amp_database`         | `DRAMP` | Reference AMP database for AMPcombi2 (`DRAMP`, `APD`, `UniRef100`).                                                                        |
| `--amp_database_dir`     | `null`  | Optional local pre-downloaded AMPcombi2 database directory.                                                                                |
| `--ripp_preset`          | `false` | Enable antiSMASH RiPP/RRE modules (`--rre`, lanthi/lasso/sacti/thiopeptides, `--asf`, `--clusterhmmer`) and extract `ripp_candidates.tsv`. |

---

## BiG-SCAPE, BiG-SLiCE & BiG-FAM Query

BiG-SCAPE and BiG-SLiCE can run independently or together in the same execution.

| Parameter                       | Default       | Description                                                                                         |
| ------------------------------- | ------------- | --------------------------------------------------------------------------------------------------- |
| `--bigscape_run`                | `false`       | Enable BiG-SCAPE v2 GCF clustering.                                                                 |
| `--bigscape_dereplicate`        | `true`        | Run `bigscape dereplicate` before clustering.                                                       |
| `--bigscape_dereplicate_cutoff` | `0.8`         | Similarity cutoff used for BiG-SCAPE dereplication.                                                 |
| `--bigscape_pfam_path`          | `null`        | Path to `Pfam-A.hmm`. Resolved automatically from `--antismash_db` if absent.                       |
| `--bigscape_gcf_cutoffs`        | `0.3 0.5 0.7` | Space-separated list of GCF similarity cutoffs.                                                     |
| `--bigscape_mix`                | `true`        | Analyse all BGC classes together in a mixed network (`--mix`).                                      |
| `--bigscape_classify`           | `none`        | Classification mode: `none` or `classify`.                                                          |
| `--bigscape_include_singletons` | `true`        | Include singleton BGCs in output.                                                                   |
| `--bigscape_zip_output`         | `false`       | Archive BiG-SCAPE output as `bigscape.tar.gz`.                                                      |
| `--bigscape_extra_args`         | `''`          | Extra CLI arguments passed to BiG-SCAPE.                                                            |
| `--bigslice_run`                | `false`       | Enable BiG-SLiCE v2 clustering (automatically incorporates GTDB-Tk taxonomy).                       |
| `--bigslice_db`                 | `null`        | Path to `bigslice-models`. Auto-downloaded when absent.                                             |
| `--bigslice_extra_args`         | `''`          | Extra CLI arguments passed to BiG-SLiCE.                                                            |
| `--bigfam_db`                   | `null`        | Path to pre-computed BiG-FAM/BiG-SLiCE DB for `bigslice --query` novelty scoring.                   |
| `--bigfam_n_ranks`              | `1`           | Top N GCF hits returned per BGC in query mode.                                                      |
| `--bigfam_novelty_threshold`    | `900`         | Distance `d`-score above which a queried BGC is flagged as putative novel.                          |
| `--bigfam_no_normalize`         | `true`        | Pass `--no_normalize_feature` for BiG-FAM v1 DBs (`false` + threshold `~0.4` for BiG-SLiCE v2 DBs). |

---

## Comparative genomics & visualization (Tier 3)

| Parameter                | Default | Description                                                                                |
| ------------------------ | ------- | ------------------------------------------------------------------------------------------ |
| `--clinker_run`          | `false` | Generate interactive clinker HTML synteny plots per multi-member GCF.                      |
| `--clinker_max_bgcs`     | `50`    | Maximum BGCs included per GCF plot.                                                        |
| `--clinker_extra_args`   | `''`    | Extra CLI arguments passed to clinker.                                                     |
| `--cblaster_run`         | `false` | Search representative/orphan BGC clusters via cblaster.                                    |
| `--cblaster_mode`        | `local` | Search mode: `local` (offline cohort or `--cblaster_db`) or `remote` (NCBI).               |
| `--cblaster_db`          | `null`  | Pre-built local cblaster `.dmnd` path (paired with `.sqlite3`). Built from cohort if null. |
| `--cblaster_max_queries` | `50`    | Maximum representative/orphan BGCs searched by cblaster.                                   |
| `--cblaster_extra_args`  | `''`    | Extra CLI arguments passed to cblaster.                                                    |

---

## Metagenomic & metatranscriptomic abundance (Tier 4)

| Parameter                  | Default   | Description                                                              |
| -------------------------- | --------- | ------------------------------------------------------------------------ |
| `--abundance_mapper`       | `bowtie2` | Read aligner for BGC/GCF quantification: `bowtie2` or `minimap2`.        |
| `--abundance_min_identity` | `95.0`    | Minimum read percent identity passed to CoverM contig quantification.    |
| `--minimap2_extra_args`    | `-ax sr`  | Extra arguments passed to `minimap2` when `--abundance_mapper minimap2`. |

---

## Tabulation & consensus summary

| Parameter                 | Default | Description                                                                                         |
| ------------------------- | ------- | --------------------------------------------------------------------------------------------------- |
| `--run_tabulation`        | `true`  | Emit `all_regions.tsv`, `region_counts.tsv`, `consensus_loci.tsv`, and `genome_qc_bgc_summary.tsv`. |
| `--count_per_contig`      | `false` | Emit one row per contig in `region_counts.tsv` instead of per assembly.                             |
| `--split_hybrids`         | `false` | Count each product class in hybrid BGCs individually.                                               |
| `--consensus_min_overlap` | `0.25`  | Minimum pairwise interval Jaccard overlap required to merge calls across predictors.                |

---

## Pipeline internals

| Parameter           | Default | Description                            |
| ------------------- | ------- | -------------------------------------- |
| `--monochrome_logs` | `false` | Disable colour in Nextflow log output. |
| `--validate_params` | `true`  | Validate CLI parameters via nf-schema. |
| `--help`            | `false` | Print parameter help and exit.         |
| `--version`         | `false` | Print pipeline version and exit.       |
