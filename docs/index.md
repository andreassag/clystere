# clystere

**clystere** is a Nextflow pipeline for end-to-end biosynthetic gene cluster (BGC) and antimicrobial peptide (AMP)
discovery, assembly QC, taxonomic annotation, gene cluster family (GCF) clustering, comparative genomics, and
metagenomic/metatranscriptomic abundance quantification.

---

## Overview

```mermaid
flowchart TD
    IN["Samplesheet (--input)<br/>GenBank or raw FASTA"] --> PRE["1. GENOME_PREPROCESSING<br/>GenBank->FASTA, skani/dRep, Pyrodigal/Bakta"]
    PRE --> QCTAX["2. GENOME_QC_TAXONOMY<br/>QUAST, CheckM2, GTDB-Tk"]
    PRE --> ANN["3. BGC_ANNOTATION<br/>antiSMASH 8 (+RiPP preset), GECCO, deepBGC, SanntiS, Macrel+AMPcombi"]
    QCTAX -- "taxonomy.tsv" --> CLUST
    ANN --> CLUST["4. BGC_CLUSTERING<br/>comBGC (4 predictors), BiG-SCAPE v2, BiG-SLiCE v2, BiG-FAM Query"]
    CLUST --> COMP["5. BGC_COMPARATIVE<br/>clinker synteny & cblaster search"]
    READS["Reads (--reads_input / fastq_1,2)"] --> ABUND["6. BGC_ABUNDANCE<br/>Bowtie2/Minimap2 + Samtools + CoverM"]
    CLUST --> ABUND
    ANN --> SUM["7. SUMMARY & CONSENSUS<br/>all_regions, region_counts, ripp_candidates, consensus_loci, genome_qc_bgc_summary"]
    CLUST --> SUM
    QCTAX --> SUM
    ABUND --> MQC["8. MULTIQC"]
    COMP --> MQC
    SUM --> MQC
```

---

## Features

- **Tier 1 (Input Flexibility, MAG QC & Taxonomy)**: Automatic Pyrodigal/Bakta pre-annotation on raw FASTA genomes,
  GenBank-to-FASTA extraction, optional QUAST, CheckM2, and GTDB-Tk feeding BiG-SLiCE taxonomy and summaries.
- **Tier 2 (Expanded BGC & Peptide Discovery)**: Parallel antiSMASH 8, GECCO, deepBGC, optional SanntiS (4th predictor),
  Macrel + AMPcombi2 antimicrobial peptides, and `--ripp_preset` with `ripp_candidates.tsv`.
- **Tier 3 (Comparative Genomics & Novelty Scoring)**: Co-executable BiG-SCAPE v2 + BiG-SLiCE v2, BiG-FAM `--query`
  novelty scoring, clinker interactive HTML synteny plots, cblaster cohort/remote neighborhood search, and cross-tool
  consensus locus summary.
- **Tier 4 (Ecological Context)**: Pre-annotation genome dereplication (`skani` / `dRep`) and metagenomic /
  metatranscriptomic BGC & GCF read-mapping matrices (`Bowtie2` / `Minimap2` + `SAMtools` + `CoverM`).

---

## Quick start

```bash
nextflow run andreassag/clystere \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile docker
```

See [Installation](installation.md) and [Usage](usage.md) for full details.
