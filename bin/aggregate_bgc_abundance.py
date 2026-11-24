#!/usr/bin/env python3
"""Aggregate CoverM contig outputs into BGC and GCF abundance/expression matrices."""

import argparse
from collections import defaultdict
import csv
from pathlib import Path
import re


def normalize_coverm_sample(file_path: Path) -> str:
    """Derive sample_id from a CoverM depth TSV filename."""
    name = file_path.name
    for suffix in (".depth.tsv", ".coverm.tsv", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return file_path.stem


def load_assay_manifest(manifest_path: Path | None) -> dict[str, str]:
    """Load sample_id -> assay_type ('metagenome' | 'metatranscriptome') mapping."""
    assays: dict[str, str] = {}
    if manifest_path is None or not manifest_path.is_file() or manifest_path.stat().st_size == 0:
        return assays
    with manifest_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sid = (row.get("sample_id") or row.get("sample") or "").strip()
            assay = (row.get("assay_type") or "metagenome").strip().lower()
            if sid:
                assays[sid] = assay or "metagenome"
    return assays


def load_catalog_metadata(metadata_path: Path | None) -> dict[str, dict[str, str]]:
    """Load bgc_id -> {sample_id, product, gcf_id, length_bp} metadata."""
    meta: dict[str, dict[str, str]] = {}
    if metadata_path is None or not metadata_path.is_file() or metadata_path.stat().st_size == 0:
        return meta
    with metadata_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            bgc_id = (row.get("bgc_id") or "").strip()
            if bgc_id:
                meta[bgc_id] = {
                    "sample_id": (row.get("sample_id") or "NA").strip(),
                    "product": (row.get("product") or "Other").strip(),
                    "gcf_id": (row.get("gcf_id") or "NA").strip(),
                    "length_bp": (row.get("length_bp") or "0").strip(),
                }
    return meta


def load_gcf_override(gcf_path: Path | None) -> dict[str, str]:
    """Load bgc_stem -> gcf_id mapping from gcf_membership_summary.tsv."""
    mapping: dict[str, str] = {}
    if gcf_path is None or not gcf_path.is_file() or gcf_path.stat().st_size == 0:
        return mapping
    with gcf_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            bgc_file = (row.get("bgc_file") or "").strip()
            gcf_id = (row.get("gcf_id") or "").strip()
            if bgc_file and gcf_id:
                stem = Path(bgc_file).stem
                mapping.setdefault(stem, gcf_id)
                mapping.setdefault(bgc_file, gcf_id)
    return mapping


def _classify_coverm_col(header: str) -> str | None:
    lower = header.strip().lower()
    if re.search(r"\btpm$", lower):
        return "tpm"
    if re.search(r"\brpkm$", lower):
        return "rpkm"
    if re.search(r"\bmean$", lower):
        return "mean"
    if "covered fraction" in lower or "covered_fraction" in lower:
        return "covered_fraction"
    if "read count" in lower or "read_count" in lower or re.search(r"\bcount$", lower):
        return "read_count"
    return None


def parse_coverm_file(coverm_path: Path) -> tuple[str, dict[str, dict[str, float]]]:
    """Parse one CoverM contig output TSV into (sample_id, bgc_id -> metric dict)."""
    sample_id = normalize_coverm_sample(coverm_path)
    metrics_by_bgc: dict[str, dict[str, float]] = {}
    if not coverm_path.is_file() or coverm_path.stat().st_size == 0:
        return sample_id, metrics_by_bgc

    with coverm_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return sample_id, metrics_by_bgc
        contig_col = reader.fieldnames[0]
        col_map: dict[str, str] = {}
        for col in reader.fieldnames[1:]:
            metric_key = _classify_coverm_col(col)
            if metric_key is not None and metric_key not in col_map:
                col_map[metric_key] = col

        for row in reader:
            raw_contig = (row.get(contig_col) or "").strip()
            if not raw_contig:
                continue
            bgc_id = raw_contig.split()[0]
            values: dict[str, float] = {}
            for metric_key, col_name in col_map.items():
                try:
                    values[metric_key] = float(row.get(col_name, 0.0) or 0.0)
                except ValueError:
                    values[metric_key] = 0.0
            metrics_by_bgc[bgc_id] = values
    return sample_id, metrics_by_bgc


def write_matrix(path: Path, row_key: str, samples: list[str], data: dict[str, dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_key, *samples])
        for item_id in sorted(data):
            writer.writerow([item_id, *[f"{data[item_id].get(s, 0.0):.4f}" for s in samples]])


def aggregate_abundance(
    coverm_files: list[Path],
    outdir: Path,
    catalog_metadata: Path | None = None,
    gcf_membership: Path | None = None,
    assay_manifest: Path | None = None,
    min_covered_fraction: float = 0.5,
) -> None:
    """Aggregate per-sample CoverM outputs into BGC/GCF matrices and activity summaries."""
    meta = load_catalog_metadata(catalog_metadata)
    gcf_override = load_gcf_override(gcf_membership)
    assays = load_assay_manifest(assay_manifest)

    tpm_mat: dict[str, dict[str, float]] = defaultdict(dict)
    rpkm_mat: dict[str, dict[str, float]] = defaultdict(dict)
    cov_mat: dict[str, dict[str, float]] = defaultdict(dict)
    mean_mat: dict[str, dict[str, float]] = defaultdict(dict)
    samples: list[str] = []

    for cpath in sorted(coverm_files):
        sample_id, parsed = parse_coverm_file(cpath)
        if sample_id not in samples:
            samples.append(sample_id)
        for bgc_id, vals in parsed.items():
            tpm_mat[bgc_id][sample_id] = vals.get("tpm", 0.0)
            rpkm_mat[bgc_id][sample_id] = vals.get("rpkm", 0.0)
            cov_mat[bgc_id][sample_id] = vals.get("covered_fraction", 0.0)
            mean_mat[bgc_id][sample_id] = vals.get("mean", 0.0)

    # Ensure every catalog BGC appears even if zero reads mapped
    for bgc_id in meta:
        tpm_mat.setdefault(bgc_id, {})
        rpkm_mat.setdefault(bgc_id, {})
        cov_mat.setdefault(bgc_id, {})
        mean_mat.setdefault(bgc_id, {})

    samples.sort()
    write_matrix(outdir / "bgc_abundance_tpm.tsv", "bgc_id", samples, tpm_mat)
    write_matrix(outdir / "bgc_abundance_rpkm.tsv", "bgc_id", samples, rpkm_mat)
    write_matrix(outdir / "bgc_covered_fraction.tsv", "bgc_id", samples, cov_mat)

    # Aggregate by GCF
    gcf_tpm: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    gcf_rpkm: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for bgc_id in sorted(tpm_mat):
        gcf_id = gcf_override.get(bgc_id) or meta.get(bgc_id, {}).get("gcf_id", "NA")
        if not gcf_id or gcf_id == "NA":
            gcf_id = f"singleton_{bgc_id}"
        for s in samples:
            gcf_tpm[gcf_id][s] += tpm_mat[bgc_id].get(s, 0.0)
            gcf_rpkm[gcf_id][s] += rpkm_mat[bgc_id].get(s, 0.0)

    write_matrix(outdir / "gcf_abundance_tpm.tsv", "gcf_id", samples, gcf_tpm)
    write_matrix(outdir / "gcf_abundance_rpkm.tsv", "gcf_id", samples, gcf_rpkm)

    # DNA vs RNA ecological activity summary
    dna_samples = [s for s in samples if assays.get(s, "metagenome") != "metatranscriptome"]
    rna_samples = [s for s in samples if assays.get(s, "metagenome") == "metatranscriptome"]

    activity_path = outdir / "bgc_activity_summary.tsv"
    activity_cols = [
        "bgc_id",
        "gcf_id",
        "product",
        "prevalence_samples",
        "mean_dna_tpm",
        "mean_rna_tpm",
        "rna_dna_tpm_ratio",
    ]
    with activity_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=activity_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for bgc_id in sorted(tpm_mat):
            gcf_id = gcf_override.get(bgc_id) or meta.get(bgc_id, {}).get("gcf_id", "NA")
            product = meta.get(bgc_id, {}).get("product", "Other")
            prevalence = sum(1 for s in samples if cov_mat[bgc_id].get(s, 0.0) >= min_covered_fraction)
            mean_dna = sum(tpm_mat[bgc_id].get(s, 0.0) for s in dna_samples) / len(dna_samples) if dna_samples else None
            mean_rna = sum(tpm_mat[bgc_id].get(s, 0.0) for s in rna_samples) / len(rna_samples) if rna_samples else None
            if mean_dna is not None and mean_rna is not None and mean_dna > 0:
                ratio = f"{mean_rna / mean_dna:.4f}"
            else:
                ratio = "NA"
            writer.writerow(
                {
                    "bgc_id": bgc_id,
                    "gcf_id": gcf_id,
                    "product": product,
                    "prevalence_samples": prevalence,
                    "mean_dna_tpm": f"{mean_dna:.4f}" if mean_dna is not None else "NA",
                    "mean_rna_tpm": f"{mean_rna:.4f}" if mean_rna is not None else "NA",
                    "rna_dna_tpm_ratio": ratio,
                }
            )

    # MultiQC per-sample summary table
    mqc_path = outdir / "bgc_abundance_mqc.tsv"
    mqc_cols = ["sample_id", "assay_type", "detected_bgcs", "total_tpm", "mean_depth"]
    with mqc_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=mqc_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for s in samples:
            detected = sum(1 for bgc_id in cov_mat if cov_mat[bgc_id].get(s, 0.0) >= min_covered_fraction)
            total_tpm = sum(tpm_mat[bgc_id].get(s, 0.0) for bgc_id in tpm_mat)
            depths = [mean_mat[bgc_id].get(s, 0.0) for bgc_id in mean_mat]
            mean_depth = sum(depths) / len(depths) if depths else 0.0
            writer.writerow(
                {
                    "sample_id": s,
                    "assay_type": assays.get(s, "metagenome"),
                    "detected_bgcs": detected,
                    "total_tpm": f"{total_tpm:.2f}",
                    "mean_depth": f"{mean_depth:.3f}",
                }
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Aggregate CoverM contig outputs into BGC/GCF abundance matrices.")
    parser.add_argument("--coverm", type=Path, nargs="*", default=[], help="CoverM *.depth.tsv output files.")
    parser.add_argument("--metadata", type=Path, default=None, help="Optional bgc_catalog_metadata.tsv.")
    parser.add_argument("--gcf-membership", type=Path, default=None, help="Optional gcf_membership_summary.tsv.")
    parser.add_argument("--manifest", type=Path, default=None, help="Optional TSV mapping sample_id -> assay_type.")
    parser.add_argument(
        "--min-covered-fraction",
        type=float,
        default=0.5,
        help="Minimum covered fraction to call a BGC detected in a sample (default: 0.5).",
    )
    parser.add_argument("-o", "--outdir", type=Path, default=Path(), help="Output directory.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    aggregate_abundance(
        coverm_files=args.coverm,
        outdir=args.outdir,
        catalog_metadata=args.metadata,
        gcf_membership=args.gcf_membership,
        assay_manifest=args.manifest,
        min_covered_fraction=args.min_covered_fraction,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
