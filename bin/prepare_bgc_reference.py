#!/usr/bin/env python3
"""Build a non-redundant nucleotide BGC catalog FASTA and metadata table from region GBKs."""

import argparse
import csv
import hashlib
from pathlib import Path
import re

from Bio import SeqIO

METADATA_COLS = ["bgc_id", "sample_id", "source_file", "length_bp", "product", "gcf_id"]


def infer_sample_id(filename: str) -> str:
    stem = Path(filename).stem
    for sep in (".sanntis.", ".deepbgc.", ".gecco.", ".region", "_cluster_"):
        if sep in stem:
            return stem.split(sep, 1)[0]
    match = re.match(r"^([^.]+)", stem)
    return match.group(1) if match else stem


def extract_product(record) -> str:
    for feat in record.features:
        if feat.type in {"region", "protocluster", "cluster"}:
            products = feat.qualifiers.get("product") or feat.qualifiers.get("product_class") or []
            if products:
                return ";".join(products)
    return "Other"


def load_gcf_map(gcf_membership: Path | None) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if gcf_membership is None or not gcf_membership.is_file() or gcf_membership.stat().st_size == 0:
        return mapping
    with gcf_membership.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            bgc_file = (row.get("bgc_file") or "").strip()
            gcf_id = (row.get("gcf_id") or "").strip()
            if bgc_file and gcf_id and bgc_file not in mapping:
                mapping[bgc_file] = gcf_id
    return mapping


def _walk_gbks(directory: Path) -> list[Path]:
    """Recursively find .gbk/.gbff/.gb files, following symlinked directories."""
    found: list[Path] = []
    for entry in sorted(directory.iterdir()):
        if entry.name.startswith("."):
            continue
        if entry.is_dir():
            found.extend(_walk_gbks(entry))
        elif entry.is_file() and entry.suffix in {".gbk", ".gbff", ".gb"}:
            found.append(entry)
    return found


def collect_gbks(paths: list[Path]) -> list[Path]:
    found: list[Path] = []
    for p in sorted(paths):
        if p.is_dir():
            found.extend(_walk_gbks(p))
        elif p.is_file() and p.suffix in {".gbk", ".gbff", ".gb"}:
            found.append(p)
    return found


def build_catalog(
    gbk_inputs: list[Path],
    fasta_out: Path,
    metadata_out: Path,
    derep_dir: Path | None = None,
    gcf_membership: Path | None = None,
) -> int:
    """Extract deduplicated BGC nucleotide sequences into a FASTA catalog + metadata TSV."""
    source_gbks: list[Path] = []
    if derep_dir is not None and derep_dir.is_dir():
        rep_dir = derep_dir / "representative_clusters"
        target = rep_dir if rep_dir.is_dir() else derep_dir
        source_gbks = collect_gbks([target])
    if not source_gbks:
        source_gbks = collect_gbks(gbk_inputs)

    gcf_map = load_gcf_map(gcf_membership)
    seen_hashes: set[str] = set()
    used_ids: set[str] = set()
    meta_rows: list[dict[str, str | int]] = []

    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    metadata_out.parent.mkdir(parents=True, exist_ok=True)

    with fasta_out.open("w", encoding="utf-8") as fa_handle:
        for gbk in source_gbks:
            if gbk.stat().st_size == 0:
                continue
            try:
                records = list(SeqIO.parse(gbk, "genbank"))
            except Exception:
                continue
            sample_id = infer_sample_id(gbk.name)
            gcf_id = gcf_map.get(gbk.name, "NA")
            for idx, rec in enumerate(records, start=1):
                try:
                    seq_str = str(rec.seq).upper().strip()
                except Exception:
                    continue
                if not seq_str:
                    continue
                seq_hash = hashlib.sha256(seq_str.encode("utf-8")).hexdigest()
                if seq_hash in seen_hashes:
                    continue
                seen_hashes.add(seq_hash)

                base_id = gbk.stem if len(records) == 1 else f"{gbk.stem}_rec{idx}"
                safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", base_id)
                candidate = safe_id
                suffix_n = 2
                while candidate in used_ids:
                    candidate = f"{safe_id}_{suffix_n}"
                    suffix_n += 1
                used_ids.add(candidate)

                product = extract_product(rec)
                fa_handle.write(f">{candidate} sample={sample_id} gcf={gcf_id} product={product}\n")
                for pos in range(0, len(seq_str), 80):
                    fa_handle.write(seq_str[pos : pos + 80] + "\n")

                meta_rows.append(
                    {
                        "bgc_id": candidate,
                        "sample_id": sample_id,
                        "source_file": gbk.name,
                        "length_bp": len(seq_str),
                        "product": product,
                        "gcf_id": gcf_id,
                    }
                )

    with metadata_out.open("w", encoding="utf-8", newline="") as meta_handle:
        writer = csv.DictWriter(meta_handle, fieldnames=METADATA_COLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(meta_rows)

    return len(meta_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build non-redundant BGC catalog FASTA from region GBK files.")
    parser.add_argument("--gbks", type=Path, nargs="*", default=[], help="Region GBK files or directories.")
    parser.add_argument("--derep-dir", type=Path, default=None, help="Optional BiG-SCAPE dereplicate output dir.")
    parser.add_argument("--gcf-membership", type=Path, default=None, help="Optional gcf_membership_summary.tsv.")
    parser.add_argument("--out-fasta", type=Path, default=Path("bgc_catalog.fna"), help="Output catalog FASTA path.")
    parser.add_argument(
        "--out-metadata",
        type=Path,
        default=Path("bgc_catalog_metadata.tsv"),
        help="Output catalog metadata TSV path.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    build_catalog(
        gbk_inputs=args.gbks,
        fasta_out=args.out_fasta,
        metadata_out=args.out_metadata,
        derep_dir=args.derep_dir,
        gcf_membership=args.gcf_membership,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
