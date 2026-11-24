#!/usr/bin/env python3
"""Prepare plain-text Macrel inputs and a locus-tag-matched GenBank file for AMPcombi2."""

import argparse
from collections import defaultdict
import gzip
from pathlib import Path
import shutil
from typing import TextIO

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature


def _open_text_read(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def unpack_file(src: Path, dst: Path) -> None:
    """Decompress a .gz file or copy a plain file to dst."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as fin, dst.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(src, dst)


def unpack_faa_strip_stops(src: Path, dst: Path) -> None:
    """Decompress a protein FASTA and strip trailing stop-codon '*' characters for ProtParam compatibility."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with _open_text_read(src) as fin, dst.open("w", encoding="utf-8", newline="\n") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                fout.write(line.replace("*", ""))


def load_antismash_products(gbk_path: Path | None) -> dict[str, list[tuple[int, int, str]]]:
    """Index CDS (start0, end0, product) intervals per contig from an optional GenBank file."""
    products: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    if gbk_path is None or not gbk_path.is_file() or gbk_path.stat().st_size == 0:
        return products
    with _open_text_read(gbk_path) as handle:
        for rec in SeqIO.parse(handle, "genbank"):
            entries: list[tuple[int, int, str]] = []
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                prod = (feat.qualifiers.get("product") or [""])[0].strip()
                if not prod:
                    continue
                entries.append((int(feat.location.start), int(feat.location.end), prod))
            if entries:
                products[rec.id].extend(entries)
                if rec.name != rec.id:
                    products[rec.name].extend(entries)
    return products


def _match_product(
    contig_products: list[tuple[int, int, str]],
    start0: int,
    end0: int,
) -> str:
    """Find best overlapping annotated product for a Macrel ORF interval."""
    best_overlap = 0
    best_prod = "hypothetical protein"
    orf_len = max(1, end0 - start0)
    for p_start, p_end, prod in contig_products:
        overlap = min(end0, p_end) - max(start0, p_start)
        if overlap > best_overlap:
            denom = max(orf_len, p_end - p_start, 1)
            if overlap / denom >= 0.5:
                best_overlap = overlap
                best_prod = prod
    return best_prod


def _parse_prodigal_header(description: str) -> tuple[int, int, int] | None:
    """Extract 0-based (start0, end0, strand) from a Prodigal/Pyrodigal FAA header."""
    parts = [part.strip() for part in description.split("#")]
    if len(parts) < 4:
        return None
    try:
        start1 = int(parts[1])
        end1 = int(parts[2])
        strand_raw = int(parts[3])
    except ValueError:
        return None
    if end1 < start1:
        start1, end1 = end1, start1
    strand = -1 if strand_raw < 0 else 1
    return max(0, start1 - 1), end1, strand


def build_macrel_genbank(
    genome_fasta: Path,
    orfs_faa: Path,
    out_gbk: Path,
    antismash_gbk: Path | None = None,
) -> int:
    """Build a GenBank file whose CDS locus_tags match Macrel's ORF Access identifiers."""
    with _open_text_read(genome_fasta) as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    rec_by_id: dict[str, object] = {}
    for rec in records:
        rec.annotations["molecule_type"] = "DNA"
        rec.name = rec.id[:16]
        rec.features = []
        rec_by_id[rec.id] = rec

    products_by_contig = load_antismash_products(antismash_gbk)

    cds_added = 0
    if orfs_faa.is_file() and orfs_faa.stat().st_size > 0:
        with _open_text_read(orfs_faa) as handle:
            for orf in SeqIO.parse(handle, "fasta"):
                parsed = _parse_prodigal_header(orf.description)
                if parsed is None:
                    continue
                start0, end0, strand = parsed
                orf_id = orf.id
                contig_guess = orf_id.rsplit("_", 1)[0] if "_" in orf_id else orf_id
                rec = rec_by_id.get(contig_guess)
                if rec is None:
                    # Fallback to matching longest contig ID prefix followed by '_'
                    for cid in sorted(rec_by_id, key=len, reverse=True):
                        if orf_id.startswith(f"{cid}_"):
                            rec = rec_by_id[cid]
                            break
                if rec is None:
                    continue

                seq_len = len(rec.seq)
                clamped_start = min(max(0, start0), seq_len)
                clamped_end = min(max(clamped_start, end0), seq_len)
                if clamped_end <= clamped_start:
                    continue

                product = _match_product(products_by_contig.get(rec.id, []), clamped_start, clamped_end)
                translation = str(orf.seq).strip().rstrip("*")
                feat = SeqFeature(
                    FeatureLocation(clamped_start, clamped_end, strand=strand),
                    type="CDS",
                    qualifiers={
                        "locus_tag": [orf_id],
                        "product": [product],
                        "translation": [translation],
                    },
                )
                rec.features.append(feat)
                cds_added += 1

    out_gbk.parent.mkdir(parents=True, exist_ok=True)
    if cds_added == 0 and antismash_gbk is not None and antismash_gbk.is_file() and antismash_gbk.stat().st_size > 0:
        shutil.copyfile(antismash_gbk, out_gbk)
    else:
        for rec in records:
            rec.features.sort(key=lambda f: (int(f.location.start), int(f.location.end)))
        SeqIO.write(records, out_gbk, "genbank")
    return cds_added


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Unpack Macrel prediction/FAA outputs and build a locus-tag-matched GenBank for AMPcombi2."
    )
    parser.add_argument("--prediction", type=Path, required=True, help="Macrel *.prediction(.gz) file.")
    parser.add_argument("--orfs-faa", type=Path, required=True, help="Macrel *.all_orfs.faa(.gz) file.")
    parser.add_argument("--fasta", type=Path, required=True, help="Genome nucleotide FASTA (.fna/.fa[.gz]).")
    parser.add_argument(
        "--antismash-gbk", type=Path, default=None, help="Optional antiSMASH GenBank for product transfer."
    )
    parser.add_argument(
        "--out-prediction", type=Path, required=True, help="Output plain-text *.macrel.prediction path."
    )
    parser.add_argument("--out-faa", type=Path, required=True, help="Output plain-text *.faa path.")
    parser.add_argument("--out-gbk", type=Path, required=True, help="Output locus-tag-matched *.gbk path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    unpack_file(args.prediction, args.out_prediction)
    unpack_faa_strip_stops(args.orfs_faa, args.out_faa)
    build_macrel_genbank(
        genome_fasta=args.fasta,
        orfs_faa=args.out_faa,
        out_gbk=args.out_gbk,
        antismash_gbk=args.antismash_gbk,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
