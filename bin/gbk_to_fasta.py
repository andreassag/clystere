#!/usr/bin/env python3
"""Convert GenBank files to FASTA (and optional FAA/GFF3) or normalize FASTA inputs."""

import argparse
from contextlib import ExitStack
import gzip
from pathlib import Path
from typing import TextIO

from Bio import SeqIO


def _open_text_read(path: Path) -> TextIO:
    """Open a plain or gzip-compressed text file for reading."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _open_text_write(path: Path) -> TextIO:
    """Open a plain or gzip-compressed text file for writing."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("w", encoding="utf-8")


def convert_genbank(
    input_path: Path,
    fasta_out: Path,
    faa_out: Path | None = None,
    gff_out: Path | None = None,
) -> tuple[int, int]:
    """Convert a GenBank record file into FASTA and optional FAA/GFF3 outputs.

    Returns:
        Tuple of (nucleotide_records_written, cds_proteins_written).
    """
    n_records = 0
    n_cds = 0

    with ExitStack() as stack:
        in_handle = stack.enter_context(_open_text_read(input_path))
        fa_handle = stack.enter_context(_open_text_write(fasta_out))
        faa_handle = stack.enter_context(_open_text_write(faa_out)) if faa_out else None
        gff_handle = stack.enter_context(_open_text_write(gff_out)) if gff_out else None

        if gff_handle is not None:
            gff_handle.write("##gff-version 3\n")

        for record in SeqIO.parse(in_handle, "genbank"):
            try:
                seq_str = str(record.seq)
            except Exception:
                continue
            if not seq_str:
                continue

            SeqIO.write(record, fa_handle, "fasta")
            n_records += 1

            if gff_handle is not None:
                gff_handle.write(f"##sequence-region {record.id} 1 {len(seq_str)}\n")

            cds_index = 0
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                cds_index += 1
                locus_tag = (
                    feature.qualifiers.get("locus_tag", [None])[0]
                    or feature.qualifiers.get("protein_id", [None])[0]
                    or f"{record.id}_CDS_{cds_index}"
                )
                product = feature.qualifiers.get("product", ["hypothetical protein"])[0]

                if faa_handle is not None:
                    translations = feature.qualifiers.get("translation", [])
                    protein_seq = translations[0].strip() if translations else ""
                    if not protein_seq:
                        try:
                            protein_seq = str(feature.extract(record.seq).translate(to_stop=True))
                        except Exception:
                            protein_seq = ""
                    if protein_seq:
                        faa_handle.write(f">{locus_tag} {product}\n{protein_seq}\n")
                        n_cds += 1

                if gff_handle is not None:
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "-" if feature.location.strand == -1 else "+"
                    safe_product = product.replace(";", "%3B").replace("=", "%3D")
                    attrs = f"ID={locus_tag};locus_tag={locus_tag};product={safe_product}"
                    gff_handle.write(f"{record.id}\tgbk_to_fasta\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attrs}\n")

    return n_records, n_cds


def normalize_fasta(input_path: Path, fasta_out: Path) -> int:
    """Copy/decompress an input FASTA file into a clean uncompressed (or .gz) FASTA."""
    n_records = 0
    with _open_text_read(input_path) as in_handle, _open_text_write(fasta_out) as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            if len(record.seq) == 0:
                continue
            SeqIO.write(record, out_handle, "fasta")
            n_records += 1
    return n_records


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert GenBank to FASTA (+ optional FAA/GFF3) or normalize FASTA.")
    parser.add_argument(
        "input", type=Path, help="Input GenBank (.gbk/.gbff[.gz]) or FASTA (.fa/.fna/.fasta[.gz]) file."
    )
    parser.add_argument("output", type=Path, help="Output nucleotide FASTA file (.fna or .fna.gz).")
    parser.add_argument("--faa", type=Path, default=None, help="Optional output protein FASTA path.")
    parser.add_argument("--gff", type=Path, default=None, help="Optional output GFF3 path.")
    parser.add_argument(
        "--normalize-fasta",
        action="store_true",
        help="Treat input as FASTA and write normalized FASTA output.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.normalize_fasta:
        normalize_fasta(args.input, args.output)
    else:
        convert_genbank(args.input, args.output, faa_out=args.faa, gff_out=args.gff)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
