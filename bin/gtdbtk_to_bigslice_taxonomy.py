#!/usr/bin/env python3
"""Convert GTDB-Tk summary TSV(s) into BiG-SLiCE taxonomy.tsv format."""

import argparse
import csv
from pathlib import Path

FASTA_SUFFIXES = (".fasta.gz", ".fna.gz", ".fa.gz", ".gbff.gz", ".gbk.gz", ".fasta", ".fna", ".fa", ".gbff", ".gbk")
RANK_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")


def normalize_genome_id(raw_id: str) -> str:
    """Strip common genome file extensions from a GTDB-Tk user_genome identifier."""
    name = Path(raw_id.strip()).name
    lower = name.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def parse_gtdb_classification(classification: str) -> list[str]:
    """Parse a 7-rank GTDB lineage string into [Domain, Phylum, Class, Order, Family, Genus, Species]."""
    ranks = [""] * 7
    if not classification:
        return ranks
    parts = [part.strip() for part in classification.split(";")]
    for idx, prefix in enumerate(RANK_PREFIXES):
        if idx < len(parts):
            token = parts[idx]
            if token.startswith(prefix):
                token = token[len(prefix) :]
            ranks[idx] = token.strip()
    return ranks


def convert_gtdbtk_summaries(
    summary_paths: list[Path],
    output_path: Path,
    dataset_folder: str = "",
) -> int:
    """Write a BiG-SLiCE 9-column taxonomy.tsv file from GTDB-Tk summary TSV(s)."""
    seen: dict[str, list[str]] = {}
    for summary_path in sorted(summary_paths):
        if not summary_path.is_file() or summary_path.stat().st_size == 0:
            continue
        with summary_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                user_genome = row.get("user_genome", "").strip()
                if not user_genome:
                    continue
                sample_id = normalize_genome_id(user_genome)
                ranks = parse_gtdb_classification(row.get("classification", ""))
                organism = ranks[6] or sample_id
                prefix = f"{dataset_folder}/{sample_id}." if dataset_folder else f"{sample_id}/"
                seen[sample_id] = [prefix, *ranks, organism]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# path_startswith\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOrganism\n")
        for sample_id in sorted(seen):
            handle.write("\t".join(seen[sample_id]) + "\n")
    return len(seen)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert GTDB-Tk summary TSV(s) to BiG-SLiCE taxonomy.tsv.")
    parser.add_argument("summaries", type=Path, nargs="*", help="GTDB-Tk *.summary.tsv file(s).")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output BiG-SLiCE taxonomy.tsv path.")
    parser.add_argument(
        "--dataset-folder",
        type=str,
        default="",
        help="Optional shared BiG-SLiCE dataset subfolder name (default: per-sample <sample_id>/ folders).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    convert_gtdbtk_summaries(args.summaries, args.output, dataset_folder=args.dataset_folder)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
