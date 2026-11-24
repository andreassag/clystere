#!/usr/bin/env python3
"""Convert CheckM2 quality_report.tsv to dRep --genomeInfo CSV format."""

import argparse
import csv
from pathlib import Path

from select_derep_genomes import load_checkm2_quality, stem_sample_id


def convert_checkm2_to_drep_csv(
    checkm2_tsv: Path,
    genome_paths: list[Path],
    output_csv: Path,
) -> int:
    """Write genome,completeness,contamination CSV matching staged dRep FASTA filenames.

    Args:
        checkm2_tsv: Path to CheckM2 ``quality_report.tsv``.
        genome_paths: Genome FASTA paths passed to dRep.
        output_csv: Destination CSV path for ``dRep --genomeInfo``.

    Returns:
        Number of genome rows written.
    """
    quality = load_checkm2_quality(checkm2_tsv)
    rows: list[dict[str, str]] = []

    if genome_paths:
        missing: list[str] = []
        for gpath in sorted(genome_paths, key=lambda p: p.name):
            sid = stem_sample_id(gpath.name)
            if sid not in quality:
                if quality:
                    missing.append(gpath.name)
                continue
            comp, cont = quality[sid]
            rows.append(
                {
                    "genome": gpath.name,
                    "completeness": f"{comp:.2f}",
                    "contamination": f"{cont:.2f}",
                }
            )
        if missing:
            raise ValueError(f"CheckM2 report '{checkm2_tsv}' is missing quality entries for: {', '.join(missing)}")
    else:
        for sid, (comp, cont) in sorted(quality.items()):
            rows.append(
                {
                    "genome": f"{sid}.fna",
                    "completeness": f"{comp:.2f}",
                    "contamination": f"{cont:.2f}",
                }
            )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["genome", "completeness", "contamination"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


def build_parser() -> argparse.ArgumentParser:
    """Create CLI argument parser."""
    parser = argparse.ArgumentParser(description="Convert CheckM2 report TSV to dRep genomeInfo CSV.")
    parser.add_argument("--checkm2", type=Path, required=True, help="Path to CheckM2 quality_report.tsv.")
    parser.add_argument("--genomes", type=Path, nargs="*", default=[], help="Staged genome FASTA paths for dRep.")
    parser.add_argument("--out", type=Path, default=Path("drep_genome_info.csv"), help="Output dRep genomeInfo CSV.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entrypoint for converting CheckM2 TSV to dRep genomeInfo CSV."""
    args = build_parser().parse_args(argv)
    convert_checkm2_to_drep_csv(args.checkm2, args.genomes, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
