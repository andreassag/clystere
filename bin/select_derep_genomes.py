#!/usr/bin/env python3
"""Select non-redundant representative genomes from skani ANI comparisons."""

import argparse
import csv
import gzip
import math
from pathlib import Path

FASTA_SUFFIXES = (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa")


def stem_sample_id(path_or_name: str | Path) -> str:
    """Normalize a genome filename to its sample identifier."""
    name = Path(str(path_or_name)).name
    lower = name.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def compute_fasta_stats(fasta_path: Path) -> tuple[int, int, int]:
    """Compute (total_length, n50, num_contigs) for a plain or gzipped FASTA."""
    lengths: list[int] = []
    current = 0
    opener = gzip.open if fasta_path.suffix == ".gz" else Path.open
    with opener(fasta_path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                    current = 0
            else:
                current += len(line)
        if current > 0:
            lengths.append(current)

    if not lengths:
        return 0, 0, 0

    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2.0
    running = 0
    n50 = lengths[-1]
    for length in lengths:
        running += length
        if running >= half:
            n50 = length
            break
    return total, n50, len(lengths)


def load_checkm2_quality(checkm2_path: Path | None) -> dict[str, tuple[float, float]]:
    """Load (completeness, contamination) keyed by normalized sample ID."""
    quality: dict[str, tuple[float, float]] = {}
    if checkm2_path is None or not checkm2_path.is_file() or checkm2_path.stat().st_size == 0:
        return quality
    with checkm2_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            name = row.get("Name") or row.get("Bin Id") or ""
            if not name:
                continue
            try:
                comp = float(row["Completeness"])
                cont = float(row["Contamination"])
            except (KeyError, TypeError, ValueError):
                continue
            quality[stem_sample_id(name)] = (comp, cont)
    return quality


def parse_skani_ani(skani_path: Path) -> list[tuple[str, str, float]]:
    """Parse pairwise (sample_a, sample_b, ani) edges from skani triangle output.

    Supports both `--sparse` edge-list TSV and dense lower-triangular PHYLIP-like output.
    """
    if not skani_path.is_file() or skani_path.stat().st_size == 0:
        return []

    lines = [line.rstrip("\n") for line in skani_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        return []

    edges: list[tuple[str, str, float]] = []
    header_parts = lines[0].split("\t")

    # Sparse edge-list format: Ref_file\tQuery_file\tANI\t...
    if len(header_parts) >= 3 and "Ref_file" in header_parts[0] and "ANI" in header_parts[2]:
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                ani = float(parts[2])
            except ValueError:
                continue
            edges.append((stem_sample_id(parts[0]), stem_sample_id(parts[1]), ani))
        return edges

    # Dense lower-triangular matrix format: first line is N, then "<file>\t<ani1>\t<ani2>..."
    start_idx = 1 if lines[0].strip().isdigit() else 0
    names: list[str] = []
    for row_idx, line in enumerate(lines[start_idx:]):
        parts = line.split("\t")
        sample = stem_sample_id(parts[0])
        names.append(sample)
        for col_idx, val in enumerate(parts[1:]):
            if col_idx >= row_idx:
                break
            try:
                ani = float(val)
            except ValueError:
                continue
            edges.append((names[col_idx], sample, ani))
    return edges


def cluster_and_select(
    genome_paths: list[Path],
    edges: list[tuple[str, str, float]],
    ani_threshold: float = 99.0,
    checkm2_quality: dict[str, tuple[float, float]] | None = None,
) -> tuple[list[str], list[dict[str, str | int | float]]]:
    """Build ANI connected components and select the highest-scoring representative per cluster."""
    quality = checkm2_quality or {}
    stats: dict[str, tuple[int, int, int]] = {}
    for gpath in sorted(genome_paths):
        sid = stem_sample_id(gpath)
        stats[sid] = compute_fasta_stats(gpath)

    # Ensure any samples present in edges are registered even if FASTA list was empty (e.g., stub)
    for sample_a, sample_b, _ani in edges:
        stats.setdefault(sample_a, (0, 0, 0))
        stats.setdefault(sample_b, (0, 0, 0))

    adj: dict[str, dict[str, float]] = {sid: {} for sid in stats}
    for sample_a, sample_b, ani in edges:
        if sample_a == sample_b:
            continue
        if ani >= ani_threshold:
            prev = adj[sample_a].get(sample_b, 0.0)
            best = max(prev, ani)
            adj[sample_a][sample_b] = best
            adj[sample_b][sample_a] = best

    visited: set[str] = set()
    components: list[list[str]] = []
    for sid in sorted(stats):
        if sid in visited:
            continue
        stack = [sid]
        comp: list[str] = []
        visited.add(sid)
        while stack:
            node = stack.pop()
            comp.append(node)
            for nxt in sorted(adj[node]):
                if nxt not in visited:
                    visited.add(nxt)
                    stack.append(nxt)
        components.append(sorted(comp))

    def score_genome(sid: str) -> tuple[float, int, str]:
        total_len, n50, _contigs = stats[sid]
        n50_term = 0.5 * math.log10(max(n50, 1))
        if sid in quality:
            completeness, contamination = quality[sid]
            score = completeness - 5.0 * contamination + n50_term
        else:
            score = n50_term + 0.1 * math.log10(max(total_len, 1))
        return score, total_len, sid

    representatives: list[str] = []
    rows: list[dict[str, str | int | float]] = []

    for cluster_idx, members in enumerate(components, start=1):
        cluster_id = f"DEREP_{cluster_idx:04d}"
        scored = [(score_genome(m), m) for m in members]
        # Highest score, then highest total length, then lexicographically smallest sample_id
        scored.sort(key=lambda item: (-item[0][0], -item[0][1], item[0][2]))
        rep_id = scored[0][1]
        representatives.append(rep_id)

        for (score, total_len, _), sid in scored:
            _, n50, _ = stats[sid]
            comp_val = f"{quality[sid][0]:.2f}" if sid in quality else "NA"
            cont_val = f"{quality[sid][1]:.2f}" if sid in quality else "NA"
            ani_to_rep = 100.0 if sid == rep_id else adj[sid].get(rep_id, ani_threshold)
            rows.append(
                {
                    "cluster_id": cluster_id,
                    "sample_id": sid,
                    "is_representative": "true" if sid == rep_id else "false",
                    "score": round(score, 4),
                    "completeness": comp_val,
                    "contamination": cont_val,
                    "n50": n50,
                    "total_length": total_len,
                    "ani_to_rep": round(ani_to_rep, 3),
                }
            )

    return sorted(representatives), rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Select representative dereplicated genomes from skani ANI output.")
    parser.add_argument("--skani", type=Path, required=True, help="Path to skani triangle TSV output.")
    parser.add_argument("--genomes", type=Path, nargs="*", default=[], help="All input genome FASTA paths.")
    parser.add_argument("--ani", type=float, default=99.0, help="ANI clustering threshold (default: 99.0).")
    parser.add_argument("--checkm2", type=Path, default=None, help="Optional CheckM2 quality_report.tsv.")
    parser.add_argument(
        "--out-ids",
        type=Path,
        default=Path("dereplicated_genome_ids.txt"),
        help="Output file listing winning representative sample IDs.",
    )
    parser.add_argument(
        "--out-clusters",
        type=Path,
        default=Path("dereplication_clusters.tsv"),
        help="Output TSV describing dereplication clusters and scores.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    edges = parse_skani_ani(args.skani)
    quality = load_checkm2_quality(args.checkm2)
    reps, rows = cluster_and_select(args.genomes, edges, ani_threshold=args.ani, checkm2_quality=quality)

    args.out_ids.parent.mkdir(parents=True, exist_ok=True)
    args.out_ids.write_text("\n".join(reps) + ("\n" if reps else ""), encoding="utf-8")

    args.out_clusters.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "cluster_id",
        "sample_id",
        "is_representative",
        "score",
        "completeness",
        "contamination",
        "n50",
        "total_length",
        "ani_to_rep",
    ]
    with args.out_clusters.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
