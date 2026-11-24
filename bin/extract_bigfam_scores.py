#!/usr/bin/env python3
"""Extract BiG-FAM / BiG-SLiCE query novelty scores from a BiG-SLiCE report SQLite DB."""

import argparse
import csv
from pathlib import Path
import sqlite3

OUTPUT_COLS = [
    "bgc_id",
    "bgc_filename",
    "gcf_id",
    "distance_d",
    "rank",
    "membership_type",
    "is_putative_novel",
]


def classify_membership(distance_d: float, novelty_threshold: float) -> tuple[str, bool]:
    """Classify a BiG-SLiCE distance score into membership type and putative novelty flag."""
    if distance_d <= novelty_threshold * 0.5:
        return "core", False
    if distance_d <= novelty_threshold:
        return "putative", False
    return "orphan_novel", True


def _walk_dbs(directory: Path) -> list[Path]:
    found: list[Path] = []
    for entry in sorted(directory.iterdir()):
        if entry.name.startswith("."):
            continue
        if entry.is_dir():
            found.extend(_walk_dbs(entry))
        elif entry.is_file() and entry.name == "data.db" and entry.stat().st_size > 0:
            found.append(entry)
    return found


def find_report_dbs(reports_dir: Path) -> list[Path]:
    """Locate BiG-SLiCE query report SQLite databases under a reports directory."""
    if reports_dir.is_file() and reports_dir.suffix == ".db":
        return [reports_dir]
    if not reports_dir.is_dir():
        return []
    return _walk_dbs(reports_dir)


def extract_scores_from_db(db_path: Path, novelty_threshold: float = 900.0) -> list[dict[str, str | int | float]]:
    """Query a BiG-SLiCE query report data.db and return novelty score rows."""
    rows: list[dict[str, str | int | float]] = []
    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        tables = {row[0] for row in cur.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "bgc" not in tables or "gcf_membership" not in tables:
            return rows

        query = """
            SELECT b.name, b.orig_filename, gm.gcf_id, gm.membership_value, gm.rank
            FROM bgc b
            LEFT JOIN gcf_membership gm ON b.id = gm.bgc_id
            ORDER BY b.name ASC, gm.rank ASC
        """
        for bgc_name, orig_filename, gcf_id, membership_value, rank in cur.execute(query):
            raw_name = Path(str(bgc_name)).name if "/" in str(bgc_name) else str(bgc_name)
            clean_id = raw_name[:-4] if raw_name.lower().endswith(".gbk") else raw_name
            clean_fname = Path(str(orig_filename or bgc_name)).name
            if gcf_id is None or membership_value is None:
                rows.append(
                    {
                        "bgc_id": clean_id,
                        "bgc_filename": clean_fname,
                        "gcf_id": "NA",
                        "distance_d": "NA",
                        "rank": "NA",
                        "membership_type": "unassigned",
                        "is_putative_novel": "true",
                    }
                )
            else:
                dist = float(membership_value)
                mtype, is_novel = classify_membership(dist, novelty_threshold)
                rows.append(
                    {
                        "bgc_id": clean_id,
                        "bgc_filename": clean_fname,
                        "gcf_id": f"GCF_{int(gcf_id)}",
                        "distance_d": round(dist, 4),
                        "rank": int(rank) if rank is not None else 0,
                        "membership_type": mtype,
                        "is_putative_novel": "true" if is_novel else "false",
                    }
                )
    finally:
        conn.close()
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract BiG-FAM query novelty scores to TSV.")
    parser.add_argument("reports_path", type=Path, help="Path to BiG-SLiCE reports directory or report data.db.")
    parser.add_argument("output", type=Path, help="Output TSV path (bigfam_novelty_scores.tsv).")
    parser.add_argument(
        "--novelty-threshold",
        type=float,
        default=900.0,
        help="Distance threshold above which a BGC is flagged as putative novel (default: 900.0).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    all_rows: list[dict[str, str | int | float]] = []
    for db_path in find_report_dbs(args.reports_path):
        all_rows.extend(extract_scores_from_db(db_path, novelty_threshold=args.novelty_threshold))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
