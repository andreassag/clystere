#!/usr/bin/env python3
"""Group region GBKs by GCF for clinker synteny plotting and select cblaster queries."""

import argparse
from collections import defaultdict
import csv
from pathlib import Path
import re
import shutil
import sqlite3

MEMBERSHIP_COLS = [
    "bgc_file",
    "sample_id",
    "tool",
    "gcf_id",
    "gcf_size",
    "is_representative",
    "cutoff",
]


def infer_sample_id(filename: str) -> str:
    """Infer sample_id from a staged region GBK filename."""
    stem = Path(filename).stem
    # Match patterns like sample1.region001, sample1.deepbgc.region001, sample1.sanntis.region001, sample1.gecco.x
    for sep in (".sanntis.", ".deepbgc.", ".gecco.", ".region", "_cluster_"):
        if sep in stem:
            return stem.split(sep, 1)[0]
    match = re.match(r"^([^.]+)", stem)
    return match.group(1) if match else stem


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


def index_gbk_files(gbk_paths: list[Path]) -> tuple[dict[str, Path], dict[str, Path]]:
    """Build filename -> Path and stem -> Path lookup indices for staged GBK files."""
    by_name: dict[str, Path] = {}
    by_stem: dict[str, Path] = {}
    for gbk in sorted(gbk_paths):
        if gbk.is_dir():
            for child in _walk_gbks(gbk):
                by_name[child.name] = child
                by_stem[child.stem] = child
        elif gbk.is_file() and gbk.suffix in {".gbk", ".gbff", ".gb"}:
            by_name[gbk.name] = gbk
            by_stem[gbk.stem] = gbk
    return by_name, by_stem


def resolve_gbk(ref_name: str, by_name: dict[str, Path], by_stem: dict[str, Path]) -> Path | None:
    """Resolve a database path or BGC identifier to a staged GBK file."""
    base = Path(ref_name).name
    if base in by_name:
        return by_name[base]
    if base in by_stem:
        return by_stem[base]
    stem = Path(base).stem
    if stem in by_stem:
        return by_stem[stem]
    return None


def parse_bigscape_db(
    db_path: Path,
    by_name: dict[str, Path],
    by_stem: dict[str, Path],
    target_cutoff: float | None = None,
) -> tuple[dict[str, list[Path]], set[str], list[dict[str, str | int]]]:
    """Extract GCF groupings and membership rows from a BiG-SCAPE v2 SQLite database."""
    gcf_members: dict[str, list[Path]] = defaultdict(list)
    center_names: set[str] = set()
    rows: list[dict[str, str | int]] = []

    if not db_path.is_file() or db_path.stat().st_size == 0:
        return gcf_members, center_names, rows

    conn = sqlite3.connect(str(db_path))
    try:
        cur = conn.cursor()
        tables = {r[0] for r in cur.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if not {"gbk", "bgc_record", "family", "bgc_record_family"}.issubset(tables):
            return gcf_members, center_names, rows

        cutoffs = [r[0] for r in cur.execute("SELECT DISTINCT cutoff FROM family ORDER BY cutoff ASC")]
        if not cutoffs:
            return gcf_members, center_names, rows

        selected_cutoff = target_cutoff if target_cutoff in cutoffs else cutoffs[0]

        query = """
            SELECT f.id, f.center_id, f.cutoff, r.id, g.path
            FROM family f
            JOIN bgc_record_family brf ON f.id = brf.family_id
            JOIN bgc_record r ON brf.record_id = r.id
            JOIN gbk g ON r.gbk_id = g.id
            WHERE f.cutoff = ?
            ORDER BY f.id ASC, r.id ASC
        """
        raw_groups: dict[int, list[tuple[int, int, float, str, Path | None]]] = defaultdict(list)
        for fam_id, center_id, cutoff, rec_id, gbk_path in cur.execute(query, (selected_cutoff,)):
            matched = resolve_gbk(gbk_path, by_name, by_stem)
            raw_groups[fam_id].append((center_id, rec_id, cutoff, gbk_path, matched))

        for fam_id, entries in sorted(raw_groups.items()):
            gcf_label = f"bigscape_GCF_{fam_id:04d}"
            # Put center first so it's always kept if capped
            entries.sort(key=lambda e: (0 if e[1] == e[0] else 1, Path(e[3]).name))
            seen_files: set[str] = set()
            matched_paths: list[Path] = []
            for center_id, rec_id, cutoff, gbk_path, matched in entries:
                fname = matched.name if matched is not None else Path(gbk_path).name
                if fname in seen_files:
                    continue
                seen_files.add(fname)
                is_rep = rec_id == center_id
                if is_rep:
                    center_names.add(fname)
                if matched is not None:
                    matched_paths.append(matched)
                rows.append(
                    {
                        "bgc_file": fname,
                        "sample_id": infer_sample_id(fname),
                        "tool": "BiG-SCAPE",
                        "gcf_id": gcf_label,
                        "gcf_size": len(entries),
                        "is_representative": "true" if is_rep else "false",
                        "cutoff": str(cutoff),
                    }
                )
            if matched_paths:
                gcf_members[gcf_label] = matched_paths
    finally:
        conn.close()

    return gcf_members, center_names, rows


def parse_bigslice_tsv_dir(
    tsv_dir: Path,
    by_name: dict[str, Path],
    by_stem: dict[str, Path],
) -> tuple[dict[str, list[Path]], set[str], list[dict[str, str | int]]]:
    """Extract GCF groupings and membership rows from BiG-SLiCE tsv_export directory."""
    gcf_members: dict[str, list[Path]] = defaultdict(list)
    rep_names: set[str] = set()
    rows: list[dict[str, str | int]] = []

    if not tsv_dir.is_dir():
        return gcf_members, rep_names, rows

    membership_file = tsv_dir / "gcf_membership.tsv"
    bgc_meta_file = tsv_dir / "bgc_metadata.tsv"
    if not bgc_meta_file.is_file():
        bgc_meta_file = tsv_dir / "bgcs.tsv"

    if not membership_file.is_file() or membership_file.stat().st_size == 0:
        return gcf_members, rep_names, rows

    bgc_id_to_file: dict[str, str] = {}
    if bgc_meta_file.is_file() and bgc_meta_file.stat().st_size > 0:
        with bgc_meta_file.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = list(reader.fieldnames or [])
            index_col = fieldnames[0] if fieldnames and fieldnames[0] == "" else None
            for row in reader:
                bid = str(
                    row.get("bgc_id") or row.get("id") or (row.get(index_col) if index_col is not None else "") or ""
                ).strip()
                fname = (
                    row.get("orig_filename") or row.get("bgc") or row.get("bgc_name") or row.get("name") or ""
                ).strip()
                if bid and fname:
                    bgc_id_to_file[bid] = fname

    raw_groups: dict[str, list[tuple[float, str, Path | None]]] = defaultdict(list)
    with membership_file.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        # Real BiG-SLiCE 2.0.2 `--export-tsv` writes wide pivot format: `bgc_id\trun#1[\trun#2...]`
        run_cols = [col for col in fieldnames if re.match(r"^run#\d+$", col.strip())]
        if run_cols and "gcf_id" not in fieldnames:
            run_cols.sort(key=lambda c: int(c.split("#", 1)[1]))
            selected_run_col = run_cols[-1]
            for row in reader:
                bgc_raw = str(row.get("bgc_id") or row.get("") or "").strip()
                gcf_raw = str(row.get(selected_run_col) or "").strip()
                if not bgc_raw or not gcf_raw or gcf_raw in {"-1", "NA", "nan"}:
                    continue
                ref_name = bgc_id_to_file.get(bgc_raw, bgc_raw)
                matched = resolve_gbk(ref_name, by_name, by_stem)
                raw_groups[gcf_raw].append((0.0, ref_name, matched))
        else:
            for row in reader:
                rank = str(row.get("rank", "0")).strip()
                if rank not in {"", "0", "1"}:
                    continue
                gcf_raw = str(row.get("gcf_id", "")).strip()
                bgc_raw = str(row.get("bgc_id", "")).strip()
                if not gcf_raw or not bgc_raw or gcf_raw == "-1":
                    continue
                ref_name = bgc_id_to_file.get(bgc_raw, bgc_raw)
                try:
                    dist = float(row.get("membership_value") or 0.0)
                except ValueError:
                    dist = 0.0
                matched = resolve_gbk(ref_name, by_name, by_stem)
                raw_groups[gcf_raw].append((dist, ref_name, matched))

    for gcf_raw, entries in sorted(raw_groups.items(), key=lambda kv: kv[0]):
        gcf_label = f"bigslice_GCF_{int(gcf_raw):04d}" if gcf_raw.isdigit() else f"bigslice_{gcf_raw}"
        # Closest to GCF centroid first
        entries.sort(key=lambda e: (e[0], Path(e[1]).name))
        seen_files: set[str] = set()
        matched_paths: list[Path] = []
        for idx, (_dist, ref_name, matched) in enumerate(entries):
            fname = matched.name if matched is not None else Path(ref_name).name
            if fname in seen_files:
                continue
            seen_files.add(fname)
            is_rep = idx == 0
            if is_rep:
                rep_names.add(fname)
            if matched is not None:
                matched_paths.append(matched)
            rows.append(
                {
                    "bgc_file": fname,
                    "sample_id": infer_sample_id(fname),
                    "tool": "BiG-SLiCE",
                    "gcf_id": gcf_label,
                    "gcf_size": len(entries),
                    "is_representative": "true" if is_rep else "false",
                    "cutoff": "NA",
                }
            )
        if matched_paths:
            gcf_members[gcf_label] = matched_paths

    return gcf_members, rep_names, rows


def group_gcfs(
    gbk_paths: list[Path],
    outdir: Path,
    bigscape_db: Path | None = None,
    bigslice_tsv_dir: Path | None = None,
    clinker_max_bgcs: int = 50,
    cblaster_max_queries: int = 50,
    bigscape_cutoff: float | None = None,
) -> tuple[int, int]:
    """Organize multi-member GCF folders for clinker and representative queries for cblaster."""
    by_name, by_stem = index_gbk_files(gbk_paths)
    gcf_dir = outdir / "gcf_groups"
    cblaster_dir = outdir / "cblaster_queries"
    gcf_dir.mkdir(parents=True, exist_ok=True)
    cblaster_dir.mkdir(parents=True, exist_ok=True)

    all_groups: dict[str, list[Path]] = {}
    all_reps: set[str] = set()
    membership_rows: list[dict[str, str | int]] = []

    if bigscape_db is not None:
        bs_groups, bs_reps, bs_rows = parse_bigscape_db(bigscape_db, by_name, by_stem, target_cutoff=bigscape_cutoff)
        all_groups.update(bs_groups)
        all_reps.update(bs_reps)
        membership_rows.extend(bs_rows)

    if bigslice_tsv_dir is not None:
        sl_groups, sl_reps, sl_rows = parse_bigslice_tsv_dir(bigslice_tsv_dir, by_name, by_stem)
        all_groups.update(sl_groups)
        all_reps.update(sl_reps)
        membership_rows.extend(sl_rows)

    # Emit multi-member GCF directories for clinker
    multi_gcf_count = 0
    for gcf_label, members in sorted(all_groups.items()):
        non_empty = [p for p in members if p.stat().st_size > 0]
        if len(non_empty) < 2:
            continue
        multi_gcf_count += 1
        target_dir = gcf_dir / gcf_label
        target_dir.mkdir(parents=True, exist_ok=True)
        for member in non_empty[:clinker_max_bgcs]:
            shutil.copy2(member, target_dir / member.name)

    # Select cblaster query GBKs: GCF representatives + unclustered singletons (or all GBKs if unclustered)
    clustered_files = {row["bgc_file"] for row in membership_rows}
    query_candidates: list[Path] = []
    for fname, fpath in sorted(by_name.items()):
        if fpath.stat().st_size == 0:
            continue
        if not membership_rows or fname in all_reps or fname not in clustered_files:
            query_candidates.append(fpath)

    for qpath in query_candidates[:cblaster_max_queries]:
        shutil.copy2(qpath, cblaster_dir / qpath.name)

    summary_path = outdir / "gcf_membership_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MEMBERSHIP_COLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(membership_rows)

    return multi_gcf_count, min(len(query_candidates), cblaster_max_queries)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Group BGC region GBKs by GCF for clinker and cblaster.")
    parser.add_argument("--gbks", type=Path, nargs="*", default=[], help="Region GBK files or directories.")
    parser.add_argument("--bigscape-db", type=Path, default=None, help="Optional BiG-SCAPE v2 bigscape.db path.")
    parser.add_argument("--bigslice-tsv", type=Path, default=None, help="Optional BiG-SLiCE tsv_export directory.")
    parser.add_argument("--bigscape-cutoff", type=float, default=None, help="Specific BiG-SCAPE cutoff to select.")
    parser.add_argument("--clinker-max-bgcs", type=int, default=50, help="Max BGCs per clinker GCF directory.")
    parser.add_argument("--cblaster-max-queries", type=int, default=50, help="Max query GBKs for cblaster.")
    parser.add_argument("-o", "--outdir", type=Path, default=Path(), help="Output directory.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    group_gcfs(
        gbk_paths=args.gbks,
        outdir=args.outdir,
        bigscape_db=args.bigscape_db,
        bigslice_tsv_dir=args.bigslice_tsv,
        clinker_max_bgcs=args.clinker_max_bgcs,
        cblaster_max_queries=args.cblaster_max_queries,
        bigscape_cutoff=args.bigscape_cutoff,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
