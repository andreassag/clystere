#!/usr/bin/env python3
"""Extract RiPP precursor candidates and RRE-Finder hits from antiSMASH 8 JSON outputs."""

import argparse
import csv
from pathlib import Path
import re

from antismash_io import list_antismash_json_files, load_antismash_json, run_threaded

OUTPUT_COLS = [
    "sample_id",
    "record_id",
    "region",
    "ripp_subclass",
    "locus_tag",
    "rre_bitscore",
    "leader_sequence",
    "core_sequence",
    "predicted_molecular_weight",
    "tailoring_modifications",
]

# Average residue masses (Da) + H2O for molecular weight estimation when absent from JSON
AA_WEIGHTS = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}
WATER_MASS = 18.01528

RIPP_MODULES = (
    ("antismash.modules.lanthipeptides", "lanthipeptide"),
    ("antismash.modules.lassopeptides", "lassopeptide"),
    ("antismash.modules.sactipeptides", "sactipeptide"),
    ("antismash.modules.thiopeptides", "thiopeptide"),
)


def parse_location_bounds(location_str: str) -> tuple[int, int] | None:
    """Extract (start, end) integer bounds from a BioPython-style location string."""
    coords = re.findall(r"\d+", str(location_str))
    if len(coords) < 2:
        return None
    return int(coords[0]), int(coords[-1])


def estimate_molecular_weight(core_seq: str) -> str:
    """Estimate average molecular weight (Da) of a core peptide sequence."""
    seq = re.sub(r"[^A-Za-z]", "", core_seq).upper()
    if not seq:
        return ""
    total = WATER_MASS + sum(AA_WEIGHTS.get(aa, 110.0) for aa in seq)
    return f"{total:.2f}"


def _find_enclosing_region(
    bounds: tuple[int, int] | None,
    regions: list[tuple[int, int, str]],
) -> str:
    if bounds is None:
        return ""
    start, end = bounds
    for reg_start, reg_end, reg_num in regions:
        if start <= reg_end and end >= reg_start:
            return reg_num
    return ""


def _extract_modifications(motif: dict) -> str:
    mods: list[str] = []
    detailed = motif.get("detailed_info") or {}
    if isinstance(detailed, dict):
        for key in ("predicted_class", "class", "bridges", "alternative_weights", "monoisotopic_mass", "tailoring"):
            val = detailed.get(key)
            if val:
                joined = ",".join(str(x) for x in val) if isinstance(val, list) else str(val)
                mods.append(f"{key}={joined}")
    tail = motif.get("tail")
    if tail:
        mods.append(f"tail={tail}")
    return ";".join(mods)


def _is_positive_weight(val: object) -> bool:
    if val in (None, ""):
        return False
    try:
        return float(str(val)) > 0.0
    except ValueError:
        return True


def _extract_weight(motif: dict, core_seq: str) -> str:
    for key in ("molecular_weight", "monoisotopic_mass", "predicted_molecular_weight"):
        val = motif.get(key)
        if _is_positive_weight(val):
            return str(val)
    detailed = motif.get("detailed_info") or {}
    if isinstance(detailed, dict):
        for key in ("molecular_weight", "monoisotopic_mass"):
            val = detailed.get(key)
            if isinstance(val, list) and val and _is_positive_weight(val[0]):
                return str(val[0])
            if not isinstance(val, list) and _is_positive_weight(val):
                return str(val)
    return estimate_molecular_weight(core_seq)


def parse_antismash_ripp_json(json_path: Path) -> list[dict[str, str]]:
    """Parse one antiSMASH 8 JSON file and return RiPP/RRE candidate rows."""
    data = load_antismash_json(json_path)
    sample_id = json_path.stem
    rows: list[dict[str, str]] = []

    for record in data.get("records", []):
        record_id = record.get("name") or record.get("id") or ""
        features = record.get("features", [])
        modules = record.get("modules") or {}

        # Map region intervals and CDS locus_tag -> coordinates
        regions: list[tuple[int, int, str]] = []
        cds_bounds: dict[str, tuple[int, int]] = {}
        for feat in features:
            ftype = feat.get("type")
            bounds = parse_location_bounds(feat.get("location", ""))
            quals = feat.get("qualifiers") or {}
            if ftype == "region" and bounds is not None:
                reg_num = str(quals.get("region_number", [""])[0])
                regions.append((bounds[0], bounds[1], reg_num))
            elif ftype == "CDS" and bounds is not None:
                locus = quals.get("locus_tag", [""])[0]
                if locus:
                    cds_bounds[locus] = bounds

        # Extract RRE-Finder bitscores per locus_tag
        rre_scores: dict[str, tuple[float, str]] = {}
        rre_mod = modules.get("antismash.modules.rrefinder") or {}
        hits_by_cds = rre_mod.get("hits_by_cds") or {}
        for locus_tag, hits in hits_by_cds.items():
            best_score = -1.0
            best_domain = "RRE"
            for hit in hits or []:
                score = float(hit.get("bitscore", hit.get("score", 0.0)))
                if score > best_score:
                    best_score = score
                    best_domain = str(hit.get("domain", "RRE"))
            if best_score >= 0.0:
                rre_scores[locus_tag] = (best_score, best_domain)

        seen_loci: set[tuple[str, str, str]] = set()

        # Extract RiPP class-specific motifs from antiSMASH 8 modules
        for mod_key, default_subclass in RIPP_MODULES:
            mod_data = modules.get(mod_key) or {}
            motifs_field = mod_data.get("motifs")
            motif_items: list[tuple[str, dict]] = []
            if isinstance(motifs_field, dict):
                for locus, motif_list in motifs_field.items():
                    for motif in motif_list or []:
                        motif_items.append((locus, motif))
            elif isinstance(motifs_field, list):
                for motif in motifs_field:
                    motif_items.append((str(motif.get("locus_tag", "")), motif))

            for locus_key, motif in motif_items:
                locus = str(motif.get("locus_tag") or locus_key)
                bounds = parse_location_bounds(motif.get("location", "")) or cds_bounds.get(locus)
                region_num = _find_enclosing_region(bounds, regions)
                leader = str(motif.get("leader") or "")
                core = str(motif.get("core") or "")
                rre_val = f"{rre_scores[locus][0]:.1f}" if locus in rre_scores else ""
                subclass = str(motif.get("peptide_subclass") or default_subclass)
                sig = (locus, subclass, core)
                if sig in seen_loci:
                    continue
                seen_loci.add(sig)
                rows.append(
                    {
                        "sample_id": sample_id,
                        "record_id": record_id,
                        "region": region_num,
                        "ripp_subclass": subclass,
                        "locus_tag": locus,
                        "rre_bitscore": rre_val,
                        "leader_sequence": leader,
                        "core_sequence": core,
                        "predicted_molecular_weight": _extract_weight(motif, core),
                        "tailoring_modifications": _extract_modifications(motif),
                    }
                )

        # Fallback: capture RiPP CDS_motif features carrying core_sequence / leader_sequence
        for feat in features:
            if feat.get("type") != "CDS_motif":
                continue
            quals = feat.get("qualifiers") or {}
            core = str(quals.get("core_sequence", [""])[0])
            leader = str(quals.get("leader_sequence", [""])[0])
            if not core and not leader:
                continue
            locus = str(quals.get("locus_tag", [""])[0])
            subclass = str(quals.get("prepeptide", quals.get("tool", ["RiPP"]))[0])
            sig = (locus, subclass, core)
            if sig in seen_loci:
                continue
            seen_loci.add(sig)
            bounds = parse_location_bounds(feat.get("location", "")) or cds_bounds.get(locus)
            region_num = _find_enclosing_region(bounds, regions)
            rre_val = f"{rre_scores[locus][0]:.1f}" if locus in rre_scores else ""
            mw = str(quals.get("molecular_weight", quals.get("monoisotopic_mass", [""]))[0])
            if not mw:
                mw = estimate_molecular_weight(core)
            rows.append(
                {
                    "sample_id": sample_id,
                    "record_id": record_id,
                    "region": region_num,
                    "ripp_subclass": subclass,
                    "locus_tag": locus,
                    "rre_bitscore": rre_val,
                    "leader_sequence": leader,
                    "core_sequence": core,
                    "predicted_molecular_weight": mw,
                    "tailoring_modifications": ";".join(quals.get("note", [])),
                }
            )

        # Emit standalone RRE-Finder hits that were not attached to a specific prepeptide motif
        reported_loci = {row["locus_tag"] for row in rows if row["record_id"] == record_id}
        for locus, (bitscore, domain) in sorted(rre_scores.items()):
            if locus in reported_loci:
                continue
            bounds = cds_bounds.get(locus)
            region_num = _find_enclosing_region(bounds, regions)
            rows.append(
                {
                    "sample_id": sample_id,
                    "record_id": record_id,
                    "region": region_num,
                    "ripp_subclass": f"RRE:{domain}",
                    "locus_tag": locus,
                    "rre_bitscore": f"{bitscore:.1f}",
                    "leader_sequence": "",
                    "core_sequence": "",
                    "predicted_molecular_weight": "",
                    "tailoring_modifications": "",
                }
            )

    return rows


def main(directory: Path, output: Path, threads: int | None = None) -> int:
    json_files = list_antismash_json_files(directory)
    all_rows: list[dict[str, str]] = []
    for parsed in run_threaded(json_files, parse_antismash_ripp_json, threads):
        all_rows.extend(parsed)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract RiPP and RRE-Finder candidates from antiSMASH 8 JSON outputs."
    )
    parser.add_argument("directory", type=Path, help="Directory containing antiSMASH JSON outputs.")
    parser.add_argument("output", type=Path, help="Output TSV path (ripp_candidates.tsv).")
    parser.add_argument("--threads", type=int, default=None, help="Number of worker threads.")
    return parser


if __name__ == "__main__":
    parsed_args = build_parser().parse_args()
    raise SystemExit(main(parsed_args.directory, parsed_args.output, parsed_args.threads))
