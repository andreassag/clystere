#!/usr/bin/env python3
"""Merge multi-predictor BGC calls into consensus loci and join QC, taxonomy, and novelty metrics."""

import argparse
from collections import defaultdict
import contextlib
import csv
from pathlib import Path
import re

FASTA_SUFFIXES = (".fasta.gz", ".fna.gz", ".fa.gz", ".gbff.gz", ".gbk.gz", ".fasta", ".fna", ".fa", ".gbff", ".gbk")

CONSENSUS_COLS = [
    "locus_id",
    "sample_id",
    "contig_id",
    "start",
    "end",
    "length",
    "tools",
    "tool_count",
    "support_tier",
    "mean_jaccard",
    "product_classes",
    "complete_status",
    "mibig_ids",
    "knownclusterblast_hit",
    "knownclusterblast_similarity",
    "bigscape_gcf",
    "bigslice_gcf",
    "bigfam_gcf",
    "bigfam_distance_d",
    "is_putative_novel",
]

GENOME_SUMMARY_COLS = [
    "sample_id",
    "gtdb_classification",
    "completeness",
    "contamination",
    "assembly_length_bp",
    "num_contigs",
    "n50",
    "antismash_bgcs",
    "gecco_bgcs",
    "deepbgc_bgcs",
    "sanntis_bgcs",
    "consensus_loci",
    "high_support_loci",
    "medium_support_loci",
    "single_tool_loci",
    "putative_novel_bgcs",
    "bgcs_per_mbp",
    "complete_bgc_fraction",
    "contig_edge_rate",
    "amp_count",
]


def normalize_sample_id(raw_id: str) -> str:
    """Normalize genome/sample names across QUAST, CheckM2, GTDB-Tk, and comBGC."""
    name = Path(raw_id.strip()).name
    lower = name.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def interval_jaccard(start_a: int, end_a: int, start_b: int, end_b: int) -> float:
    """Compute 1-based inclusive coordinate Jaccard similarity between two intervals."""
    inter = max(0, min(end_a, end_b) - max(start_a, start_b) + 1)
    if inter == 0:
        return 0.0
    union = max(end_a, end_b) - min(start_a, start_b) + 1
    return inter / union if union > 0 else 0.0


def support_tier(tool_count: int) -> str:
    """Map number of distinct agreeing predictors to a human-readable support tier."""
    if tool_count >= 3:
        return "High (3-4 tools)"
    if tool_count == 2:
        return "Medium (2 tools)"
    return "Single-tool (1 tool)"


def load_combgc_summaries(paths: list[Path]) -> list[dict[str, str | int]]:
    """Load and normalize rows from one or more per-sample combgc_summary.tsv files."""
    calls: list[dict[str, str | int]] = []
    for path in sorted(paths):
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample_raw = (row.get("sample_id") or "").strip()
                contig_id = (row.get("contig_id") or "").strip()
                tool = (row.get("Prediction_tool") or "").strip()
                if not sample_raw or not contig_id or not tool:
                    continue
                try:
                    start = int(float(row.get("BGC_start", 0)))
                    end = int(float(row.get("BGC_end", 0)))
                except (TypeError, ValueError):
                    continue
                if end < start:
                    start, end = end, start
                calls.append(
                    {
                        "sample_id": normalize_sample_id(sample_raw),
                        "contig_id": contig_id,
                        "tool": tool,
                        "product": (row.get("Product_class") or "Other").strip(),
                        "complete": (row.get("BGC_complete") or "NA").strip(),
                        "start": start,
                        "end": end,
                        "mibig_id": (row.get("MIBiG_ID") or "NA").strip(),
                    }
                )
    return calls


def parse_bgc_filename(filename: str) -> tuple[str, str, int | None]:
    """Parse (sample_id, predictor_tool, region_index) from a comBGC/clustering region GBK filename."""
    stem = Path(filename.strip()).stem
    if not stem:
        return "", "", None

    sanntis_m = re.match(r"^(.+?)\.sanntis\.region0*(\d+)$", stem)
    if sanntis_m:
        return normalize_sample_id(sanntis_m.group(1)), "SanntiS", int(sanntis_m.group(2))

    deepbgc_m = re.match(r"^(.+?)\.deepbgc\.region0*(\d+)$", stem)
    if deepbgc_m:
        return normalize_sample_id(deepbgc_m.group(1)), "deepBGC", int(deepbgc_m.group(2))

    gecco_m = re.match(r"^(.+?)\.gecco\..*?(?:region|_cluster_)0*(\d+)$", stem)
    if gecco_m:
        return normalize_sample_id(gecco_m.group(1)), "GECCO", int(gecco_m.group(2))

    gecco_legacy_m = re.match(r"^(.+?)_cluster_0*(\d+)$", stem)
    if gecco_legacy_m:
        return normalize_sample_id(gecco_legacy_m.group(1)), "GECCO", int(gecco_legacy_m.group(2))

    as_m = re.match(r"^(.+?)\.region0*(\d+)$", stem)
    if as_m:
        prefix = as_m.group(1)
        return normalize_sample_id(prefix), "antiSMASH", int(as_m.group(2))

    sid = re.split(r"\.(?:sanntis|deepbgc|gecco|region)|_cluster_", stem)[0]
    return normalize_sample_id(sid), "", None


def _parse_region_index(region_str: str) -> int | None:
    """Parse integer region number from antiSMASH region labels like '1', '2', or '1.2'."""
    text = str(region_str or "").strip()
    if not text:
        return None
    tail = text.rsplit(".", 1)[-1]
    return int(tail) if tail.isdigit() else None


def load_all_regions(path: Path | None) -> dict[tuple[str, str], list[dict[str, str | int]]]:
    """Load antiSMASH all_regions.tsv indexed by (sample_id, record_id)."""
    by_contig: dict[tuple[str, str], list[dict[str, str | int]]] = defaultdict(list)
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return by_contig
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sid = normalize_sample_id(row.get("file") or "")
            rec = (row.get("record_id") or "").strip()
            if not sid or not rec:
                continue
            try:
                start = int(row.get("start", 0))
                end = int(row.get("end", 0))
            except ValueError:
                continue
            by_contig[(sid, rec)].append(
                {
                    "start": start,
                    "end": end,
                    "region": (row.get("region") or "").strip(),
                    "contig_edge": (row.get("contig_edge") or "").strip(),
                    "kcb_hit": (row.get("knownclusterblast_hit") or "").strip(),
                    "kcb_sim": (row.get("knownclusterblast_similarity") or "").strip(),
                }
            )
    return by_contig


def load_quast(paths: list[Path]) -> dict[str, dict[str, str]]:
    """Parse QUAST report.tsv (transposed metrics x samples table) into per-sample dict."""
    metrics: dict[str, dict[str, str]] = defaultdict(dict)
    for path in sorted(paths):
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open("r", encoding="utf-8") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        if len(rows) < 2:
            continue
        header = rows[0]
        sample_cols = [(idx, normalize_sample_id(col)) for idx, col in enumerate(header[1:], start=1) if col.strip()]
        metric_map = {
            "Total length": "assembly_length_bp",
            "# contigs": "num_contigs",
            "N50": "n50",
        }
        for row in rows[1:]:
            if not row:
                continue
            key = metric_map.get(row[0].strip())
            if key is None:
                continue
            for col_idx, sid in sample_cols:
                if col_idx < len(row):
                    metrics[sid][key] = row[col_idx].strip()
    return metrics


def load_checkm2(paths: list[Path]) -> dict[str, dict[str, str]]:
    """Parse CheckM2 quality_report.tsv into per-sample dict."""
    metrics: dict[str, dict[str, str]] = {}
    for path in sorted(paths):
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                name = row.get("Name") or row.get("Bin Id") or ""
                if not name:
                    continue
                sid = normalize_sample_id(name)
                metrics[sid] = {
                    "completeness": (row.get("Completeness") or "NA").strip(),
                    "contamination": (row.get("Contamination") or "NA").strip(),
                    "assembly_length_bp": (row.get("Genome_Size") or "").strip(),
                    "num_contigs": (row.get("Total_Contigs") or "").strip(),
                    "n50": (row.get("Contig_N50") or "").strip(),
                }
    return metrics


def load_gtdbtk(paths: list[Path]) -> dict[str, str]:
    """Parse GTDB-Tk *.summary.tsv file(s) into sample_id -> classification."""
    taxonomy: dict[str, str] = {}
    for path in sorted(paths):
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                user_genome = (row.get("user_genome") or "").strip()
                classification = (row.get("classification") or "").strip()
                if user_genome and classification:
                    taxonomy[normalize_sample_id(user_genome)] = classification
    return taxonomy


def load_gcf_membership(
    path: Path | None,
) -> tuple[dict[str, dict[str, list[str]]], dict[tuple[str, str, int], dict[str, list[str]]]]:
    """Load gcf_membership_summary.tsv into sample-level and locus-level GCF mappings."""
    by_sample: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    by_key: dict[tuple[str, str, int], dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return by_sample, by_key
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sid = normalize_sample_id(row.get("sample_id") or "")
            tool = (row.get("tool") or "").strip()
            gcf = (row.get("gcf_id") or "").strip()
            bgc_file = (row.get("bgc_file") or "").strip()
            if sid and tool and gcf and gcf not in by_sample[sid][tool]:
                by_sample[sid][tool].append(gcf)
            if bgc_file and tool and gcf:
                parsed_sid, pred_tool, reg_idx = parse_bgc_filename(bgc_file)
                key_sid = sid or parsed_sid
                if key_sid and pred_tool and reg_idx is not None:
                    if gcf not in by_key[(key_sid, pred_tool, reg_idx)][tool]:
                        by_key[(key_sid, pred_tool, reg_idx)][tool].append(gcf)
                    # Also index under raw parsed prefix in case filename was <sample>_<contig>.region001
                    if (
                        parsed_sid
                        and parsed_sid != key_sid
                        and gcf not in by_key[(parsed_sid, pred_tool, reg_idx)][tool]
                    ):
                        by_key[(parsed_sid, pred_tool, reg_idx)][tool].append(gcf)
    return by_sample, by_key


def load_bigfam_scores(
    path: Path | None,
) -> tuple[dict[str, list[dict[str, str]]], dict[tuple[str, str, int], list[dict[str, str]]]]:
    """Load bigfam_novelty_scores.tsv into (by_sample, by_(sample, predictor, region_idx))."""
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_key: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return by_sample, by_key
    seen_bgc: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rank = (row.get("rank") or "0").strip()
            if rank not in {"", "0", "1", "NA"}:
                continue
            fname = (row.get("bgc_filename") or row.get("bgc_id") or "").strip()
            if not fname:
                continue
            bgc_key = Path(fname).name
            # Keep only best (first) rank per BGC file to avoid duplicate counting
            if bgc_key in seen_bgc:
                continue
            seen_bgc.add(bgc_key)
            parsed_sid, pred_tool, reg_idx = parse_bgc_filename(fname)
            # Derive base sample_id (before _contig if applicable) for sample-level aggregation
            stem = Path(fname).stem
            sid = normalize_sample_id(re.split(r"\.(?:sanntis|deepbgc|gecco|region)|_cluster_", stem)[0])
            by_sample[sid].append(row)
            if pred_tool and reg_idx is not None:
                by_key[(sid, pred_tool, reg_idx)].append(row)
                if parsed_sid and parsed_sid != sid:
                    by_key[(parsed_sid, pred_tool, reg_idx)].append(row)
    return by_sample, by_key


def load_amp_summary(path: Path | None) -> dict[str, int]:
    """Count AMP candidates per sample from AMPcombi summary TSV."""
    counts: dict[str, int] = defaultdict(int)
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return counts
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = (row.get("sample") or row.get("sample_id") or row.get("Sample") or "").strip()
            if sample:
                counts[normalize_sample_id(sample)] += 1
    return counts


def merge_consensus_loci(
    calls: list[dict[str, str | int]],
    all_regions: dict[tuple[str, str], list[dict[str, str | int]]],
    gcf_by_sample: dict[str, dict[str, list[str]]],
    bigfam_by_sample: dict[str, list[dict[str, str]]],
    min_overlap: float = 0.25,
    gcf_by_key: dict[tuple[str, str, int], dict[str, list[str]]] | None = None,
    bigfam_by_key: dict[tuple[str, str, int], list[dict[str, str]]] | None = None,
) -> list[dict[str, str | int | float]]:
    """Cluster overlapping BGC calls on the same (sample_id, contig_id) via connected components."""
    gcf_key_map = gcf_by_key or {}
    bigfam_key_map = bigfam_by_key or {}

    # Assign per-(sample, tool) and per-(sample, tool, contig) 1-based ordinal indices
    tool_ord_counter: dict[tuple[str, str], int] = defaultdict(int)
    contig_ord_counter: dict[tuple[str, str, str], int] = defaultdict(int)
    enriched_calls: list[dict[str, str | int]] = []
    for call in sorted(
        calls, key=lambda c: (str(c["sample_id"]), str(c["contig_id"]), int(c["start"]), int(c["end"]), str(c["tool"]))
    ):
        sid = str(call["sample_id"])
        cid = str(call["contig_id"])
        tool = str(call["tool"])
        tool_ord_counter[(sid, tool)] += 1
        contig_ord_counter[(sid, tool, cid)] += 1
        entry = dict(call)
        entry["_tool_ord"] = tool_ord_counter[(sid, tool)]
        entry["_contig_ord"] = contig_ord_counter[(sid, tool, cid)]
        enriched_calls.append(entry)

    grouped: dict[tuple[str, str], list[dict[str, str | int]]] = defaultdict(list)
    for call in enriched_calls:
        grouped[(str(call["sample_id"]), str(call["contig_id"]))].append(call)

    # First pass: build consensus components and collect exact locus-level GCF / BiG-FAM matches
    loci: list[dict[str, str | int | float]] = []
    locus_members_meta: list[tuple[str, list[tuple[str, str, int]]]] = []
    sample_has_exact_gcf: dict[str, set[str]] = defaultdict(set)
    sample_has_exact_bigfam: set[str] = set()
    locus_counter = 0

    for (sample_id, contig_id), contig_calls in sorted(grouped.items()):
        contig_calls.sort(key=lambda c: (int(c["start"]), int(c["end"]), str(c["tool"])))
        n = len(contig_calls)
        adj: list[list[tuple[int, float]]] = [[] for _ in range(n)]
        for i in range(n):
            sa, ea = int(contig_calls[i]["start"]), int(contig_calls[i]["end"])
            for j in range(i + 1, n):
                sb, eb = int(contig_calls[j]["start"]), int(contig_calls[j]["end"])
                jacc = interval_jaccard(sa, ea, sb, eb)
                if jacc >= min_overlap:
                    adj[i].append((j, jacc))
                    adj[j].append((i, jacc))

        visited = [False] * n
        for i in range(n):
            if visited[i]:
                continue
            stack = [i]
            visited[i] = True
            comp_indices: list[int] = []
            while stack:
                node = stack.pop()
                comp_indices.append(node)
                for nxt, _jacc in adj[node]:
                    if not visited[nxt]:
                        visited[nxt] = True
                        stack.append(nxt)
            comp_indices.sort()

            members = [contig_calls[idx] for idx in comp_indices]
            locus_counter += 1
            locus_id = f"{sample_id}_LOCUS_{locus_counter:04d}"
            start = min(int(m["start"]) for m in members)
            end = max(int(m["end"]) for m in members)
            tools = sorted({str(m["tool"]) for m in members})
            products = sorted({str(m["product"]) for m in members if str(m["product"]) not in {"", "NA"}})
            mibigs = sorted(
                {token for m in members for token in str(m["mibig_id"]).split(";") if token and token != "NA"}
            )

            # Mean pairwise Jaccard across members (1.0 for singleton)
            if len(members) == 1:
                mean_jacc = 1.0
            else:
                j_vals: list[float] = []
                for a_idx in range(len(members)):
                    for b_idx in range(a_idx + 1, len(members)):
                        j_vals.append(
                            interval_jaccard(
                                int(members[a_idx]["start"]),
                                int(members[a_idx]["end"]),
                                int(members[b_idx]["start"]),
                                int(members[b_idx]["end"]),
                            )
                        )
                mean_jacc = sum(j_vals) / len(j_vals) if j_vals else 1.0

            completes = {str(m["complete"]) for m in members if str(m["complete"]) in {"Yes", "No"}}
            if "No" in completes:
                complete_status = "Partial"
            elif "Yes" in completes:
                complete_status = "Complete"
            else:
                complete_status = "NA"

            # Match antiSMASH KnownClusterBlast hit and exact region index by coordinate overlap
            kcb_hit = "NA"
            kcb_sim = "NA"
            matched_as_reg_num: int | None = None
            for reg in all_regions.get((sample_id, contig_id), []):
                if interval_jaccard(start, end, int(reg["start"]), int(reg["end"])) > 0:
                    if matched_as_reg_num is None:
                        matched_as_reg_num = _parse_region_index(str(reg.get("region", "")))
                    if reg["kcb_hit"]:
                        kcb_hit = str(reg["kcb_hit"])
                        kcb_sim = str(reg["kcb_sim"]) or "NA"
                    if complete_status == "NA" and reg["contig_edge"]:
                        complete_status = "Partial" if str(reg["contig_edge"]).lower() == "true" else "Complete"
                    break

            # Build candidate region keys for exact locus-level GCF / BiG-FAM lookup
            cand_keys: list[tuple[str, str, int]] = []
            for m in members:
                m_tool = str(m["tool"])
                m_ord = int(m["_tool_ord"])
                c_ord = int(m["_contig_ord"])
                if m_tool == "antiSMASH":
                    if matched_as_reg_num is not None:
                        cand_keys.append((sample_id, "antiSMASH", matched_as_reg_num))
                        cand_keys.append((f"{sample_id}_{contig_id}", "antiSMASH", matched_as_reg_num))
                    cand_keys.append((sample_id, "antiSMASH", c_ord))
                    cand_keys.append((f"{sample_id}_{contig_id}", "antiSMASH", c_ord))
                    if m_ord != c_ord:
                        cand_keys.append((sample_id, "antiSMASH", m_ord))
                else:
                    cand_keys.append((sample_id, m_tool, m_ord))
                    cand_keys.append((f"{sample_id}_{contig_id}", m_tool, c_ord))
                    if c_ord != m_ord:
                        cand_keys.append((sample_id, m_tool, c_ord))

            # Deduplicate candidate keys preserving order
            uniq_keys = list(dict.fromkeys(cand_keys))
            for key in uniq_keys:
                if key in gcf_key_map:
                    for clust_tool in gcf_key_map[key]:
                        sample_has_exact_gcf[sample_id].add(clust_tool)
                if key in bigfam_key_map:
                    sample_has_exact_bigfam.add(sample_id)

            locus_members_meta.append((sample_id, uniq_keys))
            loci.append(
                {
                    "locus_id": locus_id,
                    "sample_id": sample_id,
                    "contig_id": contig_id,
                    "start": start,
                    "end": end,
                    "length": max(0, end - start + 1),
                    "tools": ";".join(tools),
                    "tool_count": len(tools),
                    "support_tier": support_tier(len(tools)),
                    "mean_jaccard": round(mean_jacc, 4),
                    "product_classes": ";".join(products) if products else "Other",
                    "complete_status": complete_status,
                    "mibig_ids": ";".join(mibigs) if mibigs else "NA",
                    "knownclusterblast_hit": kcb_hit,
                    "knownclusterblast_similarity": kcb_sim,
                    "bigscape_gcf": "NA",
                    "bigslice_gcf": "NA",
                    "bigfam_gcf": "NA",
                    "bigfam_distance_d": "NA",
                    "is_putative_novel": "NA",
                }
            )

    # Second pass: attach locus-level GCF and BiG-FAM metrics (falling back to sample-level only if no locus matched)
    for locus, (sample_id, uniq_keys) in zip(loci, locus_members_meta, strict=True):
        for clust_tool, col_name in (("BiG-SCAPE", "bigscape_gcf"), ("BiG-SLiCE", "bigslice_gcf")):
            if clust_tool in sample_has_exact_gcf.get(sample_id, set()):
                matched_gcfs: list[str] = []
                for key in uniq_keys:
                    for gcf in gcf_key_map.get(key, {}).get(clust_tool, []):
                        if gcf not in matched_gcfs:
                            matched_gcfs.append(gcf)
                locus[col_name] = ";".join(matched_gcfs) if matched_gcfs else "NA"
            else:
                locus[col_name] = ";".join(gcf_by_sample.get(sample_id, {}).get(clust_tool, [])) or "NA"

        if sample_id in sample_has_exact_bigfam:
            bf_entries: list[dict[str, str]] = []
            seen_bf: set[str] = set()
            for key in uniq_keys:
                for entry in bigfam_key_map.get(key, []):
                    sig = entry.get("bgc_filename") or entry.get("bgc_id") or ""
                    if sig not in seen_bf:
                        seen_bf.add(sig)
                        bf_entries.append(entry)
        else:
            bf_entries = bigfam_by_sample.get(sample_id, [])

        if bf_entries:
            bf_gcfs = sorted({e.get("gcf_id", "NA") for e in bf_entries if e.get("gcf_id")})
            locus["bigfam_gcf"] = ";".join(bf_gcfs) if bf_gcfs else "NA"
            numeric_dists = []
            for e in bf_entries:
                with contextlib.suppress(KeyError, ValueError):
                    numeric_dists.append(float(e["distance_d"]))
            if numeric_dists:
                locus["bigfam_distance_d"] = f"{max(numeric_dists):.2f}"
            locus["is_putative_novel"] = (
                "true" if any(e.get("is_putative_novel") == "true" for e in bf_entries) else "false"
            )

    return loci


def build_genome_summary(
    calls: list[dict[str, str | int]],
    loci: list[dict[str, str | int | float]],
    quast: dict[str, dict[str, str]],
    checkm2: dict[str, dict[str, str]],
    gtdbtk: dict[str, str],
    bigfam_by_sample: dict[str, list[dict[str, str]]],
    amp_counts: dict[str, int],
) -> list[dict[str, str | int | float]]:
    """Aggregate QC, taxonomy, predictor counts, consensus tiers, and density metrics per genome."""
    all_samples = sorted(
        {str(c["sample_id"]) for c in calls}
        | set(quast.keys())
        | set(checkm2.keys())
        | set(gtdbtk.keys())
        | set(amp_counts.keys())
    )

    calls_by_sample: dict[str, list[dict[str, str | int]]] = defaultdict(list)
    for call in calls:
        calls_by_sample[str(call["sample_id"])].append(call)

    loci_by_sample: dict[str, list[dict[str, str | int | float]]] = defaultdict(list)
    for locus in loci:
        loci_by_sample[str(locus["sample_id"])].append(locus)

    rows: list[dict[str, str | int | float]] = []
    for sid in all_samples:
        scalls = calls_by_sample[sid]
        sloci = loci_by_sample[sid]
        q_info = quast.get(sid, {})
        c_info = checkm2.get(sid, {})

        assembly_len_str = q_info.get("assembly_length_bp") or c_info.get("assembly_length_bp") or "NA"
        num_contigs_str = q_info.get("num_contigs") or c_info.get("num_contigs") or "NA"
        n50_str = q_info.get("n50") or c_info.get("n50") or "NA"

        tool_counts: dict[str, int] = defaultdict(int)
        for c in scalls:
            tool_counts[str(c["tool"])] += 1

        high_loci = sum(1 for loc in sloci if int(loc["tool_count"]) >= 3)
        med_loci = sum(1 for loc in sloci if int(loc["tool_count"]) == 2)
        single_loci = sum(1 for loc in sloci if int(loc["tool_count"]) == 1)

        novel_count = sum(1 for e in bigfam_by_sample.get(sid, []) if e.get("is_putative_novel") == "true")

        try:
            assembly_len = float(assembly_len_str)
            bgcs_per_mbp: str | float = (
                round(len(sloci) / (assembly_len / 1_000_000.0), 3) if assembly_len > 0 else "NA"
            )
        except ValueError:
            bgcs_per_mbp = "NA"

        known_completeness = [loc for loc in sloci if loc["complete_status"] in {"Complete", "Partial"}]
        if known_completeness:
            complete_n = sum(1 for loc in known_completeness if loc["complete_status"] == "Complete")
            partial_n = sum(1 for loc in known_completeness if loc["complete_status"] == "Partial")
            complete_frac: str | float = round(complete_n / len(known_completeness), 4)
            edge_rate: str | float = round(partial_n / len(known_completeness), 4)
        else:
            complete_frac = "NA"
            edge_rate = "NA"

        rows.append(
            {
                "sample_id": sid,
                "gtdb_classification": gtdbtk.get(sid, "NA"),
                "completeness": c_info.get("completeness", "NA"),
                "contamination": c_info.get("contamination", "NA"),
                "assembly_length_bp": assembly_len_str,
                "num_contigs": num_contigs_str,
                "n50": n50_str,
                "antismash_bgcs": tool_counts["antiSMASH"],
                "gecco_bgcs": tool_counts["GECCO"],
                "deepbgc_bgcs": tool_counts["deepBGC"],
                "sanntis_bgcs": tool_counts["SanntiS"],
                "consensus_loci": len(sloci),
                "high_support_loci": high_loci,
                "medium_support_loci": med_loci,
                "single_tool_loci": single_loci,
                "putative_novel_bgcs": novel_count,
                "bgcs_per_mbp": bgcs_per_mbp,
                "complete_bgc_fraction": complete_frac,
                "contig_edge_rate": edge_rate,
                "amp_count": amp_counts.get(sid, 0),
            }
        )
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str | int | float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_multiqc_files(
    outdir: Path,
    genome_rows: list[dict[str, str | int | float]],
    has_amp: bool = False,
) -> None:
    """Emit MultiQC custom-content TSV tables."""
    support_rows = [
        {
            "sample_id": row["sample_id"],
            "High_3_to_4_tools": row["high_support_loci"],
            "Medium_2_tools": row["medium_support_loci"],
            "Single_tool": row["single_tool_loci"],
        }
        for row in genome_rows
    ]
    write_tsv(
        outdir / "consensus_support_mqc.tsv",
        ["sample_id", "High_3_to_4_tools", "Medium_2_tools", "Single_tool"],
        support_rows,
    )

    overview_cols = [
        "sample_id",
        "gtdb_classification",
        "completeness",
        "contamination",
        "n50",
        "antismash_bgcs",
        "gecco_bgcs",
        "deepbgc_bgcs",
        "sanntis_bgcs",
        "consensus_loci",
        "putative_novel_bgcs",
        "bgcs_per_mbp",
        "contig_edge_rate",
    ]
    write_tsv(outdir / "genome_overview_mqc.tsv", overview_cols, genome_rows)

    if has_amp:
        amp_rows = [{"sample_id": row["sample_id"], "amp_candidates": row["amp_count"]} for row in genome_rows]
        write_tsv(outdir / "amp_summary_mqc.tsv", ["sample_id", "amp_candidates"], amp_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build cross-tool consensus BGC loci and genome QC/taxonomy summary.")
    parser.add_argument("--combgc", type=Path, nargs="*", default=[], help="Per-sample combgc_summary.tsv file(s).")
    parser.add_argument("--all-regions", type=Path, default=None, help="Optional antiSMASH all_regions.tsv.")
    parser.add_argument("--gcf-membership", type=Path, default=None, help="Optional gcf_membership_summary.tsv.")
    parser.add_argument("--bigfam-scores", type=Path, default=None, help="Optional bigfam_novelty_scores.tsv.")
    parser.add_argument("--quast", type=Path, nargs="*", default=[], help="Optional QUAST report.tsv file(s).")
    parser.add_argument(
        "--checkm2", type=Path, nargs="*", default=[], help="Optional CheckM2 quality_report.tsv file(s)."
    )
    parser.add_argument("--gtdbtk", type=Path, nargs="*", default=[], help="Optional GTDB-Tk *.summary.tsv file(s).")
    parser.add_argument("--amp-summary", type=Path, default=None, help="Optional AMPcombi summary TSV.")
    parser.add_argument(
        "--min-overlap",
        type=float,
        default=0.25,
        help="Minimum interval Jaccard overlap to merge cross-tool calls (default: 0.25).",
    )
    parser.add_argument("-o", "--outdir", type=Path, default=Path(), help="Output directory.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    calls = load_combgc_summaries(args.combgc)
    all_regions = load_all_regions(args.all_regions)
    gcf_by_sample, gcf_by_key = load_gcf_membership(args.gcf_membership)
    bigfam_by_sample, bigfam_by_key = load_bigfam_scores(args.bigfam_scores)
    quast = load_quast(args.quast)
    checkm2 = load_checkm2(args.checkm2)
    gtdbtk = load_gtdbtk(args.gtdbtk)
    amp_counts = load_amp_summary(args.amp_summary)

    loci = merge_consensus_loci(
        calls=calls,
        all_regions=all_regions,
        gcf_by_sample=gcf_by_sample,
        bigfam_by_sample=bigfam_by_sample,
        min_overlap=args.min_overlap,
        gcf_by_key=gcf_by_key,
        bigfam_by_key=bigfam_by_key,
    )
    genome_rows = build_genome_summary(
        calls=calls,
        loci=loci,
        quast=quast,
        checkm2=checkm2,
        gtdbtk=gtdbtk,
        bigfam_by_sample=bigfam_by_sample,
        amp_counts=amp_counts,
    )

    write_tsv(args.outdir / "consensus_loci.tsv", CONSENSUS_COLS, loci)
    write_tsv(args.outdir / "genome_qc_bgc_summary.tsv", GENOME_SUMMARY_COLS, genome_rows)
    write_multiqc_files(args.outdir, genome_rows, has_amp=args.amp_summary is not None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
