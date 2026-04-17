#!/usr/bin/env python3
## Given a bunch of antismash results, tabulate BGC regions
#
## Usage:
#   $ python antismash/tabulate_regions.py -h
#   usage: tabulate_regions.py [-h] [--knownclusters] [--threads THREADS] directory output
#
#   Given a bunch of antismash results, tabulate BGC regions
#
#   positional arguments:
#     directory            Directory containing antiSMASH directories
#     output               Desired path/to/filename for the output TSV
#
#   options:
#     -h, --help           show this help message and exit
#     --knownclusters      Include KnownClusterBlast columns in the output.
#                          Should be set when antiSMASH was run with --cb-knownclusters.
#     --threads THREADS    Number of threads to use for parallel processing. Defaults to number of CPUs.

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
from pathlib import Path
import re


def extract_top_known_cluster_hit(known_cluster_blast_results: list | None, region_index: int) -> dict[str, str]:
    """Extract the top KnownClusterBlast hit for a biosynthetic region.

    Args:
        known_cluster_blast_results: The KnownClusterBlast results list for all
            regions in the record, or None if the module was not run.
        region_index: Index of the region within the record's region list.

    Returns:
        Dictionary with knownclusterblast_hit, knownclusterblast_accession, and
        knownclusterblast_similarity keys. All values are empty strings when no
        significant hit is found.
    """
    no_hit = {"knownclusterblast_hit": "", "knownclusterblast_accession": "", "knownclusterblast_similarity": ""}

    if not known_cluster_blast_results:
        return no_hit

    ranked_hits = known_cluster_blast_results[region_index]["ranking"]
    if not ranked_hits:
        return no_hit

    percent_similarity = ranked_hits[0][1]["similarity"]
    if percent_similarity <= 15:
        return no_hit

    match percent_similarity:
        case s if s > 75:
            similarity_category = "high"
        case s if s > 50:
            similarity_category = "medium"
        case _:
            similarity_category = "low"

    return {
        "knownclusterblast_hit": ranked_hits[0][0]["description"],
        "knownclusterblast_accession": ranked_hits[0][0]["accession"],
        "knownclusterblast_similarity": similarity_category,
    }


def extract_known_cluster_blast_results(genomic_record: dict) -> list | None:
    """Safely extract KnownClusterBlast results from an antiSMASH genomic record.

    Navigates the nested module results structure to retrieve per-region
    KnownClusterBlast rankings.

    Args:
        genomic_record: A dictionary representing one antiSMASH sequence record,
            containing a nested ``modules`` key with clusterblast results.

    Returns:
        A list of per-region KnownClusterBlast result dictionaries if found,
        or None if the clusterblast module was not run or the key path is absent.
    """
    try:
        return genomic_record["modules"]["antismash.modules.clusterblast"]["knowncluster"]["results"]
    except (KeyError, TypeError, AttributeError):
        return None


def parse_antismash_json(json_path: Path) -> list[dict[str, str]]:
    """Parse an antiSMASH JSON output file and extract per-region annotation data.

    Args:
        json_path: Path to the antiSMASH JSON file for one genome.

    Returns:
        List of row dictionaries, one per biosynthetic region, containing
        genomic coordinates, product class, and optionally KnownClusterBlast hits.
    """
    with json_path.open() as json_file:
        antismash_data = json.load(json_file)

    region_rows: list[dict[str, str]] = []

    for genomic_record in antismash_data["records"]:
        if not genomic_record["areas"]:
            continue

        region_features = [feat for feat in genomic_record["features"] if feat["type"] == "region"]
        known_cluster_blast_results = extract_known_cluster_blast_results(genomic_record)

        for region_index, region_feature in enumerate(region_features):
            start, end = re.findall(r"\d+", region_feature["location"])
            region_qualifiers = region_feature["qualifiers"]

            region_row = {
                "file": json_path.stem,
                "record_id": genomic_record["name"],
                "region": region_qualifiers["region_number"][0],
                "start": start,
                "end": end,
                "contig_edge": region_qualifiers["contig_edge"][0],
                "product": " / ".join(region_qualifiers["product"]),
                "record_desc": genomic_record["description"],
            } | extract_top_known_cluster_hit(known_cluster_blast_results, region_index)

            region_rows.append(region_row)

    return region_rows


def main(directory: Path, output: Path, knownclusters: bool = False, threads: int | None = None):
    region_rows: list[dict[str, str]] = []

    antismash_json_files = list(directory.glob("*/*.json"))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        parse_futures = {
            executor.submit(parse_antismash_json, json_path): json_path for json_path in antismash_json_files
        }

        for parse_future in as_completed(parse_futures):
            region_rows.extend(parse_future.result())

    output_columns = [
        "file",
        "record_id",
        "region",
        "start",
        "end",
        "contig_edge",
        "product",
    ]
    if knownclusters:
        output_columns += [
            "knownclusterblast_hit",
            "knownclusterblast_accession",
            "knownclusterblast_similarity",
        ]
    output_columns.append("record_desc")

    with output.open("w") as output_file:
        tsv_writer = csv.DictWriter(output_file, fieldnames=output_columns, delimiter="\t", extrasaction="ignore")
        tsv_writer.writeheader()
        tsv_writer.writerows(region_rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Given a bunch of antismash results, tabulate BGC regions")

    parser.add_argument("directory", type=Path, help="Directory containing antiSMASH directories")
    parser.add_argument("output", type=Path, help="Desired path/to/filename for the output TSV")
    parser.add_argument(
        "--knownclusters",
        action="store_true",
        help=(
            "Include KnownClusterBlast columns in the output. "
            "Should be set when antiSMASH was run with --cb-knownclusters."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads to use for parallel processing. Defaults to number of CPUs.",
    )

    args = parser.parse_args()

    main(args.directory, args.output, args.knownclusters, args.threads)
