#!/usr/bin/env python3
"""Count BGC regions from antiSMASH JSON outputs."""

import argparse
from collections import Counter
import csv
from pathlib import Path

from antismash_io import list_antismash_json_files, load_antismash_json, run_threaded


def parse_antismash_json(json_path: Path) -> tuple[str, dict[str, list[list[str]]], dict[str, str]]:
    """Parse an antiSMASH JSON output file and extract per-contig BGC product data.

    Args:
        json_path: Path to the antiSMASH JSON file for one genome.

    Returns:
        A tuple containing:
            - str: The original input file name recorded in the JSON.
            - dict[str, list[list[str]]]: BGC product lists per biosynthetic area,
                keyed by sequence record (contig) name.
            - dict[str, str]: Sequence descriptions keyed by record name.
    """
    antismash_data = load_antismash_json(json_path)

    products_by_contig = {
        record["name"]: [biosynthetic_area["products"] for biosynthetic_area in record["areas"]]
        for record in antismash_data["records"]
    }

    contig_descriptions = {record["name"]: record.get("description", "") for record in antismash_data["records"]}

    return antismash_data["input_file"], products_by_contig, contig_descriptions


def tabulate_bgc_counts(
    genome_products: dict[str, dict[str, list[list[str]]]],
    contig_descriptions: dict[str, dict[str, str]],
    per_contig: bool = False,
    split_hybrids: bool = False,
) -> list[dict[str, str | int]]:
    """Tabulate BGC product-type counts from parsed antiSMASH data.

    Converts a nested mapping of genome → contig → biosynthetic region products
    into a flat list of count dictionaries suitable for TSV output.

    Args:
        genome_products: Nested mapping of genome name → contig name → list of
            BGC areas, where each area is a list of product type strings.
        contig_descriptions: Nested mapping of genome name → contig name →
            sequence description string.
        per_contig: If True, emit one row per contig. If False, aggregate counts
            across all contigs for each genome assembly. Defaults to False.
        split_hybrids: If True, count each product type in a hybrid BGC
            individually. If False, count multi-product regions as a single
            "hybrid" category. Defaults to False.

    Returns:
        A list of row dictionaries, each containing BGC product-type counts,
        a ``record`` identifier, a ``total_count``, and a ``description``.
    """
    count_rows = []

    for genome_name, contig_products in genome_products.items():
        for contig_name, bgc_areas in contig_products.items():
            bgc_type_counts: Counter[str] = Counter()
            for bgc_products in bgc_areas:
                if len(bgc_products) > 1 and not split_hybrids:
                    bgc_type_counts["hybrid"] += 1
                else:
                    bgc_type_counts.update(bgc_products)

            if per_contig:
                count_rows.append(
                    {
                        **bgc_type_counts,
                        "record": f"{genome_name}|{contig_name}",
                        "total_count": len(bgc_areas),
                        "description": contig_descriptions[genome_name][contig_name],
                    }
                )

        if not per_contig:
            # Aggregate BGC type counts across all contigs for this genome assembly
            genome_bgc_counts: Counter[str] = Counter()
            for bgc_areas in contig_products.values():
                for bgc_products in bgc_areas:
                    if len(bgc_products) > 1 and not split_hybrids:
                        genome_bgc_counts["hybrid"] += 1
                    else:
                        genome_bgc_counts.update(bgc_products)

            representative_contig = next(iter(contig_products))
            contig_count = len(contig_products)
            plural_suffix = "s" if contig_count > 1 else ""
            genome_description = (
                f"{contig_descriptions[genome_name][representative_contig]} "
                f"[{contig_count} total record{plural_suffix}]"
            )
            count_rows.append(
                {
                    **genome_bgc_counts,
                    "record": genome_name,
                    "total_count": sum(len(areas) for areas in contig_products.values()),
                    "description": genome_description,
                }
            )

    return count_rows


def main(
    directory: Path,
    output_path: Path,
    per_contig: bool = False,
    split_hybrids: bool = False,
    threads: int | None = None,
):
    """Main function to process antiSMASH results and output BGC counts."""
    genome_products: dict[str, dict[str, list[list[str]]]] = {}
    contig_descriptions: dict[str, dict[str, str]] = {}

    antismash_json_files = list_antismash_json_files(directory)
    for genome_name, contig_products, genome_contig_descriptions in run_threaded(
        antismash_json_files,
        parse_antismash_json,
        threads,
    ):
        genome_products[genome_name] = contig_products
        contig_descriptions[genome_name] = genome_contig_descriptions

    count_rows = tabulate_bgc_counts(genome_products, contig_descriptions, per_contig, split_hybrids)

    # Collect all BGC product types observed across all genomes
    all_bgc_product_types = {k for row in count_rows for k in row} - {
        "record",
        "total_count",
        "hybrid",
        "description",
    }

    output_columns = ["record", "total_count", *sorted(all_bgc_product_types)]
    if not split_hybrids:
        output_columns.append("hybrid")
    output_columns.append("description")

    with output_path.open("w") as output_file:
        tsv_writer = csv.DictWriter(output_file, fieldnames=output_columns, delimiter="\t", restval=0)
        tsv_writer.writeheader()
        tsv_writer.writerows(count_rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Given a bunch of antismash results, count the BGC regions")

    parser.add_argument("directory", type=Path, help="Directory containing antiSMASH directories")
    parser.add_argument("output", type=Path, help="Desired path+name for the output TSV")
    parser.add_argument(
        "--by_contig",
        action="store_true",
        help="Count regions per each individual contig rather than per assembly",
    )
    parser.add_argument(
        "--split_hybrids",
        action="store_true",
        help=(
            "Count each hybrid region multiple times, once for each "
            "constituent BGC class. The total_count column is unaffected."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads to use for parallel processing. Defaults to number of CPUs.",
    )

    args = parser.parse_args()

    main(args.directory, args.output, args.by_contig, args.split_hybrids, args.threads)
