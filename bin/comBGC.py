#!/usr/bin/env python3
"""Aggregate BGC predictions from antiSMASH, deepBGC, and GECCO into combgc_summary.tsv.

Written by Jasmin Frangenberg and released under the MIT license.

MIT License

Copyright (c) 2023 Jasmin Frangenberg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import argparse
from pathlib import Path
import re
import sys

from Bio import SeqIO
import pandas as pd

TOOL_VERSION = "0.6.3"

WELCOME = f"""\
                ........................
                    * comBGC v.{TOOL_VERSION} *
                ........................
    This tool aggregates the results of BGC prediction tools:
                antiSMASH, deepBGC, and GECCO
    For detailed usage documentation please refer
    to https://nf-co.re/funcscan
    ........................................................."""

INTERNAL_COLS = [
    "Sample_ID",
    "Prediction_tool",
    "Contig_ID",
    "Product_class",
    "BGC_probability",
    "BGC_complete",
    "BGC_start",
    "BGC_end",
    "BGC_length",
    "CDS_ID",
    "CDS_count",
    "PFAM_domains",
    "MIBiG_ID",
    "InterPro_ID",
]

OUTPUT_COLS = [
    "sample_id",
    "contig_id",
    "Prediction_tool",
    "Product_class",
    "BGC_probability",
    "BGC_complete",
    "BGC_start",
    "BGC_end",
    "BGC_length",
    "BGC_region_contig_ids",
    "CDS_count",
    "PFAM_domains",
    "MIBiG_ID",
    "InterPro_ID",
]


def prepare_multisample_input_antismash(antismash_dir: str | Path) -> list[list[str]]:
    """Prepare input path groups from a multi-sample antiSMASH output directory."""
    sample_paths: list[list[str]] = []
    for index_file in sorted(Path(antismash_dir).rglob("index.html")):
        sample_dir = index_file.parent
        sample = sample_dir.name
        gbk_path = sample_dir / f"{sample}.gbk"
        kcb_path = sample_dir / "knownclusterblast"
        if kcb_path.exists():
            sample_paths.append([str(gbk_path), str(kcb_path)])
        else:
            sample_paths.append([str(gbk_path)])
    return sample_paths


def parse_knownclusterblast(kcb_file_path: str | Path) -> list[str]:
    """Extract MIBiG IDs from a KnownClusterBlast TXT file."""
    with Path(kcb_file_path).open() as kcb_file:
        in_hits = False
        mibig_ids: list[str] = []

        for line in kcb_file:
            if line.strip() == "Significant hits:" and not in_hits:
                in_hits = True
            elif line == "\n" and in_hits:
                break
            elif line.strip() != "Significant hits:" and in_hits:
                match = re.search(r"(BGC\d+)", line)
                if match:
                    mibig_ids.append(match.group(1))
    return mibig_ids


def _antismash_row(
    sample_id: str,
    contig_id: str,
    product_class: list[str] | str,
    bgc_complete: str,
    bgc_start: int | str,
    bgc_end: int | str,
    bgc_length: int | str,
    cds_ids: list[str],
    cds_count: int,
    pfam_domains: list[str],
    mibig_id: str,
) -> dict[str, str | int]:
    return {
        "Sample_ID": sample_id,
        "Prediction_tool": "antiSMASH",
        "Contig_ID": contig_id,
        "Product_class": ";".join(product_class) if isinstance(product_class, list) else product_class,
        "BGC_probability": "NA",
        "BGC_complete": bgc_complete,
        "BGC_start": bgc_start,
        "BGC_end": bgc_end,
        "BGC_length": bgc_length,
        "CDS_ID": ";".join(cds_ids) if cds_ids else "NA",
        "CDS_count": cds_count,
        "PFAM_domains": ";".join(pfam_domains) if pfam_domains else "NA",
        "MIBiG_ID": mibig_id,
        "InterPro_ID": "NA",
    }


def antismash_workflow(antismash_paths: list[str], verbose: bool = False) -> pd.DataFrame:
    """Create a DataFrame with aggregated antiSMASH BGC predictions."""
    kcb_path: Path | None = None
    gbk_path: Path | None = None
    for raw_path in antismash_paths:
        candidate = Path(raw_path)
        if candidate.name == "knownclusterblast" or "knownclusterblast" in candidate.parts:
            kcb_path = candidate
        else:
            gbk_path = candidate

    if gbk_path is None or not gbk_path.is_file():
        return pd.DataFrame(columns=INTERNAL_COLS)

    kcb_files: set[str] = set()
    if kcb_path is not None and kcb_path.is_dir():
        kcb_files = {p.name for p in kcb_path.iterdir() if p.is_file() and p.name.endswith(".txt")}

    sample_id = gbk_path.stem
    if verbose:
        print(f"\nParsing antiSMASH file(s): {sample_id}\n... ", end="")

    rows: list[dict[str, str | int]] = []
    with gbk_path.open() as gbk:
        for record in SeqIO.parse(gbk, "genbank"):
            cluster_num = 1
            contig_id = record.id
            product_class: list[str] = []
            bgc_complete = ""
            bgc_start: int | str = ""
            bgc_end: int | str = ""
            bgc_length: int | str = ""
            pfam_domains: list[str] = []
            cds_ids: list[str] = []
            cds_count = 0
            mibig_id = "NA"
            in_cluster = False

            for feature in record.features:
                if feature.type == "protocluster":
                    if in_cluster:
                        rows.append(
                            _antismash_row(
                                sample_id,
                                contig_id,
                                product_class,
                                bgc_complete,
                                bgc_start,
                                bgc_end,
                                bgc_length,
                                cds_ids,
                                cds_count,
                                pfam_domains,
                                mibig_id,
                            )
                        )
                        cds_ids = []
                        cds_count = 0
                        pfam_domains = []
                        mibig_id = "NA"

                    in_cluster = True
                    product_class = [p[0].upper() + p[1:] if p else p for p in feature.qualifiers.get("product", [])]

                    contig_edge = feature.qualifiers.get("contig_edge", [])
                    if contig_edge == ["True"]:
                        bgc_complete = "No"
                    elif contig_edge == ["False"]:
                        bgc_complete = "Yes"
                    else:
                        bgc_complete = "NA"

                    bgc_start = int(feature.location.start) + 1
                    bgc_end = int(feature.location.end)
                    bgc_length = int(feature.location.end) - int(feature.location.start)

                    if kcb_files and kcb_path is not None:
                        kcb_file = f"{record.id}_c{cluster_num}.txt"
                        if kcb_file in kcb_files:
                            mibig_joined = ";".join(parse_knownclusterblast(kcb_path / kcb_file))
                            if mibig_joined:
                                mibig_id = mibig_joined
                            cluster_num += 1

                elif (
                    feature.type == "CDS"
                    and "translation" in feature.qualifiers
                    and in_cluster
                    and isinstance(bgc_end, int)
                    and int(feature.location.end) <= bgc_end
                ):
                    if "locus_tag" in feature.qualifiers:
                        cds_ids.append(feature.qualifiers["locus_tag"][0])
                    cds_count += 1
                    for pfam_domain in feature.qualifiers.get("sec_met_domain", []):
                        match = re.search(r"(.+) \(E-value", pfam_domain)
                        pfam_domains.append(match.group(1) if match else pfam_domain)

            if in_cluster:
                rows.append(
                    _antismash_row(
                        sample_id,
                        contig_id,
                        product_class,
                        bgc_complete,
                        bgc_start,
                        bgc_end,
                        bgc_length,
                        cds_ids,
                        cds_count,
                        pfam_domains,
                        mibig_id,
                    )
                )

    if verbose:
        print("Done.")
    return pd.DataFrame(rows, columns=INTERNAL_COLS)


def deepbgc_workflow(deepbgc_path: str | Path, verbose: bool = False) -> pd.DataFrame:
    """Create a DataFrame with aggregated deepBGC predictions."""
    if verbose:
        print("\nParsing deepBGC file\n... ", end="")

    deepbgc_map_dict = {
        "sequence_id": "Contig_ID",
        "nucl_start": "BGC_start",
        "nucl_end": "BGC_end",
        "nucl_length": "BGC_length",
        "num_proteins": "CDS_count",
        "deepbgc_score": "BGC_probability",
        "product_class": "Product_class",
        "protein_ids": "CDS_ID",
        "pfam_ids": "PFAM_domains",
    }
    deepbgc_unused_cols = [
        "detector_version",
        "detector_label",
        "bgc_candidate_id",
        "num_domains",
        "num_bio_domains",
        "product_activity",
        "antibacterial",
        "cytotoxic",
        "inhibitor",
        "antifungal",
        "Alkaloid",
        "NRP",
        "Other",
        "Polyketide",
        "RiPP",
        "Saccharide",
        "Terpene",
        "bio_pfam_ids",
    ]

    path = Path(deepbgc_path)
    sample = path.name.rsplit(".bgc", 1)[0]
    deepbgc_df = pd.read_csv(path, sep="\t")
    if deepbgc_df.empty:
        if verbose:
            print("Done.")
        return pd.DataFrame(columns=INTERNAL_COLS)

    deepbgc_df = deepbgc_df.drop(deepbgc_unused_cols, axis=1, errors="ignore").rename(columns=deepbgc_map_dict)
    deepbgc_df["Sample_ID"] = sample
    deepbgc_df["Prediction_tool"] = "deepBGC"
    deepbgc_df["BGC_complete"] = "NA"
    deepbgc_df["MIBiG_ID"] = "NA"
    deepbgc_df["InterPro_ID"] = "NA"

    for col in INTERNAL_COLS:
        if col not in deepbgc_df.columns:
            deepbgc_df[col] = "NA"

    if verbose:
        print("Done.")
    return deepbgc_df[INTERNAL_COLS]


def get_interpro_id(gbk_path: str | Path) -> str:
    """Retrieve InterPro IDs from a GECCO GBK file."""
    ip_ids: list[str] = []
    id_pattern = r'InterPro:(.*)"'
    with Path(gbk_path).open() as gbk:
        for line in gbk:
            if "InterPro:" in line:
                match = re.search(id_pattern, line)
                if match:
                    ip_ids.append(match.group(1))
    return ";".join(ip_ids) if ip_ids else "NA"


def gecco_workflow(gecco_paths: list[str], verbose: bool = False) -> pd.DataFrame:
    """Create a DataFrame with aggregated GECCO predictions."""
    if verbose:
        print("\nParsing GECCO files\n... ", end="")

    map_dict = {
        "sequence_id": "Contig_ID",
        "bgc_id": "cluster_id",
        "type": "Product_class",
        "average_p": "BGC_probability",
        "start": "BGC_start",
        "end": "BGC_end",
        "domains": "PFAM_domains",
        "proteins": "CDS_ID",
    }
    unused_cols = [
        "max_p",
        "alkaloid_probability",
        "polyketide_probability",
        "ripp_probability",
        "saccharide_probability",
        "terpene_probability",
        "nrp_probability",
    ]

    tsv_path: Path | None = None
    gbk_paths: list[Path] = []
    for raw_path in gecco_paths:
        candidate = Path(raw_path)
        if candidate.suffix == ".tsv":
            tsv_path = candidate
        else:
            gbk_paths.append(candidate)

    if tsv_path is None:
        return pd.DataFrame(columns=INTERNAL_COLS)

    sample = tsv_path.name.split(".")[0]
    gecco_df = pd.read_csv(tsv_path, sep="\t")
    if gecco_df.empty:
        if verbose:
            print("Done.")
        return pd.DataFrame(columns=INTERNAL_COLS)

    gecco_df = gecco_df.drop(unused_cols, axis=1, errors="ignore").rename(columns=map_dict)
    gecco_df["Sample_ID"] = sample
    gecco_df["BGC_length"] = gecco_df["BGC_end"] - gecco_df["BGC_start"]
    gecco_df["CDS_count"] = [
        0 if pd.isna(val) or not str(val).strip() else len(str(val).split(";")) for val in gecco_df["CDS_ID"]
    ]
    gecco_df["Prediction_tool"] = "GECCO"
    gecco_df["InterPro_ID"] = "NA"

    for gbk_path in gbk_paths:
        bgc_id = re.sub(r"\.region0*(\d+)$", r"_cluster_\1", gbk_path.stem)
        if "cluster_id" in gecco_df.columns:
            gecco_df.loc[gecco_df["cluster_id"] == bgc_id, "InterPro_ID"] = get_interpro_id(gbk_path)

    gecco_df["BGC_complete"] = "NA"
    gecco_df["MIBiG_ID"] = "NA"
    gecco_df["PFAM_domains"] = gecco_df["PFAM_domains"].fillna("NA")

    for col in INTERNAL_COLS:
        if col not in gecco_df.columns:
            gecco_df[col] = "NA"

    if verbose:
        print("Done.")
    return gecco_df[INTERNAL_COLS]


def build_parser() -> argparse.ArgumentParser:
    """Create the comBGC CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="comBGC",
        formatter_class=argparse.RawTextHelpFormatter,
        description=WELCOME,
        add_help=True,
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="PATH(s)",
        dest="input",
        nargs="*",
        help="""path(s) to the required output file(s) of antiSMASH, DeepBGC and/or GECCO
these can be:
- antiSMASH: <sample name>.gbk and (optional) knownclusterblast/ directory
- DeepBGC:   <sample name>.bgc.tsv
- GECCO:     <sample name>.clusters.tsv
Note: Please provide files from a single sample only. If you would like to
summarize multiple samples, please see the --antismash_multiple_samples flag.""",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        metavar="PATH",
        dest="outdir",
        nargs="?",
        help="directory for comBGC output. Default: current directory",
        type=str,
        default=".",
    )
    parser.add_argument(
        "-a",
        "--antismash_multiple_samples",
        metavar="PATH",
        dest="antismash_multiple_samples",
        nargs="?",
        help="""directory of antiSMASH output. Should contain subfolders (one per
sample). Can only be used if --input is not specified.""",
        type=str,
    )
    parser.add_argument("-vv", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-v", "--version", help="show version number and exit", action="store_true")
    return parser


def classify_inputs(input_paths: list[str] | None) -> tuple[list[str], str, list[str]]:
    """Categorize input paths into antiSMASH, deepBGC, and GECCO lists."""
    input_antismash: list[str] = []
    input_deepbgc = ""
    input_gecco: list[str] = []

    if not input_paths:
        return input_antismash, input_deepbgc, input_gecco

    for raw_path in input_paths:
        path_obj = Path(raw_path)
        if raw_path.endswith(".gbk"):
            with path_obj.open() as infile:
                for line in infile:
                    if "##GECCO-Data-START##" in line:
                        input_gecco.append(raw_path)
                        break
                    if "##antiSMASH-Data-START##" in line and not re.search(r"region\d\d\d\.gbk$", raw_path):
                        input_antismash.append(raw_path)
                        break
        elif raw_path.endswith("bgc.tsv"):
            input_deepbgc = raw_path
        elif raw_path.endswith("clusters.tsv"):
            input_gecco.append(raw_path)
        elif raw_path.rstrip("/").endswith("knownclusterblast") or (
            path_obj.is_dir() and path_obj.name == "knownclusterblast"
        ):
            input_antismash.append(raw_path)

    return input_antismash, input_deepbgc, input_gecco


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for comBGC aggregation."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(f"comBGC {TOOL_VERSION}")
        return 0

    if args.input and args.antismash_multiple_samples:
        sys.exit(
            "The flags --input and --antismash_multiple_samples are mutually exclusive.\n"
            "Please use only one of them (or see --help for how to use)."
        )

    input_antismash, input_deepbgc, input_gecco = classify_inputs(args.input)
    dir_antismash = args.antismash_multiple_samples

    if not (input_antismash or input_gecco or input_deepbgc or dir_antismash):
        sys.exit("Please specify at least one input file (i.e. output from antismash, deepbgc, or gecco) or see --help")

    tools_provided: dict[str, list[str] | str] = {}
    if input_antismash:
        tools_provided["antiSMASH"] = input_antismash
    elif dir_antismash:
        tools_provided["antiSMASH"] = dir_antismash
    if input_deepbgc:
        tools_provided["deepBGC"] = input_deepbgc
    if input_gecco:
        tools_provided["GECCO"] = input_gecco

    if args.verbose:
        print(WELCOME)
        print("\nYou provided input for: " + ", ".join(tools_provided))

    frames: list[pd.DataFrame] = []
    for tool, tool_input in tools_provided.items():
        if tool == "antiSMASH":
            if dir_antismash:
                for sample_paths in prepare_multisample_input_antismash(dir_antismash):
                    frames.append(antismash_workflow(sample_paths, verbose=args.verbose))
            else:
                assert isinstance(tool_input, list)
                frames.append(antismash_workflow(tool_input, verbose=args.verbose))
        elif tool == "deepBGC":
            assert isinstance(tool_input, str)
            frames.append(deepbgc_workflow(tool_input, verbose=args.verbose))
        elif tool == "GECCO":
            assert isinstance(tool_input, list)
            frames.append(gecco_workflow(tool_input, verbose=args.verbose))

    non_empty_frames = [df for df in frames if not df.empty]
    if non_empty_frames:
        summary_all = pd.concat(non_empty_frames, ignore_index=True)
        summary_all = summary_all.sort_values(
            by=["Sample_ID", "Contig_ID", "BGC_start", "BGC_length", "Prediction_tool"],
            axis=0,
        )
        summary_all = summary_all.iloc[:, [0, 2, 1] + list(range(3, len(summary_all.columns)))]
        summary_all = summary_all.rename(
            columns={"Sample_ID": "sample_id", "Contig_ID": "contig_id", "CDS_ID": "BGC_region_contig_ids"},
        )
    else:
        summary_all = pd.DataFrame(columns=OUTPUT_COLS)

    outdir_path = Path(args.outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    summary_file = outdir_path / "combgc_summary.tsv"
    summary_all.to_csv(summary_file, sep="\t", index=False)
    print(f"Your BGC summary file is: {summary_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
