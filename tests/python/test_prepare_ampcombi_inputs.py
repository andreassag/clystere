"""Unit tests for bin/prepare_ampcombi_inputs.py."""

import gzip

from Bio import SeqIO
import prepare_ampcombi_inputs as pai


def test_prepare_ampcombi_inputs_unpacks_and_builds_matched_gbk(tmp_path):
    pred_gz = tmp_path / "sample1.prediction.gz"
    with gzip.open(pred_gz, "wt", encoding="utf-8") as fh:
        fh.write(
            "# Prediction from macrel v1.4.0\n"
            "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability\n"
            "ctg1_1\tMKKLLVVAA\tADP\t0.91\tNHEMO\t0.10\n"
        )

    faa_gz = tmp_path / "sample1.all_orfs.faa.gz"
    with gzip.open(faa_gz, "wt", encoding="utf-8") as fh:
        fh.write(">ctg1_1 # 1 # 30 # 1 # ID=1_1\nMKKLLVVAA*\n>ctg1_2 # 40 # 90 # -1 # ID=1_2\nMAAAACCCGGGTT*\n")

    fasta = tmp_path / "sample1.fna"
    fasta.write_text(">ctg1\n" + ("ATG" * 40) + "\n")

    as_gbk = tmp_path / "antismash.gbk"
    as_gbk.write_text(
        "LOCUS       ctg1                     120 bp    DNA     linear   BCT 01-JAN-2026\n"
        "FEATURES             Location/Qualifiers\n"
        "     CDS             complement(40..90)\n"
        '                     /locus_tag="AS_0002"\n'
        '                     /product="ABC transporter ATP-binding protein"\n'
        '                     /translation="MAAAACCCGGGTT"\n'
        "ORIGIN\n"
        "        1 " + ("atg" * 20) + "\n"
        "       61 " + ("atg" * 20) + "\n"
        "//\n"
    )

    out_pred = tmp_path / "sample1.macrel.prediction"
    out_faa = tmp_path / "sample1.faa"
    out_gbk = tmp_path / "sample1.gbk"
    rc = pai.main(
        [
            "--prediction",
            str(pred_gz),
            "--orfs-faa",
            str(faa_gz),
            "--fasta",
            str(fasta),
            "--antismash-gbk",
            str(as_gbk),
            "--out-prediction",
            str(out_pred),
            "--out-faa",
            str(out_faa),
            "--out-gbk",
            str(out_gbk),
        ]
    )
    assert rc == 0
    assert "ctg1_1" in out_pred.read_text()
    faa_text = out_faa.read_text()
    assert ">ctg1_1" in faa_text
    assert "*" not in faa_text

    records = list(SeqIO.parse(out_gbk, "genbank"))
    assert len(records) == 1
    cds_by_tag = {f.qualifiers["locus_tag"][0]: f for f in records[0].features if f.type == "CDS"}
    assert set(cds_by_tag) == {"ctg1_1", "ctg1_2"}
    assert cds_by_tag["ctg1_1"].qualifiers["translation"][0] == "MKKLLVVAA"
    # Overlapping antiSMASH transporter annotation should be preserved on ctg1_2
    assert "transporter" in cds_by_tag["ctg1_2"].qualifiers["product"][0].lower()
