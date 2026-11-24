"""Unit tests for bin/gbk_to_fasta.py."""

import gzip
from pathlib import Path

from Bio import SeqIO
import gbk_to_fasta

FIXTURES = Path(__file__).parents[1] / "data" / "combgc"
GENOMES = Path(__file__).parents[1] / "data" / "genomes"


def test_convert_plain_genbank_to_fasta_faa_gff(tmp_path):
    fa_out = tmp_path / "sample1.fna"
    faa_out = tmp_path / "sample1.faa"
    gff_out = tmp_path / "sample1.gff3"

    n_rec, n_cds = gbk_to_fasta.convert_genbank(FIXTURES / "sample1.gbk", fa_out, faa_out=faa_out, gff_out=gff_out)
    assert n_rec == 1
    assert n_cds >= 1

    records = list(SeqIO.parse(fa_out, "fasta"))
    assert len(records) == 1
    assert len(records[0].seq) > 0

    proteins = list(SeqIO.parse(faa_out, "fasta"))
    assert len(proteins) == n_cds
    assert "##gff-version 3" in gff_out.read_text()
    assert "\tCDS\t" in gff_out.read_text()


def test_convert_gzipped_genbank_to_gzipped_fasta(tmp_path):
    fa_gz = tmp_path / "ecoli.fna.gz"
    n_rec, _ = gbk_to_fasta.convert_genbank(GENOMES / "GCA_000005845_mini.gbff.gz", fa_gz)
    assert n_rec >= 1
    with gzip.open(fa_gz, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    assert len(records) == n_rec
    assert len(records[0].seq) > 1000


def test_normalize_gzipped_fasta(tmp_path):
    src_gz = tmp_path / "input.fasta.gz"
    with gzip.open(src_gz, "wt") as fh:
        fh.write(">ctg1 desc\nATGCGTACGT\n>empty\n\n")
    dst = tmp_path / "output.fna"
    rc = gbk_to_fasta.main([str(src_gz), str(dst), "--normalize-fasta"])
    assert rc == 0
    records = list(SeqIO.parse(dst, "fasta"))
    assert len(records) == 1
    assert records[0].id == "ctg1"
    assert str(records[0].seq) == "ATGCGTACGT"
