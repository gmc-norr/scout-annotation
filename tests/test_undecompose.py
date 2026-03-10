from pathlib import Path
import pysam
import pytest
import sys

SCRIPT_DIR = Path(__file__).parent / "../workflow/scripts"
assert SCRIPT_DIR.exists()
sys.path.insert(0, str(SCRIPT_DIR.resolve()))

from undecompose_vcf import parse_old_clumped, undecompose, expected_decomposed_snvs


@pytest.fixture
def vcf_header():
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.add_line("##contig=<ID=chr1,length=1000000>")
    header.add_line(
        '##INFO=<ID=OLD_CLUMPED,Number=1,Type=String,Description="Original clumped variant">'
    )
    header.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
    header.add_line(
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations">'
    )
    header.add_line(
        '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Calling methods">'
    )
    header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line(
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">'
    )
    header.add_sample("S1")
    return header


@pytest.fixture
def make_record(vcf_header):
    def _make_record(
        *,
        chrom="chr1",
        pos=100,
        ref="A",
        alt="C",
        old_clumped=None,
        dp=10,
        gt=(0, 1),
        ad=(5, 5),
    ):
        rec = vcf_header.new_record(
            contig=chrom,
            start=pos - 1,
            stop=(pos - 1) + len(ref),
            alleles=(ref, alt),
        )
        rec.filter.add("PASS")
        rec.info["DP"] = dp

        if old_clumped is not None:
            rec.info["OLD_CLUMPED"] = old_clumped

        rec.samples["S1"]["GT"] = gt
        rec.samples["S1"]["AD"] = ad
        return rec

    return _make_record


def write_vcf(path, header, records):
    with pysam.VariantFile(path, "w", header=header) as out_vcf:
        for rec in records:
            out_vcf.write(rec)


def test_parse_old_clumped_valid_and_invalid_values():
    assert parse_old_clumped("chr1:1111111:AA/CC") == ("chr1", 1111111, "AA", "CC")
    assert parse_old_clumped("chrX:42:A/T") == ("chrX", 42, "A", "T")

    assert parse_old_clumped(None) is None
    assert parse_old_clumped("chr1-1111111-AA/CC") is None
    assert parse_old_clumped("chr1:abc:AA/CC") is None
    assert parse_old_clumped("chr1:100:A>C") is None


def test_undecompose_merges_same_old_clumped_and_preserves_first_record_fields(
    tmp_path, vcf_header, make_record
):
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"

    old = "chr1:100:AAA/CAC"

    # Same OLD_CLUMPED group: should merge into one AAA->CAC record
    rec1 = make_record(
        pos=100,
        ref="A",
        alt="C",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )
    rec2 = make_record(
        pos=102,
        ref="A",
        alt="C",
        old_clumped=old,
        dp=999,  # deliberately different: merged record should still use rec1 values
        gt=(1, 1),
        ad=(0, 20),
    )

    # Second decompose, right after the first one
    old = "chr1:103:TA/GC"
    rec3 = make_record(
        pos=103,
        ref="T",
        alt="G",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )
    rec4 = make_record(
        pos=104,
        ref="A",
        alt="C",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )

    # No OLD_CLUMPED: should pass straight through unchanged
    passthrough = make_record(
        pos=105,
        ref="G",
        alt="T",
        old_clumped=None,
        dp=25,
        gt=(1, 1),
        ad=(0, 25),
    )
    old = "chr1:200:TATT/GTTC"
    # Third decomposed, right after normal SNV
    rec5 = make_record(
        pos=200,
        ref="T",
        alt="G",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )
    rec6 = make_record(
        pos=201,
        ref="A",
        alt="T",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )
    rec7 = make_record(
        pos=203,
        ref="T",
        alt="C",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )

    old = "chr3:300:ACGT/GGGC"
    # Fourth decomposed, missing SNVs
    rec8 = make_record(
        pos=300,
        ref="A",
        alt="G",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )
    rec9 = make_record(
        pos=303,
        ref="T",
        alt="C",
        old_clumped=old,
        dp=17,
        gt=(0, 1),
        ad=(9, 8),
    )

    write_vcf(
        input_vcf,
        vcf_header,
        [rec1, rec2, rec3, rec4, passthrough, rec5, rec6, rec7, rec8, rec9],
    )

    undecompose(str(input_vcf), str(output_vcf))

    with pysam.VariantFile(output_vcf) as vcf:
        out = list(vcf)

    assert len(out) == 6

    merged1 = out[0]
    merged2 = out[1]
    merged3 = out[3]
    unmerged1 = out[4]
    unmerged2 = out[5]
    assert merged1.contig == "chr1"
    assert merged1.pos == 100
    assert tuple(merged1.alleles) == ("AAA", "CAC")

    # INFO copied from first record, except OLD_CLUMPED and CSQ
    assert merged1.info["DP"] == 17
    assert "OLD_CLUMPED" not in merged1.info

    # FORMAT copied from first record
    assert merged1.samples["S1"]["GT"] == (0, 1)
    assert tuple(merged1.samples["S1"]["AD"]) == (9, 8)

    # Test second merged
    assert merged2.contig == "chr1"
    assert merged2.pos == 103
    assert tuple(merged2.alleles) == ("TA", "GC")

    # Test third merged
    assert merged3.contig == "chr1"
    assert merged3.pos == 200
    assert tuple(merged3.alleles) == ("TATT", "GTTC")

    unchanged = out[2]
    assert unchanged.pos == 105
    assert tuple(unchanged.alleles) == ("G", "T")
    assert unchanged.info["DP"] == 25
    assert unchanged.samples["S1"]["GT"] == (1, 1)

    assert unmerged1.pos == 300
    assert tuple(unmerged1.alleles) == ("A", "G")
    assert unmerged1.info["OLD_CLUMPED"] == "chr3:300:ACGT/GGGC"

    assert unmerged2.pos == 303
    assert tuple(unmerged2.alleles) == ("T", "C")
    assert unmerged2.info["OLD_CLUMPED"] == "chr3:300:ACGT/GGGC"


def test_undecompose_flushes_when_old_clumped_group_changes(
    tmp_path, vcf_header, make_record
):
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"

    old1 = "chr1:100:AA/CC"
    old2 = "chr1:200:GG/TT"

    records = [
        make_record(pos=100, ref="A", alt="C", old_clumped=old1),
        make_record(pos=101, ref="A", alt="C", old_clumped=old1),
        make_record(pos=200, ref="G", alt="T", old_clumped=old2),
        make_record(pos=201, ref="G", alt="T", old_clumped=old2),
    ]

    write_vcf(input_vcf, vcf_header, records)
    undecompose(str(input_vcf), str(output_vcf))

    with pysam.VariantFile(output_vcf) as vcf:
        out = list(vcf)

    assert len(out) == 2
    assert out[0].pos == 100
    assert tuple(out[0].alleles) == ("AA", "CC")
    assert out[1].pos == 200
    assert tuple(out[1].alleles) == ("GG", "TT")


def test_expected_decomposed_snvs_returns_only_changed_positions():
    old = ("chr1", 1111, "ACGT", "GGGC")

    expected = expected_decomposed_snvs(old)

    assert expected == [
        ("chr1", 1111, "A", "G"),
        ("chr1", 1112, "C", "G"),
        ("chr1", 1114, "T", "C"),
    ]
