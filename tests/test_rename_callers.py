from pathlib import Path
import pysam
import pytest
import sys

SCRIPT_DIR = Path(__file__).parent / "../workflow/scripts"
assert SCRIPT_DIR.exists()
sys.path.insert(0, str(SCRIPT_DIR.resolve()))

from rename_callers import rename_callers


def write_vcf(path, header, records):
    with pysam.VariantFile(path, "w", header=header) as out_vcf:
        for rec in records:
            out_vcf.write(rec)


@pytest.fixture
def vcf_header():
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.add_line("##contig=<ID=chr1,length=1000000>")
    header.add_line("##contig=<ID=chr2,length=1000000>")
    header.add_line("##contig=<ID=chr3,length=1000000>")
    header.add_line(
        '##INFO=<ID=FOUND_IN,Number=.,Type=String,Description="Individual caller support">'
    )
    header.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
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
        found_in=["gatk_mutect2", "vardict"],
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

        rec.info["FOUND_IN"] = found_in

        rec.samples["S1"]["GT"] = gt
        rec.samples["S1"]["AD"] = ad
        return rec

    return _make_record


def test_rename_callers(tmp_path, vcf_header, make_record):
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"

    rec1 = make_record(found_in=["gatk_mutect2", "vardict"])
    rec2 = make_record(found_in=["vardict"])
    rec3 = make_record(found_in=["gatk_mutect2"])
    rec4 = make_record(found_in=["other", "vardict", "gatk_mutect2"])

    write_vcf(
        input_vcf,
        vcf_header,
        [rec1, rec2, rec3, rec4],
    )

    rename_callers(str(input_vcf), str(output_vcf), {"gatk_mutect2": "mutect"})

    with pysam.VariantFile(output_vcf) as vcf:
        out = list(vcf)

    assert len(out) == 4

    renamed1 = out[0]
    renamed2 = out[1]
    renamed3 = out[2]
    renamed4 = out[3]

    assert renamed1.info["FOUND_IN"] == ("mutect", "vardict")
    assert renamed2.info["FOUND_IN"] == ("vardict",)
    assert renamed3.info["FOUND_IN"] == ("mutect",)
    assert renamed4.info["FOUND_IN"] == ("other", "vardict", "mutect")

    rename_callers(
        str(input_vcf),
        str(output_vcf),
        {"gatk_mutect2": "mutect", "vardict": "VarDict"},
    )

    with pysam.VariantFile(output_vcf) as vcf:
        out = list(vcf)

    renamed1 = out[0]
    renamed2 = out[1]
    renamed3 = out[2]
    renamed4 = out[3]

    assert len(out) == 4

    assert renamed1.info["FOUND_IN"] == ("mutect", "VarDict")
    assert renamed2.info["FOUND_IN"] == ("VarDict",)
    assert renamed3.info["FOUND_IN"] == ("mutect",)
    assert renamed4.info["FOUND_IN"] == ("other", "VarDict", "mutect")
