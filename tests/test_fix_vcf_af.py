from pathlib import Path
import pysam
import pytest
import sys

SCRIPT_DIR = Path(__file__).parent / "../workflow/scripts"
print(SCRIPT_DIR)
assert SCRIPT_DIR.exists()
sys.path.insert(0, str(SCRIPT_DIR.resolve()))

from fix_vcf_af import fix_vcf_af, write_variants


@pytest.fixture(scope="function")
def hd832_vcf():
    vcf_filename = "tests/integration/data/HD832_chr7_twist-solid-0.1.5-alpha.vcf"
    return pysam.VariantFile(vcf_filename)


@pytest.fixture(scope="function")
def af_ao_vcf():
    vcf_filename = "tests/integration/data/af_ao.vcf"
    return pysam.VariantFile(vcf_filename)


def test_fix_vcf_af(hd832_vcf, af_ao_vcf):
    variants = fix_vcf_af(hd832_vcf)
    assert len(variants) == 162
    for v in variants:
        assert "AO" in v.format
        assert "AF" in v.format
        assert "DP" in v.format

    variants = fix_vcf_af(af_ao_vcf)
    assert len(variants) == 3
    for v in variants:
        assert "AO" in v.format
        assert "AF" in v.format
        assert "DP" in v.format

    print(variants[0])

    assert variants[0].samples[0]["AO"] == 248
    assert variants[0].samples[0]["DP"] == 1175
    assert variants[0].samples[0]["AF"] == pytest.approx(0.2110638298)

    assert variants[1].samples[0]["AO"] is None
    assert variants[1].samples[0]["DP"] == 1175
    assert variants[1].samples[0]["AF"] is None

    assert variants[2].samples[0]["AO"] is 248
    assert variants[2].samples[0]["DP"] is None
    assert variants[2].samples[0]["AF"] is None


def test_write_vcf(hd832_vcf, tmpdir):
    variants = fix_vcf_af(hd832_vcf)
    out_vcf = tmpdir / "out.vcf"

    write_variants(out_vcf, variants, hd832_vcf.header)
    assert out_vcf.exists

    n_variants = 0
    with open(out_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_variants += 1

    assert n_variants == 162
