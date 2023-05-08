from pathlib import Path
import pysam
import pytest
import sys

SCRIPT_DIR = Path(__file__).parent / "../workflow/scripts"
print(SCRIPT_DIR)
assert SCRIPT_DIR.exists()
sys.path.insert(0, str(SCRIPT_DIR.resolve()))

from fix_vcf_af import fix_vcf_af, write_variants

@pytest.fixture(scope="session")
def vcf():
    return Path("tests/integration/data/HD832_chr7_twist-solid-0.1.5-alpha.vcf")


def test_fix_vcf_af(vcf):
    variantfile = pysam.VariantFile(vcf)
    variants = fix_vcf_af(variantfile)

    assert len(variants) == 162
    for v in variants:
        assert "AO" in v.format
        assert "AF" in v.format
        assert "DP" in v.format


def test_write_vcf(vcf, tmpdir):
    variantfile = pysam.VariantFile(vcf)
    variants = fix_vcf_af(variantfile)
    out_vcf = tmpdir / "out.vcf"

    write_variants(out_vcf, variants, variantfile.header)
    assert out_vcf.exists
    
    n_variants = 0
    with open(out_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_variants += 1

    assert n_variants == 162
