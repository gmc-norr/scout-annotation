from pathlib import Path
import pytest

from scout_annotation.path import WildcardPath


@pytest.fixture(scope="module")
def twist_solid_results(tmp_path_factory):
    samples = ["sample1", "sample2", "sample3"]
    types = ["T", "N", "T"]
    root = tmp_path_factory.mktemp("twist_solid")
    results = root / "results" / "dna"
    results.mkdir(parents=True)
    bam_dna = root / "bam_dna"
    bam_dna.mkdir()

    alt_bam_dir = bam_dna / "mutect2_indel_bam"
    alt_bam_dir.mkdir()

    for sample, sample_type in zip(samples, types):
        (bam_dna / f"{sample}_{sample_type}.bam").touch()
        (alt_bam_dir / f"{sample}_{sample_type}.bam").touch()
        vcf = results / f"{sample}_{sample_type}" / "vcf"
        vcf.mkdir(parents=True)
        (vcf / f"{sample}_{sample_type}.filtered.exons_only.vcf.gz").touch()

    return root


@pytest.mark.parametrize(
    "input, expected",
    [
        (
            "results/{sample}_{type}.txt",
            ["sample", "type"],
        ),
        (
            "results/this_is_a_prefix_{sample}_{type}.txt",
            ["sample", "type"],
        ),
        (
            "results/{sample}_{type}/this_is_a_prefix_{sample}_{type}.txt",
            ["sample", "type"],
        ),
    ],
)
def test_extract_wildcards(input, expected, tmp_path_factory):
    wcp = WildcardPath(Path("/path/to/analysis", input))
    assert wcp.wildcards == expected


def test_expand_bam(twist_solid_results):
    wcp = WildcardPath(twist_solid_results / "bam_dna" / "{sample}_{type}.bam")
    bam_paths = wcp.expand()
    assert len(bam_paths) == 3
    wildcards = [x[1] for x in bam_paths]
    assert {"sample": "sample1", "type": "T"} in wildcards
    assert {"sample": "sample2", "type": "N"} in wildcards
    assert {"sample": "sample3", "type": "T"} in wildcards


def test_expand_vcf(twist_solid_results):
    wcp = WildcardPath(
        twist_solid_results
        / "results"
        / "dna"
        / "{sample}_{type}"
        / "vcf"
        / "{sample}_{type}.filtered.exons_only.vcf.gz"
    )
    vcf_paths = wcp.expand()
    assert len(vcf_paths) == 3
    wildcards = [x[1] for x in vcf_paths]
    assert {"sample": "sample1", "type": "T"} in wildcards
    assert {"sample": "sample2", "type": "N"} in wildcards
    assert {"sample": "sample3", "type": "T"} in wildcards


def test_expand_missing(twist_solid_results):
    wcp = WildcardPath(twist_solid_results / "reports" / "{sample}_{type}.cnv.html")
    report_paths = wcp.expand()
    assert len(report_paths) == 0


def test_no_wildcards(twist_solid_results):
    with pytest.raises(ValueError):
        WildcardPath(twist_solid_results / "bam_dna" / "sample1_T.bam")


@pytest.mark.parametrize(
    "input, expected",
    [
        ("/path/to/some/wildcard/prefix_{sample}.txt", "/path/to/some/wildcard"),
        ("/path/to/some/wildcard/{sample}_{type}.txt", "/path/to/some/wildcard"),
    ],
)
def test_determined_root(input, expected):
    wcp = WildcardPath(input)
    assert wcp.determined_root() == Path(expected)


@pytest.mark.parametrize(
    "input",
    ["results/dna/{type}_{sample}/vcf/{sample}_{type}.filtered.exons_only.vcf.gz"],
)
def test_mismatching_wildcards(input, twist_solid_results):
    wcp = WildcardPath(twist_solid_results / input)
    assert len(wcp.expand()) == 0
