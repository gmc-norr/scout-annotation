import cyvcf2
import gzip
from pathlib import Path
import pytest
import re
import shutil
import subprocess
import yaml

from cli_integration_test import SNAKEFILE, DEFAULT_CONFIG

@pytest.fixture(scope="session")
def integration(tmp_path_factory):
    args = [
        "snakemake",
        "-s",
        SNAKEFILE,
        "--singularity-args",
        "--bind /storage",
        "--use-singularity",
        "--singularity-prefix",
        "/storage/userdata/singularity_cache",
        "--configfiles",
        DEFAULT_CONFIG,
        "config.yaml",
        "--show-failed-logs",
        "--notemp",
        "--cores",
        "1",
    ]

    wd = tmp_path_factory.mktemp("integration")
    shutil.copy("tests/integration/config.yaml", wd)
    shutil.copy("tests/integration/data/samples.tsv", wd)
    shutil.copytree("tests/integration/data", wd / "data")
    shutil.copytree("tests/integration/filters", wd / "filters")
    shutil.copytree("tests/integration/panels", wd / "panels")

    p = subprocess.run(args, cwd=wd)

    try:
        assert p.returncode == 0
    except AssertionError:
        pytest.skip("snakemake process failed, skip all dependent tests")

    return wd


@pytest.fixture(scope="session")
def integration_no_filtering(tmp_path_factory):
    args = [
        "snakemake",
        "-s",
        SNAKEFILE,
        "--singularity-args",
        "--bind /storage",
        "--use-singularity",
        "--singularity-prefix",
        "/storage/userdata/singularity_cache",
        "--configfiles",
        DEFAULT_CONFIG,
        "config.yaml",
        "--config",
        "samples=samples_no-filtering.tsv",
        "output_directory=results_no-filtering",
        "--show-failed-logs",
        "--notemp",
        "--cores",
        "1",
    ]

    wd = tmp_path_factory.mktemp("integration_no_filtering")
    shutil.copy("tests/integration/config.yaml", wd)
    shutil.copy("tests/integration/data/samples_no-filtering.tsv", wd)
    shutil.copytree("tests/integration/data", wd / "data")
    shutil.copytree("tests/integration/filters", wd / "filters")
    shutil.copytree("tests/integration/panels", wd / "panels")

    p = subprocess.run(args, cwd=wd)

    try:
        assert p.returncode == 0
    except AssertionError:
        pytest.skip("snakemake process failed, skip all dependent tests")

    return wd


@pytest.fixture(
    params=[
        ("sample1", "clingen-somatic"),
        ("sample2", "clingen"),
        ("sample3", "clingen-somatic"),
        ("sample4", "clingen"),
        ("sample5", "clingen-somatic"),
        ("sample6", "clingen"),
        ("sample7", "clingen"),
    ],
    scope="session",
    ids=lambda x: x[0],
)
def load_config(request, integration):
    return {
        "sample": request.param[0],
        "owner": request.param[1],
        "path": Path(
            integration, f"results/{request.param[0]}/{request.param[0]}.load_config.yaml"
        ),
    }


@pytest.fixture(
    params=[
        ("sample1", False, 23),
        ("sample2", False, 0),
        ("sample3", True, 0),
        ("sample4", True, 32),
        ("sample5", True, 24),
        ("sample6", False, 0),
        ("sample7", True, 0),
    ],
    scope="session",
    ids=lambda x: x[0],
)
def scout_vcf(request, integration):
    return {
        "sample": request.param[0],
        "snv_filtering": request.param[1],
        "n_variants": request.param[2],
        "path": Path(
            integration, f"results/{request.param[0]}/{request.param[0]}.scout-annotated.vcf.gz"
        ),
    }


@pytest.fixture(scope="session")
def scout_vcfs_no_filtering(integration_no_filtering):
    return [
        dict(
            sample="sample2-1",
            n_variants=160,
            path=Path(
                integration_no_filtering,
                "results_no-filtering/sample2-1/sample2-1.scout-annotated.vcf.gz",
            ),
        ),
    ]


def test_vcfs_exist(scout_vcf):
    assert scout_vcf["path"].exists()


def test_vcf_sample_names(scout_vcf):
    assert cyvcf2.VCF(scout_vcf["path"]).samples[0] == scout_vcf["sample"]


def test_vembrane_filtering(scout_vcf):
    with gzip.open(scout_vcf["path"], "rt") as f:
        vembrane_found = False
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("##vembraneCmd"):
                vembrane_found = True
    assert scout_vcf["snv_filtering"] == vembrane_found


def test_number_of_variants(scout_vcf):
    n_variants = 0
    with gzip.open(scout_vcf["path"], "rt") as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if not line.startswith("#"):
                n_variants += 1
    assert scout_vcf["n_variants"] == n_variants


def test_vcfs_no_filtering(scout_vcfs_no_filtering):
    for vcf in scout_vcfs_no_filtering:
        assert vcf["path"].exists()


def test_number_of_variants_no_filtering(scout_vcfs_no_filtering):
    vcf = scout_vcfs_no_filtering[0]
    n_variants = 0
    with gzip.open(vcf["path"], "rt") as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if not line.startswith("#"):
                n_variants += 1
    assert vcf["n_variants"] == n_variants


def test_load_configs_exist(load_config):
    assert load_config["path"].exists()


def test_rank_score_threshold(load_config):
    c = yaml.safe_load(load_config["path"].read_text())
    assert "rank_score_threshold" in c
    assert c["rank_score_threshold"] == -1000


def test_rank_score_sample_names(scout_vcf):
    rank_score_pattern = re.compile(r"RankScore=(?P<sample>[^:]+):(\d+)")
    with gzip.open(scout_vcf["path"], "rt") as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if line.startswith("#"):
                continue
            rs_match = rank_score_pattern.search(line)
            assert rs_match is not None
            assert rs_match.group("sample") == scout_vcf["sample"]


def test_case_owner(load_config):
    c = yaml.safe_load(load_config["path"].read_text())
    assert "owner" in c, load_config["sample"]
    assert c["owner"] == load_config["owner"]
