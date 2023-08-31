import cyvcf2
import gzip
from pathlib import Path
import pytest
import re
import subprocess
import yaml

SNAKEFILE = Path("./scout_annotation/workflow/Snakefile").resolve()
CONFIG = Path("./scout_annotation/default_config/config.yaml").resolve()

@pytest.fixture(scope="session")
def integration():
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
        CONFIG,
        "config.yaml",
        "--show-failed-logs",
        "--cores",
        "1",
    ]

    subprocess.run(args, cwd=Path(Path(__file__).parent, "integration"))

@pytest.fixture(scope="session")
def integration_no_filtering():
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
        CONFIG,
        "config.yaml",
        "--config",
        "samples=samples_no-filtering.tsv",
        "output_directory=results_no-filtering",
        "--show-failed-logs",
        "--cores",
        "1",
    ]

    subprocess.run(args, cwd=Path(Path(__file__).parent, "integration"))

@pytest.fixture(scope="session")
def load_configs(integration):
    return [
        dict(
            sample="sample1",
            owner="clingen-somatic",
            path=Path("tests/integration/results/sample1/sample1.load_config.yaml"),
        ),
        dict(
            sample="sample2",
            owner="clingen",
            path=Path("tests/integration/results/sample2/sample2.load_config.yaml"),
        ),
        dict(
            sample="sample3",
            owner="clingen-somatic",
            path=Path("tests/integration/results/sample3/sample3.load_config.yaml"),
        ),
        dict(
            sample="sample4",
            owner="clingen",
            path=Path("tests/integration/results/sample4/sample4.load_config.yaml"),
        ),
        dict(
            sample="sample5",
            owner="clingen-somatic",
            path=Path("tests/integration/results/sample5/sample5.load_config.yaml"),
        ),
        dict(
            sample="sample6",
            owner="clingen",
            path=Path("tests/integration/results/sample6/sample6.load_config.yaml"),
        ),
        dict(
            sample="sample7",
            owner="clingen",
            path=Path("tests/integration/results/sample7/sample7.load_config.yaml"),
        ),
    ]

@pytest.fixture(scope="session")
def scout_vcfs(integration):
    return [
        dict(
            sample="sample1",
            snv_filtering=False,
            n_variants=23,
            path=Path("tests/integration/results/sample1/sample1.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample2",
            snv_filtering=False,
            n_variants=0,
            path=Path("tests/integration/results/sample2/sample2.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample3",
            snv_filtering=True,
            n_variants=0,
            path=Path("tests/integration/results/sample3/sample3.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample4",
            snv_filtering=True,
            n_variants=32,
            path=Path("tests/integration/results/sample4/sample4.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample5",
            snv_filtering=True,
            n_variants=24,
            path=Path("tests/integration/results/sample5/sample5.scout-annotated.vcf.gz"),
        ),
        dict(
            # This sample has SNV filtering enabled, but it should be empty
            # after the panel filtering.
            sample="sample6",
            snv_filtering=False,
            n_variants=0,
            path=Path("tests/integration/results/sample6/sample6.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample7",
            snv_filtering=True,
            n_variants=0,
            path=Path("tests/integration/results/sample7/sample7.scout-annotated.vcf.gz"),
        ),
    ]

@pytest.fixture(scope="session")
def scout_vcfs_no_filtering(integration_no_filtering):
    return [
        dict(
            sample="sample2-1",
            n_variants=160,
            path=Path("tests/integration/results_no-filtering/sample2-1/sample2-1.scout-annotated.vcf.gz"),
        ),
    ]

def test_vcfs_exist(scout_vcfs):
    for vcf in scout_vcfs:
        assert vcf["path"].exists()

def test_vcf_sample_names(scout_vcfs):
    for vcf in scout_vcfs:
        assert cyvcf2.VCF(vcf["path"]).samples[0] == vcf["sample"], vcf["sample"]

def test_vembrane_filtering(scout_vcfs):
    for vcf in scout_vcfs:
        with gzip.open(vcf["path"], "rt") as f:
            vembrane_found = False
            for line in f:
                if not line.startswith("#"):
                    break
                if line.startswith("##vembraneCmd"):
                    vembrane_found = True
        assert vcf["snv_filtering"] == vembrane_found

def test_number_of_variants(scout_vcfs):
    for vcf in scout_vcfs:
        n_variants = 0
        with gzip.open(vcf["path"], "rt") as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                if not line.startswith("#"):
                    n_variants += 1
        assert vcf["n_variants"] == n_variants, vcf["sample"]

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

def test_load_configs_exist(load_configs):
    for config in load_configs:
        assert config["path"].exists()

def test_rank_score_threshold(load_configs):
    for config in load_configs:
        c = yaml.safe_load(config["path"].read_text())
        assert "rank_score_threshold" in c, config["sample"]
        assert c["rank_score_threshold"] == -1000, config["sample"]

def test_rank_score_sample_names(scout_vcfs):
    rank_score_pattern = re.compile(r"RankScore=(?P<sample>[^:]+):(\d+)")
    for vcf in scout_vcfs:
        with gzip.open(vcf["path"], "rt") as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                if line.startswith("#"):
                    continue
                rs_match = rank_score_pattern.search(line)
                assert rs_match is not None, vcf["sample"]
                assert rs_match.group("sample") == vcf["sample"], vcf["sample"]

def test_case_owner(load_configs):
    for config in load_configs:
        c = yaml.safe_load(config["path"].read_text())
        assert "owner" in c, config["sample"]
        assert c["owner"] == config["owner"], config["sample"]
