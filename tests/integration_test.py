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
        "../config.yaml",
        "--show-failed-logs",
        "--cores",
        "1",
    ]

    wd = Path(Path(__file__).parent, "integration", "test")

    return subprocess.run(args, cwd=wd), wd


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
        "../config.yaml",
        "--config",
        "samples=samples_no-filtering.tsv",
        "output_directory=results_no-filtering",
        "--show-failed-logs",
        "--cores",
        "1",
    ]

    wd = Path(Path(__file__).parent, "integration", "no_filtering")

    return subprocess.run(args, cwd=wd), wd


@pytest.fixture(scope="session")
def snakemake_trio():
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
        "../config.yaml",
        "--config",
        "samples=samples_trio.tsv",
        "output_directory=results_trio",
        "--show-failed-logs",
        "--cores",
        "1",
    ]

    wd = Path(Path(__file__).parent, "integration", "trio")

    return subprocess.run(args, cwd=wd), wd


@pytest.fixture(scope="session")
def cli_single_from_outside(tmp_path_factory):
    args = [
        "python",
        "-m",
        "scout_annotation",
        "--use-apptainer",
        "--apptainer-args",
        "--bind /storage",
        "--apptainer-prefix",
        "/storage/userdata/singularity_cache",
        "--cores",
        "1",
        "single",
        "--snv-filter",
        "rare_disease",
        "-o",
        "cli_single_results_from_outside",
        Path("tests/integration/data/HD832_chr7_twist-solid-0.1.5-alpha.vcf").resolve(),
    ]

    wd = tmp_path_factory.mktemp("cli_workdir")

    return subprocess.run(args, cwd=wd), wd


def test_cli_single_from_outside(cli_single_from_outside):
    assert cli_single_from_outside[0].returncode == 0
    results_dir = Path(cli_single_from_outside[1])
    assert results_dir.exists()
    assert results_dir.is_dir()


@pytest.fixture(scope="session")
def cli_single():
    args = [
        "python",
        "-m",
        "scout_annotation",
        "--use-apptainer",
        "--apptainer-args",
        "--bind /storage",
        "--cores",
        "1",
        "single",
        "-o",
        "cli_single_results",
        "../data/HD832_chr7_twist-solid-0.1.5-alpha.vcf",
    ]

    wd = Path(Path(__file__).parent, "integration", "cli_single")

    return subprocess.run(args, cwd=wd), wd


def test_cli_single(cli_single):
    assert cli_single[0].returncode == 0
    results_dir = cli_single[1]
    assert results_dir.exists()
    assert results_dir.is_dir()


@pytest.fixture(scope="session")
def cli_batch():
    args = [
        "python",
        "-m",
        "scout_annotation",
        "--use-apptainer",
        "--apptainer-args",
        "--bind /storage",
        "--cores",
        "1",
        "batch",
        "-o",
        "cli_batch_results",
        "../batch_data/",
    ]

    wd = Path(Path(__file__).parent, "integration", "cli_batch")

    return subprocess.run(args, cwd=wd), wd


def test_cli_batch(cli_batch):
    assert cli_batch[0].returncode == 0
    results_dir = cli_batch[1]
    assert results_dir.exists()
    assert results_dir.is_dir()


@pytest.fixture(scope="session")
def load_configs(integration):
    return [
        dict(
            sample="sample1",
            owner="clingen-somatic",
            path=Path(integration[1], "results/sample1/sample1.load_config.yaml"),
        ),
        dict(
            sample="sample2",
            owner="clingen",
            path=Path(integration[1], "results/sample2/sample2.load_config.yaml"),
        ),
        dict(
            sample="sample3",
            owner="clingen-somatic",
            path=Path(integration[1], "results/sample3/sample3.load_config.yaml"),
        ),
        dict(
            sample="sample4",
            owner="clingen",
            path=Path(integration[1], "results/sample4/sample4.load_config.yaml"),
        ),
        dict(
            sample="sample5",
            owner="clingen-somatic",
            path=Path(integration[1], "results/sample5/sample5.load_config.yaml"),
        ),
        dict(
            sample="sample6",
            owner="clingen",
            path=Path(integration[1], "results/sample6/sample6.load_config.yaml"),
        ),
        dict(
            sample="sample7",
            owner="clingen",
            path=Path(integration[1], "results/sample7/sample7.load_config.yaml"),
        ),
    ]


@pytest.fixture(scope="session")
def scout_vcfs(integration):
    return [
        dict(
            sample="sample1",
            snv_filtering=False,
            n_variants=23,
            path=Path(integration[1], "results/sample1/sample1.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample2",
            snv_filtering=False,
            n_variants=0,
            path=Path(integration[1], "results/sample2/sample2.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample3",
            snv_filtering=True,
            n_variants=0,
            path=Path(integration[1], "results/sample3/sample3.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample4",
            snv_filtering=True,
            n_variants=32,
            path=Path(integration[1], "results/sample4/sample4.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample5",
            snv_filtering=True,
            n_variants=24,
            path=Path(integration[1], "results/sample5/sample5.scout-annotated.vcf.gz"),
        ),
        dict(
            # This sample has SNV filtering enabled, but it should be empty
            # after the panel filtering.
            sample="sample6",
            snv_filtering=False,
            n_variants=0,
            path=Path(integration[1], "results/sample6/sample6.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample7",
            snv_filtering=True,
            n_variants=0,
            path=Path(integration[1], "results/sample7/sample7.scout-annotated.vcf.gz"),
        ),
    ]


@pytest.fixture(scope="session")
def scout_vcfs_no_filtering(integration_no_filtering):
    return [
        dict(
            sample="sample2-1",
            n_variants=160,
            path=Path(
                integration_no_filtering[1],
                "results_no-filtering/sample2-1/sample2-1.scout-annotated.vcf.gz",
            ),
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


def test_trio_exit_code(snakemake_trio):
    assert snakemake_trio[0].returncode == 0


def test_trio_vcf_samples(snakemake_trio):
    vcf = Path(snakemake_trio[1], "results_trio/f1234/f1234.scout-annotated.vcf.gz")
    assert vcf.exists()
    vcf_samples = cyvcf2.VCF(vcf).samples
    assert len(vcf_samples) == 3
    assert "mother" in vcf_samples
    assert "father" in vcf_samples
    assert "child" in vcf_samples


def test_trio_config_samples(snakemake_trio):
    config_file = Path(snakemake_trio[1], "results_trio/f1234/f1234.load_config.yaml")
    assert config_file.exists()
    with open(config_file) as f:
        load_config = yaml.safe_load(f)
    assert "samples" in load_config
    assert len(load_config["samples"]) == 3
    load_config_sample_ids = [s["sample_id"] for s in load_config["samples"]]
    assert "mother" in load_config_sample_ids
    assert "father" in load_config_sample_ids
    assert "child" in load_config_sample_ids
