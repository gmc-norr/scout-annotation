import cyvcf2
import gzip
from pathlib import Path
import pytest
import re
import shutil
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
        "--notemp",
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
        "--notemp",
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
        "--notemp",
        "--cores",
        "1",
    ]

    wd = Path(Path(__file__).parent, "integration", "trio")

    return subprocess.run(args, cwd=wd), wd


@pytest.fixture(scope="session")
def snakemake_trio_config(snakemake_trio):
    config_file = Path(
        snakemake_trio[1], "results_trio/ceph1463/ceph1463.load_config.yaml"
    )
    assert config_file.exists()
    with open(config_file) as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="session")
def cli_trio():
    args = [
        "scout-annotation",
        "--use-apptainer",
        "--apptainer-prefix",
        "/storage/userdata/singularity_cache",
        "--apptainer-args",
        "--bind /storage",
        "trio",
        "--seq-type",
        "wes",
        "../data/NA12877.vcf",
        "../data/NA12878.vcf",
        "../data/NA12879.vcf",
        "../data/ceph1463_trio.ped",
    ]

    wd = Path(Path(__file__).parent, "integration", "cli_trio")
    return subprocess.run(args, cwd=wd), wd


def test_cli_trio(cli_trio):
    assert cli_trio[0].returncode == 0


@pytest.fixture(scope="session")
def cli_single(tmp_path_factory):
    vcf = Path("tests/integration/data/HD832_chr7_twist-solid-0.1.5-alpha.vcf")
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
        str(vcf.name),
    ]

    wd = tmp_path_factory.mktemp("cli_single")
    shutil.copy(str(vcf), wd)

    return subprocess.run(args, cwd=wd), wd


def test_cli_single(cli_single):
    assert cli_single[0].returncode == 0
    results_dir = Path(cli_single[1])
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
        "--bam-dir",
        "../batch_data/bam",
        "--notemp",
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


def test_cli_batch_load_config(cli_batch):
    config_path = Path(cli_batch[1], "cli_batch_results/HD832/HD832.load_config.yaml")
    bam_path = Path(cli_batch[1], "cli_batch_results/HD832/HD832.bam")
    with open(config_path) as f:
        load_config = yaml.safe_load(f)

    assert load_config["samples"][0]["alignment_path"] == "HD832.bam"
    assert bam_path.exists()


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
    vcf = Path(
        snakemake_trio[1], "results_trio/ceph1463/ceph1463.scout-annotated.vcf.gz"
    )
    assert vcf.exists()
    vcf_samples = cyvcf2.VCF(vcf).samples
    assert len(vcf_samples) == 3
    assert "NA12877" in vcf_samples
    assert "NA12878" in vcf_samples
    assert "NA12879" in vcf_samples


def test_trio_config_samples(snakemake_trio):
    config_file = Path(
        snakemake_trio[1], "results_trio/ceph1463/ceph1463.load_config.yaml"
    )
    assert config_file.exists()
    with open(config_file) as f:
        load_config = yaml.safe_load(f)
    assert "samples" in load_config
    assert len(load_config["samples"]) == 3
    load_config_sample_ids = [s["sample_id"] for s in load_config["samples"]]
    assert "NA12877" in load_config_sample_ids
    assert "NA12878" in load_config_sample_ids
    assert "NA12879" in load_config_sample_ids


def test_trio_peddy(snakemake_trio):
    results_dir = Path(snakemake_trio[1], "results_trio/ceph1463")
    peddy_het_check = results_dir / "ceph1463.peddy.het_check.csv"
    peddy_ped_check = results_dir / "ceph1463.peddy.ped_check.csv"
    peddy_sex_check = results_dir / "ceph1463.peddy.sex_check.csv"
    peddy_html = results_dir / "ceph1463.peddy.html"
    peddy_ped = results_dir / "ceph1463.peddy.ped"
    madeline2_pedigree = results_dir / "ceph1463.pedigree.svg"

    assert peddy_het_check.exists()
    assert peddy_ped_check.exists()
    assert peddy_sex_check.exists()
    assert peddy_html.exists()
    assert peddy_ped.exists()
    assert madeline2_pedigree.exists()


def test_trio_config_pedigree(snakemake_trio_config):
    assert "madeline" in snakemake_trio_config
    assert snakemake_trio_config["madeline"] == "ceph1463.pedigree.svg"
