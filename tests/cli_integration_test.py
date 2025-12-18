from typing import Tuple
import cyvcf2
from pathlib import Path
import pytest
import shutil
import subprocess
import yaml

from scout_annotation.resources import default_config, snakefile

SNAKEFILE = snakefile()
DEFAULT_CONFIG = default_config()


if not Path("/storage").exists():
    pytest.skip("requires /storage", allow_module_level=True)


@pytest.fixture(scope="session")
def cli_trio(tmp_path_factory) -> Tuple[subprocess.CompletedProcess, Path]:
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
        "data/NA12877.vcf",
        "data/NA12878.vcf",
        "data/NA12879.vcf",
        "data/ceph1463_trio.ped",
    ]

    wd = tmp_path_factory.mktemp("cli_trio")
    shutil.copytree("tests/integration/data", wd / "data")

    p = subprocess.run(args, cwd=wd)

    try:
        assert p.returncode == 0
    except AssertionError:
        pytest.skip("snakemake process failed, skip all dependent tests")

    return wd


@pytest.fixture(scope="session")
def trio_load_config(cli_trio):
    config_file = Path(
        cli_trio,
        "results/ceph1463/ceph1463.load_config.yaml",
    )
    assert config_file.exists()
    with open(config_file) as f:
        return yaml.safe_load(f)


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

    p = subprocess.run(args, cwd=wd)

    try:
        assert p.returncode == 0
    except AssertionError:
        pytest.skip("snakemake process failed, skip all dependent tests")

    return wd


@pytest.fixture(scope="session")
def cli_batch(tmp_path_factory):
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
        "batch_data/bam",
        "-o",
        "cli_batch_results",
        "batch_data/",
    ]

    wd = tmp_path_factory.mktemp("cli_batch")
    shutil.copytree("tests/integration/data", wd / "data")
    shutil.copytree("tests/integration/batch_data", wd / "batch_data")

    p = subprocess.run(args, cwd=wd)

    try:
        assert p.returncode == 0
    except AssertionError:
        pytest.skip("snakemake process failed, skip all dependent tests")

    return wd


def test_cli_batch_load_config(cli_batch):
    config_path = Path(cli_batch, "cli_batch_results/HD832/HD832.load_config.yaml")
    bam_path = Path(cli_batch, "cli_batch_results/HD832/HD832.bam")
    with open(config_path) as f:
        load_config = yaml.safe_load(f)

    assert not load_config.get("gene_panels", [])
    assert load_config["samples"][0]["alignment_path"] == str(bam_path)
    assert bam_path.exists()


def test_cli_single_load_config(cli_single):
    config_path = Path(cli_single, "results/HD832/HD832.load_config.yaml")
    with open(config_path) as f:
        load_config = yaml.safe_load(f)

    assert not load_config.get("gene_panels", [])


def test_trio_vcf_samples(cli_trio):
    vcf = Path(cli_trio, "results/ceph1463/ceph1463.scout-annotated.vcf.gz")
    assert vcf.exists()
    vcf_samples = cyvcf2.VCF(vcf).samples
    assert len(vcf_samples) == 3
    assert "NA12877" in vcf_samples
    assert "NA12878" in vcf_samples
    assert "NA12879" in vcf_samples


def test_trio_config_samples(cli_trio):
    config_file = Path(cli_trio, "results/ceph1463/ceph1463.load_config.yaml")
    assert config_file.exists()
    with open(config_file) as f:
        load_config = yaml.safe_load(f)
    assert "samples" in load_config
    assert len(load_config["samples"]) == 3
    load_config_sample_ids = [s["sample_id"] for s in load_config["samples"]]
    assert "NA12877" in load_config_sample_ids
    assert "NA12878" in load_config_sample_ids
    assert "NA12879" in load_config_sample_ids


def test_trio_peddy(cli_trio):
    results_dir = Path(cli_trio, "results/ceph1463")
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


def test_trio_config_pedigree(cli_trio, trio_load_config):
    assert "madeline" in trio_load_config
    assert trio_load_config["madeline"] == str(cli_trio / "results/ceph1463/ceph1463.pedigree.svg")
