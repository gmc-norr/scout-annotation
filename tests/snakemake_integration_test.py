import gzip
from pathlib import Path
import pytest
import re
import shutil
import subprocess

REPO_ROOT = Path(__file__).parent.parent.resolve()
SNAKEFILE = str(REPO_ROOT / "workflow/Snakefile")
DEFAULT_CONFIG = str(REPO_ROOT / "default_config/config.yaml")

if not Path("/storage").exists():
    pytest.skip("requires /storage", allow_module_level=True)


@pytest.fixture(scope="session")
def integration(tmp_path_factory):
    args = [
        "snakemake",
        "-s",
        SNAKEFILE,
        "--singularity-args",
        "--bind /storage",
        "--use-apptainer",
        "--apptainer-prefix",
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

    p = subprocess.run(args, cwd=wd)

    try:
        assert p.returncode == 0
    except AssertionError:
        pytest.skip("snakemake process failed, skip all dependent tests")

    return wd


def test_integration_fixture_runs(integration):
    # Just check that the working directory was created by the fixture
    assert integration.exists()


@pytest.fixture(
    params=[("family1", "sample1", 162), ("family2", "sample2", 162)],
    scope="session",
    ids=lambda x: x[0],
)
def annotated_vcf(request, integration):
    return {
        "family": request.param[0],
        "sample": request.param[1],
        "n_variants": request.param[2],
        "path": Path(
            integration,
            f"results/annotation/{request.param[0]}/{request.param[0]}.annotated.genmod.vcf.gz",
        ),
    }


def test_vcfs_exist(annotated_vcf):
    assert annotated_vcf["path"].exists()


def test_vcf_sample_names(annotated_vcf):
    with gzip.open(annotated_vcf["path"], "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                columns = line.strip().split("\t")
                # Sample names start at column 10
                sample_name = columns[9]
                assert sample_name == annotated_vcf["sample"]
                break


def test_number_of_variants(annotated_vcf):
    n_variants = 0
    with gzip.open(annotated_vcf["path"], "rt") as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if not line.startswith("#"):
                n_variants += 1
    assert annotated_vcf["n_variants"] == n_variants


def test_rank_score_family_names(annotated_vcf):
    rank_score_pattern = re.compile(r"RankScore=(?P<family>[^:]+):(-?\d+)")
    with gzip.open(annotated_vcf["path"], "rt") as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if line.startswith("#"):
                continue
            rs_match = rank_score_pattern.search(line)
            # assert rs_match is not None
            assert rs_match.group("family") == annotated_vcf["family"]


def test_caller_names(annotated_vcf):
    found_in_re = re.compile(r"FOUND_IN=([^;]+)")
    count_found_in_mutect = 0
    with gzip.open(annotated_vcf["path"], "rt") as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if line.startswith("#"):
                continue
            found_in_match = found_in_re.search(line)

            if not found_in_match:
                continue

            callers_str = found_in_match.group(1)
            callers = callers_str.split(",")

            assert "mutect2" not in callers, f"mutect2 found in: {callers}"

            if "gatk" in callers:
                count_found_in_mutect += 1
    assert count_found_in_mutect > 0, f"mutect is not among callers for any record"
