import gzip
from pathlib import Path
import pytest
import subprocess

@pytest.fixture(scope="session")
def integration():
    args = [
        "snakemake",
        "-s", "../../workflow/Snakefile",
        "--singularity-args",
        "--bind /storage",
        "--use-singularity",
        "--singularity-prefix",
        "/storage/userdata/singularity_cache",
        "--configfiles",
        "../../config/config.yaml",
        "config.yaml",
        "--cores",
        "1",
    ]

    subprocess.run(args, cwd=Path(Path(__file__).parent, "integration"))

@pytest.fixture(scope="session")
def integration_no_filtering():
    args = [    
        "snakemake",
        "-s", "../../workflow/Snakefile",
        "--singularity-args",
        "--bind /storage",
        "--use-singularity",
        "--singularity-prefix",
        "/storage/userdata/singularity_cache",
        "--configfiles",
        "../../config/config.yaml",
        "config.yaml",
        "--config",
        "samples=samples_no-filtering.tsv",
        "output_directory=results_no-filtering",
        "--cores",
        "1",
    ]

    subprocess.run(args, cwd=Path(Path(__file__).parent, "integration"))

@pytest.fixture(scope="session")
def scout_vcfs(integration):
    return [
        dict(
            sample="sample1",
            empty=False,
            path=Path("tests/integration/results/sample1/sample1.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample2",
            empty=True,
            path=Path("tests/integration/results/sample2/sample2.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample3",
            empty=True,
            path=Path("tests/integration/results/sample3/sample3.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample4",
            empty=False,
            path=Path("tests/integration/results/sample4/sample4.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample5",
            empty=False,
            path=Path("tests/integration/results/sample5/sample5.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample6",
            empty=True,
            path=Path("tests/integration/results/sample6/sample6.scout-annotated.vcf.gz"),
        ),
        dict(
            sample="sample7",
            empty=True,
            path=Path("tests/integration/results/sample7/sample7.scout-annotated.vcf.gz"),
        ),
    ]

@pytest.fixture(scope="session")
def scout_vcfs_no_filtering(integration_no_filtering):
    return [
        dict(
            sample="sample1",
            empty=False,
            n_variants=162,
            path=Path("tests/integration/results/sample1/sample1.scout-annotated.vcf.gz"),
        ),
    ]

def test_vcfs_exist(scout_vcfs):
    for vcf in scout_vcfs:
        assert vcf["path"].exists()

def test_number_of_variants(scout_vcfs):
    for vcf in scout_vcfs:
        n_variants = 0
        with gzip.open(vcf["path"], "rt") as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                if not line.startswith("#"):
                    n_variants += 1
        assert vcf["empty"] == (n_variants == 0), vcf["sample"]

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
