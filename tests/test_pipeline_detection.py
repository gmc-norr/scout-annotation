import json
import pytest
import yaml

from scout_annotation import pipeline_utils


@pytest.fixture(
    params=[("Twist_Solid", "v0.23.0", True), ("poppy", "v1.0.1", False)],
    ids=lambda x: f"{x[0]}_{x[1]}",
)
def snakemake_dir(tmp_path_factory, request):
    pipeline = request.param[0]
    version = request.param[1]
    has_general_report = request.param[2]
    d = tmp_path_factory.mktemp("snakemake_pipeline")
    (d / ".snakemake").mkdir()
    version_d = d / "results" / "versions"
    software_d = version_d / f"software_{pipeline}"
    software_d.mkdir(parents=True)

    # general report json
    if has_general_report:
        reports_d = d / "general_json_report"
        reports_d.mkdir()

        with open(reports_d / "sample1.general.json", "w") as f:
            json.dump({"pipeline": {"name": pipeline, "version": version}}, f)

    # make some other yaml files
    (version_d / f"config_{pipeline}.yaml").touch()
    (software_d / "softwares_mqc_versions.yaml").touch()

    # make a yaml file that resembles the one we're looking for
    decoy_yaml = software_d / "decoy.yaml"
    with open(decoy_yaml, "w") as f:
        yaml.safe_dump(
            {
                "not_a_pipeline": "not_a_version",
            },
            f,
        )

    version_f = software_d / f"{pipeline}__{version}_mqc_versions.yaml"
    with open(version_f, "w") as f:
        version_dict = {}
        version_dict[pipeline] = version
        yaml.safe_dump(version_dict, f)

    return d, (pipeline, version)


@pytest.fixture(
    params=[("nf-core/raredisease", "2.4.0")],
    ids=lambda x: f"{x[0]}_{x[1]}",
)
def nextflow_dir(tmp_path_factory, request):
    pipeline = request.param[0]
    version = request.param[1]
    d = tmp_path_factory.mktemp("nextflow_pipeline")
    (d / ".nextflow").mkdir()
    logfile = d / ".nextflow.log"
    with open(logfile, "w") as f:
        f.writelines(
            [
                "------------------------------------------------------\n",
                "                                        ,--./,-.\n",
                "        ___     __   __   __   ___     /,-._.--~'\n",
                "  |\\ | |__  __ /  ` /  \\ |__) |__         }  {\n",
                "  | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,\n",
                "                                        `._,._,'\n",
                "  nf-core/raredisease 2.4.0\n",
                "------------------------------------------------------\n",
            ]
        )
    return d, (pipeline, version)


@pytest.fixture()
def invalid_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("invalid_pipeline")
    return d


def test_snakemake_detection(snakemake_dir):
    d = pipeline_utils.detect_pipeline(snakemake_dir[0])
    assert d == snakemake_dir[1]


def test_nextflow_detection(nextflow_dir):
    d = pipeline_utils.detect_pipeline(nextflow_dir[0])
    assert d == nextflow_dir[1]


def test_invalid_pipeline_detection(invalid_dir):
    with pytest.raises(ValueError):
        pipeline_utils.detect_pipeline(invalid_dir)


@pytest.mark.parametrize(
    "input,expected", [("v0.23.0", True), ("0.1.0", True), ("version", False), ("1.0", False)]
)
def test_is_version(input, expected):
    assert pipeline_utils.is_version(input) == expected
