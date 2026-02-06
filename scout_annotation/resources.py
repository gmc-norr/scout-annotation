import importlib.resources
import pathlib
from typing import Dict
import yaml

def _package_dir(path=None):
    package_path = importlib.resources.files("scout_annotation")

    if path is None:
        return package_path

    return package_path / path


def default_config() -> pathlib.Path:
    return _package_dir("default_config/config.yaml")


def pipeline_files() -> Dict[str, Dict[str, str]]:
    definition_file = _package_dir("pipelines/snakemake.yaml")
    with open(definition_file) as f:
        return yaml.safe_load(f)


def snakefile() -> pathlib.Path:
    return _package_dir("workflow/Snakefile")
