import importlib.resources
import pathlib

def _package_dir(path=None):
    package_path = importlib.resources.files("scout_annotation")._paths[0]

    if path is None:
        return package_path

    return package_path / path


def default_config() -> pathlib.Path:
    return _package_dir("default_config/config.yaml")


def default_resources() -> pathlib.Path:
    return _package_dir("default_config/resources.yaml")


def snakefile() -> pathlib.Path:
    return _package_dir("workflow/Snakefile")

def snakefile_download() -> pathlib.Path:
    return _package_dir("workflow/Snakefile_download")
