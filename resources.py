import importlib.resources
import pathlib

def _package_dir(path=None):
    package_path = importlib.resources.files("scout_annotation")

    if path is None:
        return package_path

    return package_path / path


def default_config() -> pathlib.Path:
    return _package_dir("default_config/config.yaml")


def snakefile() -> pathlib.Path:
    return _package_dir("workflow/Snakefile")
