from pathlib import Path

from scout_annotation.resources import default_config, default_resources, snakefile


def test_default_config_path():
    config_path = default_config()
    assert isinstance(config_path, Path)
    assert config_path.name == "config.yaml"
    assert config_path.exists()


def test_default_resource_path():
    resource_path = default_resources()
    assert isinstance(resource_path, Path)
    assert resource_path.name == "resources.yaml"
    assert resource_path.exists()


def test_snakefile_path():
    snakefile_path = snakefile()
    assert isinstance(snakefile_path, Path)
    assert snakefile_path.name == "Snakefile"
    assert snakefile_path.exists()
