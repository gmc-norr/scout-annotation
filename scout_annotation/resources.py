import importlib.resources
import pathlib

def snakefile() -> pathlib.Path:
    package_dir = importlib.resources.files(__name__.split(".")[0])
    snakefile_path = package_dir.joinpath("workflow/Snakefile")
    with importlib.resources.as_file(snakefile_path) as sfp:
        return pathlib.Path(sfp)

