import importlib.resources
import pathlib


def get_panels():
    package_dir = importlib.resources.files(__name__.split(".")[0])
    panel_dir = package_dir.joinpath("data/panels")
    with importlib.resources.as_file(panel_dir) as pd:
        panel_paths = pathlib.Path(pd).glob("*.tsv")

    panel_dict = {}
    for p in panel_paths:
        with open(p) as f:
            n_lines = sum(1 for line in f if not line.startswith("#"))
            panel_dict[p.stem] = (n_lines, p)
    return panel_dict
