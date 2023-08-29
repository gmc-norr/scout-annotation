import pathlib


def get_panels():
    # TODO: fetch panel directory from config
    panel_paths = pathlib.Path(__file__).parent.glob("panels/*.tsv")
    panel_dict = {}
    for p in panel_paths:
        with open(p) as f:
            n_lines = sum(1 for line in f if not line.startswith("#"))
            panel_dict[p.stem] = (n_lines, p)
    return panel_dict
