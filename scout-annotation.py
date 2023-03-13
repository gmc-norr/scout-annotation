#!/usr/bin/env python

import click
import hashlib
import pathlib
import subprocess
import sys
import time
from typing import Dict, List

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


class Config:
    def __init__(self, config, resources):
        self.config = config
        self.resources = resources


def get_panels():
    panel_paths = pathlib.Path(__file__).parent.glob("panels/*.tsv")
    panel_dict = {}
    for p in panel_paths:
        with open(p) as f:
            n_lines = sum(1 for line in f if not line.startswith("#"))
            panel_dict[p.stem] = (n_lines, p)
    return panel_dict


def write_samples(samples: List[Dict], directory: pathlib.Path):
    if not directory.is_dir():
        directory.mkdir()

    sample_md5 = hashlib.md5(
        ",".join(d["sample"] for d in samples).encode("utf8")
    ).hexdigest()[:8]
    filename = pathlib.Path(
        directory, f"{time.strftime('%Y%m%d')}-{sample_md5}_samples.txt"
    )

    # Header
    cols = ("sample", "sex", "type", "track", "vcf", "ped", "panels")
    with open(filename, "w") as f:
        print("\t".join(cols), file=f)
        for s in samples:
            print("\t".join(str(s[c]) for c in cols), file=f)

    return filename


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("-c", "--config", help="config file", default="config/config.yaml")
@click.option("-r", "--resources", help="resources file", default=None)
@click.pass_context
def cli(ctx, config, resources):
    ctx.obj = Config(config, resources)


@cli.command(
    epilog="""
    Gene panels that are added to a case should exist in Scout. To see
    panels available to the pipeline, run `scout-annotation.py panels`."""
)
@click.argument(
    "vcf",
    type=click.Path(
        path_type=pathlib.Path,
        exists=True,
        dir_okay=False,
        file_okay=True,
        resolve_path=True,
    ),
)
@click.option("-n", "--name", help="sample name")
@click.option(
    "--profile",
    help="snakemake profile to use",
    type=click.Path(
        path_type=pathlib.Path,
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
    ),
)
@click.option(
    "-t",
    "--track",
    help="sample track",
    type=click.Choice(["cancer", "rare_disease"]),
    default="rare_disease",
)
@click.option(
    "--samples-dir",
    help="path to directory where to store sample files",
    type=click.Path(
        path_type=pathlib.Path, dir_okay=True, file_okay=False, resolve_path=True
    ),
    default=pathlib.Path("./sample_files"),
)
@click.option(
    "-s",
    "--seq-type",
    help="type of sequencing",
    type=click.Choice(["panel", "wes", "wgs"]),
    default="panel",
)
@click.option(
    "--sex",
    help="sex of the patient",
    type=click.Choice(["unknown", "male", "female"]),
    default="unknown",
)
@click.option(
    "-b",
    "--bam-file",
    help="path to alignment BAM file",
    type=click.Path(
        path_type=pathlib.Path, dir_okay=False, file_okay=True, resolve_path=True
    ),
)
@click.option(
    "-p",
    "--panel",
    help="gene panel to filter by, can be passed multiple times",
    multiple=True,
)
@click.pass_obj
def single(config, vcf, profile, name, track, samples_dir, seq_type, sex, bam_file, panel):
    """Annotate a single sample."""

    if name is None:
        name = vcf.stem.split("_")[0]

    sample = {
        "sample": name,
        "sex": sex,
        "type": seq_type,
        "track": track,
        "vcf": vcf,
        "bam": bam_file,
        "panels": ",".join(panel),
        "ped": "",
    }

    gene_panels = get_panels()
    for p in panel:
        if p not in gene_panels:
            print(f"error: panel not found: {p}", file=sys.stderr)
            exit(1)
    print(f"annotate {vcf} with panels: {panel}")

    samples_file = write_samples([sample], samples_dir)

    print(config.config, config.resources)

    args = [
        "snakemake",
        "-s",
        pathlib.Path(pathlib.Path(__file__).parent, "workflow/Snakefile"),
        "--configfile",
        config.config,
        "--config",
        f"samples={samples_file}",
    ]
    if config.resources is not None:
        args.append(f"resources={config.resources}")
    if profile is not None:
        args.append("--profile")
        args.append(profile)
    subprocess.Popen(args).communicate()


@cli.command()
@click.pass_context
def panels(ctx):
    """List available gene panels"""
    print("#name\tn_genes\tpath")
    for panel_name, (n_genes, p) in get_panels().items():
        print(f"{panel_name}\t{n_genes}\t{p}")


if __name__ == "__main__":
    cli()
