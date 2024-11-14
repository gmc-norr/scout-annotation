import click
import pathlib
import subprocess
import sys

from scout_annotation.resources import default_config, default_resources, snakefile, snakefile_download


@click.command()
@click.option(
    "--dryrun",
    help="perform a snakemake dryrun, no results are produced",
    is_flag=True,
)
@click.option(
    "--verbose",
    help="print verbose messages from snakemake",
    is_flag=True,
)
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
@click.pass_obj
def download(
    config,
    profile,
    verbose,
    dryrun,
):

    args = [
        "snakemake",
        "-s",
        snakefile_download(),
        "--rerun-incomplete",
        "--cores",
        str(config.cores),
        "--configfiles",
        default_config(),
    ]

    if config.config is not None:
        args.append(config.config)

    if profile is not None:
        args.append("--profile")
        args.append(profile)
    if dryrun:
        args.append("--dryrun")
    if verbose:
        args.append("-p")
        args.append("--verbose")
    
    args.extend(
        [
            "--config",
            f"mode=download",
        ]
    )

    p = subprocess.run(args)

    sys.exit(p.returncode)
