# import logging
import subprocess
import sys
import click
from pathlib import Path
from pydantic import ValidationError

from scout_annotation import pipeline_definition, pipeline_utils
import scout_annotation.parsers as parsers
from scout_annotation.samples import write_samples
from scout_annotation import resources


@click.command()
@click.argument(
    "input_dir",
    metavar="dir",
    type=click.Path(
        path_type=Path,
        exists=True,
        dir_okay=True,
        file_okay=False,
    ),
)
@click.option(
    "--pipeline-definitions",
    "pipeline_definitions",
    help="path to configuration of the autodetection",
    type=click.Path(path_type=Path, dir_okay=False, file_okay=True, exists=True),
    default=Path("~/.config/scout-annotation/pipeline_definitions.yaml"),
)
@click.option(
    "-n",
    "--pipeline-name",
    "pipeline_name",
    help="name of pipeline that generated the results",
)
@click.option(
    "-v",
    "--pipeline-version",
    "pipeline_version",
    help="version of pipeline that generated the results",
)
@click.option(
    "--out-dir",
    "-o",
    help="directory to write output files to",
    type=click.Path(path_type=Path, dir_okay=True, file_okay=False),
)
@click.option(
    "--profile",
    help="snakemake profile to use",
    type=click.Path(
        path_type=Path,
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
    ),
)
@click.option(
    "--dryrun",
    help="perform a snakemake dryrun, no results are produced",
    is_flag=True,
)
@click.option(
    "--notemp",
    help="do not output files marked as temporary",
    is_flag=True,
)
@click.option(
    "-s",
    "--seq-type",
    help="type of sequencing",
    type=click.Choice(["panel", "wes", "wgs"]),
    default="panel",
)
@click.option(
    "-p",
    "--panel",
    help="gene panel to filter by, can be passed multiple times",
    multiple=True,
)
@click.option(
    "--snv-filter",
    help="SNV filter to apply",
)
@click.option(
    "--owner",
    help="Scout institiute ID that should own the case",
    default="clingen",
    show_default=True,
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
    type=click.Path(path_type=Path, dir_okay=True, file_okay=False, resolve_path=True),
    default=Path("./sample_files"),
)
@click.pass_obj
def auto(
    config,
    out_dir: Path,
    dryrun: bool,
    notemp: bool,
    seq_type: str,
    panel: list[str],
    owner: str,
    profile: Path,
    snv_filter: str,
    pipeline_definitions: Path,
    pipeline_name: str,
    pipeline_version: str,
    input_dir: Path,
    track: str,
    samples_dir: Path,
):
    logger = config.logger.getChild("auto")
    if pipeline_name is None or pipeline_version is None:
        logger.info("trying to autodetect pipeline from output")
        try:
            pipeline_name, pipeline_version = pipeline_utils.detect_pipeline(input_dir)
        except ValueError as e:
            logger.error(e)
            raise click.Abort
        logger.info(f"detected pipeline: name={pipeline_name}, version={pipeline_version}")

    logger.info(f"annotating results from {pipeline_name} {pipeline_version}")

    try:
        pipeline_defs = pipeline_definition.read(pipeline_definitions)
    except ValidationError as e:
        for err in e.errors():
            config.logger.error(
                f"validation error: {'.'.join(map(str, err['loc']))}: {err['msg']}"
            )
        sys.exit(1)

    pdef = pipeline_defs.find(pipeline_name, pipeline_version)
    if pdef is None:
        logger.error("unsupported pipeline")
        raise click.Abort

    try:
        samples = pdef.resolve_paths(input_dir)
    except ValueError as e:
        config.logger.error(e)
        sys.exit(1)

    if len(samples) == 0:
        config.logger.error("no samples found, aborting")
        sys.exit(1)

    samples_file = write_samples(samples.values(), samples_dir)

    args = [
        "snakemake",
        "-s",
        resources.snakefile(),
        "--rerun-incomplete",
        "--verbose",
        "--debug-dag",
        # "--cores",
        # str(config.cores),
        "--configfiles",
        resources.default_config(),
    ]

    if config.config is not None:
        args.append(config.config)
    if config.use_apptainer:
        args.append("--use-singularity")
    if config.apptainer_args is not None:
        args.append("--singularity-args")
        args.append(config.apptainer_args)
    if config.apptainer_prefix is not None:
        args.append("--singularity-prefix")
        args.append(config.apptainer_prefix)

    if profile is not None:
        args.append("--profile")
        args.append(profile)
    if dryrun:
        args.append("--dryrun")
    if notemp:
        args.append("--notemp")

    args.extend(
        [
            "--config",
            f"samples={samples_file}",
        ]
    )

    if out_dir is not None:
        args.append(f"output_directory={out_dir}")

    config.logger.debug(f"command line: {" ".join(map(str, args))}")
    p = subprocess.run(args)

    sys.exit(p.returncode)
