import click
import pathlib
import subprocess
import sys

from scout_annotation.resources import default_config, default_resources, snakefile
from scout_annotation.panels import get_panels
from scout_annotation.samples import write_samples


@click.command()
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
@click.option(
    "--out-dir",
    "-o",
    help="directory to write output files to",
    type=click.Path(path_type=pathlib.Path, dir_okay=True, file_okay=False)
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
@click.pass_obj
def single(
    config,
    vcf,
    out_dir,
    profile,
    dryrun,
    notemp,
    name,
    track,
    samples_dir,
    seq_type,
    sex,
    bam_file,
    panel,
    snv_filter,
    owner,
):
    """Annotate a single sample."""

    if name is None:
        name = vcf.stem.split("_")[0]

    sample = {
        "sample": name,
        "owner": owner,
        "sex": sex,
        "type": seq_type,
        "track": track,
        "vcf": vcf,
        "bam": bam_file if bam_file is not None else "",
        "panels": ",".join(panel),
        "ped": "",
    }

    if snv_filter is not None:
        sample["filtering"] = snv_filter

    gene_panels = get_panels()
    for p in panel:
        if p not in gene_panels:
            print(f"error: panel not found: {p}", file=sys.stderr)
            exit(1)

    samples_file = write_samples([sample], samples_dir)

    args = [
        "snakemake",
        "-s",
        snakefile(),
        "--rerun-incomplete",
        "--cores",
        str(config.cores),
        "--configfiles",
        default_config(),
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

    args.extend([
        "--config",
        f"samples={samples_file}",
    ])

    if config.resources is not None:
        args.append(f"resources={config.resources}")
    else:
        args.append(f"resources={default_resources()}")
    if out_dir is not None:
        args.append(f"output_directory={out_dir}")

    subprocess.Popen(args).communicate()
