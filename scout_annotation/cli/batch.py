import click
import pathlib
import subprocess
import sys

from scout_annotation.resources import default_config, default_resources, snakefile
from scout_annotation.panels import get_panels
from scout_annotation.samples import write_samples


@click.command()
@click.argument(
    "vcf-dir", type=click.Path(path_type=pathlib.Path, dir_okay=True, file_okay=False)
)
@click.option(
    "--out-dir",
    "-o",
    help="directory to write output files to",
    type=click.Path(path_type=pathlib.Path, dir_okay=True, file_okay=False)
)
@click.option(
    "--bam-dir",
    help="directory to look for BAM files in (default: vcf-dir)",
    type=click.Path(path_type=pathlib.Path, dir_okay=True, file_okay=False),
)
@click.option(
    "--ped-dir",
    help="directory to look for PED files in (default: vcf-dir)",
    type=click.Path(path_type=pathlib.Path, dir_okay=True, file_okay=False),
)
@click.option(
    "--sep",
    help="character separating sample name from the rest of the file name",
    default="_",
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
@click.pass_obj
def batch(
    config,
    vcf_dir,
    out_dir,
    bam_dir,
    ped_dir,
    sep,
    track,
    samples_dir,
    profile,
    dryrun,
    notemp,
    seq_type,
    panel,
    snv_filter,
    owner,
):
    """Annotate a batch of samples."""
    if bam_dir is None:
        bam_dir = vcf_dir

    if ped_dir is None:
        ped_dir = vcf_dir

    gene_panels = get_panels()
    for p in panel:
        if p not in gene_panels:
            print(f"error: panel not found: {p}", file=sys.stderr)
            exit(1)

    vcf_files = list(vcf_dir.glob("*.vcf"))
    vcf_files.extend(vcf_dir.glob("*.vcf.gz"))
    sample_names = set()
    samples = []

    for vf in vcf_files:
        try:
            sample_split = vf.name.split(sep)
        except ValueError as ve:
            print(f"error: {ve}", file=sys.stderr)
            exit(1)
        if len(sample_split) < 2:
            print(
                f'error: separator "{sep}" not found in filename: {vf.name}',
                file=sys.stderr,
            )
            exit(1)
        sample_name = sample_split[0]
        if sample_name in sample_names:
            print(
                f"error: duplicated sample name found, check input files: {sample_name}",
                file=sys.stderr,
            )
            exit(1)
        sample_names.add(sample_name)

        bam_filename = list(bam_dir.glob(f"{sample_name}*.bam"))
        if len(bam_filename) == 0:
            bam_filename = ""
        elif len(bam_filename) > 1:
            print(
                f"error: found more than one possible bam file for {sample_name}",
                file=sys.stderr,
            )
            exit(1)
        else:
            bam_filename = bam_filename[0]

        ped_filename = list(bam_dir.glob(f"{sample_name}*.ped"))
        if len(ped_filename) == 0:
            ped_filename = ""
        elif len(ped_filename) > 1:
            print(
                f"error: found more than one possible ped file for {sample_name}",
                file=sys.stderr,
            )
            exit(1)
        else:
            ped_filename = ped_filename[0]

        sample = {
            "sample": sample_name,
            "owner": owner,
            "sex": "unknown",
            "type": seq_type,
            "track": track,
            "vcf": vf,
            "bam": bam_filename,
            "panels": ",".join(panel),
            "ped": ped_filename,
        }

        if snv_filter is not None:
            sample["filtering"] = snv_filter

        samples.append(sample)

    samples_file = write_samples(samples, samples_dir)

    args = [
        "snakemake",
        "-s",
        snakefile(),
        "--rerun-incomplete",
        "--configfiles",
        default_config(),
        "--config",
        f"samples={samples_file}",
    ]
    if config.config is not None:
        args.insert(6, config.config)
    if config.resources is not None:
        args.append(f"resources={config.resources}")
    else:
        args.append(f"resources={default_resources()}")
    if out_dir is not None:
        args.append(f"output_directory={out_dir}")
    if profile is not None:
        args.append("--profile")
        args.append(profile)
    if dryrun:
        args.append("--dryrun")
    if notemp:
        args.append("--notemp")

    subprocess.Popen(args).communicate()

