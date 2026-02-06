import click
import cyvcf2
import pathlib
import subprocess
import sys

from scout_annotation.resources import default_config, snakefile
from scout_annotation.panels import get_panels
from scout_annotation.samples import write_samples


def get_ped_sex(s):
    match s:
        case 1 | "1":
            return "male"
        case 2 | "2":
            return "female"
        case _:
            return "unknown"


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
    nargs=3,
)
@click.argument(
    "ped",
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
    type=click.Path(path_type=pathlib.Path, dir_okay=True, file_okay=False),
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
    "-b",
    "--bam-files",
    help="path to alignment BAM files, one for each sample",
    type=click.Path(
        path_type=pathlib.Path, dir_okay=False, file_okay=True, resolve_path=True
    ),
    default=[None, None, None],
    nargs=3,
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
def trio(
    config,
    vcf,
    ped,
    out_dir,
    profile,
    dryrun,
    notemp,
    track,
    samples_dir,
    seq_type,
    bam_files,
    panel,
    snv_filter,
    owner,
):
    """Annotate a trio of samples

    Supply three VCF files, one each for the mother, father and child."""
    config.logger.info(f"annotating trio")
    config.logger.info(f"fetching sample information from PED file: {ped}")

    family = None
    samples = []
    sexes = []

    with open(ped) as f:
        for line in f:
            try:
                f, s, _, _, sex, _ = line.strip().split()
            except ValueError:
                config.logger.error(f"invalid PED file: {ped}")
                raise click.Abort()
            if family is not None and f != family:
                config.logger.error(f"more than one family in PED file: {ped}")
                raise click.Abort()
            family = f
            samples.append(s)
            sexes.append(get_ped_sex(sex))
    if len(samples) != 3:
        config.logger.error(f"expected 3 samples in PED, found {len(samples)}")
        raise click.Abort()

    for s, sex, vcf_path in zip(samples, sexes, vcf):
        config.logger.debug(f"{s} ({sex}): {vcf_path}")

    gene_panels = get_panels()
    for p in panel:
        if p not in gene_panels:
            config.logger.error(f"panel not found: {p}")
            raise click.Abort()

    sample_vcf = {}
    for vcf_path in vcf:
        vcf_samples = cyvcf2.VCF(vcf_path).samples
        if len(vcf_samples) != 1:
            config.logger.error(
                f"expected 1 sample in vcf, found {len(vcf_samples)}: {vcf_path}"
            )
            raise click.Abort()
        sample_vcf[vcf_samples[0]] = vcf_path

    sample_rows = []
    for sname, ssex, sbam in zip(samples, sexes, bam_files):
        config.logger.debug(f"generating sample row for {sname}")
        try:
            svcf = sample_vcf[sname]
        except KeyError:
            config.logger.error(f"sample {sname} not found in any VCF")
            raise click.Abort()
        sample_rows.append(
            {
                "sample": sname,
                "family": family,
                "owner": owner,
                "sex": ssex,
                "type": seq_type,
                "track": track,
                "filtering": snv_filter if snv_filter is not None else "",
                "panels": ",".join(panel),
                "vcf": svcf,
                "ped": ped,
                "bam": sbam if sbam is not None else "",
            }
        )

    config.logger.debug(f"generated sample rows for {len(sample_rows)} samples")

    samples_file = write_samples(sample_rows, samples_dir)

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

    args.extend(
        [
            "--config",
            f"samples={samples_file}",
        ]
    )

    if out_dir is not None:
        args.append(f"output_directory={out_dir}")

    config.logger.info("starting snakemake")

    p = subprocess.run(args)

    sys.exit(p.returncode)
