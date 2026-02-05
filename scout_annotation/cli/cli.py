import click
import logging

from scout_annotation.cli.batch import batch
from scout_annotation.cli.panels import panels
from scout_annotation.cli.single import single
from scout_annotation.cli.trio import trio

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


class Config:
    def __init__(
        self,
        config,
        cores,
        use_apptainer,
        apptainer_args,
        apptainer_prefix,
        loglevel,
    ):
        self.config = config
        self.cores = cores
        self.use_apptainer = use_apptainer
        self.apptainer_args = apptainer_args
        self.apptainer_prefix = apptainer_prefix

        self.logger = logging.getLogger()
        self.logger.setLevel(loglevel)
        ch = logging.StreamHandler()
        logformatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        ch.setFormatter(logformatter)
        self.logger.addHandler(ch)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("-c", "--config", help="config file used for overwriting defaults")
@click.option(
    "--cores", help="number of cores available for snakemake", type=click.INT, default=1
)
@click.option(
    "--use-apptainer",
    "--use-singularity",
    help="use apptainer as executor",
    is_flag=True,
)
@click.option("--apptainer-args", "--singularity-args", help="arguments for apptainer")
@click.option(
    "--apptainer-prefix",
    "--singularity-prefix",
    help="path to cached apptainer containers",
)
@click.option(
    "--loglevel",
    help="set logging level",
    default="WARNING",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
)
@click.version_option(prog_name="scout_annotation")
@click.pass_context
def cli(
    ctx,
    config,
    cores,
    use_apptainer,
    apptainer_args,
    apptainer_prefix,
    loglevel,
):
    ctx.obj = Config(
        config,
        cores,
        use_apptainer,
        apptainer_args,
        apptainer_prefix,
        loglevel,
    )
    ctx.obj.logger.info("starting annotation")


cli.add_command(batch)
cli.add_command(panels)
cli.add_command(single)
cli.add_command(trio)
