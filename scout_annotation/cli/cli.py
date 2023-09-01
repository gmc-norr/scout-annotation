import click

from scout_annotation.cli.batch import batch
from scout_annotation.cli.panels import panels
from scout_annotation.cli.single import single

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


class Config:
    def __init__(self, config, resources, cores, use_apptainer, apptainer_args, apptainer_prefix):
        self.config = config
        self.resources = resources
        self.cores = cores
        self.use_apptainer = use_apptainer
        self.apptainer_args = apptainer_args
        self.apptainer_prefix = apptainer_prefix


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("-c", "--config", help="config file used for overwriting defaults")
@click.option("-r", "--resources", help="resources file for overwriting defaults")
@click.option("--cores", help="number of cores available for snakemake", type=click.INT, default=1)
@click.option(
    "--use-apptainer", "--use-singularity",
    help="use apptainer as executor",
    is_flag=True
)
@click.option(
    "--apptainer-args", "--singularity-args",
    help="arguments for apptainer"
)
@click.option(
    "--apptainer-prefix", "--singularity-prefix",
    help="path to store cached apptainer containers"
)
@click.version_option(prog_name="scout_annotation")
@click.pass_context
def cli(ctx, config, resources, cores, use_apptainer, apptainer_args, apptainer_prefix):
    ctx.obj = Config(config, resources, cores, use_apptainer, apptainer_args, apptainer_prefix)

cli.add_command(batch)
cli.add_command(panels)
cli.add_command(single)
