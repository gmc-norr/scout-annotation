import click

from scout_annotation.cli.batch import batch
from scout_annotation.cli.panels import panels
from scout_annotation.cli.single import single

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


class Config:
    def __init__(self, config, resources):
        self.config = config
        self.resources = resources


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("-c", "--config", help="config file", default="config/config.yaml")
@click.option("-r", "--resources", help="resources file", default=None)
@click.version_option(prog_name="scout_annotation")
@click.pass_context
def cli(ctx, config, resources):
    ctx.obj = Config(config, resources)

cli.add_command(batch)
cli.add_command(panels)
cli.add_command(single)
