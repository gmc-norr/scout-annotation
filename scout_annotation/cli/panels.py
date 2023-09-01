import click

from scout_annotation.panels import get_panels


@click.command()
def panels():
    """List available gene panels"""
    print("#name\tn_genes\tpath")
    for panel_name, (n_genes, p) in get_panels().items():
        print(f"{panel_name}\t{n_genes}\t{p}")
