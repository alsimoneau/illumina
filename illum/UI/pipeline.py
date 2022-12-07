#!/usr/bin/env python3

import click

import illum


# pipelines
@click.group()
def pipe():
    """Execution pipelines."""
    pass


# pipe pts
@pipe.command(name="pts")
@click.argument("inventory", type=click.Path(exists=True))
@click.option(
    "--topography",
    default="TOPOGRAPHY",
    type=click.Path(exists=True),
    show_default=True,
    help="Topography input file or folder name.",
)
def pipe_pts(inventory, topography):
    """VIIRS pipeline."""
    illum.domain()
    illum.warp(topography, "topography")
    illum.input_pts(inventory)
    click.echo("Done.")


# pipe viirs
@pipe.command(name="viirs")
@click.argument("inventory", type=click.Path(exists=True))
@click.option(
    "--viirs",
    default="VIIRS",
    type=click.Path(exists=True),
    show_default=True,
    help="VIIRS input file or folder name.",
)
@click.option(
    "--topography",
    default="TOPOGRAPHY",
    type=click.Path(exists=True),
    show_default=True,
    help="Topography input file or folder name.",
)
def pipe_viirs(inventory, viirs, topography):
    """VIIRS pipeline."""
    illum.domain()
    illum.warp(viirs, "viirs")
    illum.warp(topography, "topography")
    illum.warp("hydropolys.zip", "water")
    illum.input_viirs(inventory)
    click.echo("Done.")
