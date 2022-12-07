#!/usr/bin/env python3

import importlib
import os

import click

import illum


# main
@click.group()
@click.version_option(illum.__version__, prog_name="Illumina model")
def main():
    r"""Illumina model.

    See 'illum COMMAND --help' to read about specific subcommand.
    The user's guide is available at:
    http://dome.obsand.org:2080/wiki/pmwiki.php/Prof/IlluminaGuidev22

    Below is detailled the typical use case:

    1. init -- Prep the experiment folder.

    2. domain -- Define the simulation domain once the desired values
    are set in `domain_params.in`.

    3. warp -- Process the satellite imagery.

    4a. inputs -- Prepare the execution folder once the desired values
    are set in `inputs_params.in`.

    4b. alternate -- (Optionnal) Prepare conversion scenarios to be executed.

    At this point, the `Inputs` folder can be moved to a compute cluster.

    5. batches -- Prepare the model to run in parallel.

    At this point, the produced batches files can be executed.

    6a. failed -- (Optionnal) Identify failed execution that need to be rerun.

    6b. extract -- Extract the results from the executions.

    At this point, the analysis of the data is left to the user.
    We recommend Python and provide some functions to assist in that regard
    in the `MultiScaleData` module.


    One may want to convert the custom HDF5 format used by Illumina to a more
    standard format for use with GIS programs.
    This can be acheived with `convert`.
    """
    pass  # Entry point


# help
@main.command()
@click.argument("subcommand")
@click.pass_context
def help(ctx, subcommand):
    "Display detailled usage information of a command."
    ctx.info_name = subcommand
    subcommand_obj = main.get_command(ctx, subcommand)
    if subcommand_obj is None:
        click.echo("I don't know that command.")
    else:
        click.echo(subcommand_obj.get_help(ctx))


# alternate
@main.command()
@click.argument("name")
@click.option(
    "-z",
    "--zones",
    type=click.Path(exists=True),
    help="New zones inventory filename.",
)
@click.option(
    "-l",
    "--lights",
    type=click.Path(exists=True),
    help="New discrete lights inventory filename.",
)
def alternate(name, zones, lights):
    """Generates an alternate scenario at constant lumen.

    This scenatio will be based on the content of the `Inputs` folder and
    will be placed in a folder named `Inputs_NAME`.
    """
    illum.alternate(name, zones, lights)
    click.echo("Done.")


# batches
@main.command()
@click.argument("input_path", type=click.Path(exists=True), default=".")
@click.argument("batch_name", required=False)
@click.option(
    "-c",
    "--compact",
    is_flag=True,
    help=(
        "If given, will chain similar executions. Reduces the overall number"
        " of runs at the cost of longuer individual executions."
    ),
)
@click.option(
    "-N",
    "--batch_size",
    type=int,
    default=300,
    show_default=True,
    help="Number of runs per produced batch file.",
)
def batches(input_path, compact, batch_size, batch_name=None):
    """Makes the execution batches.

    INPUT_PATH is the path to the folder containing the inputs.

    BATCH_NAME is an optional name for the produced batch files.
    It overwrites the one defined in 'inputs_params.in' is given.
    """
    illum.batches(input_path, compact, batch_size, batch_name)
    click.echo("Done.")


# domain
@main.command()
def domain():
    """Defines the simulation domain.

    Reads domain definition parameters from 'domain_params.in'.

    Outputs the definition in 'domain.ini'.
    """
    illum.domain()
    click.echo("Done.")


# extract
@main.command()
@click.argument("exec_dir", default=".", type=click.Path(exists=True))
@click.option(
    "-c",
    "--contrib",
    is_flag=True,
    help="If present, extract contribution maps.",
)
@click.option(
    "-p",
    "--params",
    multiple=True,
    nargs=2,
    help=(
        "Parameter name,value pair to extract. Can be provided more than once."
    ),
)
@click.option(
    "-f",
    "--full",
    is_flag=True,
    help="If present, will extract all available outputs.",
)
@click.option(
    "-x",
    "--profile",
    is_flag=True,
    help="If present, extract profile along line of sight.",
)
def extract(exec_dir, contrib, params, full, profile):
    """Extract Illumina outputs.

    Will walk the EXEC_DIR to locate and extract illumina outputs.
    May fail if some execution failed, so one should validate the completude
    of the runs with 'illum failed' before using this.

    If not given, EXEC_DIR will default to the current directory.
    """
    illum.extract(exec_dir, contrib, params, full, profile)


# failed
@main.command()
@click.option(
    "-e",
    "--executable",
    is_flag=True,
    help="If given, returns the executable code to rerun failed executions.",
)
def failed(executable):
    "Find failed ILLUMINA executions."
    illum.failed(executable)


# init
@main.command()
def init():
    """Initialize an execution folder."""
    illum.init()
    click.echo("Done.")


# input
@main.group()
def input():
    """Prepares the executions inputs."""
    pass


# input iss
@input.command()
def iss():
    """Processes ISS images from cities at night."""
    illum.input_iss()
    click.echo("Done.")


# input pts
@input.command()
@click.argument("inventory", type=click.Path(exists=True))
@click.option(
    "-p",
    "--params",
    default="inputs_params.in",
    type=click.Path(exists=True),
    show_default=True,
    help="Path to the input parameters file.",
)
@click.option(
    "-o",
    "--output",
    default="Inputs",
    show_default=True,
    help="Folder name for the created input files.",
)
def pts(inventory, params, output):
    """Processes discrete light sources inventories."""
    illum.input_pts(inventory, params, output)
    click.echo("Done.")


# input viirs
@input.command()
@click.argument("inventory", type=click.Path(exists=True))
@click.option(
    "-p",
    "--params",
    default="inputs_params.in",
    type=click.Path(exists=True),
    show_default=True,
    help="Path to the input parameters file.",
)
@click.option(
    "-o",
    "--output",
    default="Inputs",
    show_default=True,
    help="Folder name for the created input files.",
)
def viirs(inventory, params, output):
    """Processes VIIRS data."""
    illum.input_viirs(inventory, params, output)
    click.echo("Done.")


# roads
@main.command()
@click.option(
    "-r",
    "--resolution",
    default=1.0,
    show_default=True,
    help="Output grid resolution [m]",
)
@click.option(
    "-x", "--xsize", default=10000.0, show_default=True, help="Tile x size [m]"
)
@click.option(
    "-y", "--ysize", default=10000.0, show_default=True, help="Tile y size [m]"
)
@click.option(
    "-b",
    "--buffer",
    default=200.0,
    show_default=True,
    help="Compute buffer [m]",
)
@click.option(
    "--crs",
    type=int,
    help=(
        "EPSG code of the projection to use. Determined automatically if not"
        " provided."
    ),
)
def roads(domain, resolution, xsize, ysize, buffer, crs):
    "Performs road analysis on the specified domain."
    illum.roads(domain, resolution, xsize, ysize, buffer, crs)
    click.echo("Done.")


# run
@main.command()
@click.argument(
    "inputs",
    nargs=-1,
    default=["Inputs"],
    type=click.Path(exists=True, file_okay=False),
)
def run(inputs):
    """Executes the illumina model for the specified inputs folder."""
    illum.run(inputs)
    click.echo("Done.")


# warp
@main.command()
@click.argument("target", type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    help="Base name of the output file. Based on TARGET when ommited.",
)
@click.option(
    "-d",
    "--domain",
    default="domain.parr",
    help="Reference polar array or parameter file.",
)
@click.option(
    "-f",
    "--field",
    help=(
        "Field name to burn for polygons. Inverse burns all-touched when"
        " ommited."
    ),
)
def warp(target, output, domain, field):
    """Warps georeferenced data.

    Will warp rasters and burn polygons.

    'target' can be a file or a folder containing multiple files at different
    scales or folders of tiles.
    """
    if output is None:
        output = os.path.splitext(os.path.basename(target))[0]
    illum.warp(target, output, domain, field)
    click.echo(f"Writting to '{output}.parr'.")


pipes = importlib.import_module(".UI.pipeline", package="illum")
main.add_command(pipes.pipe)
