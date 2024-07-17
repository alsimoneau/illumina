#!/usr/bin/env python3

import click

import illum


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
@main.command
@click.argument("name")
@click.option(
    "-z", "--zones", type=click.Path(exists=True), help="New zones inventory filename."
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


# batches
@main.command
@click.argument("input_path", type=click.Path(exists=True), default=".")
@click.argument("batch_name", required=False)
@click.option(
    "-c",
    "--compact",
    is_flag=True,
    help="If given, will chain similar executions. Reduces the overall number "
    "of runs at the cost of longuer individual executions.",
)
@click.option(
    "-N",
    "--batch_size",
    type=int,
    default=300,
    show_default=True,
    help="Number of runs per produced batch file.",
)
@click.option(
    "-s",
    "--scheduler",
    type=click.Choice(["parallel", "sequential", "slurm"]),
    default="sequential",
    help="Job scheduler",
)
def batches(input_path, compact, batch_size, scheduler, batch_name=None):
    """Makes the execution batches.

    INPUT_PATH is the path to the folder containing the inputs.

    BATCH_NAME is an optional name for the produced batch files.
    It overwrites the one defined in 'inputs_params.in' is given.
    """
    illum.batches(input_path, compact, batch_size, scheduler, batch_name)


# convert
@main.command
@click.option("--vector/--raster", "-v/-r", default=True, help="Output type.")
@click.option("-log", is_flag=True, help="Logarithmic scale (base 10)")
@click.option("-area", is_flag=True, help="Normalized by pixel area (in km²)")
@click.argument("filename", type=click.Path(exists=True))
@click.argument("outname")
def convert(filename, outname, vector, log, area):
    """Convert an Illumina HDF file to a georeferenced format.

    Converts FILENAME to OUTNAME.EXT where ext is defined based on the output
    format. The output format is either vector (GeoJSON) or raster (Tiff).
    """
    illum.convert(filename, outname, vector, log, area)


# domain
@main.command
def domain():
    """Defines the simulation domain.

    Reads domain definition parameters from 'domain_params.in'.

    Outputs the definition in 'domain.ini'.
    """
    illum.domain()


# extract
@main.command
@click.argument("exec_dir", default=".", type=click.Path(exists=True))
@click.option(
    "-c", "--contrib", is_flag=True, help="If present, extract contribution maps."
)
@click.option(
    "-p",
    "--params",
    multiple=True,
    nargs=2,
    help="Parameter name,value pair to extract." " Can be provided more than once.",
)
@click.option(
    "-f", "--full", is_flag=True, help="If present, will extract all available outputs."
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
@main.command
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
@main.command
def init():
    """Initialize an execution folder."""
    illum.init()


# inputs
@main.command
def inputs():
    """Prepares the executions inputs."""
    illum.inputs()


# warp
@main.command
@click.argument("output_name", required=False)
@click.argument("infiles", required=False, nargs=-1)
def warp(output_name=None, infiles=[]):
    """Warps the satellite imagery.

    Warps the satellite imagery based on the domain defined in
    'domain.ini'.

    \b
    Requires the folowing data:
        Unzipped SRTM data in a folder named 'SRTM'.
        If used, VIIRS data in a volder named 'VIIRS-DNB'.
        If VIIRS data is used, the 'hydropolys.zip' file.

    Can also be used on specific files, in wich case an output name and a list
    of files to process must be given (the use of bash wildcards is encouraged)
    """
    illum.warp(output_name, infiles)
