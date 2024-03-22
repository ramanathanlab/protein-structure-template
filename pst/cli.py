"""CLI for pst."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from pst.output import OutputType

app = typer.Typer(
    add_completion=False,
    pretty_exceptions_show_locals=False,
    help="PST - a program for generating natural language descriptions of protein structures.",
)

dssp_app = typer.Typer(
    add_completion=False,
    pretty_exceptions_show_locals=False,
    help="Module for running DSSP and processing its output.",
)
tmalign_app = typer.Typer(
    add_completion=False,
    pretty_exceptions_show_locals=False,
    help="Module for running TM-Align and formatting its output.",
)
stats_app = typer.Typer(
    add_completion=False,
    pretty_exceptions_show_locals=False,
    help="Module for generating simple statistics on protein structures.",
)

app.add_typer(dssp_app, name="dssp")
app.add_typer(tmalign_app, name="tmalign")
app.add_typer(stats_app, name="stats")


def get_residue_fields() -> list[str]:
    from pst.dssp import Residue

    return list(Residue.model_fields.keys())


# DSSP CLI
@dssp_app.command()
def process_dssp(
    dssp_dir: Path = typer.Option(
        ...,
        "--dssp-dir",
        help="Path to a directory containing DSSP files, will recursively glob",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Path to output location (if none, will output to stdout)",
    ),
    glob_pattern: str = typer.Option(
        "**/*.dssp",
        "--glob-pattern",
        "-g",
        help="Pattern to search for in the `dssp-dir`",
    ),
    output_format: OutputType = typer.Option(
        OutputType.stdout,
        "--format",
        "-f",
        help="Output format, if none specified will print to stdout",
    ),
    fields: list[str] = typer.Option(
        get_residue_fields(),
        "--fields",
        help=f"Fields to dump to file, available are {get_residue_fields()}",
    ),
    num_cpus: int = typer.Option(
        1, "--num-cpus", "-n", help="Number of processors to use for parsing output"
    ),
    test_choices: str = typer.Option(
        "stdout",
        "--test-choices",
        help="Test the choices feature of typer",
    ),
) -> None:
    from pst.dssp import process_dssp

    process_dssp(dssp_dir, output, glob_pattern, output_format, fields, num_cpus)


# TMAlign CLI
@tmalign_app.command()
def run_tmalign(
    tmalign_path: Path = typer.Option(
        "TMalign", "--tm-align", help="Path to the TMAlign executable"
    ),
    input_dir: Path = typer.Option(
        ...,
        help="Directory of PDB/CIF files, will peform pairwise comparison on all files in the directory",
    ),
    glob_pattern: str = typer.Option(
        "**/*.cif",
        "--glob-pattern",
        "-g",
        help="Pattern to search for in the input dir. NOTE: cannot process *.gz files, must be unzipped",
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", help="Path to output location (can be file)"
    ),
    output_format: OutputType = typer.Option(
        OutputType.stdout,
        "--format",
        "-f",
        help="Output format, defaults to stdout if not specified",
    ),
    num_cpus: int = typer.Option(
        1, "--num-cpus", "-n", help="Number of processors to use for parsing output"
    ),
) -> None:
    from pst.tmalign import run_tmalign

    run_tmalign(tmalign_path, input_dir, glob_pattern, output, output_format, num_cpus)


# Stats CLI
@stats_app.command()
def generate_stats(
    input_dir: Path = typer.Option(
        ...,
        help="Directory of PDB/CIF files, will peform pairwise comparison on all files in the directory",
    ),
    glob_pattern: str = typer.Option(
        "**/*.cif",
        "--glob-pattern",
        "-g",
        help="Pattern to search for in the input dir. NOTE: cannot process *.gz files, must be unzipped",
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", help="Path to output location (can be file)"
    ),
    output_format: Optional[OutputType] = typer.Option(
        None,
        "--format",
        "-f",
        help="Output format, defaults to stdout if not specified",
    ),
    fields: list[str] = typer.Option(
        get_residue_fields(),
        "--fields",
        help=f"Fields to dump to file, available are {get_residue_fields()}",
    ),
    num_cpus: int = typer.Option(
        1, "--num-cpus", "-n", help="Number of processors to use for parsing output"
    ),
) -> None:
    from pst.stats import structure_statistics

    structure_statistics(
        input_dir, glob_pattern, output, output_format, fields, num_cpus
    )


def main() -> None:
    app()


if __name__ == "__main__":
    main()
