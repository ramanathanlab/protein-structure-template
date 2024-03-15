from __future__ import annotations
from pathlib import Path
import subprocess
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from typing import Any
import json
import re
import itertools


def run_tmalign(file1: Path, file2: Path) -> tuple[float, float]:
    """Runs tm-align on two pdb files and returns the TM-scores, the first normalized to the first input and the second normalized to the second input

    Parameters
    ----------
    file1 : Path
        filepath to first file
    file2 : Path
        Filepath to second file

    Returns
    -------
    tuple[float, float]
        TM-Scores, first normalized to first input, second normalized to second input
    """
    cmd = f"{str(args.tm_align)} {str(file1)} {str(file2)}"
    res = subprocess.run(cmd.split(), capture_output=True)
    # Score 1 is normalized by pdb 2 score 2 is normalized by pdb 1
    scores = pattern.findall(res.stdout.decode("utf-8"))

    return float(scores[0]), float(scores[1])


def output_json(parsed_data: dict[str, Any], output_file: Path) -> None:
    """Output JSON representation of PDB structure general statistics

    Parameters
    ----------
    parsed_data : dict[str, Any]
        Parsed structure data, keys are structure names, values are the statistics, all json serializable
    output_file : Path
        Path to save overall output file
    """

    with open(output_file, "w") as f:
        json.dump(parsed_data, f)


def output_template(parsed_data: dict[str, Any], output_file: Path) -> None:
    pass


if __name__ == "__main__":
    pattern = re.compile(r"^TM-score= ([+-]?[0-9]*[.]?[0-9]+)", flags=re.MULTILINE)

    parser = ArgumentParser("Structure statistics from PDBx/mmCIF file")
    parser.add_argument(
        "--tm-align",
        default="TMalign",
        type=Path,
        help="Path to the tm-align executable",
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Get header of structure file. If input is a directory will get all structures in directory. If input is a file will get header of that file",
    )
    parser.add_argument(
        "--glob-pattern",
        "-g",
        default="**/*.cif",
        help="Pattern to search for in the input dir. NOTE: cannor process *.gz files, must be unzipped",
    )
    parser.add_argument(
        "--output", type=Path, help="Path to output location (can be file)"
    )
    parser.add_argument(
        "--format",
        "-f",
        choices=["json", "template"],
        help="Output format, defaults to stdout if not specified",
    )
    parser.add_argument(
        "--num-cpus",
        "-n",
        type=int,
        default=1,
        help="Number of processors to use for parsing output",
    )

    args = parser.parse_args()

    file_input: Path = args.input
    if file_input.is_file():
        raise ValueError("Input must be a directory")
    elif file_input.is_dir():
        files = list(file_input.glob(args.glob_pattern))
    file_pairs = list(itertools.combinations(files, 2))

    data = {}
    # TODO: Same as dssp, think about streaming to an output format here
    with ProcessPoolExecutor(max_workers=args.num_cpus) as executor:
        # Here's the crucial change: Use a lambda to unpack the file pairs
        futures = [executor.submit(run_tmalign, *file_pair) for file_pair in file_pairs]

        for future, file_pair in zip(futures, file_pairs):
            res = future.result()
            struct_name = "--".join(
                str(fp.stem) for fp in file_pair
            )  # Changed how struct_name is formed
            data[struct_name] = res

    match args.format:
        case "json":
            output_json(data, args.output)

        case None:
            for fname, stats in data.items():
                print(f"{fname}: {stats}")
