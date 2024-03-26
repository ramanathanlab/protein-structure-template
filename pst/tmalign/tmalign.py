from __future__ import annotations
from pathlib import Path
import subprocess
from concurrent.futures import ProcessPoolExecutor
from typing import Any
import json
import re
import itertools
from functools import partial

from pst.output import OutputType

__all__ = ["run_tmalign"]


def align(
    file1: Path, file2: Path, tmalign_path: Path, score_pattern: re.Pattern
) -> tuple[float, float]:
    """Runs tm-align on two pdb files and returns the TM-scores, the first normalized to the first input and the second normalized to the second input

    Parameters
    ----------
    file1 : Path
        filepath to first file
    file2 : Path
        Filepath to second file
    tmalign_path : Path
        Path to the tm-align executable
    score_pattern : re.Pattern
        Compiled regex pattern to extract TM-scores from tm-align output

    Returns
    -------
    tuple[float, float]
        TM-Scores, first normalized to first input, second normalized to second input
    """
    cmd = f"{str(tmalign_path)} {str(file1)} {str(file2)}"
    res = subprocess.run(cmd.split(), capture_output=True)
    # Score 1 is normalized by pdb 2 score 2 is normalized by pdb 1
    scores = score_pattern.findall(res.stdout.decode("utf-8"))

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


# Entry points
def run_tmalign(
    tm_align: Path,
    input_dir: Path,
    glob_pattern: str,
    output: Path | None = None,
    output_format: str | None = None,
    num_cpus: int = 1,
) -> None:
    score_pattern = re.compile(
        r"^TM-score= ([+-]?[0-9]*[.]?[0-9]+)", flags=re.MULTILINE
    )

    if input_dir.is_file():
        raise ValueError("Input must be a directory")
    elif input_dir.is_dir():
        files = list(input_dir.glob(glob_pattern))
    file_pairs = list(itertools.combinations(files, 2))

    align_prefilled = partial(align, tmalign_path=tm_align, score_pattern=score_pattern)

    data = {}
    # TODO: Same as dssp, think about streaming to an output format here
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        # Here's the crucial change: Use a lambda to unpack the file pairs
        futures = [
            executor.submit(align_prefilled, *file_pair) for file_pair in file_pairs
        ]

        for future, file_pair in zip(futures, file_pairs):
            res = future.result()
            struct_name = "--".join(
                str(fp.stem) for fp in file_pair
            )  # Changed how struct_name is formed
            data[struct_name] = res

    match output_format:
        case OutputType.json:
            output_json(data, output)

        case OutputType.stdout:
            for fname, stats in data.items():
                print(f"{fname}: {stats}")
