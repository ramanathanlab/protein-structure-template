from __future__ import annotations
from pathlib import Path
import gzip
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from typing import Any
import json

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Structure import Structure
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils import seq1
from Bio.PDB.SASA import ShrakeRupley


def get_structure(fp: Path) -> Structure:
    """Get a structure from a PDBx/mmCIF file

    Parameters
    ----------
    fp : Path
        Location of structure file

    Returns
    -------
    Structure
        Biopython Structure object representing the structure
    """

    struct_id = str(fp).rstrip("".join(fp.suffixes))
    if ".pdb" in fp.suffixes:
        parser = PDBParser(QUIET=True)
    elif ".cif" in fp.suffixes:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"File {fp} type not supported")

    if ".gz" in fp.suffixes:
        file = gzip.open(fp, "rt")
    else:
        file = open(fp, "r")

    structure = parser.get_structure(struct_id, file)

    file.close()

    return structure


def get_structure_statistics(fp: Path) -> dict[str, str]:
    """Get statistics from a structure file

    Parameters
    ----------
    fp : Path
        Location of structure file

    Returns
    -------
    dict[str, str]
        Dictionary of statistics from the structure
    """

    structure = get_structure(fp)
    amino_acids = seq1("".join([res.resname for res in structure.get_residues()]))

    structure_output = {}
    # Simple statistics
    structure_output["num_models"] = len(structure)
    structure_output["num_chains"] = len([chain for chain in structure.get_chains()])
    structure_output["num_residues"] = len([res for res in structure.get_residues()])
    structure_output["num_atoms"] = len([atom for atom in structure.get_atoms()])
    structure_output["amino_acid_count  "] = dict(Counter(amino_acids))

    # Isoelectic point + various related info
    ip = IsoelectricPoint(amino_acids)
    structure_output["isoelectric_point"] = ip.pi()
    # structure_output['charged_aa_content'] = ip.charged_aas_content

    # Solvent accessible surface area
    sr = ShrakeRupley()
    sr.compute(structure, level="S")
    structure_output["sasa"] = structure.sasa

    return structure_output


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
    parser = ArgumentParser("Structure statistics from PDBx/mmCIF file")
    parser.add_argument(
        "--input",
        type=Path,
        help="Get header of structure file. If input is a directory will get all structures in directory. If input is a file will get header of that file",
    )
    parser.add_argument(
        "--glob-pattern",
        "-g",
        default="**/*.cif*",
        help="Pattern to search for in the `dssp-dir`",
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
        files = [file_input]
    elif file_input.is_dir():
        files = list(file_input.glob(args.glob_pattern))

    data = {}
    # TODO: Same as dssp, think about streaming to an output format here
    with ProcessPoolExecutor(max_workers=args.num_cpus) as executor:
        for fp, res in zip(files, executor.map(get_structure_statistics, files)):
            struct_name = str(fp).rstrip("".join(fp.suffixes))
            data[struct_name] = res

    match args.format:
        case "json":
            output_json(data, args.output)

        case None:
            for fname, stats in data.items():
                print(f"{fname}: {stats}")
