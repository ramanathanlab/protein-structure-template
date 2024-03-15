from pathlib import Path
import gzip
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Structure import Structure

"""
Notes:
-----------

- This probably isn't worth pursuing more as most prediction methods don't populate the header with anything useful.
  Here is an example from a ESM output:

    python header.py --input "/vol/structure/bvbrc/tmp/esmfold/fig|1000646.7.peg.1|KM576_sSgp1.pdb"
    >>> /vol/structure/bvbrc/tmp/esmfold/fig|100064: {'name': '', 'head': '', 'idcode': '', 'deposition_date': '1909-01-08', 'release_date': '1909-01-08', 'structure_method': 'unknown', 'resolution': None, 'structure_reference': [], 'journal_reference': '', 'author': '', 'compound': {'1': {'misc': ''}}, 'source': {'1': {'misc': ''}}, 'has_missing_residues': False, 'missing_residues': []}

  And from a alphafold output:
    python header.py --input /vol/structure/alphafold/taxa/1000002/AF-F8U1Q0-F1-model_v4.cif.gz
    >>> /vol/structure/alphafold/taxa/1000002/AF-F8U1Q0-F1-model_v4: {'name': '', 'head': '', 'idcode': '', 'deposition_date': '2022-06-01', 'structure_method': '', 'resolution': None}


So at this point its probably only useful to serve as an example of loading/basic parsing of a structure file coming from BVBRC group.
"""


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


if __name__ == "__main__":
    parser = ArgumentParser("Structure header from PDBx/mmCIF file")
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
        default="json",
        help="Output format, defaults to json",
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
    # TODO: Think about streaming to an output format here
    with ProcessPoolExecutor(max_workers=args.num_cpus) as executor:
        for fp, res in zip(files, executor.map(get_structure, files)):
            struct_name = str(fp).rstrip("".join(fp.suffixes))
            data[struct_name] = res

    for k, v in data.items():
        print(f"{k}: {v.header}")
