from __future__ import annotations
from pathlib import Path
from typing import Any, List, TypeVar, Union
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from concurrent.futures import ProcessPoolExecutor
import json

from pydantic import BaseModel, model_validator, field_validator

T = TypeVar('T', bound='BaseModel')
PathLike = Union[str, Path]


class Residue(BaseModel):
    """Container for a residue as it is output from DSSP"""
    number: int
    """Global index (regardless of chain)"""
    index: int
    """Resiude index in chain"""
    chain: str
    """Chain index (letter indexed)"""
    residue: str
    """Residue one letter code"""
    ss: str
    """Secondary structure single letter code (7 designations)"""
    acc: int
    """Solvent accessability (unsure of scaling here)"""
    h_bonds: Union[List[float], str]
    """List of hbonds, currently unparsed because they are nasty"""
    tco: float
    """cosine of angle between C=O of residue i and C=O of residue i-1"""
    kappa: float
    """virtual bond angle (bend angle) defined by the three C-alpha atoms of residues i-2,i,i+2."""
    alpha: float
    """virtual torsion angle (dihedral angle) defined by the four C-alpha atoms of residues i-1,i,i+1,i+2"""
    phi: float
    """IUPAC peptide backbone torsion angles"""
    psi: float
    """IUPAC peptide backbone torsion angles"""
    x: float
    """X coordinate in 3D space"""
    y: float
    """Y coordinate in 3D space"""
    z: float
    """Z coordinate in 3D space"""

    @model_validator(mode="before")
    def strip_whitespace(cls, values: dict[str, Any]) -> dict[str, Any]:
        stripped_values = {}
        for k, v in values.items():
            if isinstance(v, str):
                stripped_values[k] = v.strip()
            else:
                stripped_values[k] = v
        return stripped_values
    
    @field_validator("ss")
    @classmethod
    def validate_ss(cls, v: Any) -> Any:
        if isinstance(v, str):
            if len(v.strip()) == 0:
                return "-"

        return v

class Protein(BaseModel):
    residues: List[Residue]

    @classmethod
    def from_dssp(cls: type[T], path: PathLike) -> T:
        # How/where to grab fields from the dssp files, can view at
        # https://web.archive.org/web/20230305032352/https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
        fields = {
            'number': slice(0, 5),          # Global index (regardless of chain)
            'index': slice(5, 10),          # index in chain
            'chain': slice(10, 12),         # Chain index (letter based)
            'residue': slice(12, 16),       # Single residue letter # code
            'ss': slice(16, 17),            # Secondary structure single letter code
            'acc': slice(35, 38),           # Solvent accessability
            'h_bonds': slice(38, 85),       # string list of hbonds
            'tco': slice(85, 91),           # cosine of angle between C=O of residue i and C=O of residue i-1
            'kappa': slice(91, 97),         # virtual bond angle (bend angle) defined by the three C-alpha atoms of residues i-2,i,i+2.
            'alpha': slice(97, 103),        # virtual torsion angle (dihedral angle) defined by the four C-alpha atoms of residues i-1,i,i+1,i+2
            'phi': slice(103, 109),         # IUPAC peptide backbone torsion angles
            'psi': slice(109, 115),         # IUPAC peptide backbone torsion angles
            'x': slice(115, 122),           # X coordinate
            'y': slice(122, 129),           # Y coordinate
            'z': slice(129, 136),           # Z coordinate
        }

        residues = []
        past_header = False
        with open(path, 'r') as file:
            for line in file.readlines():
                # Need to iterate past a large/variable length header
                if not past_header:
                    if '#' in line.split()[0]:  # need to check if its first elem
                        past_header = True
                    continue

                residue = Residue(**{f: line[e] for f, e in fields.items()})

                # A '!' specifies a break in the chain or invalid residue
                if "!" in residue.residue:
                    continue

                residues.append(residue)

        return cls(residues=residues)


def output_json(parsed_data: dict[str, Protein], output_file: Path) -> None:
    """Output single jsonl file with key being the file it came from (without suffix) and the values being a dictionary of lists of included fields

    *Will not work on especially large datasets (no incremental writing at this point, and json is inherently limiting)*

    Parameters
    ----------
    parsed_data : dict[str, Protein]
        Parsed DSSP data per structure, keyed by the structure filename
    output_file : Path
        output json file path
    """
    output_data = {}
    for filename, prot in parsed_data.items():
        localout = {}
        for k in args.fields:
            localout[k] = [getattr(elem, k) for elem in prot.residues]

        output_data[filename] = localout

    with open(output_file, 'w') as f:
        json.dump(output_data, f)


# TODO
def output_template(parsed_data: dict[str, Protein], output_file: Path) -> None:
    ...


if __name__ == "__main__":
    parser = ArgumentParser("DSSP Secondary Structure Parser", formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("--dssp-dir", type=Path, required=True, help="Path to a directory containing DSSP files, will recursively glob")
    parser.add_argument("--glob-pattern", '-g', default="**/*.dssp", help="Pattern to search for in the `dssp-dir`")
    parser.add_argument("--output", type=Path, help="Path to output location (can be file)")
    parser.add_argument("--format", "-f", choices=['json', 'template'], default='json', help="Output format, defaults to json")
    parser.add_argument("--fields", choices=list(Residue.model_fields.keys()), action="append", help=f"Fields to dump to file, available are {list(Residue.model_fields.keys())}")
    parser.add_argument("--num-cpus", "-n", type=int, default=1, help="Number of processors to use for parsing output")

    args = parser.parse_args()

    dssp_dir = args.dssp_dir
    files = list(dssp_dir.glob(args.glob_pattern))

    data = {}

    # Because appending and defaults don't mix well, manually set output fields if they are not given
    if args.fields is None:
        args.fields = list(Residue.model_fields())

    # TODO: Think about streaming to an output format here
    with ProcessPoolExecutor(max_workers=args.num_cpus) as executor:
        for file, res in zip(files, executor.map(Protein.from_dssp, files)):
            data[file.name.replace('.dssp', '')] = res

    match args.format:
        case "json":
            output_json(data, args.output)
