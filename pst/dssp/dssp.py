"""DSSP module, hooks into the DSSP executable and parses its output."""

from __future__ import annotations

import json
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any
from typing import Optional
from typing import TypeVar
from typing import Union

from pydantic import BaseModel
from pydantic import field_validator
from pydantic import model_validator

from pst.output import OutputType

T = TypeVar("T", bound="BaseModel")
PathLike = Union[str, Path]

__all__ = ["Protein", "Residue", "process_dssp", "run_and_parse_dssp"]


# TODO: believe it or not the dataclass is cumbersome... figure out how to cleanup with respect to the necessary output
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
    """Solvent accessibility (unsure of scaling here)"""
    h_bonds: Union[list[float], str]
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
    residues: list[Residue]

    @classmethod
    def from_str(cls: type[T], data: str) -> T:
        return cls.parse_data(data)

    @classmethod
    def from_dssp(cls: type[T], path: PathLike) -> T:
        data = path.read_text()
        return cls.parse_data(data)

    @classmethod
    def parse_data(cls: type[T], data: str) -> T:
        # How/where to grab fields from the dssp files, can view at
        # https://web.archive.org/web/20230305032352/https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
        fields = {
            "number": slice(0, 5),  # Global index (regardless of chain)
            "index": slice(5, 10),  # index in chain
            "chain": slice(10, 12),  # Chain index (letter based)
            "residue": slice(12, 16),  # Single residue letter # code
            "ss": slice(16, 17),  # Secondary structure single letter code
            "acc": slice(35, 38),  # Solvent accessibility
            "h_bonds": slice(38, 85),  # string list of hbonds
            "tco": slice(
                85,
                91,
            ),  # cosine of angle between C=O of residue i and C=O of residue i-1
            "kappa": slice(
                91,
                97,
            ),  # virtual bond angle (bend angle) defined by the three C-alpha atoms of residues i-2,i,i+2.
            "alpha": slice(
                97,
                103,
            ),  # virtual torsion angle (dihedral angle) defined by the four C-alpha atoms of residues i-1,i,i+1,i+2
            "phi": slice(103, 109),  # IUPAC peptide backbone torsion angles
            "psi": slice(109, 115),  # IUPAC peptide backbone torsion angles
            "x": slice(115, 122),  # X coordinate
            "y": slice(122, 129),  # Y coordinate
            "z": slice(129, 136),  # Z coordinate
        }

        residues = []
        past_header = False
        for line in data.strip().split("\n"):
            # Need to iterate past a large/variable length header
            if not past_header:
                if "#" in line.split()[0]:  # need to check if its first elem
                    past_header = True
                continue

            residue = Residue(**{f: line[e] for f, e in fields.items()})

            # A '!' specifies a break in the chain or invalid residue
            if "!" in residue.residue:
                continue

            residues.append(residue)

        return cls(residues=residues)


def output_json(
    parsed_data: dict[str, Protein],
    output_file: Path,
    fields: list[str],
) -> None:
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
        for k in fields:
            localout[k] = [getattr(elem, k) for elem in prot.residues]

        output_data[filename] = localout

    with open(output_file, "w") as f:
        json.dump(output_data, f)


# TODO
def output_template(parsed_data: dict[str, Protein], output_file: Path) -> None: ...


def process_dssp(
    dssp_dir: Path,
    output: Path,
    glob_pattern: str,
    output_format: str | None = None,
    fields: list[str] | None = None,
    num_cpus: int = 1,
) -> None:
    files = list(dssp_dir.glob(glob_pattern))

    data = {}

    # Because appending and defaults don't mix well, manually set output fields if they are not given
    if fields is None:
        fields = list(Residue.model_fields)

    # TODO: Think about streaming to an output format here
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        for file, res in zip(files, executor.map(Protein.from_dssp, files)):
            data[file.name.replace(".dssp", "")] = res

    match output_format:
        case OutputType.json:
            assert output is not None, "Output file must be specified for json output"
            output_json(data, output, fields)

        case OutputType.stdout:
            for filename, prot in data.items():
                localout = {}
                for k in fields:
                    localout[k] = " ".join(
                        [str(getattr(elem, k)) for elem in prot.residues],
                    )
                print(f"{filename}:")
                for k, v in localout.items():
                    print(f"\t{k}: {v}")
                print("\n")


def run_and_parse_dssp(
    pdb_file: Path,
    dssp_path: Path,
    fields: Optional[list[str]] = None,
) -> dict[str, Any]:
    """Run DSSP and parse the output into a dataclass

    Parameters
    ----------
    pdb_file : Path
        Path to the pdb/mmcif file
    dssp_path : Path
        Path to the DSSP executable

    Returns
    -------
    str
        The contents of the output
    """
    cmd = f"{dssp_path} -i {pdb_file}"
    res = subprocess.run(cmd.split(), capture_output=True, check=False)

    assert res.returncode == 0, f"Error running DSSP: {res.stderr}"

    dssp_data = res.stdout.decode("utf-8")
    protein = Protein.from_str(dssp_data)

    if fields is None:
        fields = list(Residue.model_fields.keys())

    structure_data = {}
    for k in fields:
        structure_data[k] = [getattr(elem, k) for elem in protein.residues]

    # Hard coded calculations for now but we need this for the `context`
    aa_len = len(protein.residues)
    helix_code = ["H", "G", "I"]
    beta_code = ["E", "B"]
    coil_code = ["T", "S", " "]
    alpha_percent = (
        sum([1 for ss in structure_data["ss"] if ss in helix_code]) / aa_len
    ) * 100
    beta_percent = (
        sum([1 for ss in structure_data["ss"] if ss in beta_code]) / aa_len
    ) * 100
    coil_percent = (
        sum([1 for ss in structure_data["ss"] if ss in coil_code]) / aa_len
    ) * 100

    structure_data["alpha_percent"] = alpha_percent
    structure_data["beta_percent"] = beta_percent
    structure_data["coil_percent"] = coil_percent

    most_common = [
        ("alpha helix", alpha_percent),
        ("beta sheet", beta_percent),
        ("coil", coil_percent),
    ]
    most_common = sorted(most_common, key=lambda x: x[1], reverse=True)
    structure_data["ordered_ss"] = most_common

    return structure_data


if __name__ == "__main__":
    pdb_path = Path("/home/khippe/small-structs/1000001/AF-F8U1P7-F1-model_v4.cif")

    dssp_path = Path("/home/khippe/miniforge3/envs/dssp/bin/mkdssp")

    print(run_and_parse_dssp(pdb_path, dssp_path))
