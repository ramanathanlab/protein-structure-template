from __future__ import annotations

import gzip
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

from pst.output import OutputType

__all__ = ['structure_statistics', 'get_structure_statistics']

AA_NAME_CONVERSION = {
    'A': 'Alanine',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic Acid',
    'C': 'Cysteine',
    'E': 'Glutamic Acid',
    'Q': 'Glutamine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'S': 'Serine',
    'T': 'Threonine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine',
    'V': 'Valine',
}


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
    struct_id = str(fp).rstrip(''.join(fp.suffixes))
    if '.pdb' in fp.suffixes:
        parser = PDBParser(QUIET=True)
    elif '.cif' in fp.suffixes:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f'File {fp} type not supported')

    if '.gz' in fp.suffixes:
        file = gzip.open(fp, 'rt')
    else:
        file = open(fp)

    structure = parser.get_structure(struct_id, file)

    file.close()

    return structure


# TODO: might be good to parameterize some of the fields (e.g top N)
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
    amino_acids = seq1(
        ''.join([res.resname for res in structure.get_residues()]),
    )

    structure_output = {}
    # Simple statistics
    structure_output['num_models'] = len(structure)
    structure_output['num_chains'] = len(
        [chain for chain in structure.get_chains()],
    )
    structure_output['num_residues'] = len(
        [res for res in structure.get_residues()],
    )
    structure_output['num_atoms'] = len(
        [atom for atom in structure.get_atoms()],
    )
    structure_output['amino_acid_count'] = dict(Counter(amino_acids))
    # top 3 most common amino acids and their total percentage
    total_aa = structure_output['num_residues']
    common_residues = sorted(
        structure_output['amino_acid_count'].items(),
        key=lambda x: x[1],
        reverse=True,
    )
    structure_output['top3_residues'] = [
        AA_NAME_CONVERSION[res[0]] for res in common_residues[:3]
    ]
    structure_output['top3_percent'] = (
        sum(res[1] for res in common_residues[:3]) / total_aa
    ) * 100

    # Isoelectic point + various related info
    ip = IsoelectricPoint(amino_acids)
    structure_output['isoelectric_point'] = ip.pi()

    # Solvent accessible surface area
    sr = ShrakeRupley()
    sr.compute(structure, level='S')
    sr.compute(structure, level='R')
    structure_output['sasa'] = structure.sasa

    return structure_output


def output_json(
    parsed_data: dict[str, Any],
    output_file: Path,
    fields: list[str],
) -> None:
    """Output JSON representation of PDB structure general statistics

    Parameters
    ----------
    parsed_data : dict[str, Any]
        Parsed structure data, keys are structure names, values are the statistics, all json serializable
    output_file : Path
        Path to save overall output file
    """
    output_data = {k: v for k, v in parsed_data.items() if k in fields}
    with open(output_file, 'w') as f:
        json.dump(output_data, f)


def output_template(parsed_data: dict[str, Any], output_file: Path) -> None:
    pass


def structure_statistics(
    input_dir: Path,
    glob_pattern: str,
    output: Path | None = None,
    output_format: str | None = None,
    fields: list[str] | None = None,
    num_cpus: int = 1,
) -> None:
    if input_dir.is_file():
        files = [input_dir]
    elif input_dir.is_dir():
        files = list(input_dir.glob(glob_pattern))

    data = {}
    # TODO: Same as dssp, think about streaming to an output format here
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        for fp, res in zip(
            files,
            executor.map(get_structure_statistics, files),
        ):
            struct_name = str(fp).rstrip(''.join(fp.suffixes))
            data[struct_name] = res

    match output_format:
        case OutputType.json:
            assert (
                output is not None
            ), 'Output file must be specified for json output'
            output_json(data, output, fields)

        case OutputType.stdout:
            for fname, stats in data.items():
                print(f'{Path(fname).stem}: {stats}')
