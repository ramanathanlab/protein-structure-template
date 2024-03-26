from __future__ import annotations

from pathlib import Path
from random import choice
from random import seed
from typing import Any
from typing import Optional

from jinja2 import Environment
from jinja2 import FileSystemLoader

from pst.dssp import run_and_parse_dssp
from pst.stats import get_structure_statistics

seed(42)


# TODO: might need module specific kwargs? this header could get bulky if not
def generate_description(
    structure_fp: Path,
    dssp_path: Path,
    dssp_fields: Optional[list[str]],
) -> dict[str, Any]:
    """Generate a description of a structure

    Parameters
    ----------
    structure : Path
        Location of structure file

    Returns
    -------
    str
        Description of the structure
    """
    single_structure_env = (
        Path(__file__).parent.parent / 'templates' / 'single'
    )
    # choose random template from directory of single templates
    available_templates = [
        e.name for e in list(single_structure_env.glob('*'))
    ]
    template_name = choice(available_templates)
    # TODO: add logging instead of printing
    print('Using template:', template_name)

    env = Environment(loader=FileSystemLoader(single_structure_env))
    template = env.get_template(template_name)

    gene_id = structure_fp.stem
    # TODO: fix this hack for Alphafold DB
    gene_id = gene_id.replace('AF-', '')
    gene_id = gene_id.replace('-F1-model_v4', '')

    context = {'gene_id': gene_id}

    simple_stats = get_structure_statistics(structure_fp)
    context.update(simple_stats)

    if dssp_fields is None:
        # Default to just secondary structure
        dssp_fields = ['ss']

    ss_stats = run_and_parse_dssp(structure_fp, dssp_path, dssp_fields)
    context.update(ss_stats)

    return {'text': template.render(context), 'context': context}


if __name__ == '__main__':
    struct_path = Path(
        '/home/khippe/small-structs/1000001/AF-F8U1P7-F1-model_v4.cif',
    )
    dssp_path = Path('/home/khippe/miniforge3/envs/dssp/bin/mkdssp')
    dssp_fields = ['ss']

    text = generate_description(struct_path, dssp_path, dssp_fields)

    print(text)
