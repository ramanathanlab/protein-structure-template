# protein-structure-template
Set of tools for writing templates for structured data extraction from protein structures.


# Tools

## DSSP

__Installation__

```
conda install -c salilab dssp
```

If you get errors about `libboost_threads*` you can install the package via conda with 

```
conda install anaconda::libboost==1.73.0
```

__Usage__

```
mkdssp -i [pdb/cif] -o [output.dssp]
```


# Goals

- Single protein/gene
    - [ ] DSSP for secondary structure
    - [ ] PDB header extraction - syntactically meaningful metadata extraction
    - [ ] point cloud stats(??, spread?)
    - [ ] structure bio-chemical stats (hydrophobicity, isoelectric points, solvent accessibility)
    - [ ] binding domain (this might be hard)
- Two protein/genes (includes all above)
	- [ ] TM-align for pairwise structural similarity
	- [ ] ??? (some tool, leaving in for edit) for identifying poorly folded regions
- Many protein/genes (includes all above)
    - [ ] ??? ideally this would be a MSA for structure that shows the relationship between them all but this is a full scale project I think