# protein-structure-template
Set of tools for writing templates for structured data extraction from protein structures.


# Templating

__Available Fields__
|  Field  | Type  | Description            |
| :-----: | :---: | :--------------------- |
| `$name` | `str` | identifyer for protein |
| `$residue_dist`| `dict[str, float]` | distribution of amino acid occurences. Ordered dict (high to lowest occurences)
| `$secondary_structure` | `dict[str, float]`| secondary structure        distribution, keys are helix, sheet, coil, and values are percentage of overall structure, sorted by highest occurring (first key) |
|`$isoelectric` | `float` | single value isoelectric point value, pH value for neutral charge |
| `$sasa` | `float` | solvent accessible surface area |
| | | |
| `$TMAlign` | :rotating_light:TODO:rotating_light: |


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

__Output__
Taken directly from the (now defunct) [DSSP webpage](https://web.archive.org/web/20230214164302/https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html)

```
HEADER    HYDROLASE   (SERINE PROTEINASE)         17-MAY-76   1EST
...
  240  1  4  4  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS,
                  NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .
 10891.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)
  162 67.5   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  ; PER 100 RESIDUES
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES; PER 100 RESIDUES
   84 35.0   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES; PER 100 RESIDUES
...
   26 10.8   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+2)
   30 12.5   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+3)
   10  4.2   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4)
...
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N
    2   17   V  B 3   +A  182   0A   8  180,-2.5 180,-1.9   1,-0.2 134,-0.1
                                   ...Next two lines wrapped as a pair...
                                    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
                                  -0.776 360.0   8.1 -84.5 125.5  -14.7   34.4   34.8
                                   ...Next two lines wrapped as a pair...
                                               CHAIN AUTHCHAIN
                                                   A         A
....;....1....;....2....;....3....;....4....;....5....;....6....;....7..
    .-- sequential resnumber, including chain breaks as extra residues
    |    .-- original PDB resname, not nec. sequential, may contain letters
    |    | .-- one-letter chain ID, if any
    |    | | .-- amino acid sequence in one letter code
    |    | | |  .-- secondary structure summary based on columns 19-38
    |    | | |  | xxxxxxxxxxxxxxxxxxxx recommend columns for secstruc details
    |    | | |  | .-- 3-turns/helix
    |    | | |  | |.-- 4-turns/helix
    |    | | |  | ||.-- 5-turns/helix
    |    | | |  | |||.-- geometrical bend
    |    | | |  | ||||.-- chirality
    |    | | |  | |||||.-- beta bridge label
    |    | | |  | ||||||.-- beta bridge label
    |    | | |  | |||||||   .-- beta bridge partner resnum
    |    | | |  | |||||||   |   .-- beta bridge partner resnum
    |    | | |  | |||||||   |   |.-- beta sheet label
    |    | | |  | |||||||   |   ||   .-- solvent accessibility
    |    | | |  | |||||||   |   ||   |
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC
    |    | | |  | |||||||   |   ||   |
   35   47 A I  E     +     0   0    2
   36   48 A R  E >  S- K   0  39C  97
   37   49 A Q  T 3  S+     0   0   86
   38   50 A N  T 3  S+     0   0   34
   39   51 A W  E <   -KL  36  98C   6
   ```



## TMAlign

__Installation__

1. Get [TMAlign.cpp](https://zhanggroup.org/TM-align/TMalign.cpp) from [ZhangLab](https://zhanggroup.org/TM-align)
2. compile with: `g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp`, this will serve as an input to the tmalign.py file. _Note:_ some machines may not support the static flag, feel free to remove. 

_Note:_ this is available on `ash.cels.anl.gov` at `/home/khippe/github/tm-align/TMalign`

__Usage__

```
./TMAlign struct1.pdb struct2.pdb
```

_Note:_ It appears to work on mmCIF files, but does not advertise support for them.

__Output__
Defaults to output to standard out. Look for the line with `TMAlign = [float]` and will show you two scores, the first is a alignment normalized to the length of the first structure, the second is the alignment score normalized to the second structure.



# Goals

- Single protein/gene
	- [X] DSSP for secondary structure
	- [X] PDB header extraction - syntactically meaningful metadata extraction __Update:__ If we are only using computational predictions the headers will not be meaningful
	- [ ] point cloud stats(??, spread?)
	- [ ] structural motifs (TIM barrell, Greek key)
	- [X] structure bio-chemical stats (aa-distribution, isoelectric points, solvent accessibility, etc) __Update:__ could do more stats, but basic ones are implemented
		- [ ] Radius of Gyration
		- [ ] number/type of inter-residue contacts
	- [ ] binding domain, if applicable (this might be hard)
	- [ ] (some tool, leaving in for edit) for identifying poorly folded regions
- Two protein/genes (includes all above)
	- [x] TM-align for pairwise structural similarity
- Many protein/genes (includes all above)
	- [ ] Ideally this would be a MSA for structure that shows the relationship between them all but this is a full scale project I think


__TODO__

- [ ] generalized argument parser
- [ ] generalized output formatter
