# PROPTIMUS PREOPTIMER

## Description

PROPTIMUS PREOPTIMER is a parameterizable tool for setting bond lengths and bond angles in protein structures. 
In the `preoptimize` mode, the bond lengths and bond angles are set by a short AMBER minimization with the equilibrium values taken from the parameters. 
Optionally, only the positions of hydrogen atoms can be preoptimized while all other atoms are kept in place. 
In the `parameterize` mode, the parameters are calculated from a set of protein structures. 

The `parameters` directory currently contains ready-to-use parameters derived from structures optimized with the [GFN-FF](https://onlinelibrary.wiley.com/doi/10.1002/anie.202004239) force field (`parameters/GFN-FF_parameters.json`).

## Getting Started
Before running PROPTIMUS PREOPTIMER for the first time, you need to create and activate a Python virtual environment where the script will always be executed. Python 3.11 or newer is required. You can set up the virtual environment in two steps:

**1. Clone GitHub repository:**

```
git clone https://github.com/sb-ncbr/proptimus_preoptimer
cd proptimus_preoptimer
```

**2. Set up a virtual environment**

* Create a virtual environment named `.proptimus_preoptimer_venv`:

```
python3 -m venv .proptimus_preoptimer_venv
```

* Activate the virtual environment:

```
source .proptimus_preoptimer_venv/bin/activate
```

* Install the dependencies:

```
pip install -r requirements.txt
```

Once the virtual environment is created and activated, no further installation is required, and PROPTIMUS PREOPTIMER is ready to run.

## Executing the calculation

### Required arguments

`--mode {preoptimize,parameterize}`          Mode `preoptimize` sets bond lengths and bond angles according to the parameters. Mode `parameterize` creates these parameters from a set of protein structures.

`--parameters FILE`              JSON file with parameters. In mode `preoptimize`, the parameters are read from it (available parameters are in the `parameters` directory). In mode `parameterize`, the calculated parameters are saved to it; the file must not exist.

### Arguments for `--mode preoptimize`

`--input_PDB_file FILE`          PDB file with structure, which should be preoptimized (required).

`--output_PDB_file FILE`         PDB file to store preoptimized structure. The file must not exist (required).

`--only_hydrogens`               Preoptimize only hydrogen atoms, positions of all other atoms are kept.

### Arguments for `--mode parameterize`

`--proteins_directory DIR`       Directory with protein structures (`*.pdb` files) used to calculate the parameters (required).

`--workers N`                    Number of worker processes (default: number of CPUs).

## Preoptimization

In the `preoptimize` mode, PROPTIMUS PREOPTIMER sets the bond lengths and bond angles of a protein structure according to the parameters. It is intended as a fast preparation step: it adjusts local geometry before a more expensive optimization, for example after adding hydrogen atoms or building missing atoms. The input structure must be protonated and use standard PDB atom names (see [Parameterization](#parameterization)). The output file is the input file with the new coordinates (all other records are kept, except ANISOU). Only the first model of the PDB file is preoptimized and saved; of alternate locations of an atom, only the first one is kept.

The preoptimization is a short AMBER minimization with [OpenMM](https://openmm.org/) (AMBER19 force field in vacuum), in which the equilibrium bond lengths and bond angles of the AMBER force field are replaced by the valence parameters (see [Valence parameters](#valence-parameters)) and their force constants are multiplied by 5 (`VALENCE_FORCE_SCALE` in `src/amber.py`). It optimizes all atoms except the CA atoms, with the heavy atoms harmonically restrained to their original positions (2·10<sup>3</sup> kJ/mol/nm<sup>2</sup>, `AMBER_RESTRAINT_FORCE_CONSTANT` in `src/protein.py`; they move by ~0.1–0.2 Å, at most ~1–3 Å at the ends of charged side chains, so the valence geometry is set but the conformation is kept), or only the hydrogen atoms if `--only_hydrogens` is used, and runs until the RMS force drops below 10 kJ/mol/nm. As all bonds and angles of the topology are included (including those of hydrogens), it also sets the peptide bonds, the closures of the rings and the disulfide bridges, and it resolves overlapping atoms and steric clashes. The minimized structure is the result of the preoptimization; the dihedral angles are given by the AMBER torsions and nonbonded interactions.

With `--only_hydrogens`, only hydrogen atoms are moved. All heavy atoms stay at their original positions and the hydrogens are optimized by the AMBER minimization described above. This mode is useful when the heavy atoms come from an experiment or from a previous optimization and only the hydrogens need to be adjusted.

### How to run

Preoptimization of a structure with the GFN-FF parameters from the `parameters` directory:

```bash
python preoptimer.py --mode preoptimize --input_PDB_file structure.pdb --output_PDB_file structure_preoptimized.pdb --parameters parameters/GFN-FF_parameters.json
```

Preoptimization of hydrogen atoms only with the GFN-FF parameters (all other atoms are kept in place):

```bash
python preoptimer.py --mode preoptimize --input_PDB_file structure.pdb --output_PDB_file structure_preoptimized.pdb --parameters parameters/GFN-FF_parameters.json --only_hydrogens
```

## Parameterization

In the `parameterize` mode, PROPTIMUS PREOPTIMER measures all bonds and bond angles of the AMBER topology of a set of protein structures. Their averages give the valence parameters (see [Valence parameters](#valence-parameters)), which are used by the AMBER minimization of the `preoptimize` mode.

### Input structures

The parameters reproduce the geometry of the structures they were calculated from, so the input set has to be consistent:

* All structures in `--proteins_directory` must be optimized **with the same method** (e.g. the same force field or the same quantum chemical method). Mixing structures optimized by different methods produces parameters that correspond to none of them.
* The structures must be protonated (hydrogen atoms are parameterized as well) and use standard PDB atom names (`HB2`, `HB3`, `HD21`, ...).
* All `*.pdb` files in the directory are used (the extension is case-insensitive). Structures which cannot be processed (e.g. residues without an AMBER template, such as ligands) are skipped and listed after the calculation.
* Keys with fewer than 20 values (`MIN_VALENCE_SAMPLES` in `src/parameterization.py`) are not stored and the less specific keys of the same terms are used instead, so rare residue variants (e.g. deprotonated tyrosine or terminal residues) get specific parameters only if they occur often enough.

### How to run

Calculation of new parameters from a set of structures using 8 worker processes:

```bash
python preoptimer.py --mode parameterize --proteins_directory my_structures --parameters my_parameters.json --workers 8
```

The parameters are saved to `my_parameters.json` and can be used for preoptimization:

```bash
python preoptimer.py --mode preoptimize --input_PDB_file structure.pdb --output_PDB_file structure_preoptimized.pdb --parameters my_parameters.json
```

### Valence parameters

The valence parameters contain the equilibrium values of all bonds (`bonds`, Å) and bond angles (`angles`, °) of the AMBER topology, including the bonds and angles between residues (peptide bonds, disulfide bridges). Each value is stored as `[mean, standard deviation, number of values]` under a key composed of the atoms of the term, `<residue><offset>:<atom>` separated by `|`:

```json
"VALENCE": {
    "bonds": {
        "ALA:C|PRO+1:N": [1.3901, 0.0067, 507],
        "*+1:N|ALA:C": [1.3777, 0.0054, 12879],
        "CYX:SG|CYXx:SG": [2.0606, 0.0144, 155]
    },
    "angles": {
        "ALA:C|ALA:CA|ALA:N": [110.2403, 1.1940, 12879]
    }
}
```

The residue is the name of the AMBER residue template (it distinguishes the termini, e.g. `NMET`, and the protonation states, e.g. `ASH`, `HIE`, `CYX`; residues without an AMBER template get names like `TYR-HH`), the offset is relative to the home residue of the term (the residue of the central atom of an angle, the first residue of a bond): empty for the home residue, `+1`/`-1` for the sequence neighbours and `x` for other residues. Each term is stored under keys of three levels: with all residues specified (`ALA:C|PRO+1:N`), with the other residues replaced by `*` (`*+1:N|ALA:C`), and in addition with the home residue given by its residue name regardless of the template (`*+1:N|~ALA:C`). The preoptimization uses the most specific key available. Keys with too few values are not stored (see [Input structures](#input-structures)), and terms without any key keep their AMBER values.

## Testing and use cases

The results of testing and examples of use cases are available in the [GitHub wiki](https://github.com/sb-ncbr/proptimus_preoptimer/wiki) of this repository.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
