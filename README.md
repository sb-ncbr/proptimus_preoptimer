# PROPTIMUS PREOPTIMER

## Description

PROPTIMUS PREOPTIMER is a parameterizable tool for setting bond lengths, bond angles and dihedral angles in protein structures. 
In the `preoptimize` mode, the internal coordinates of the atoms are set according to the parameters. 
Optionally, only the positions of hydrogen atoms can be preoptimized while all other atoms are kept in place. 
In the `parameterize` mode, the parameters are calculated from a set of protein structures. 

The `parameters` directory currently contains parameters derived from structures optimised with the [GFN-FF](https://onlinelibrary.wiley.com/doi/10.1002/anie.202004239) force field.

## Getting Started
Before running PROPTIMUS PREOPTIMER for the first time, you need to create and activate a Python virtual environment where the script will always be executed. Python 3.11 or newer is required. You can set up the virtual environment in two steps:

**1. Clone GitHub repository:**

```
git clone https://github.com/sb-ncbr/proptimus_preoptimer
cd proptimus_preoptimer
```

**2. Set up a virtual environment**

* Create a virtual environment named `.venv`:

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

`--mode {preoptimize,parameterize}`          Mode `preoptimize` sets bond lengths, bond angles and dihedral angles according to the parameters. Mode `parameterize` creates these parameters from a set of protein structures.

### Arguments for `--mode preoptimize`

`--input_PDB_file FILE`          PDB file with structure, which should be preoptimised (required).

`--output_PDB_file FILE`         PDB file to store preoptimized structure. The file must not exist (required).

`--parameters FILE`              JSON file with parameters. Available parameters are in the `parameters` directory (required).

`--only_hydrogens`               Preoptimize only hydrogen atoms, positions of all other atoms are kept.

`--noopt`                        Skip the short AMBER optimization (OpenMM) that is run before the preoptimization.

### Arguments for `--mode parameterize`

`--proteins_directory DIR`       Directory with protein structures (`*.pdb` files) used to calculate the parameters (required).

`--data_dir DIR`                 Directory to store the results. The directory must not exist. The calculated parameters are saved to `<data_dir>/parameters.json` and HTML graphs of all atom types of all residues to the `<data_dir>/graphs` subdirectory (required).

`--workers N`                    Number of worker processes (default: number of CPUs).

`--open_browser`                 Open the HTML graphs in the default web browser when they are ready.

## Preoptimization

In the `preoptimize` mode, PROPTIMUS PREOPTIMER sets the bond lengths, bond angles and dihedral angles of a protein structure according to the parameters. It is intended as a fast preparation step: it adjusts local geometry before a more expensive optimization, for example after adding hydrogen atoms or building missing atoms. The input structure must be protonated and use standard PDB atom names, and the parameters must contain all atom types present in the structure (see [Parameterization](#parameterization)).

Ready-to-use parameters derived from structures optimized with the [GFN-FF](https://onlinelibrary.wiley.com/doi/10.1002/anie.202004239) force field are available in the `parameters` directory as `parameters/GFN-FF_parameters.json`.

Before the preoptimization itself, a short AMBER optimization is run with [OpenMM](https://openmm.org/) (AMBER19 force field with GBn2 implicit water). It optimizes all atoms except the CA atoms, or only the hydrogen atoms if `--only_hydrogens` is used. The minimization runs until the RMS force drops below 10 kJ/mol/nm. This step resolves overlapping atoms and steric clashes, which the internal coordinates alone cannot fix. It can be skipped with `--noopt`.

With `--only_hydrogens`, only hydrogen atoms are placed. All heavy atoms stay at their original positions. Each hydrogen is placed from its ancestors with the parameterized bond length and bond angle, and with the closest dihedral cluster center (or its current dihedral angle if the dihedral is not parameterized). This mode is useful when the heavy atoms come from an experiment or from a previous optimization and only the hydrogens need to be adjusted.

### How to run

Preoptimization of a structure with the GFN-FF parameters from the `parameters` directory:

```bash
python preoptimer.py --mode preoptimize --input_PDB_file structure.pdb --output_PDB_file structure_preoptimized.pdb --parameters parameters/GFN-FF_parameters.json
```

Preoptimization of hydrogen atoms only with the GFN-FF parameters (all other atoms are kept in place):

```bash
python preoptimer.py --mode preoptimize --input_PDB_file structure.pdb --output_PDB_file structure_preoptimized.pdb --parameters parameters/GFN-FF_parameters.json --only_hydrogens
```

Preoptimization without the preceding AMBER optimization:

```bash
python preoptimer.py --mode preoptimize --input_PDB_file structure.pdb --output_PDB_file structure_preoptimized.pdb --parameters parameters/GFN-FF_parameters.json --noopt
```

## Parameterization

In the `parameterize` mode, PROPTIMUS PREOPTIMER measures the bond lengths, bond angles and dihedral angles of all atoms in a set of protein structures and clusters them per atom type. The resulting cluster centers form the parameters used by the `preoptimize` mode.

### Input structures

The parameters reproduce the geometry of the structures they were calculated from, so the input set has to be consistent:

* All structures in `--proteins_directory` must be optimized **with the same method** (e.g. the same force field or the same quantum chemical method). Mixing structures optimized by different methods produces parameters that correspond to none of them.
* The structures must be protonated (hydrogen atoms are parameterized as well) and use standard PDB atom names (`HB2`, `HB3`, `HD21`, ...).
* All `*.pdb` files in the directory are used.
* Every atom type must occur at least 100 times in the whole set (`MIN_SAMPLES` in `src/parameterization.py`). Otherwise, the calculation ends with an error listing the atom types with too few values, and more structures have to be added. Keep in mind that rare atom types (e.g. deprotonated tyrosine or protonated histidine) need enough occurrences too.

### How to run

Calculation of new parameters from a set of structures using 8 worker processes:

```bash
python preoptimer.py --mode parameterize --proteins_directory my_structures --data_dir my_parameterization --workers 8
```

The parameters are then saved to `my_parameterization/parameters.json` and can be used for preoptimization:

```bash
python preoptimer.py --mode preoptimize --input_PDB_file structure.pdb --output_PDB_file structure_preoptimized.pdb --parameters my_parameterization/parameters.json
```

### Results

All results are stored in the directory given by `--data_dir` (it must not exist, it is created by the program):

```
<data_dir>/
├── parameters.json        calculated parameters
└── graphs/
    ├── index.html         links to all graphs, grouped by residue
    ├── <RES>_<TYPE>.html  one graph per atom type, e.g. ALA_HB.html, LYS_NZ_+1.html
    └── plotly.min.js      shared by all graphs
```

* `parameters.json` can be used directly as the `--parameters` file of the `preoptimize` mode. For every residue and atom type, it contains the cluster centers of `bond_length` (Å), `bond_angle` (°) and `dihedral_angle` (°), e.g.:

  ```json
  "ALA": {
      "HB": {
          "bond_length": [1.0898],
          "bond_angle": [110.8533],
          "dihedral_angle": [61.9776, -178.5177, -57.2619]
      }
  }
  ```

* The HTML graphs show, for every atom type, the distributions of bond lengths and bond angles (histogram with a rug plot) and of dihedral angles (polar histogram). The calculated cluster centers are drawn as red lines. Hovering over the rug plot shows which structure, residue and atom a value belongs to, which is useful for checking the input data. Use `--open_browser` to open `index.html` automatically.
* An outlier analysis (atoms with a bond length or bond angle far from the rest of the distribution) is printed to the standard output at the end of the calculation.

### Custom atom types

Atom types are assigned in the constructor of the `Protein` class (`src/protein.py`) by a sequence of classification functions:

```python
self.atom_types = {}  # atom -> atom type, only for atoms which are processed
self._assign_basic_atomic_types()
self._assign_histidine_types()
self._assign_formal_charges()
```

* `_assign_basic_atomic_types` assigns the basic atom type defined in `src/data.py` (e.g. `HB` for `HB2` and `HB3` of alanine).
* `_assign_histidine_types` and `_assign_formal_charges` refine the basic types by appending a suffix separated by `/`, e.g. `CE1/HID` for the tautomer of histidine or `NZ/+1` for protonated lysine.

To distinguish other atom types (e.g. atoms in a disulfide bridge or a specific environment), write another classification function and call it in `__init__` after `_assign_basic_atomic_types`. The function modifies `self.atom_types`, typically by appending a new suffix using the `_add_suffix(residue, atom_names, suffix)` helper:

```python
def _assign_my_types(self):
    for residue in self.residues:
        if <condition>:
            self._add_suffix(residue, ("CB", "HB2", "HB3"), "MY")
```

Nothing else has to be changed. The part before the first `/` (the basic type) is used to look up the ancestors and the number of clusters in `src/data.py`, while the full type (including suffixes) is used as the key in the parameters. Each new type is therefore parameterized separately and gets its own graph. Note that the parameters have to be recalculated afterwards, because the `preoptimize` mode needs parameters for all atom types assigned to the structure.

### Metrics and number of clusters

The definition of the basic atom types is in `src/data.py`. Each atom type is created by the `_type` function:

```python
"HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
```

* `atoms` – names of the atoms belonging to this atom type.
* `ancestors` – the three atoms `a2 a3 a4` that define the measured metrics: the bond length `atom-a2`, the bond angle `atom-a2-a3` and the dihedral angle `atom-a2-a3-a4`. A suffix `/-1` or `/+1` marks an atom of the previous or the next residue, and `x` is replaced by the third character of the atom name (e.g. `CDx` means `CD1` for `HD11`).
* `n_clusters` – number of clusters for `(bond length, bond angle, dihedral angle)`. `0` means that the metric is not parameterized: it is not stored in the parameters and its current value is kept during preoptimization. The number of dihedral clusters corresponds to the number of preferred conformations (e.g. `3` for a methyl group, `2` for the two hydrogens of a CH2 group). During preoptimization, the dihedral angle is set to the cluster center closest to its current value, while the bond length and bond angle are always set to the first cluster center, so they should be kept at `1`.
* `custom` – atoms with a single chemical bond, which are placed directly from their ancestors (outside of the Biopython internal coordinates).
* `reverse_order` – the bond angle is stored as `atom:a2:a3` instead of `a3:a2:atom` in the Biopython internal coordinates.

The backbone atom types (`HA`, `O`, `C`, `N`, `H`) are shared by all residues and are defined in the `_residue` function. Atoms listed in `unprocessed_atoms` (`CA`) are never typed and are kept in place.

The atoms of the terminal residues are typed as atoms of two pseudo-residues shared by all residues, because each protein has only one N- and one C-terminus: `NTERM` (`N` and the hydrogens `H`, `H1`, `H2`, `H3` of the first residue) and `CTERM` (`C`, `O` and `OXT` of the last residue). They have their own parameters and graphs (e.g. `NTERM_H.html`) and are placed after all other atoms during preoptimization. On an N-terminal proline, one of the three positions of the NH3+ hydrogens is taken by `CD`, so its two hydrogens are placed to the remaining tetrahedral positions (dihedral angle of `CD` ± 120°).

## Testing and use cases

The results of testing and examples of use cases are available in the [GitHub wiki](https://github.com/sb-ncbr/proptimus_preoptimer/wiki) of this repository.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
