"""Atom type definitions of the 20 standard amino acids.

Every atom type is described by
    atoms          names of the atoms of this type,
    ancestors      the three atoms defining the bond length (atom-a2), bond angle (atom-a2-a3)
                   and dihedral angle (atom-a2-a3-a4). "/-1" and "/+1" mark an atom in the previous /
                   next residue, "x" is replaced by the third character of the atom name
                   (e.g. "CDx" for HD1 -> CD1),
    n_clusters     number of clusters of (bond length, bond angle, dihedral angle); 0 = not parameterized,
    custom         atoms with a single chemical bond, placed using the closest cluster center of the dihedral,
    reverse order  the bond angle is stored as atom:a2:a3 instead of a3:a2:atom.

The atoms of the N- and C-terminus are typed as atoms of the pseudo-residues N_TERMINUS and C_TERMINUS,
shared by all residues, because each protein has only one N- and one C-terminal residue.
"""
from collections import defaultdict

unprocessed_atoms = ["CA"]
N_TERMINUS = "NTERM"
C_TERMINUS = "CTERM"


def _type(atoms, ancestors, n_clusters, custom=False, reverse_order=False):
    atom_type = {"atoms": set(atoms.split()), "ancestors": tuple(ancestors.split()), "n_clusters": n_clusters}
    if custom:
        atom_type["ancestors_type"] = "custom"
    if reverse_order:
        atom_type["reverse order"] = True
    return atom_type


def _residue(side_chain, ha="HA", ha_clusters=1, amide_h=True):
    """Backbone atom types shared by all residues, plus the side chain atom types."""
    backbone = {
        "HA": _type(ha, "CA C N", (1, 1, ha_clusters), custom=True, reverse_order=True),
        "O": _type("O", "C CA N/+1", (1, 1, 1), custom=True),
        "C": _type("C", "CA N/+1 CA/+1", (1, 1, 0), custom=True),
        "N": _type("N", "CA C/-1 CA/-1", (1, 1, 0), custom=True),
    }
    if amide_h:
        backbone["H"] = _type("H", "N CA C/-1", (1, 1, 0), reverse_order=True)
    return {**backbone, **side_chain}


atomic_types = {
    "GLY": _residue({
    }, ha="HA2 HA3", ha_clusters=2),
    "ALA": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB1 HB2 HB3", "CB CA N", (1, 1, 3)),
    }),
    "SER": _residue({
        "HG": _type("HG", "OG CB CA", (1, 1, 0)),
        "OG": _type("OG", "CB CA N", (1, 1, 0)),
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA OG", (1, 1, 2), custom=True),
    }),
    "THR": _residue({
        "HG1": _type("HG1", "OG1 CB CA", (1, 1, 0)),
        "OG1": _type("OG1", "CB CA N", (1, 1, 0)),
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "CG2": _type("CG2", "CB CA N", (1, 1, 0)),
        "HB": _type("HB", "CB CA CG2", (1, 1, 1), custom=True),
        "HG": _type("HG21 HG22 HG23", "CG2 CB CA", (1, 1, 3)),
    }),
    "VAL": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB", "CB CA CG1", (1, 1, 1), custom=True),
        "CG": _type("CG1 CG2", "CB CA N", (1, 1, 0)),
        "HG": _type("HG11 HG12 HG13 HG21 HG22 HG23", "CGx CB CA", (1, 1, 3)),
    }),
    "LEU": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "HG": _type("HG", "CG CB CD1", (1, 1, 1), custom=True),
        "CD": _type("CD1 CD2", "CG CB CA", (1, 1, 0)),
        "HD": _type("HD11 HD12 HD13 HD21 HD22 HD23", "CDx CG CB", (1, 1, 3)),
    }),
    "ILE": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB", "CB CA CG1", (1, 1, 1), custom=True),
        "CG1": _type("CG1", "CB CA N", (1, 1, 0)),
        "CG2": _type("CG2", "CB CA N", (1, 1, 0)),
        "HG1": _type("HG12 HG13", "CG1 CB CD1", (1, 1, 2), custom=True),
        "HG2": _type("HG21 HG22 HG23", "CG2 CB CA", (1, 1, 3)),
        "CD1": _type("CD1", "CG1 CB CA", (1, 1, 0)),
        "HD": _type("HD11 HD12 HD13", "CD1 CG1 CB", (1, 1, 3)),
    }),
    "PRO": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "CD": _type("CD", "CG CB CA", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "HG": _type("HG2 HG3", "CG CB CD", (1, 1, 2), custom=True),
        "HD": _type("HD2 HD3", "CD CG N", (1, 1, 2), custom=True),
    }, amide_h=False),
    "ASN": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "OD1": _type("OD1", "CG CB CA", (1, 1, 0)),
        "ND2": _type("ND2", "CG CB CA", (1, 1, 0)),
        "HD": _type("HD21 HD22", "ND2 CG CB", (1, 1, 0)),
    }),
    "GLN": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "HG": _type("HG2 HG3", "CG CB CD", (1, 1, 2), custom=True),
        "CD": _type("CD", "CG CB CA", (1, 1, 0)),
        "OE1": _type("OE1", "CD CG CB", (1, 1, 0)),
        "NE2": _type("NE2", "CD CG CB", (1, 1, 0)),
        "HE": _type("HE21 HE22", "NE2 CD CG", (1, 1, 0)),
    }),
    "CYS": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA SG", (1, 1, 2), custom=True),
        "SG": _type("SG", "CB CA N", (1, 1, 0)),
        "HG": _type("HG", "SG CB CA", (1, 1, 0)),
    }),
    "MET": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "HG": _type("HG2 HG3", "CG CB SD", (1, 1, 2), custom=True),
        "SD": _type("SD", "CG CB CA", (1, 1, 0)),
        "CE": _type("CE", "SD CG CB", (1, 1, 0)),
        "HE": _type("HE1 HE2 HE3", "CE SD CG", (1, 1, 0)),
    }),
    "PHE": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "CD": _type("CD1 CD2", "CG CB CA", (1, 1, 0)),
        "HD": _type("HD1 HD2", "CDx CG CB", (1, 1, 1)),
        "CE": _type("CE1 CE2", "CDx CG CB", (1, 1, 0)),
        "HE": _type("HE1 HE2", "CEx CDx CG", (1, 1, 1)),
        "CZ": _type("CZ", "CE1 CD1 CG", (1, 1, 0)),
        "HZ": _type("HZ", "CZ CE1 CD1", (1, 1, 1)),
    }),
    "TYR": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "CD": _type("CD1 CD2", "CG CB CA", (1, 1, 0)),
        "HD": _type("HD1 HD2", "CDx CG CB", (1, 1, 1)),
        "CE": _type("CE1 CE2", "CDx CG CB", (1, 1, 0)),
        "HE": _type("HE1 HE2", "CEx CDx CG", (1, 1, 1)),
        "CZ": _type("CZ", "CE1 CD1 CG", (1, 1, 0)),
        "OH": _type("OH", "CZ CE1 CD1", (1, 1, 1)),
        "HH": _type("HH", "OH CZ CE1", (1, 1, 0)),
    }),
    "TRP": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "CD1": _type("CD1", "CG CB CA", (1, 1, 0)),
        "CD2": _type("CD2", "CG CB CA", (1, 1, 0)),
        "HD1": _type("HD1", "CD1 CG CB", (1, 1, 1)),
        "NE1": _type("NE1", "CD1 CG CB", (1, 1, 0)),
        "HE1": _type("HE1", "NE1 CD1 CG", (1, 1, 1)),
        "CE2": _type("CE2", "CD2 CG CB", (1, 1, 0)),
        "CE3": _type("CE3", "CD2 CG CB", (1, 1, 0)),
        "HE3": _type("HE3", "CE3 CD2 CG", (1, 1, 1)),
        "CZ2": _type("CZ2", "CE2 CD2 CG", (1, 1, 0)),
        "CZ3": _type("CZ3", "CH2 CZ2 CE2", (1, 1, 0)),
        "HZ2": _type("HZ2", "CZ2 CE2 CD2", (1, 1, 1)),
        "HZ3": _type("HZ3", "CZ3 CH2 CZ2", (1, 1, 1)),
        "CH2": _type("CH2", "CZ2 CE2 CD2", (1, 1, 0)),
        "HH2": _type("HH2", "CH2 CZ2 CE2", (1, 1, 1)),
    }),
    "LYS": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "HG": _type("HG2 HG3", "CG CB CD", (1, 1, 2), custom=True),
        "CD": _type("CD", "CG CB CA", (1, 1, 0)),
        "HD": _type("HD2 HD3", "CD CG CE", (1, 1, 2), custom=True),
        "CE": _type("CE", "CD CG CB", (1, 1, 0)),
        "HE": _type("HE2 HE3", "CE CD NZ", (1, 1, 2), custom=True),
        "NZ": _type("NZ", "CE CD CG", (1, 1, 0)),
        "HZ": _type("HZ1 HZ2 HZ3", "NZ CE CD", (1, 1, 0)),
    }),
    "ARG": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "HG": _type("HG2 HG3", "CG CB CD", (1, 1, 2), custom=True),
        "CD": _type("CD", "CG CB CA", (1, 1, 0)),
        "HD": _type("HD2 HD3", "CD CG NE", (1, 1, 2), custom=True),
        "CZ": _type("CZ", "NE CD CG", (1, 1, 0)),
        "NE": _type("NE", "CD CG CB", (1, 1, 0)),
        "HE": _type("HE", "NE CD CG", (1, 1, 0)),
        "NH": _type("NH1 NH2", "CZ NE CD", (1, 1, 0)),
        "HH": _type("HH11 HH12 HH21 HH22", "NHx CZ NE", (1, 1, 0)),
    }),
    "ASP": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "OD1": _type("OD1", "CG CB CA", (1, 1, 0)),
        "OD2": _type("OD2", "CG CB CA", (1, 1, 0)),
        "HD2": _type("HD2", "OD2 CG CB", (1, 1, 0)),
    }),
    "GLU": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "HG": _type("HG2 HG3", "CG CB CD", (1, 1, 2), custom=True),
        "CD": _type("CD", "CG CB CA", (1, 1, 0)),
        "OE1": _type("OE1", "CD CG CB", (1, 1, 0)),
        "OE2": _type("OE2", "CD CG CB", (1, 1, 0)),
        "HE2": _type("HE2", "OE2 CD CG", (1, 1, 0)),
    }),
    "HIS": _residue({
        "CB": _type("CB", "CA N C", (1, 1, 0)),
        "HB": _type("HB2 HB3", "CB CA CG", (1, 1, 2), custom=True),
        "CG": _type("CG", "CB CA N", (1, 1, 0)),
        "CD2": _type("CD2", "CG CB CA", (1, 1, 0)),
        "HD2": _type("HD2", "CD2 CG CB", (1, 1, 1)),
        "ND1": _type("ND1", "CG CB CA", (1, 1, 0)),
        "HD1": _type("HD1", "ND1 CG CB", (1, 1, 0)),
        "CE1": _type("CE1", "ND1 CG CB", (1, 1, 0)),
        "NE2": _type("NE2", "CD2 CG CB", (1, 1, 0)),
        "HE1": _type("HE1", "CE1 ND1 CG", (1, 1, 1)),
        "HE2": _type("HE2", "NE2 CD2 CG", (1, 1, 0)),
    }),
    # terminal atoms, they are placed directly from their ancestors (after the other atoms)
    N_TERMINUS: {
        "N": _type("N", "CA C N/+1", (1, 1, 0), custom=True),
        "H": _type("H H1 H2 H3", "N CA C", (1, 1, 3), custom=True, reverse_order=True),
    },
    C_TERMINUS: {
        "C": _type("C", "CA N C/-1", (1, 1, 0), custom=True),
        "OXT": _type("OXT", "C CA N", (1, 1, 0), custom=True),
        "O": _type("O", "C CA OXT", (1, 1, 1), custom=True),
    },
}

# atom name -> atom type, for every residue
at_types_lookup = defaultdict(dict)
for res_name, res_at_types in atomic_types.items():
    for res_at_type, res_at_type_values in res_at_types.items():
        for atom_name in res_at_type_values["atoms"]:
            at_types_lookup[res_name][atom_name] = res_at_type
