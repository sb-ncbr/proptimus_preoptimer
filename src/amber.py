"""AMBER system (OpenMM) of a protein and its valence parameters.

The bond lengths and bond angles of the AMBER system (all bonds and angles of the topology, including the peptide
bonds, the closures of the rings and the disulfide bridges) can be replaced by the valence parameters calculated
from a set of reference structures, so that the AMBER minimization reproduces the geometry of the reference method.

A valence term (bond or angle) is identified by a key composed of the tokens of its atoms "<residue><offset>:<atom>".
The residue is the AMBER template name (e.g. ALA, NMET, CYX, HIE, ASH) or the residue name, the offset is
relative to the home residue of the term (the residue of the central atom of an angle, the first residue
of a bond): empty for the home residue, +1/-1 for the sequence neighbours and x for other residues
(disulfide bridges). The tokens are ordered so that the key does not depend on the direction of the term.
Each term has keys of three levels, from the most to the least specific:
    1. all residues by their template names, e.g. ALA:C|PRO+1:N
    2. the home residue by its template name, other residues "*", e.g. *+1:N|ALA:C
    3. the home residue by its residue name marked by "~" (any template), other residues "*", e.g. *+1:N|~ALA:C
Keys of the levels which do not differ (e.g. a term within one residue on levels 1 and 2) are identical.
The most specific key with parameters is used, terms without parameters keep their AMBER values.

The AMBER minimization does not reproduce the equilibrium values exactly, the nonbonded terms and torsions shift
them. The force constants of the terms with valence parameters are therefore increased by VALENCE_FORCE_SCALE.
"""

from collections import defaultdict

import numpy as np
from openmm import HarmonicAngleForce, HarmonicBondForce, unit
from openmm.app import CutoffNonPeriodic, ForceField, PDBFile

AMBER_FORCE_FIELD = ("amber19-all.xml", "implicit/gbn2.xml")  # AMBER19 with GBn2 implicit water
# AMBER19 in vacuum, used by the AMBER minimization of the preoptimization; GBn2 makes it ~10x slower
# and the restrained minimization does not need the solvent
AMBER_VACUUM_FORCE_FIELD = ("amber19-all.xml",)
AMBER_CUTOFF = 1.0  # nm, cutoff of nonbonded interactions in the AMBER minimization
# Protonation variants of AMBER residue templates, used for residues without a template (e.g. TYR-, neutral ARG)
AMBER_VARIANTS = {"ASP": ("ASP", "ASH"), "GLU": ("GLU", "GLH"), "CYS": ("CYS", "CYM", "CYX"), "LYS": ("LYS", "LYN"),
                  "HIS": ("HIE", "HID", "HIP")}
AMBER_BACKBONE_ATOMS = {"N", "H", "H1", "H2", "H3", "CA", "HA", "C", "O", "OXT"}
VALENCE = "VALENCE"  # key of the valence parameters in the parameters file
KEY_LEVELS = 3
# The force constants of the bonds and angles with valence parameters are multiplied by it. With 5, the valence
# geometry of the minimized structures is close to the parameters, 20 brings almost no further improvement.
VALENCE_FORCE_SCALE = 5


def create_amber_system(pdb_file, force_field=AMBER_FORCE_FIELD):
    """Return the OpenMM PDBFile, the AMBER system without constraints and the template name of each residue."""
    pdb = PDBFile(str(pdb_file))
    forcefield = ForceField(*force_field)
    _add_templates_for_unmatched_residues(forcefield, pdb.topology)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                                     nonbondedCutoff=AMBER_CUTOFF * unit.nanometer, constraints=None)
    template_names = [template.name for template in forcefield.getMatchingTemplates(pdb.topology)]
    return pdb, system, template_names


def _add_templates_for_unmatched_residues(forcefield, topology):
    """AMBER has no templates for some protonation states (e.g. TYR-, neutral ARG, protonated C-terminal GLU).
    For such residues, a template is created whose atoms take the atom types and charges of the same-named
    atoms of a reference template: backbone atoms from the standard template with the terminal prefix
    (e.g. CGLU), side chain atoms from the protonation variant with the most similar side chain (e.g. GLH).
    The charges are then shifted uniformly so that the residue charge is the charge of the standard template
    minus the number of missing and plus the number of extra hydrogens. The template is named after
    the standard template and the missing and extra hydrogens, e.g. TYR-HH or CGLU+HE2."""
    templates, residues = forcefield.generateTemplatesForUnmatchedResidues(topology)
    for template, residue in zip(templates, residues):
        chain_residues = list(residue.chain.residues())
        prefix = "N" if residue is chain_residues[0] else "C" if residue is chain_residues[-1] else ""
        variants = [forcefield._templates[name] for name in AMBER_VARIANTS.get(residue.name, (residue.name,))
                    if name in forcefield._templates]
        standard = forcefield._templates.get(prefix + residue.name)
        if standard is None or not variants:
            exit(f"error AMBER optimization: no template for residue {residue.name} {residue.id}")
        side_chain = {a.name for a in template.atoms} - AMBER_BACKBONE_ATOMS
        variants.sort(key=lambda variant: len(side_chain ^ ({a.name for a in variant.atoms}
                                                             - AMBER_BACKBONE_ATOMS)))

        for atom in template.atoms:
            references = [standard] + variants if atom.name in AMBER_BACKBONE_ATOMS else variants + [standard]
            reference_atom = next((a for reference in references for a in reference.atoms
                                   if a.name == atom.name), None)
            if reference_atom is None:
                exit(f"error AMBER optimization: no atom type for atom {atom.name} "
                     f"of residue {residue.name} {residue.id}")
            atom.type = reference_atom.type
            atom.parameters = dict(reference_atom.parameters)

        standard_hydrogens = {a.name for a in standard.atoms if a.element.symbol == "H"}
        hydrogens = {a.name for a in template.atoms if a.element.symbol == "H"}
        charge = (sum(a.parameters["charge"] for a in standard.atoms)
                  - len(standard_hydrogens - hydrogens) + len(hydrogens - standard_hydrogens))
        shift = (round(charge) - sum(a.parameters["charge"] for a in template.atoms)) / len(template.atoms)
        for atom in template.atoms:
            atom.parameters["charge"] += shift
        name = (standard.name + "".join(f"-{h}" for h in sorted(standard_hydrogens - hydrogens))
                + "".join(f"+{h}" for h in sorted(hydrogens - standard_hydrogens)))
        while name in forcefield._templates:  # the same hydrogens, but different heavy atoms
            name += "'"
        template.name = name
        forcefield.registerResidueTemplate(template)


def valence_key(tokens, level):
    """Key of a valence term from the tokens (template name, residue name, offset, atom name) of its atoms,
    see the module docstring. The offset is 0 for the home residue, +1/-1 for the sequence neighbours and None
    for other residues. The key does not depend on the direction in which the atoms are listed."""
    def token(template_name, residue_name, offset, atom_name):
        if offset == 0:
            residue = template_name if level < 3 else f"~{residue_name}"
            return f"{residue}:{atom_name}"
        residue = template_name if level == 1 else "*"
        return f"{residue}{'x' if offset is None else f'{offset:+d}'}:{atom_name}"
    forward = [token(*t) for t in tokens]
    return "|".join(min(forward, forward[::-1]))


def _valence_terms(system, topology, template_names):
    """Yield (force, index, kind, key tokens, atom indices) of all bonds and angles of the AMBER system."""
    atoms = list(topology.atoms())
    chain_position = {}  # residue index -> (chain index, position in the chain)
    for chain in topology.chains():
        for position, residue in enumerate(chain.residues()):
            chain_position[residue.index] = (chain.index, position)

    def tokens(indices, home):
        result = []
        for i in indices:
            residue = atoms[i].residue
            (chain, position), (home_chain, home_position) = chain_position[residue.index], chain_position[home]
            offset = position - home_position if chain == home_chain and abs(position - home_position) <= 1 \
                else None
            result.append((template_names[residue.index], residue.name, offset, atoms[i].name))
        return result

    for force in system.getForces():
        if isinstance(force, HarmonicBondForce):
            for index in range(force.getNumBonds()):
                indices = force.getBondParameters(index)[:2]
                home = min(atoms[i].residue.index for i in indices)
                yield force, index, "bonds", tokens(indices, home), indices
        elif isinstance(force, HarmonicAngleForce):
            for index in range(force.getNumAngles()):
                indices = force.getAngleParameters(index)[:3]
                yield force, index, "angles", tokens(indices, atoms[indices[1]].residue.index), indices


def measure_valence_terms(pdb_file):
    """Bond lengths (Å) and bond angles (°) of all valence terms of the structure, as
    {"bonds"/"angles": {key: [values]}} with the keys of all levels."""
    pdb, system, template_names = create_amber_system(pdb_file)
    positions = pdb.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    values = {"bonds": defaultdict(list), "angles": defaultdict(list)}
    for _, _, kind, tokens, indices in _valence_terms(system, pdb.topology, template_names):
        if kind == "bonds":
            value = float(np.linalg.norm(positions[indices[0]] - positions[indices[1]]))
        else:
            value = angle(*(positions[i] for i in indices))
        for key in {valence_key(tokens, level) for level in range(1, KEY_LEVELS + 1)}:
            values[kind][key].append(value)
    return values


def set_valence_parameters(system, topology, template_names, valence_parameters, force_scale=VALENCE_FORCE_SCALE):
    """Replace the equilibrium bond lengths and bond angles of the system by the valence parameters
    (the most specific key with parameters is used) and multiply their force constants by force_scale.
    Returns the number of replaced and of all terms."""
    replaced = total = 0
    for force, index, kind, tokens, indices in _valence_terms(system, topology, template_names):
        total += 1
        value = next((valence_parameters[kind][key][0] for key in
                      (valence_key(tokens, level) for level in range(1, KEY_LEVELS + 1))
                      if key in valence_parameters[kind]), None)
        if value is None:
            continue
        replaced += 1
        if kind == "bonds":
            *_, k = force.getBondParameters(index)
            force.setBondParameters(index, *indices, value * unit.angstrom, k * force_scale)
        else:
            *_, k = force.getAngleParameters(index)
            force.setAngleParameters(index, *indices, np.radians(value) * unit.radian, k * force_scale)
    return replaced, total


def angle(a1, a2, a3):
    """Bond angle a1-a2-a3 in degrees."""
    v1, v2 = a1 - a2, a3 - a2
    return float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))))
