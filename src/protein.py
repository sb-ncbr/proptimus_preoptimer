import json

import numpy as np
from openmm import CustomExternalForce, LocalEnergyMinimizer, Platform, VerletIntegrator, unit
from openmm.app import PDBFile, Simulation

from . import PreoptimerError
from .amber import VALENCE, create_amber_system, set_valence_parameters

AMBER_TOLERANCE = 10  # kJ/mol/nm, RMS force at which the AMBER minimization stops
# kJ/mol/nm^2, harmonic restraints of the heavy atoms to their original positions in the AMBER minimization;
# they move by ~0.1-0.2 Å (max ~1-3 Å, charged side chain ends), which sets the valence geometry but keeps
# the conformation; weaker restraints let the charged side chains move more (from ~1e3) and form salt bridges
# in vacuum (<= 3e2)
AMBER_RESTRAINT_FORCE_CONSTANT = 2e3


class Protein:
    def __init__(self, pdb_file):
        print(f"Loading protein from {pdb_file}... ", end="", flush=True)
        with open(pdb_file, "r", encoding="utf-8") as f:
            self.pdb_lines = f.readlines()
        self.pdb = PDBFile(str(pdb_file))
        self.positions = self.pdb.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        print("ok")

    def preoptimize(self, parameters_file, only_hydrogens=False):
        """Short AMBER minimization by OpenMM in vacuum, with the equilibrium bond lengths and bond angles replaced
        by the valence parameters. Only hydrogens are moved if only_hydrogens is set, otherwise all atoms except CA,
        with the heavy atoms harmonically restrained to their original positions. Fixed atoms are made massless,
        which the minimizer keeps in place."""
        with open(parameters_file, "r", encoding="utf-8") as f:
            parameters = json.load(f)
        if VALENCE not in parameters:
            raise PreoptimerError(f"the parameters contain no valence parameters ({VALENCE}), "
                                  "calculate the parameters again")
        print("Running short AMBER optimization... ", end="", flush=True)
        topology = self.pdb.topology
        system, template_names = create_amber_system(self.pdb)
        replaced, total = set_valence_parameters(system, topology, template_names, parameters[VALENCE])
        if not only_hydrogens:
            restraints = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            restraints.addGlobalParameter("k", AMBER_RESTRAINT_FORCE_CONSTANT)
            for name in ("x0", "y0", "z0"):
                restraints.addPerParticleParameter(name)
            for omm_atom, position in zip(topology.atoms(), self.pdb.positions):
                if omm_atom.element is not None and omm_atom.element.symbol != "H":
                    restraints.addParticle(omm_atom.index, position.value_in_unit(unit.nanometer))
            system.addForce(restraints)
        for omm_atom in topology.atoms():
            if only_hydrogens:
                fixed = omm_atom.element is None or omm_atom.element.symbol != "H"
            else:
                fixed = omm_atom.name == "CA"
            if fixed:
                system.setParticleMass(omm_atom.index, 0)
        simulation = Simulation(topology, system, VerletIntegrator(0.001), Platform.getPlatformByName("CPU"))
        simulation.context.setPositions(self.pdb.positions)
        LocalEnergyMinimizer.minimize(simulation.context, tolerance=AMBER_TOLERANCE)
        positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        self.positions = positions.value_in_unit(unit.angstrom)
        print(f"ok ({replaced} of {total} bonds and angles parameterized)")

    def save(self, output_PDB_file):
        print(f"Saving preoptimized protein to {output_PDB_file}... ", end="", flush=True)
        with open(output_PDB_file, "w", encoding="utf-8") as f:
            f.writelines(self._output_lines())
        print("ok\n")

    def _output_lines(self):
        """Lines of the input PDB file with the coordinates of the atoms replaced by the current positions.
        OpenMM keeps the atoms in the order of the file, so the atom records are paired with the OpenMM atoms
        in order; the pairing is checked by the original coordinates. Of a multi-model file, only the first model
        (read by OpenMM) is written, of the alternate locations of an atom only the first one (used by OpenMM).
        ANISOU records are left out, as they do not correspond to the new positions."""
        original_positions = self.pdb.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        atoms = iter(zip(original_positions, self.positions))
        alternate_locations = set()  # atoms (name and residue columns) with an already written alternate location
        models = 0
        skip = False  # inside a model other than the first one
        lines = []
        for line in self.pdb_lines:
            record = line[:6].strip()
            if record == "MODEL":
                models += 1
                skip = models > 1
            if skip:
                skip = record != "ENDMDL"
                continue
            if record == "ANISOU":
                continue
            if record in ("ATOM", "HETATM"):
                if line[16] != " ":
                    atom_id = (line[12:16], line[17:27])
                    if atom_id in alternate_locations:
                        continue
                    alternate_locations.add(atom_id)
                    line = f"{line[:16]} {line[17:54]}  1.00{line[60:]}"
                original_position, position = next(atoms, (None, None))
                line_position = [float(line[i:i + 8]) for i in (30, 38, 46)]
                if original_position is None or not np.allclose(original_position, line_position, atol=1e-3):
                    raise PreoptimerError(f"atom record cannot be paired with OpenMM atom: {line.rstrip()}")
                line = f"{line[:30]}{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}{line[54:]}"
            lines.append(line)
        if next(atoms, None) is not None:
            raise PreoptimerError("OpenMM structure contains atoms missing in the PDB file")
        return lines
