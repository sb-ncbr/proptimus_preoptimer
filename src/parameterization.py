import json
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from sklearn.cluster import KMeans
from tqdm import tqdm

from .protein import Protein
from .data import atomic_types

MIN_SAMPLES = 100  # minimal number of values of an atom type needed for its parameterization


def _extract_protein(pdb_file):
    """Measure one protein; returns plain tuples, so the result can be sent from a worker process."""
    protein = Protein(pdb_file, output=False)
    return [(f"{pdb_file},{m.res_name}{m.res_id},{m.atom_name}", m) for m in protein.get_lengths_angles()]


def _run_kmeans(values, num_clusters, angular_data=False):
    if angular_data:  # dihedral angles: clustering on the unit circle
        radians = np.deg2rad(values)
        data = np.column_stack((np.cos(radians), np.sin(radians)))
        kmeans = KMeans(n_clusters=num_clusters, init='k-means++', n_init=10, random_state=42)
        kmeans.fit(data)
        return [float(np.rad2deg(np.arctan2(y, x))) for x, y in kmeans.cluster_centers_]

    # 1D data, bond lengths and angles
    data = np.array(values, dtype=float).reshape(-1, 1)
    kmeans = KMeans(n_clusters=num_clusters, init='k-means++', n_init='auto', random_state=42)
    kmeans.fit(data)
    return kmeans.cluster_centers_.flatten().tolist()


def _calculate_atom_type(job):
    """Parameters of one atom type; a top-level function so it can run in a worker process."""
    residue_name, at_type, at_type_values = job
    bond_length_n_clusters, bond_angle_n_clusters, dihedral_angle_n_clusters = \
        atomic_types[residue_name][at_type.split("/")[0]]["n_clusters"]
    result = {}
    if bond_length_n_clusters:
        result["bond_length"] = _run_kmeans(at_type_values["bond_length"], bond_length_n_clusters)
    if bond_angle_n_clusters:
        result["bond_angle"] = _run_kmeans(at_type_values["bond_angle"], bond_angle_n_clusters)
    if dihedral_angle_n_clusters:
        result["dihedral_angle"] = _run_kmeans(at_type_values["dihedral_angle"], dihedral_angle_n_clusters,
                                               angular_data=True)
    return residue_name, at_type, result


class Parameterization:
    def __init__(self, structures_dir, parameters_file, workers=None):
        self.structures_dir = Path(structures_dir)
        self.parameters_file = parameters_file
        self.workers = workers  # None = number of CPUs
        self.data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.parameters = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    def run(self):
        self.extract_data()
        self.check_sample_sizes()
        self.calculate_parameters()
        self.save_parameters()

    def extract_data(self):
        """Measure bond lengths, bond angles and dihedral angles of all proteins in the directory."""
        pdb_files = sorted(self.structures_dir.glob("*.pdb"))
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            for protein_data in tqdm(executor.map(_extract_protein, pdb_files), total=len(pdb_files),
                                     desc="Extracting structural data", unit="protein"):
                for atom_id, m in protein_data:
                    atom_type_data = self.data[m.res_name][m.atom_type]
                    atom_type_data["ID"].append(atom_id)
                    atom_type_data["bond_length"].append(m.bond_length)
                    atom_type_data["bond_angle"].append(m.bond_angle)
                    atom_type_data["dihedral_angle"].append(m.dihedral_angle)

    def check_sample_sizes(self, min_samples=MIN_SAMPLES):
        """Terminate the program if any atom type has fewer than min_samples values."""
        too_small = [f"{residue_name} {at_type}: {len(values['ID'])}"
                     for residue_name, residue_at_types in self.data.items()
                     for at_type, values in residue_at_types.items() if len(values["ID"]) < min_samples]
        if too_small:
            raise SystemExit(f"\nError: fewer than {min_samples} values for these atom types, "
                             "use more proteins:\n    " + "\n    ".join(too_small))

    def calculate_parameters(self):
        """Cluster the measured values of every atom type."""
        jobs = [(residue_name, at_type, at_type_values)
                for residue_name, residue_at_types in self.data.items()
                for at_type, at_type_values in residue_at_types.items()]
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            for residue_name, at_type, result in tqdm(executor.map(_calculate_atom_type, jobs), total=len(jobs),
                                                      desc="Calculating parameters", unit="atom type"):
                self.parameters[residue_name][at_type].update(result)

    def save_parameters(self):
        with open(self.parameters_file, "w", encoding="utf-8") as f:
            json.dump(self.parameters, f, ensure_ascii=False, indent=4)

    def find_outliers(self, gap_factor=0.5):
        print("\n\nOUTLIERS ANALYSIS:")
        outliers = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for residue_name, residue_atom_types in self.data.items():
            for atom_type, atom_type_values in residue_atom_types.items():
                for metric in ("bond_length", "bond_angle"):
                    if metric not in atom_type_values:
                        continue
                    values = np.array(atom_type_values[metric], dtype=float)
                    n = len(values)
                    ordered_indices = np.argsort(values)
                    sorted_values = values[ordered_indices]
                    bulk_width = np.percentile(sorted_values, 99.5) - np.percentile(sorted_values, 0.5)
                    flagged = np.zeros(n, dtype=bool)
                    for m in range(1, n // 100 + 1):
                        if sorted_values[n - m] - sorted_values[n - m - 1] > gap_factor * bulk_width:
                            flagged[n - m:] = True
                        if sorted_values[m] - sorted_values[m - 1] > gap_factor * bulk_width:
                            flagged[:m] = True
                    outliers_mask = np.zeros(n, dtype=bool)
                    outliers_mask[ordered_indices[flagged]] = True
                    if not outliers_mask.any():
                        continue
                    center = self.parameters[residue_name][atom_type][metric][0]
                    deviations = values - center
                    distances_in_md = np.abs(deviations) / np.nanmean(np.abs(deviations[~outliers_mask]))
                    for i in np.flatnonzero(outliers_mask):
                        pdb_file, residue, atom_name = atom_type_values["ID"][i].rsplit(",", 2)
                        res_id = int(residue[len(residue_name):])
                        outliers[pdb_file][(res_id, residue_name)][atom_name].append((metric, values[i], center, distances_in_md[i]))

        units = {"bond_length": "Å", "bond_angle": "°"}
        for pdb_file in sorted(outliers):
            print(f"\n{pdb_file}")
            for (res_id, residue_name), atoms in sorted(outliers[pdb_file].items()):
                print(f"    {residue_name}{res_id}")
                for atom_name, atom_outliers in sorted(atoms.items()):
                    for metric, value, center, distance_in_md in sorted(set(atom_outliers)):
                        unit = units[metric]
                        print(f"        {atom_name}: {metric} = {value:.3f}{unit}, estimated value: {center:.3f}{unit}, {distance_in_md:.1f} MD")
        print("")
