import json
import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from tqdm import tqdm

from . import PreoptimerError
from .amber import VALENCE, measure_valence_terms

MIN_VALENCE_SAMPLES = 20  # minimal number of values of a valence term key needed for its parameterization


def _measure_statistics(pdb_file):
    """Number, sum and sum of squares of the values of each valence term key of the structure, as
    ({"bonds"/"angles": {key: (n, sum, sum of squares)}}, None), or (None, error message) if the structure cannot
    be processed. Only the sums are kept instead of all values, which keeps the memory independent of the number
    of structures."""
    try:
        values = measure_valence_terms(pdb_file)
    except Exception as error:  # the structure is skipped
        return None, f"{type(error).__name__}: {error}"
    return {kind: {key: (len(key_values), math.fsum(key_values), math.fsum(x * x for x in key_values))
                   for key, key_values in kind_values.items()}
            for kind, kind_values in values.items()}, None


class Parameterization:
    def __init__(self, structures_dir, parameters_file, workers=None):
        self.structures_dir = Path(structures_dir)
        self.parameters_file = parameters_file
        self.workers = workers  # None = number of CPUs
        # kind -> key -> [n, sum, sum of squares]
        self.statistics = {"bonds": defaultdict(lambda: [0, 0.0, 0.0]), "angles": defaultdict(lambda: [0, 0.0, 0.0])}
        self.valence_parameters = {"bonds": {}, "angles": {}}

    def run(self):
        self.extract_data()
        self.calculate_valence_parameters()
        self.save_parameters()

    def extract_data(self):
        """Measure bond lengths and bond angles of all valence terms (bonds and angles of the AMBER topology)
        of all proteins in the directory. Structures which cannot be processed are skipped."""
        pdb_files = sorted(path for path in self.structures_dir.iterdir()
                           if path.is_file() and path.suffix.lower() == ".pdb")
        if not pdb_files:
            raise PreoptimerError(f"no PDB files in {self.structures_dir}")
        skipped = []
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            for pdb_file, (statistics, error) in zip(pdb_files, tqdm(
                    executor.map(_measure_statistics, pdb_files), total=len(pdb_files),
                    desc="Extracting structural data", unit="protein")):
                if error is not None:
                    skipped.append((pdb_file, error))
                    continue
                for kind, kind_statistics in statistics.items():
                    for key, key_statistics in kind_statistics.items():
                        total = self.statistics[kind][key]
                        for i, value in enumerate(key_statistics):
                            total[i] += value
        if skipped:
            print(f"Skipped {len(skipped)} of {len(pdb_files)} structures which could not be processed:")
            for pdb_file, error in skipped:
                print(f"    {pdb_file.name}: {error}")
        if len(skipped) == len(pdb_files):
            raise PreoptimerError("no structure could be processed")

    def calculate_valence_parameters(self, min_samples=MIN_VALENCE_SAMPLES):
        """Mean, standard deviation and number of values of each valence term key with at least min_samples values.
        Keys with fewer values are left out; the less specific keys of the same terms are used instead."""
        for kind, kind_statistics in self.statistics.items():
            for key, (n, total, total_squares) in sorted(kind_statistics.items()):
                if n >= min_samples:
                    mean = total / n
                    std = math.sqrt(max(total_squares / n - mean ** 2, 0))
                    self.valence_parameters[kind][key] = [mean, std, n]

    def save_parameters(self):
        with open(self.parameters_file, "w", encoding="utf-8") as f:
            json.dump({VALENCE: self.valence_parameters}, f, ensure_ascii=False, indent=4)
