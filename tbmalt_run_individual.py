import sys
from os.path import exists
from os.path import basename, splitext
import csv
from glob import glob
from time import time

import torch
from tbmalt import Geometry, OrbitalInfo
from tbmalt.physics.dftb import Dftb2
from tbmalt.physics.dftb.feeds import SkFeed, SkfOccupationFeed, HubbardFeed, RepulsiveSplineFeed
from tbmalt.common.exceptions import ConvergenceError

from ase.io.extxyz import read_xyz

# Retrieve number of molecules to run on from first argument
nmols = int(sys.argv[1])

Tensor = torch.Tensor

# This must be set until typecasting from HDF5 databases has been implemented.
torch.set_default_dtype(torch.float64)

# ============== #
# STEP 1: Inputs #
# ============== #

# 1.1: System settings
# --------------------

# Provide information about the orbitals on each atom; this is keyed by atomic
# numbers and valued by azimuthal quantum numbers like so:
#   {Z₁: [ℓᵢ, ℓⱼ, ..., ℓₙ], Z₂: [ℓᵢ, ℓⱼ, ..., ℓₙ], ...}
shell_dict = {1: [0], 6: [0, 1], 7: [0, 1], 8: [0, 1]}

# 1.2: Model settings
# -------------------
# Location at which the DFTB parameter set database is located
parameter_db_path = 'example_dftb_parameters.h5'

# Location of a file storing the properties that will be fit to.
target_path = 'target_data.json'

# ============= #
# STEP 2: Setup #
# ============= #

# 2.1: Target system specific objects
# -----------------------------------

# Construct the `Geometry` and `OrbitalInfo` objects. The former is analogous
# to the ase.Atoms object while the latter provides information about what
# orbitals are present and which atoms they belong two.
# Have to do index 0 because there's extra junk at the end of the file
xyz_paths = glob('qm9_xyz/*.xyz')[:nmols]

ase_atoms = []
mol_names = []
for path in xyz_paths:
    mol_id, ext = splitext(basename(path))
    try:
        mol = next(read_xyz(open(path), 0))
    except:
        print(f'Failed on {path}')
        continue
    # Keep it only if every atomic number is in the keys of shell_dict
    if all(atom.number in shell_dict.keys() for atom in mol):
        ase_atoms.append(mol)
        mol_names.append(mol_id)

geometries = [Geometry.from_ase_atoms([mol]) for mol in ase_atoms]
orbs = [OrbitalInfo(geometry.atomic_numbers, shell_dict, shell_resolved=False) \
        for geometry in geometries]

# 2.2: Loading of the DFTB parameters into their associated feed objects
# ----------------------------------------------------------------------

# Construct the Hamiltonian and overlap matrix feeds; but ensure that the DFTB
# parameter set database actually exists first.
if not exists(parameter_db_path):
    raise FileNotFoundError(
        f'The DFTB parameter set database "{parameter_db_path}" could '
        f'not be found, please ensure "example_01_setup.py" has been run.')

# Identify which species are present
species = shell_dict.keys()

# Load the Hamiltonian feed model
h_feed = SkFeed.from_database(parameter_db_path, species, 'hamiltonian',
                              interpolation='spline')

# Load the overlap feed model
s_feed = SkFeed.from_database(parameter_db_path, species, 'overlap',
                              interpolation='spline')

# Load the occupation feed object
o_feed = SkfOccupationFeed.from_database(parameter_db_path, species)

# Load the Hubbard-U feed object
u_feed = HubbardFeed.from_database(parameter_db_path, species)

# Load the repulsive spline feed object
r_feed = RepulsiveSplineFeed.from_database(parameter_db_path, species)

# 2.3: Construct the SCC-DFTB calculator object
# ---------------------------------------------
# calculator object.
dftb_calculator = Dftb2(h_feed, s_feed, o_feed, u_feed, r_feed, suppress_SCF_error = False)

outpath = 'tbmalt_results.csv'
with open(outpath, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['molecule', 'status', 'energy', 'repulsive_energy', 'scc_energy', 'run_time'])
    for name, geometry, orbs in zip(mol_names, geometries, orbs):
        
        try:
            dftb_calculator(geometry, orbs)
            # Run the DFTB calculation
            start_time = time()
            results = dftb_calculator(geometry, orbs)
            repulsive_energy = dftb_calculator.repulsive_energy
            scc_energy = dftb_calculator.scc_energy
            end_time = time()
            run_time = end_time - start_time
            for energy, repulsive_energy, scc_energy \
                    in zip(results, repulsive_energy, scc_energy):
                writer.writerow([name, 'success',
                    energy.item(), repulsive_energy.item(), scc_energy.item(),
                    run_time])
        except ConvergenceError as e:
            writer.writerow([name, 'ConvergenceError',
                'NA', 'NA', 'NA', 'NA'])
