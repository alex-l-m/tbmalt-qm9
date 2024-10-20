import sys
from os.path import exists
from os.path import basename, splitext
import csv
from glob import glob
from time import time

import pandas as pd
import torch
from tbmalt import Geometry, OrbitalInfo
from tbmalt.physics.dftb import Dftb2
from tbmalt.physics.dftb.feeds import SkFeed, SkfOccupationFeed, HubbardFeed, RepulsiveSplineFeed

from ase.io.extxyz import read_xyz

# Retrieve number of molecules to run on from first argument
nmols = int(sys.argv[1])
# First half is training set
training_set_size = nmols // 2

print(f'Running {nmols} molecules with TBMaLT')

# Path for the predicted values
outpath = sys.argv[2]

# Path for the loss values
losspath = sys.argv[3]

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

# Load the target property for training as a dictionary mapping the 'mol_id' column to the 'tbmalt_target' column
target_property_dict = dict()
for row in pd.read_csv('tbmalt_target.csv').itertuples():
    target_property_dict[row.mol_id] = row.tbmalt_target


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
# List of target property values
target_property_list = []
for path in xyz_paths:
    mol_id, ext = splitext(basename(path))
    try:
        mol = next(read_xyz(open(path), 0))
    except:
        print(f'Failed on {path}')
        continue
    # Keep it only if every atomic number is in the keys of shell_dict
    if all(atom.number in shell_dict.keys() for atom in mol):
        try:
            target_property_list.append(target_property_dict[mol_id])
            ase_atoms.append(mol)
            mol_names.append(mol_id)
        except KeyError:
            print(f'No target property for {mol_id}')

geometry = Geometry.from_ase_atoms(ase_atoms)
orbs = OrbitalInfo(geometry.atomic_numbers, shell_dict, shell_resolved=False)

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
# filling_temp has to be given value to enable Fermi smearing.  without that,
# not everything converges, due to oscillating charges. In addition, some
# results converge but to the wrong value
# filling_temp=0 should do it, but then I get errors during training saying that it's trying a Cholesky decomposition on a matrix that's not positive definite
dftb_calculator = Dftb2(h_feed, s_feed, o_feed, u_feed, r_feed, filling_temp = 0.01)

# Make a torch tensor of target properties
target_properties = torch.tensor(target_property_list)

# Construct machine learning object
lr = 0.0001
# Use torch's mean squared error loss
criterion = torch.nn.MSELoss()

h_var, s_var = [], []
for key in h_feed.off_sites.keys():

    # Collect spline parameters and add to optimizer
    h_feed.off_sites[key].abcd.requires_grad_(True)
    s_feed.off_sites[key].abcd.requires_grad_(True)
    h_var.append({'params': h_feed.off_sites[key].abcd, 'lr': lr})
    s_var.append({'params': s_feed.off_sites[key].abcd, 'lr': lr})

optimizer = getattr(torch.optim, 'Adam')(h_var + s_var, lr=lr)
print('Beginning training')
# Write the headers
with open(losspath, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['epoch', 'split', 'loss'])
with open(outpath, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['epoch', 'mol_id', 'energy', 'repulsive_energy', 'scc_energy', 'run_time', 'split'])
for epoch in range(20):
    start = time()
    # Run the DFTB calculation
    dftb_calculator(geometry, orbs)
    results = getattr(dftb_calculator, 'total_energy')
    repulsive_energy = dftb_calculator.repulsive_energy
    scc_energy = dftb_calculator.scc_energy

    optimizer.zero_grad()
    # Calculate loss only for the training molecules (first half)
    loss = criterion(results[:training_set_size], target_properties[:training_set_size])
    loss.backward()
    optimizer.step()
    # Also calculate test set loss, just so I can write it to the file
    test_loss = criterion(results[training_set_size:], target_properties[training_set_size:])
    end = time()
    run_time = end - start
    print(f'Epoch {epoch}: time {run_time}, loss {loss.item()}')
    with open(losspath, 'a') as f:
        writer = csv.writer(f)
        writer.writerow([epoch, 'train', loss.item()])
        writer.writerow([epoch, 'test', test_loss.item()])
    with open(outpath, 'a') as f:
        writer = csv.writer(f)
    
        for i, (name, energy, repulsive_energy, scc_energy) \
                in enumerate(zip(mol_names, results, repulsive_energy, scc_energy)):
            split = 'train' if i < training_set_size else 'test'
            writer.writerow([epoch, name, energy.item(), repulsive_energy.item(), scc_energy.item(), run_time, split])
