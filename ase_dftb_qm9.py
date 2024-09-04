import sys
import os
from os.path import exists
from os.path import basename, splitext
import csv
from glob import glob
from time import time

from typing import Any

from ase.io.extxyz import read_xyz
from ase.calculators.dftb import Dftb

# Retrieve number of molecules to run on from first argument
nmols = int(sys.argv[1])

xyz_paths = glob('/home/alexlm/databases/qm9_xyz/*.xyz')[:nmols]

ase_atoms = []
mol_names = []
for path in xyz_paths:
    mol_id, ext = splitext(basename(path))
    try:
        mol = next(read_xyz(open(path), 0))
    except:
        print(f'Failed on {path}')
        continue
    if all(atom.number in [1, 6, 7, 8] for atom in mol):
        ase_atoms.append(mol)
        mol_names.append(mol_id)

# Directory of Slater-Koster files
# The slash at the end is necessary, which I think is a bug in ASE
sk_dir = 'auorg-1-1/'
calc = Dftb(label='current',
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='p',
            Hamiltonian_MaxAngularMomentum_C='p',
            Hamiltonian_MaxAngularMomentum_N='p',
            Hamiltonian_MaxAngularMomentum_H='s',
            # Set directory of Slater-Koster files
            slako_dir=sk_dir
            )

outfile_name = 'ase_dftb_qm9.csv'
with open(outfile_name, 'w') as outfile:
    writer = csv.writer(outfile)
    header = ['mol_id', 'energy', 'run_time']
    writer.writerow(header)

    for atoms, mol_name in zip(ase_atoms, mol_names):
        atoms.calc = calc
        start_time = time()
        energy = atoms.get_potential_energy()
        end_time = time()
        run_time = end_time - start_time
        
        # Write results
        row = [mol_name, energy, run_time]
        writer.writerow(row)
