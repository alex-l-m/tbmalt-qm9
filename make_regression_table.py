from collections import Counter
from ase.io.extxyz import read_xyz
import pandas as pd

def read_properties(xyz_path):
    with open(xyz_path) as xyz_file:
        n_atoms = int(next(xyz_file))
        property_line = next(xyz_file)
    # First two will be "gdb \d+"
    properties = property_line.split()
    return properties
tbmalt_results = pd.read_csv('trained.csv')

rows = []
for row in tbmalt_results.itertuples():
    mol_id = row.mol_id
    energy = row.energy
    xyz_path = f'qm9_xyz/{mol_id}.xyz'
    ase_atoms = next(read_xyz(open(xyz_path), 0))
    formula = Counter(ase_atoms.get_chemical_symbols())
    properties = read_properties(xyz_path)
    # I want property 13, U0. They're one indexed in the documentation
    u0 = float(properties[13 - 1])
    # I don't currently use this for anything, retrieving just in case I want
    # to compare to other analyses that use it
    qm9_id = int(properties[2 - 1])
    row = {'mol_id': mol_id,
           'qm9_id': qm9_id,
           'n_H': formula['H'],
           'n_C': formula['C'],
           'n_N': formula['N'],
           'n_O': formula['O'],
           # n_F should be zero because of previous filtering, so not bothering
           # to retrieve that
           'qm9_u0': u0}
    rows.append(row)
pd.DataFrame(rows).to_csv('regression_table.csv', index = False)
