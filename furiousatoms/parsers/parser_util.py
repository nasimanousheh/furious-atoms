import numpy as np

def float_or_zero(val):
    try:
        return float(val)
    except ValueError:
        return 0.0
    
def has_bond(bonds, new_bond):
    for i in range(len(bonds)):
        if bonds[i][0] == new_bond[0] and bonds[i][1] == new_bond[1] \
        or bonds[i][1] == new_bond[0] and bonds[i][0] == new_bond[1]:
            return True
    return False

def clean_bonds(bonds, atom_count):
    if not bonds.shape[0] > 0 or not bonds.shape[1] == 2:
        return np.empty(shape=(0, 2), dtype=int)
    return bonds[((bonds[:,0] < atom_count) & (bonds[:,1] < atom_count))]