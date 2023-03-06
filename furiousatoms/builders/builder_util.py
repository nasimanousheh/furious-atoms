import numpy as np

def copy_bonds(bonds, num_unitcell_in_lx, num_unitcell_in_ly, unit_atom_count):
    copied_bonds = np.zeros(shape=(len(bonds) * num_unitcell_in_lx * num_unitcell_in_ly, 2),
            dtype=int)
    bondIndex = 0
    for u in range(num_unitcell_in_lx * num_unitcell_in_ly):
        for bond in bonds:
            offset = u * unit_atom_count
            for i in range(0, 2):
                copied_bonds[bondIndex][i] = bond[i] + offset
            bondIndex += 1
    return copied_bonds