import numpy as np
from fury import molecular as mol

###############################################################################
# Creating a `PeriodicTable()` object to obtain atomic numbers from names of
# elements
table = mol.PTable()

def get_default_molecular_info(fine_name):
    ###############################################################################
    # Creating empty lists which will be filled with atomic information as we
    # parse the pdb file.
    all_info = False
    NumberOfAtoms = 0
    atom_coords = []
    atomic_numbers = []
    atom_types = []
    model = []
    sheets = []
    helix = []
    residue_seq = []
    chain = []
    is_hetatm = []
    current_model_number = 1
    ###############################################################################
    # Parsing the pdb file for information about coordinates and atoms
    pdbfile = open(fine_name, 'r')
    pdb_lines = pdbfile.readlines()
    for line in pdb_lines:
        line = line.split()
        try:
            if line[0] == 'ATOM' or line[0] == 'HETATM':
                if line[-1] != 'H':
                    coorX, coorY, coorZ = float(line[6]), float(line[7]), \
                                        float(line[8])
                    resi = line[5]
                    current_chain = ord(line[4])
                    atom_coords += [[coorX, coorY, coorZ]]
                    residue_seq += [resi]
                    chain += [current_chain]
                    atomic_numbers += [table.atomic_number(line[-1])]
                    atom_types += [line[2]]
                    model += [current_model_number]
                    NumberOfAtoms += 1
                    if(line[0] == 'HETATM'):
                        is_hetatm += [1]
                    else:
                        is_hetatm += [0]
            if line[0] == 'SHEET':
                start_chain = ord(line[5])
                start_resi = int(line[6])
                end_chain = ord(line[8])
                end_resi = int(line[9])
                r = [start_chain, start_resi, end_chain, end_resi]
                sheets += [r]
            if line[0] == 'HELIX':
                start_chain = ord(line[4])
                start_resi = int(line[5])
                end_chain = ord(line[7])
                end_resi = int(line[8])
                r = [start_chain, start_resi, end_chain, end_resi]
                helix += [r]
        except Exception:
            continue

    atom_coords = np.array(atom_coords)
    atomic_numbers = np.array(atomic_numbers)
    atom_types = np.array(atom_types)

    if len(helix) > 0:
        all_info = True
        residue_seq = np.array(residue_seq, dtype=int)
        chain = np.array(chain)
        model = np.array(model)
        sheets = np.array(sheets)
        helix = np.array(helix)
        is_hetatm = np.array(is_hetatm)
        molecule = mol.Molecule(atomic_numbers, atom_coords, atom_types, model,
                                residue_seq, chain, sheets, helix, is_hetatm)

    else:
        all_info = False
        molecule = mol.Molecule(atomic_numbers, atom_coords, atom_types)
    try:
        mol.compute_bonding(molecule)
        mol.get_all_bond_orders(molecule)
    except Exception:
        pass

    return molecule, all_info