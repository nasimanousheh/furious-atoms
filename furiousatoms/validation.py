from furiousatoms.bond_guesser import generate_vdw_radii_map, connections_map_to_array, generate_connections_map
from math import sqrt
import numpy as np

def atom_distance(atom1, atom2):
    distance = 0
    for i in range(0, 3):
        distance += (atom1[i] - atom2[i]) ** 2
    
    return sqrt(distance)


def finalize_bonds(structure):
    should_recreate = False #we avoid mixing up the user's bond order if possible
    final_num_bonds = len(structure.bonds)
    NUM_ATOMS = len(structure.pos)
    
    '''
    In connections_map, each index holds a list of atoms bonded to.
    Example: bond [0, 2] creates map [[2], [], [0]]
    '''
    connections_map = generate_connections_map(structure)
    
    to_remove = {}

    #Delete any bond that is too long
    FUDGE_FACTOR = 1.05
    type_to_vdw_map = generate_vdw_radii_map(structure)
    for src_id in range(NUM_ATOMS):
        src_typ = structure.atom_types[src_id]
        src_vdw = type_to_vdw_map[src_typ]
        for dst_id in connections_map[src_id]:
            dst_typ = structure.atom_types[dst_id]
            dst_vdw = type_to_vdw_map[dst_typ]

            distance = atom_distance(structure.pos[src_id], structure.pos[dst_id])
            if distance > (src_vdw + dst_vdw) * FUDGE_FACTOR:
                to_remove_sub_list = to_remove.get(src_id, [])
                to_remove_sub_list.append(dst_id)
                to_remove[src_id] = to_remove_sub_list
                final_num_bonds -= 1
                should_recreate = True
    
    if not should_recreate:
        return structure.bonds

    return connections_map_to_array(connections_map, final_num_bonds, to_remove)



if __name__ == "__main__":
    from furiousatoms.molecular import MolecularStructure
    s = MolecularStructure.create_empty
    s.bonds = np.array([
        [1,2],
        [3,2],
        [2,1],
        [1,4],
        [4,1],
    ])
    s.bonds -= 1
    s.pos = np.array([
        [0,0,0],
        [0,0,0],
        [0,0,0],
        [0,20,0],
    ])
    s.atom_types = np.array([
        'C',
        'C',
        'C',
        'C'
    ])
    print(finalize_bonds(s))