#TODO: if vdw is unknown for an element (e.g. C1), assume it is carbon and use tcarbon's vdw

from furiousatoms.molecular import periodic_table #TODO move to SHM
from math import sqrt
import numpy as np
import kdtree

def generate_vdw_radii_map(structure):
    type_to_vdw_map = {}
    for typ in np.unique(structure.atom_types):
        atomic_number = periodic_table.atomic_number(typ)
        vdw = periodic_table.GetVDWRadius(atomic_number)
        type_to_vdw_map[typ] = vdw
    return type_to_vdw_map



class KDTreeItem(object):
    def __init__(self, position, id):
        self.coords = position
        self.id = id
    
    def __len__(self):
        return len(self.coords)

    def __getitem__(self, i):
        return self.coords[i]

    def __repr__(self):
        return 'Item({}, {}, {})'.format(self.coords[0], self.coords[1], self.coords[2], self.id)

# def guess_bonds(structure):
#     structure.type_to_vdw_map = structure.compute_van_der_waals_radii()

#     # box_size = [
#     #     max(structure.pos[:,0]) - min(structure.pos[:,0]), 
#     #     max(structure.pos[:,1]) - min(structure.pos[:,1]), 
#     #     max(structure.pos[:,2]) - min(structure.pos[:,2])
#     # ]
    
#     # largest_vdw = 1.89
#     # max_cutoff = 2.0 * largest_vdw
#     # min_cutoff = max_cutoff if box is not None else None
#     # min_cutoff = 1.2

#     # max_cutoff = max(type_to_vdw_map.values())
#     # min_cutoff = min(type_to_vdw_map.values())
    
#     #TODO try making point_list a numpy arr, then use tolist() during kdtree creation
#     # point_list = structure.pos.tolist()
#     atom_count = len(structure.pos)
#     point_list = []  #[[0, 0, 0, 0] * len(structure.pos)]
#     for i in range(atom_count):
#         point_list.append(Item(structure.pos[i], i))

#     tree = kdtree.create(point_list=point_list, dimensions=3)

#     bond_list = []

#     # get_pairs(tree, min_cutoff, max_cutoff)
#     get_pairs(tree, structure, bond_list)
#     print(bond_list)

#     kdtree.visualize(tree)

#     bond_arr = np.array(bond_list, dtype=int)
#     return bond_arr






def guess_bonds(structure):
    #For each and every node A in the tree
    #get its knn (4 nearest neighbors) B, C, D, E.
    #Check if the distance between A and B is < cutff (or sum of vdw).

    MIN_CUTOFF = 1.0
    K = 4
    FUDGE_FACTOR = 1.0 #0.55

    type_to_vdw_map = generate_vdw_radii_map(structure)

    atom_count = len(structure.pos)
    point_list = []  #[[0, 0, 0, 0] * len(structure.pos)]
    for i in range(atom_count):
        point_list.append(KDTreeItem(structure.pos[i], i))
    
    tree = kdtree.create(point_list=point_list, dimensions=3)
    
    bond_list = []
    for src_node in tree.inorder():
        src_typ = structure.atom_types[src_node.data.id]
        src_vdw = type_to_vdw_map[src_typ]
        for neighbor in tree.search_knn(src_node.data, K):
            dist = neighbor[1] #distance(node.data, neighbor.data)
            neighbor_id = neighbor[0].data.id
            # print("Dist between",node.data,"and",neighbor[0].data,"is",dist)

            if tree.left:
                dst_typ = structure.atom_types[neighbor[0].data.id]
                dst_vdw = type_to_vdw_map[dst_typ]
                # print("compare ",dist,((src_vdw + dst_vdw) * 0.55))
                if dist < (src_vdw + dst_vdw) * FUDGE_FACTOR and dist > MIN_CUTOFF:
                    bond_list.append([src_node.data.id, neighbor_id])

        # print("---")
    
    print(bond_list)

    # kdtree.visualize(tree)

    bond_arr = np.array(bond_list, dtype=int)
    return bond_arr


if __name__ == "__main__":
    from furiousatoms.parsers.pdb_parser import PDBParser
    parser = PDBParser()
    # structure = parser.parse("C:\\Users\\Pete\\Desktop\\furious-atoms-exercises\\p38-MAPK14.pdb")
    structure = parser.parse("C:\\Users\\Pete\\Desktop\\furious-atoms-exercises\\CB_18\\CB_18.pdb")
    
    guessed_bonds = guess_bonds(structure)