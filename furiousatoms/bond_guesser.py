from furiousatoms.molecular import periodic_table
import numpy as np
import kdtree
from furiousatoms.element_lookup import lookup_num_bonds_by_atomic_number

def generate_vdw_radii_map(structure):
    type_to_vdw_map = {}
    for typ in np.unique(structure.atom_types):        
        type_without_numbers = ''.join([i for i in typ if not i.isdigit()]) #For ions like C6 or C4, treat as C
        atomic_number = periodic_table.atomic_number(type_without_numbers)
        vdw = periodic_table.GetVDWRadius(atomic_number)
        #If VDW is unknown, use 0 and do not form bonds with that element.

        type_to_vdw_map[typ] = vdw

    return type_to_vdw_map

def generate_max_bonds_map(structure):
    max_bonds_map = {}
    for typ in np.unique(structure.atom_types):        
        type_without_numbers = ''.join([i for i in typ if not i.isdigit()]) #For ions like C6 or C4, treat as C
        atomic_number = periodic_table.atomic_number(type_without_numbers)
        max_num_bonds = lookup_num_bonds_by_atomic_number(atomic_number)
        if max_num_bonds < 1:
            max_num_bonds = 1

        max_bonds_map[typ] = max_num_bonds

    return max_bonds_map



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



def guess_bonds(structure):
    '''
    Return an array of bonds. The formula for guessing is the 
    same as in MDAnalysis and VMD:
    If (distance between 2 atoms) < (atom 1 VDW radius + atom 2 VDW radius)*FUDGE_FACTOR,
    form a bond. VDW means van der Waals.
    '''
    #For each and every node A in the tree
    #get its knn (4 nearest neighbors) B, C, D, E.
    #Check if the distance between A and B is < cutff (or sum of vdw).

    MIN_CUTOFF = 0.2
    FUDGE_FACTOR = 0.8 #0.6
    STRICT_MAX_BONDS = 4 #by law of chemistry, no atom can form >4 bonds

    type_to_vdw_map = generate_vdw_radii_map(structure)
    max_bonds_map = generate_max_bonds_map(structure)

    atom_count = len(structure.pos)
    point_list = []
    for i in range(atom_count):
        point_list.append(KDTreeItem(structure.pos[i], i))
    
    tree = kdtree.create(point_list=point_list, dimensions=3)
    
    #Sort the atoms by the valence so that the atoms that can form the least number
    #of bonds (like hydrogen) are towards the front.
    all_nodes = []
    for node in tree.inorder():
        all_nodes.append(node)
    all_nodes.sort(key=lambda node: -max_bonds_map[structure.atom_types[node.data.id]])

    connections_map = {} #example: {1: [2,3], 2: [1], 3: [1]}

    bond_list = []
    for src_node in all_nodes: #src = source
        src_id = src_node.data.id
        src_typ = structure.atom_types[src_id]
        src_vdw = type_to_vdw_map[src_typ]
        if src_id not in connections_map:
            connections_map[src_id] = []
        '''
        k is the number of neighbors to check if there is a bond with.
        A small k gives better performance, but it may not find every 
        neighbor--so we add 1 to it. It may seem arbitrary, but it works.
        '''
        k = max_bonds_map[src_typ]
        if k == 0:
            continue
        elif k < STRICT_MAX_BONDS:
            k += 1
        
        for dst_node in tree.search_knn(src_node.data, k): #dst = destination (neighbor atom)
            distance = dst_node[1]
            dst_id = dst_node[0].data.id
            existing_dst_bonds = connections_map.get(dst_id, [])
            if len(existing_dst_bonds) >= STRICT_MAX_BONDS or src_id in existing_dst_bonds:
                continue
            
            dst_typ = structure.atom_types[dst_id]
            dst_vdw = type_to_vdw_map[dst_typ]
            if distance < (src_vdw + dst_vdw) * FUDGE_FACTOR and distance > MIN_CUTOFF:
                bond_list.append([src_id, dst_id])
                #Record the connection to avoid making a duplicate bond
                connections_map.get(src_id).append(dst_id)
                if dst_id not in connections_map:
                    connections_map[dst_id] = []
                connections_map.get(dst_id).append(src_id)

    bond_arr = np.array(bond_list, dtype=int)
    return bond_arr