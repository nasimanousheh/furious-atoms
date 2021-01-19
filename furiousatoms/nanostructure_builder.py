
"""
Single-wall nanotube, Multiple-wall nanotube and graphene builder


"""
from MDAnalysis import *
import MDAnalysis.analysis.align
import numpy as np
from numpy.linalg import norm
from fractions import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
import sys
import io

thre = 1e-10
vacuum = 4
SM = SharedMemory()

"""
  The numbers (n,m) show that your tube is obtained from taking one atom of the sheet and rolling it onto
  that atom that is at located na1+ma2 away from your original atom. N is the Number of hexagons in a unit cell. a1 and a2 are lattice vectors.
  (n,m=n) gives an “armchair” tube,e.g. (5,5). (n,m=0) gives an “zig-zag” tube, e.g. (6,0). Other tubes are “chiral”, e.g. (6,2)
"""


def SWNT_builder(H_termination_SWNT, n, m, N=1, length=None, a=1.421, species=('C', 'C'), centered=False):

    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2))
    Ch = n*a1+m*a2
    T = t1*a1+t2*a2
    if length:
        N = int(np.ceil(length/np.linalg.norm(T)))
    Ch_proj, T_proj = [v/np.linalg.norm(v)**2 for v in [Ch, T]]
    basis = [np.array((0, 0)), (a1+a2)/3]
    pts = []
    for i1, i2 in product(range(0, n+t1+1), range(t2, m+1)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            if all(-thre < pt.dot(v) < 1-thre for v in [Ch_proj, T_proj]):
                for k in range(N):
                    pts.append((sp, pt+k*T))

    # Here we define the diameter of SWNT:
    SM.diameter_SWNT = np.linalg.norm(Ch)/np.pi
    # This function converts the graphene to nanotube with given diameter:
    def gr2tube(v):
        phi = 2*np.pi*v.dot(Ch_proj)
        return np.array((SM.diameter_SWNT/2*np.cos(phi),
                         SM.diameter_SWNT/2*np.sin(phi),
                         v.dot(T_proj)*np.linalg.norm(T)))
    xyz = [gr2tube(v) for _, v in pts]
    atom_types_swnt = [v for v, _ in pts]
    m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    fragments = m.to_json(scale=1)

    # Number of atoms in SWNT:
    num_atoms_swnt = len(xyz)

    # Atom coordinates:
    coord_array_swnt = np.array(xyz)
    assert coord_array_swnt.shape == (num_atoms_swnt, 3)
    swnt = MDAnalysis.Universe.empty(num_atoms_swnt, trajectory=True, n_residues=10)
    swnt.atoms.positions = coord_array_swnt

    # Bonds information connected the atoms:
    all_bonds_swnt = np.array(fragments['bonds'])
    swnt.add_TopologyAttr('name', atom_types_swnt)
    swnt.add_bonds(all_bonds_swnt)

    # If the user chooses "None", only SWNT structure without hydrogens will be returned:
    if H_termination_SWNT == 'None':
        return swnt

   ##############################################Create Hydrogens at only one end of tube##############################################

   #To hydrogenate one end of tube, we devide the number of atoms in nanotube into two:

    half_num_atoms_swnt = int(num_atoms_swnt/2)
    num_hydrogen = 0
    H_coordinaes = []

    for a in range(half_num_atoms_swnt):

        # 'a' is the atom id. We look at the bond indices to verify where the atom 'a' index is appeared:
        indices_of_a = np.where(all_bonds_swnt == a)[0]

        # If 'a' in bond indices is repeated only once, it means it has connected to one atom. So to hydrogenate the atom, we need to connect it to two hydrogen atoms.
        if len(indices_of_a) == 1:
            # Here we identify the indices of all atoms at the end of tube that are need to be hydrogenized.
            end_atom_indices = (all_bonds_swnt[indices_of_a])
            if end_atom_indices[0][0] == a:
                end_atom_indices = (all_bonds_swnt[indices_of_a])
                f_connec_to_end_atom_index = end_atom_indices[0][1]
                core_connections = all_bonds_swnt[np.where(all_bonds_swnt == f_connec_to_end_atom_index)[0]]
                Both_connected_atoms_with_a = np.setdiff1d(core_connections, [f_connec_to_end_atom_index])
                Both_connected_atoms = np.setdiff1d(Both_connected_atoms_with_a, a)
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
            if end_atom_indices[0][1] == a:
                end_atom_indices = (all_bonds_swnt[indices_of_a])
                f_connec_to_end_atom_index = end_atom_indices[0][0]
                core_connections = all_bonds_swnt[np.where(all_bonds_swnt == f_connec_to_end_atom_index)[0]]
                Both_connected_atoms_with_a = np.setdiff1d(core_connections, [f_connec_to_end_atom_index])
                Both_connected_atoms = np.setdiff1d(Both_connected_atoms_with_a, a)
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
            H_coord_1 = xyz[a] - xyz[right_connection] + xyz[f_connec_to_end_atom_index]
            H_coord_2 = xyz[a] - xyz[left_connection] + xyz[f_connec_to_end_atom_index]
            H_coordinaes.extend([H_coord_1])
            H_coordinaes.extend([H_coord_2])
            num_hydrogen = num_hydrogen + 2

        # If a in bond indices is repeated two times, it means it has connected to two other atoms. So to hydrogenate the this atom, we need to connect it to one hydrogen atom.
        if len(indices_of_a) == 2:
            end_atom_indices = (all_bonds_swnt[indices_of_a])
            if end_atom_indices[0][0] == a:
                end_atom_index = end_atom_indices[0][0]
                f_connec_to_end_atom_index = np.where(all_bonds_swnt == end_atom_index)[0]
                core_connections = all_bonds_swnt[f_connec_to_end_atom_index]
                Both_connected_atoms = np.setdiff1d(core_connections, [end_atom_index])
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
                first_vector = xyz[end_atom_index] - xyz[right_connection]
                second_vector = xyz[end_atom_index] - xyz[left_connection]
                H_coord = (first_vector + second_vector) + xyz[end_atom_index]
                H_coordinaes.extend([H_coord])
                num_hydrogen = num_hydrogen + 1

            if end_atom_indices[0][1] == a:
                end_atom_index = end_atom_indices[0][1]
                f_connec_to_end_atom_index = np.where(all_bonds_swnt == end_atom_index)[0]
                core_connections = all_bonds_swnt[f_connec_to_end_atom_index]
                Both_connected_atoms = np.setdiff1d(core_connections, [end_atom_index])
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
                first_vector = xyz[end_atom_index] - xyz[right_connection]
                second_vector = xyz[end_atom_index] - xyz[left_connection]
                H_coord = (first_vector + second_vector) + xyz[end_atom_index]
                H_coordinaes.extend([H_coord])
                num_hydrogen = num_hydrogen + 1

# Here we define the bond information between the atoms of SWNT and hydrogen, if the number of hydrogen is not zero:
    if num_hydrogen > 0:
        h = MDAnalysis.Universe.empty(num_hydrogen, trajectory=True, n_residues=10)
        coord_array_H_indice = np.array(H_coordinaes)
        assert coord_array_H_indice.shape == (num_hydrogen, 3)
        h.atoms.positions = coord_array_H_indice
        h.add_TopologyAttr('name', ['H']*num_hydrogen)
    combined_one_end = MDAnalysis.Merge(swnt.atoms, h.atoms)
    combined_one_end.add_bonds(all_bonds_swnt)
    num_hydrogen = 0
    for x in range(half_num_atoms_swnt):
        indices_of_a = np.where(all_bonds_swnt == x)
        if len(indices_of_a[0]) == 1:
            combined_one_end.add_bonds([(x, num_atoms_swnt + num_hydrogen)])
            combined_one_end.add_bonds([(x, num_atoms_swnt + num_hydrogen + 1)])
            num_hydrogen = num_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            combined_one_end.add_bonds([(x, num_atoms_swnt + num_hydrogen)])
            num_hydrogen = num_hydrogen + 1

    # If the user chooses "One end", only SWNT structure with one end hydrogenated will be returned:
    if H_termination_SWNT == 'One end':
        return combined_one_end

#########################################################
    num_hydrogen = 0
    H_coordinaes = []
    for a in range(num_atoms_swnt):
        indices_of_a = np.where(all_bonds_swnt == a)[0]
        if len(indices_of_a) == 1:
            end_atom_indices = (all_bonds_swnt[indices_of_a])
            if end_atom_indices[0][0] == a:
                end_atom_indices = (all_bonds_swnt[indices_of_a]) #[2735 2447] main= 2735
                f_connec_to_end_atom_index = end_atom_indices[0][1] # 2447
                core_connections = all_bonds_swnt[np.where(all_bonds_swnt == f_connec_to_end_atom_index)[0]]
                Both_connected_atoms_with_a = np.setdiff1d(core_connections, [f_connec_to_end_atom_index])
                Both_connected_atoms = np.setdiff1d(Both_connected_atoms_with_a, a)
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
            if end_atom_indices[0][1] == a:
                end_atom_indices = (all_bonds_swnt[indices_of_a])
                f_connec_to_end_atom_index = end_atom_indices[0][0]
                core_connections = all_bonds_swnt[np.where(all_bonds_swnt == f_connec_to_end_atom_index)[0]]
                Both_connected_atoms_with_a = np.setdiff1d(core_connections, [f_connec_to_end_atom_index])
                Both_connected_atoms = np.setdiff1d(Both_connected_atoms_with_a, a)
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
            H_coord_1 = xyz[a] - xyz[right_connection] + xyz[f_connec_to_end_atom_index]
            H_coord_2 = xyz[a] - xyz[left_connection] + xyz[f_connec_to_end_atom_index]
            H_coordinaes.extend([H_coord_1])
            H_coordinaes.extend([H_coord_2])
            num_hydrogen = num_hydrogen + 2

        if len(indices_of_a) == 2:
            end_atom_indices = (all_bonds_swnt[indices_of_a])
            if end_atom_indices[0][0] == a:
                end_atom_index = end_atom_indices[0][0]
                f_connec_to_end_atom_index = np.where(all_bonds_swnt == end_atom_index)[0]
                core_connections = all_bonds_swnt[f_connec_to_end_atom_index]
                Both_connected_atoms = np.setdiff1d(core_connections, [end_atom_index])
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
                first_vector = xyz[end_atom_index] - xyz[right_connection]
                second_vector = xyz[end_atom_index] - xyz[left_connection]
                H_coord = (first_vector + second_vector) + xyz[end_atom_index]
                H_coordinaes.extend([H_coord])
                num_hydrogen = num_hydrogen + 1

            if end_atom_indices[0][1] == a:
                end_atom_index = end_atom_indices[0][1]
                f_connec_to_end_atom_index = np.where(all_bonds_swnt == end_atom_index)[0]
                core_connections = all_bonds_swnt[f_connec_to_end_atom_index]
                Both_connected_atoms = np.setdiff1d(core_connections, [end_atom_index])
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
                first_vector = xyz[end_atom_index] - xyz[right_connection]
                second_vector = xyz[end_atom_index] - xyz[left_connection]
                H_coord = (first_vector + second_vector) + xyz[end_atom_index]
                H_coordinaes.extend([H_coord])
                num_hydrogen = num_hydrogen + 1

    if num_hydrogen > 0:
        h = MDAnalysis.Universe.empty(num_hydrogen, trajectory=True, n_residues=10)
        coord_array_H_indice = np.array(H_coordinaes)
        assert coord_array_H_indice.shape == (num_hydrogen, 3)
        h.atoms.positions = coord_array_H_indice
        h.add_TopologyAttr('name', ['H']*num_hydrogen)
    combined = MDAnalysis.Merge(swnt.atoms, h.atoms)
    combined.add_bonds(all_bonds_swnt)
    num_hydrogen = 0
    for x in range(num_atoms_swnt):
        indices_of_a = np.where(all_bonds_swnt == x)
        if len(indices_of_a[0]) == 1:
            combined.add_bonds([(x, num_atoms_swnt + num_hydrogen)])
            combined.add_bonds([(x, num_atoms_swnt + num_hydrogen + 1)])
            num_hydrogen = num_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            combined.add_bonds([(x, num_atoms_swnt + num_hydrogen)])
            num_hydrogen = num_hydrogen + 1
    # If the user chooses "Both ends", SWNT structure with both ends hydrogenated will be returned:
    if H_termination_SWNT == 'Both ends':
        return combined

def MWNT_builder(n, m, N=1, length=True, a=1.421, species=('B', 'C'), centered=False):
    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2))
    Ch = n*a1+m*a2
    T = t1*a1+t2*a2
    # if length:
    #     N = int(np.ceil(length/norm(T)))
    Ch_proj, T_proj = [v/norm(v)**2 for v in [Ch, T]]
    basis = [np.array((0, 0)), (a1+a2)/3]
    pts = []
    xyz = []
    atom_types_swnt = []
    for i1, i2 in product(range(0, n+t1+1), range(t2, m+1)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            if all(-thre < pt.dot(v) < 1-thre for v in [Ch_proj, T_proj]):
                for k in range(N):
                    pts.append((sp, pt+k*T))
    diameter = norm(Ch)/np.pi
    print(diameter)
    def gr2tube(v):
        phi = 2*np.pi*v.dot(Ch_proj)
        return np.array((diameter/2*np.cos(phi),
                         diameter/2*np.sin(phi),
                         v.dot(T_proj)*norm(T)))
    xyz = [gr2tube(v) for _, v in pts]
    atom_types_swnt = [v for v, _ in pts]
    mol_1 = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    fragments = mol_1.to_json(scale =1)
    n_atoms_swnt = len(xyz)
    coord_array_swnt = np.array(xyz)
    assert coord_array_swnt.shape == (n_atoms_swnt, 3)
    swnt = MDAnalysis.Universe.empty(n_atoms_swnt, trajectory=True, n_residues=2)
    swnt.atoms.positions = coord_array_swnt
    all_bonds_swnt = np.array(fragments['bonds'])
    swnt.add_TopologyAttr('name', atom_types_swnt)
    swnt.add_bonds(all_bonds_swnt)
    return swnt

def graphene_builder(H_termination_graphene, n, m, N=1, length=None, a=1.421, species=('C', 'C'), centered=False):

    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0,0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2,0))
    Ch = n*a1+m*a2
    T = t1*a1+t2*a2
    if length:
        N = int(np.ceil(length/np.linalg.norm(T)))
    Ch_proj, T_proj = [v/np.linalg.norm(v)**2 for v in [Ch, T]]
    basis = [np.array((0,0,0)), (a1+a2)/3]
    pts = []
    for i1, i2 in product(range(0, n+t1+1), range(t2, m+1)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            pts.append((sp, pt))
    xyz = [v for _, v in pts]
    atom_types_graphene = [v for v, _ in pts]
    m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    fragments = m.to_json(scale =1)

    # Number of atoms in graphene:
    num_atoms_graphene = len(xyz)
    # Atom coordinates:
    coord_array_graphene = np.array(xyz)
    assert coord_array_graphene.shape == (num_atoms_graphene, 3)
    graphene = MDAnalysis.Universe.empty(num_atoms_graphene, trajectory=True, n_residues=10)
    graphene.atoms.positions = coord_array_graphene
    # Bonds information connected the atoms:
    all_bonds_graphene = np.array(fragments['bonds'])
    graphene.add_TopologyAttr('name', atom_types_graphene)
    graphene.add_bonds(all_bonds_graphene)
    # If the user chooses "None", only graphene structure without hydrogens will be returned:
    if H_termination_graphene == 'None':
        return graphene
   ##############################################Create Hydrogens at the end of graphene##############################################

    num_hydrogen = 0
    H_coordinaes = []
    for a in range(num_atoms_graphene):
        # 'a' is the atom id. We look at the bond indices to verify where the atom 'a' index is appeared:
        indices_of_a = np.where(all_bonds_graphene == a)[0]
        # If 'a' in bond indices is repeated only once, it means it has connected to one atom. So to hydrogenate the atom, we need to connect it to two hydrogen atoms.
        if len(indices_of_a) == 1:
            end_atom_indices = (all_bonds_graphene[indices_of_a])
            if end_atom_indices[0][0] == a:
                end_atom_indices = (all_bonds_graphene[indices_of_a])
                f_connec_to_end_atom_index = end_atom_indices[0][1]
                core_connections = all_bonds_graphene[np.where(all_bonds_graphene == f_connec_to_end_atom_index)[0]]
                Both_connected_atoms_with_a = np.setdiff1d(core_connections, [f_connec_to_end_atom_index])
                Both_connected_atoms = np.setdiff1d(Both_connected_atoms_with_a, a)
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
            if end_atom_indices[0][1] == a:
                end_atom_indices = (all_bonds_graphene[indices_of_a])
                f_connec_to_end_atom_index = end_atom_indices[0][0]
                core_connections = all_bonds_graphene[np.where(all_bonds_graphene == f_connec_to_end_atom_index)[0]]
                Both_connected_atoms_with_a = np.setdiff1d(core_connections, [f_connec_to_end_atom_index])
                Both_connected_atoms = np.setdiff1d(Both_connected_atoms_with_a, a)
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
            H_coord_1 = xyz[a] - xyz[right_connection] + xyz[f_connec_to_end_atom_index]
            H_coord_2 = xyz[a] - xyz[left_connection] + xyz[f_connec_to_end_atom_index]
            H_coordinaes.extend([H_coord_1])
            H_coordinaes.extend([H_coord_2])
            num_hydrogen = num_hydrogen + 2

        if len(indices_of_a) == 2:
            end_atom_indices = (all_bonds_graphene[indices_of_a])
            if end_atom_indices[0][0] == a:
                end_atom_index = end_atom_indices[0][0]
                f_connec_to_end_atom_index = np.where(all_bonds_graphene == end_atom_index)[0]
                core_connections = all_bonds_graphene[f_connec_to_end_atom_index]
                Both_connected_atoms = np.setdiff1d(core_connections, [end_atom_index])
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]
                first_vector = xyz[end_atom_index] - xyz[right_connection]
                second_vector = xyz[end_atom_index] - xyz[left_connection]
                H_coord = (first_vector + second_vector) + xyz[end_atom_index]
                H_coordinaes.extend([H_coord])
                num_hydrogen = num_hydrogen + 1

            if end_atom_indices[0][1] == a:
                end_atom_index = end_atom_indices[0][1]
                f_connec_to_end_atom_index = np.where(all_bonds_graphene == end_atom_index)[0]
                core_connections = all_bonds_graphene[f_connec_to_end_atom_index]
                Both_connected_atoms = np.setdiff1d(core_connections, [end_atom_index])
                right_connection = Both_connected_atoms[0]
                left_connection = Both_connected_atoms[1]

                first_vector = xyz[end_atom_index] - xyz[right_connection]
                second_vector = xyz[end_atom_index] - xyz[left_connection]
                H_coord = (first_vector + second_vector) + xyz[end_atom_index]
                H_coordinaes.extend([H_coord])
                num_hydrogen = num_hydrogen + 1

    if num_hydrogen > 0:
        h = MDAnalysis.Universe.empty(num_hydrogen, trajectory=True, n_residues=10)
        coord_array_H_indice = np.array(H_coordinaes)
        assert coord_array_H_indice.shape == (num_hydrogen, 3)
        h.atoms.positions = coord_array_H_indice
        h.add_TopologyAttr('name', ['H']*num_hydrogen)
    graphene_with_hydrogen = MDAnalysis.Merge(graphene.atoms, h.atoms)
    graphene_with_hydrogen.add_bonds(all_bonds_graphene)
    num_hydrogen = 0
    for x in range(num_atoms_graphene):
        indices_of_a = np.where(all_bonds_graphene == x)
        if len(indices_of_a[0]) == 1:
            graphene_with_hydrogen.add_bonds([(x, num_atoms_graphene + num_hydrogen)])
            graphene_with_hydrogen.add_bonds([(x, num_atoms_graphene + num_hydrogen + 1)])
            num_hydrogen = num_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            graphene_with_hydrogen.add_bonds([(x, num_atoms_graphene + num_hydrogen)])
            num_hydrogen = num_hydrogen + 1

    # If the user chooses "All", hydrogenated graphene structure will be returned:
    if H_termination_graphene == 'All':
        return graphene_with_hydrogen
