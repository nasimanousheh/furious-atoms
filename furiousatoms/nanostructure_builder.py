#!/usr/bin/env python3
"""Nanotube generator.

For details and notation, see

http://www.photon.t.u-tokyo.ac.jp/~maruyama/kataura/chirality.html

Usage:
    ./nanotuby.py N M LENGTH [-a DIST] [--species SP] [--centered]

Options:
    -a DIST                    Atom-atom distance [default: 1.421].
    --species SP               Comma-separated pair of species to use [default: C,C].
    --centered                 Center nanotube in unit cell.
"""
from MDAnalysis import *
import MDAnalysis.analysis.align
import numpy as np
from numpy.linalg import norm
from fractions import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
from furiousatoms.forms.widget_SWNT import Ui_Form_SWNT
import math
import sys
import io

thre = 1e-10
vacuum = 4
SM = SharedMemory()

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
        N = int(np.ceil(length/norm(T)))
    Ch_proj, T_proj = [v/norm(v)**2 for v in [Ch, T]]
    basis = [np.array((0, 0)), (a1+a2)/3]
    pts = []
    for i1, i2 in product(range(0, n+t1+1), range(t2, m+1)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            if all(-thre < pt.dot(v) < 1-thre for v in [Ch_proj, T_proj]):
                for k in range(N):
                    pts.append((sp, pt+k*T))
    SM.diameter_SWNT = norm(Ch)/np.pi
    def gr2tube(v):
        phi = 2*np.pi*v.dot(Ch_proj)
        return np.array((SM.diameter_SWNT/2*np.cos(phi),
                         SM.diameter_SWNT/2*np.sin(phi),
                         v.dot(T_proj)*norm(T)))
    xyz = [gr2tube(v) for _, v in pts]
    atom_types_swnt = [v for v, _ in pts]
    m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    fragments = m.to_json(scale =1)
    lx = SM.diameter_SWNT + 1
    ly = SM.diameter_SWNT + 1
    lz = 30+SM.diameter_SWNT/2
    n_atoms_swnt = len(xyz)
    coord_array_swnt = np.array(xyz) #+ np.array((0, 0, int(lz/2)))
    assert coord_array_swnt.shape == (n_atoms_swnt, 3)
    swnt = MDAnalysis.Universe.empty(n_atoms_swnt, trajectory=True, n_residues=10)
    swnt.atoms.positions = coord_array_swnt
    all_bonds_swnt = np.array(fragments['bonds'])
    swnt.add_TopologyAttr('name', atom_types_swnt)
    swnt.add_bonds(all_bonds_swnt)
    if H_termination_SWNT == 'None':
        return swnt
    #########################################################
    half_n_atoms_swnt = int(n_atoms_swnt/2)
    n_hydrogen = 0
    H_coordinaes = []
    for a in range(half_n_atoms_swnt):
        indices_of_a = np.where(all_bonds_swnt == a)[0]
        if len(indices_of_a) == 1:
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
            n_hydrogen = n_hydrogen + 2
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
                n_hydrogen = n_hydrogen + 1

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
                n_hydrogen = n_hydrogen + 1

    if n_hydrogen > 0:
        h = MDAnalysis.Universe.empty(n_hydrogen, trajectory=True, n_residues=10)
        coord_array_H_indice = np.array(H_coordinaes)
        assert coord_array_H_indice.shape == (n_hydrogen, 3)
        h.atoms.positions = coord_array_H_indice
        h.add_TopologyAttr('name', ['H']*n_hydrogen)
    combined_one_end = MDAnalysis.Merge(swnt.atoms, h.atoms)
    combined_one_end.add_bonds(all_bonds_swnt)
    n_hydrogen = 0
    for x in range(half_n_atoms_swnt):
        indices_of_a = np.where(all_bonds_swnt == x)
        if len(indices_of_a[0]) == 1:
            combined_one_end.add_bonds([(x, n_atoms_swnt + n_hydrogen)])
            combined_one_end.add_bonds([(x, n_atoms_swnt + n_hydrogen + 1)])
            n_hydrogen = n_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            combined_one_end.add_bonds([(x, n_atoms_swnt + n_hydrogen)])
            n_hydrogen = n_hydrogen + 1

    if H_termination_SWNT == 'One end':
        return combined_one_end

#########################################################
    n_hydrogen = 0
    H_coordinaes = []
    for a in range(n_atoms_swnt):
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
            n_hydrogen = n_hydrogen + 2

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
                n_hydrogen = n_hydrogen + 1

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
                n_hydrogen = n_hydrogen + 1

    if n_hydrogen > 0:
        h = MDAnalysis.Universe.empty(n_hydrogen, trajectory=True, n_residues=10)
        coord_array_H_indice = np.array(H_coordinaes)
        assert coord_array_H_indice.shape == (n_hydrogen, 3)
        h.atoms.positions = coord_array_H_indice
        h.add_TopologyAttr('name', ['H']*n_hydrogen)
    combined = MDAnalysis.Merge(swnt.atoms, h.atoms)
    combined.add_bonds(all_bonds_swnt)
    n_hydrogen = 0
    for x in range(n_atoms_swnt):
        indices_of_a = np.where(all_bonds_swnt == x)
        if len(indices_of_a[0]) == 1:
            combined.add_bonds([(x, n_atoms_swnt + n_hydrogen)])
            combined.add_bonds([(x, n_atoms_swnt + n_hydrogen + 1)])
            n_hydrogen = n_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            combined.add_bonds([(x, n_atoms_swnt + n_hydrogen)])
            n_hydrogen = n_hydrogen + 1

    combined.trajectory.ts.dimensions[0] = lx
    combined.trajectory.ts.dimensions[1] = ly
    combined.trajectory.ts.dimensions[2] = lz
    # print(SM.H_termination_SWNT)
    if H_termination_SWNT == 'Both ends':
        return combined
    # combined.atoms.write('C:/Users/nasim/Devel/furious-atoms/nanotube_structure_both_ends.pdb')

def MWNT_nanotube_builder(n, m, N=1, length=None, a=1.421, species=('C', 'C'), centered=False):
    SM = SharedMemory()
    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2))
    Ch = n*a1+m*a2
    T = t1*a1+t2*a2
    if length:
        N = int(np.ceil(length/norm(T)))
    Ch_proj, T_proj = [v/norm(v)**2 for v in [Ch, T]]
    basis = [np.array((0, 0)), (a1+a2)/3]
    pts = []
    for i1, i2 in product(range(0, n+t1+1), range(t2, m+1)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            if all(-thre < pt.dot(v) < 1-thre for v in [Ch_proj, T_proj]):
                for k in range(N):
                    pts.append((sp, pt+k*T))
    SM.diameter_MWNT = norm(Ch)/np.pi
    def gr2tube(v):
        phi = 2*np.pi*v.dot(Ch_proj)
        return np.array((SM.diameter_SWNT/2*np.cos(phi),
                         SM.diameter_MWNT/2*np.sin(phi),
                         v.dot(T_proj)*norm(T)))
    xyz = [gr2tube(v) for _, v in pts]
    atom_types = [v for v, _ in pts]
    m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    fragments = m.to_json(scale =1)
    n_atoms = len(xyz)
    SM.pos = np.array(xyz)
    assert SM.pos.shape == (n_atoms, 3)
    u = MDAnalysis.Universe.empty(n_atoms, trajectory=True, n_residues=10)
    u.atoms.positions = SM.pos
    SM.bonds = np.array(fragments['bonds'])
    u.add_bonds(SM.bonds)
    lx = 0
    ly = 0
    lz = 0
    u.trajectory.ts.dimensions[0] = lx
    u.trajectory.ts.dimensions[1] = ly
    u.trajectory.ts.dimensions[2] = lz
    u.add_TopologyAttr('name', atom_types)
    u.atoms.write('C:/Users/nasim/Devel/furious-atoms/MWNT_nanotube_structure.pdb')

def graphene_builder(n, m, N=1, length=None, a=1.421, species=('C', 'C'), centered=False):
    # SM = SharedMemory()
    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0,0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2,0))
    Ch = n*a1+m*a2
    T = t1*a1+t2*a2
    if length:
        N = int(np.ceil(length/norm(T)))
    Ch_proj, T_proj = [v/norm(v)**2 for v in [Ch, T]]
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
    lx = 0
    ly = 0
    lz = 0
    n_atoms_graphene = len(xyz)
    coord_array_graphene = np.array(xyz)
    assert coord_array_graphene.shape == (n_atoms_graphene, 3)
    graphene = MDAnalysis.Universe.empty(n_atoms_graphene, trajectory=True, n_residues=10)
    graphene.atoms.positions = coord_array_graphene
    all_bonds_graphene = np.array(fragments['bonds'])
    graphene.trajectory.ts.dimensions[0] = lx
    graphene.trajectory.ts.dimensions[1] = ly
    graphene.trajectory.ts.dimensions[2] = lz
    graphene.add_TopologyAttr('name', atom_types_graphene)
    graphene.add_bonds(all_bonds_graphene)
    graphene.atoms.write('C:/Users/nasim/Devel/furious-atoms/graphene_structure.pdb')
#########################################################
    n_hydrogen = 0
    H_coordinaes = []
    for a in range(n_atoms_graphene):
        indices_of_a = np.where(all_bonds_graphene == a)[0]

        if len(indices_of_a) == 1:
            end_atom_indices = (all_bonds_graphene[indices_of_a])
            if end_atom_indices[0][0] == a:
                end_atom_indices = (all_bonds_graphene[indices_of_a]) #[2735 2447] main= 2735
                f_connec_to_end_atom_index = end_atom_indices[0][1] # 2447
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
            n_hydrogen = n_hydrogen + 2

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
                n_hydrogen = n_hydrogen + 1

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
                n_hydrogen = n_hydrogen + 1

    if n_hydrogen > 0:
        h = MDAnalysis.Universe.empty(n_hydrogen, trajectory=True, n_residues=10)
        coord_array_H_indice = np.array(H_coordinaes)
        assert coord_array_H_indice.shape == (n_hydrogen, 3)
        h.atoms.positions = coord_array_H_indice
        h.add_TopologyAttr('name', ['H']*n_hydrogen)
    combined = MDAnalysis.Merge(graphene.atoms, h.atoms)
    combined.add_bonds(all_bonds_graphene)
    n_hydrogen = 0
    for x in range(n_atoms_graphene):
        indices_of_a = np.where(all_bonds_graphene == x)
        if len(indices_of_a[0]) == 1:
            combined.add_bonds([(x, n_atoms_graphene + n_hydrogen)])
            combined.add_bonds([(x, n_atoms_graphene + n_hydrogen + 1)])
            n_hydrogen = n_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            combined.add_bonds([(x, n_atoms_graphene + n_hydrogen)])
            n_hydrogen = n_hydrogen + 1
###########################################
    combined.trajectory.ts.dimensions[0] = lx
    combined.trajectory.ts.dimensions[1] = ly
    combined.trajectory.ts.dimensions[2] = lz
    combined.atoms.write('C:/Users/nasim/Devel/furious-atoms/graphene_structure.pdb')
