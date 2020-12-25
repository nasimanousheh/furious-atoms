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
import sys
import io

thre = 1e-10
vacuum = 4


def nanotube_builder(n, m, N=1, length=None, a=1.421, species=('C', 'C'), centered=False):
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
    SM.diameter_tube = norm(Ch)/np.pi
    def gr2tube(v):
        phi = 2*np.pi*v.dot(Ch_proj)
        return np.array((SM.diameter_tube/2*np.cos(phi),
                         SM.diameter_tube/2*np.sin(phi),
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
    u.atoms.write('C:/Users/nasim/Devel/furious-atoms/nanotube_structure.pdb')
    if centered:
        m = m.shifted(np.array((-1, -1, 0))*SM.diameter_tube*(vacuum+1)/2)
    return Crystal([((vacuum+1)*SM.diameter_tube, 0, 0),
                    (0, (vacuum+1)*SM.diameter_tube, 0),
                    (0, 0, N*norm(T))],
                   m.atoms)

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
    SM.diameter_tube = norm(Ch)/np.pi
    def gr2tube(v):
        phi = 2*np.pi*v.dot(Ch_proj)
        return np.array((SM.diameter_tube/2*np.cos(phi),
                         SM.diameter_tube/2*np.sin(phi),
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
    if centered:
        m = m.shifted(np.array((-1, -1, 0))*SM.diameter_tube*(vacuum+1)/2)
    return Crystal([((vacuum+1)*SM.diameter_tube, 0, 0),
                    (0, (vacuum+1)*SM.diameter_tube, 0),
                    (0, 0, N*norm(T))],
                   m.atoms)

def graphene_builder(n, m, N=1, length=None, a=1.421, species=('C', 'C'), centered=False):
    SM = SharedMemory()
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
    basis = [np.array((0, 0,0)), (a1+a2)/3]
    pts = []
    for i1, i2 in product(range(0, n+t1+1), range(t2, m+1)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            if all(-thre < pt.dot(v) < 1-thre for v in [Ch_proj, T_proj]):
                for k in range(N):
                    pts.append((sp, pt+k*T))
    xyz = [v for _, v in pts]
    atom_types = [v for v, _ in pts]
    m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    lx = 0
    ly = 0
    lz = 0
    fragments = m.to_json(scale =1)
    n_atoms = len(xyz)
    SM.pos = np.array(xyz)
    assert SM.pos.shape == (n_atoms, 3)
    u = MDAnalysis.Universe.empty(n_atoms, trajectory=True, n_residues=1)
    u.trajectory.ts.dimensions[0] = lx
    u.trajectory.ts.dimensions[1] = ly
    u.trajectory.ts.dimensions[2] = lz
    u.atoms.positions = SM.pos
    SM.bonds = np.array(fragments['bonds'])
    u.add_bonds(SM.bonds)
    u.add_TopologyAttr('name', atom_types)
    u.atoms.write('C:/Users/nasim/Devel/furious-atoms/graphene_structure.pdb')
