import re
import os
import string
from MDAnalysis import *
import MDAnalysis.analysis.align
import numpy as np
from numpy.linalg import norm
from fractions import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
import sys
from furiousatoms import io
import vtk
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


def load_CC1_file(fname, debug=False):
    lammps_file = open(fname, 'r')
    lines = lammps_file.readlines()
    no_lines = len(lines)
    lammps_dix = {'index': {}}
    frames_cnt = 0
    no_atoms = int(len(lines)-1)
    for i in range(no_lines-1):
        line = lines[i]
        lammps_dix['index'][frames_cnt] = {'no_atoms': no_atoms, 'coords': [], 'bonds':[]}
        words = lines[i].split()
        no_columns = 9
        frame_as_list_text = lines[i+1: i+1 + no_atoms]
        frame_positions = np.genfromtxt(frame_as_list_text)
        lammps_dix['index'][frames_cnt]['coords'] = frame_positions
        frames_cnt += 1

    frames_cnt= 0
    no_atoms = lammps_dix['index'][frames_cnt]['no_atoms']
    pos =  lammps_dix['index'][frames_cnt]['coords'][:, 2:5]
    bonds = lammps_dix['index'][frames_cnt]['coords'][:, 6:9]
    suffix = '.pdb'
    # dir_name = 'C:/Users/nasim/OneDrive/Desktop/cnvert/fullerene.pdb'
    # out_put_pdb = 'C:/Users/nasim/OneDrive/Desktop/cnvert/fullerene.pdb'
    # out_put_pdb = os.path.join(dir_name, fname + suffix)#os.path.join(dir_name, fname + suffix)

    file_name = 'fullerene.pdb'
    outdump = open(file_name, "w")
    outdump.write("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P1\n")
    for i in range(no_atoms):
        tmp="{:6s}{:5d}  {:5s}{:2s}{:3s} {:2s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format('ATOM',1 + i, 'C','C','A','1', float("{:.3f}".format(pos[i][0])), float("{:.3f}".format(pos[i][1])),float("{:.3f}".format(pos[i][2])),1,0)
        outdump.write(tmp)
    outdump.write("TER\n")
    for j in range(no_atoms):
        tmp_bond="{:4s}{:5d} {:4d} {:4d} {:4d}\n".format('CONECT',1 + j, int(bonds[j][0]), int(bonds[j][1]), int(bonds[j][2]))
        outdump.write(tmp_bond)
    outdump.write("{:4s} {:5d} {:5d} {:5d} {:5d}\n".format('MASTER        0    0    0    0    0    0    0    0', no_atoms, 0, no_atoms, 0))
    outdump.write("END\n")


    # window = self.win.create_mdi_child()
    # window.load_file(fname=file_name)
    # window.show()
    return file_name


# def load_CC1_file(fname, debug=False):
#     lammps_file = open(fname, 'r')
#     lines = lammps_file.readlines()
#     no_lines = len(lines)
#     lammps_dix = {'index': {}}
#     frames_cnt = 0
#     no_atoms = int(len(lines)-1)
#     for i in range(no_lines-1):
#         line = lines[i]
#         lammps_dix['index'][frames_cnt] = {'no_atoms': no_atoms, 'coords': [],
#                                            'bonds': []}
#         words = lines[i].split()
#         no_columns = 9
#         frame_as_list_text = lines[i+1: i+1 + no_atoms]
#         frame_positions = np.genfromtxt(frame_as_list_text)
#         lammps_dix['index'][frames_cnt]['coords'] = frame_positions
#         frames_cnt += 1

#     frames_cnt = 0
#     pos_fullerene = lammps_dix['index'][frames_cnt]['coords'][:, 2:5]
#     bonds_fullerene = lammps_dix['index'][frames_cnt]['coords'][:, 6:9]
#     bonds_fullerene = bonds_fullerene.astype('i8') - 1
#     # bonds_fullerene_2 = lammps_dix['index'][frames_cnt]['coords'][:, 7:9]
#     # print(bonds_fullerene[:, :2].astype('i8'))
#     # bonds = lammps_dix['index'][frames_cnt]['coords'][:, 6:9]
#     bf_n = [(row[0], row[1]) for row in bonds_fullerene]
#     bf_n += [(row[1], row[2]) for row in bonds_fullerene]

#     print(bf_n)
#     n_residues = 1
#     fullerene = MDAnalysis.Universe.empty(no_atoms, trajectory=True,
#                                           n_residues=1)
#     fullerene.atoms.positions = pos_fullerene
#     fullerene.add_TopologyAttr('name', ['C']*no_atoms)
#     fullerene.add_TopologyAttr('type', ['C']*no_atoms)
#     fullerene.add_TopologyAttr('resname', ['MOL']*n_residues)
#     fullerene.add_bonds(np.array(bf_n))
#     # fullerene.atoms.indices[:] = fullerene.atoms.indices[:] + 1
#     # fullerene.add_TopologyAttr('bonds', bonds_fullerene)#[:, :2])
#     # self.win.process_universe(fullerene)
#     return fullerene



#     # out_put_pdb = 'C:/Users/nasim/OneDrive/Desktop/cnvert/fullerene.pdb'
#     # outdump = open(out_put_pdb, "w")
#     # outdump.write("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P1\n")
#     # for i in range(no_atoms):
#     #     tmp="{:6s}{:5d}  {:5s}{:2s}{:3s} {:2s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format('ATOM',1 + i, 'C','C','A','1',float("{:.3f}".format(pos[i][0])),float("{:.3f}".format(pos[i][1])),float("{:.3f}".format(pos[i][2])),1,0)
#     #     outdump.write(tmp)
#     # outdump.write("TER\n")
#     # for j in range(no_atoms):
#     #     tmp_bond="{:4s}{:5d} {:4d} {:4d} {:4d}\n".format('CONECT',1 + j, int(bonds[j][0]),int(bonds[j][1]),int(bonds[j][2]))
#     #     outdump.write(tmp_bond)
#     # outdump.write("{:4s} {:5d} {:5d} {:5d} {:5d}\n".format('MASTER        0    0    0    0    0    0    0    0',no_atoms,0,no_atoms,0))
#     # outdump.write("END\n")
#     # return out_put_pdb
