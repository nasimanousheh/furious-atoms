import numpy as np
from numpy.linalg import norm
from fractions import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
import sys
from furiousatoms import io
import vtk
import numpy as np
from fury import window, actor, utils, pick, ui, primitive
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from furiousatoms.io import create_universe, merged_universe_with_H
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from furiousatoms.structure import bbox
import MDAnalysis
import os

SM = SharedMemory()
"""
    Ui_solution class creates a widget for building box and solution
"""

class Ui_solution(QtWidgets.QMainWindow): #QWidget

    def __init__(self, app_path=None, parent=None):
        super(Ui_solution, self).__init__(parent)
        self.solution = io.load_ui_widget("solution.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.solution)
        self.setCentralWidget(self.solution)
        self.setLayout(self.v_layout)
        self.resize(222, 202)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.solution.pushButton_build_solution.clicked.connect(self.solution_builder_callback)
        self.solution.pushButton_build_solution.clicked.connect(lambda:self.close())

    def solution_builder_callback(self):
        active_window = self.win.active_mdi_child()
        SM = active_window.universe_manager
        solution_diameter = 3.16555789
        spacing_dia = 1
        spacing = spacing_dia * solution_diameter
        total_solution_inside = int(int(SM.box_lx/spacing) * int(SM.box_ly/spacing) * int(SM.box_lz/spacing))
        n_residues = total_solution_inside
        n_atoms = n_residues * 3
        resindices = np.repeat(range(n_residues), 3)
        assert len(resindices) == n_atoms
        # all solution molecules belong to 1 segment
        segindices = [0] * n_residues
        sol = MDAnalysis.Universe.empty(n_atoms,
                                n_residues=n_residues,
                                atom_resindex=resindices,
                                residue_segindex=segindices,
                                trajectory=True)
        sol.dimensions = [SM.box_lx, SM.box_ly, SM.box_lz, 0, 0, 0]
        sol.add_TopologyAttr('name', ['O', 'H1', 'H2']*n_residues)
        sol.add_TopologyAttr('type', ['O', 'H', 'H']*n_residues)
        sol.add_TopologyAttr('resname', ['SOL']*n_residues)
        sol.add_TopologyAttr('resid', list(range(1, n_residues+1)))
        sol.add_TopologyAttr('segid', ['SOL'])
        h2o = np.array([[ 0,        0,       0      ],  # oxygen
                        [ 0.95908, -0.02691, 0.03231],  # hydrogen
                        [-0.28004, -0.58767, 0.70556]]) # hydrogen

        coordinates = []
        num_molecul_in_box_lx = int(SM.box_lx / spacing)
        num_molecul_in_box_ly = int(SM.box_ly / spacing)
        num_molecul_in_box_lz = int(SM.box_lz / spacing)
        for i in range(num_molecul_in_box_lx):
            for j in range(num_molecul_in_box_ly):
                for k in range(num_molecul_in_box_lz):
                    x = (-SM.box_lx/2 + (0.5*spacing)) + (i * spacing)
                    y = (-SM.box_ly/2 + (0.5*spacing)) + (j * spacing)
                    z = (-SM.box_lz/2 + (0.5*spacing)) + (k * spacing)
                    xyz = np.array([x, y, z])
                    coordinates.extend(h2o + xyz.T)
        coord_array = np.array(coordinates)

        assert coord_array.shape == (n_atoms, 3)
        sol.atoms.positions = coord_array
        assert not hasattr(sol, 'bonds')
        bonds = []
        for o in range(0, n_atoms, 3):
            bonds.extend([(o, o+1), (o, o+2)])
        sol.add_TopologyAttr('bonds', bonds)
        combined = MDAnalysis.Merge(SM.universe.atoms, sol.atoms)
        atom_type = ' '.join(SM.unique_types.tolist())
        combined.select_atoms("same resid as (not around 5 name C)")
        print('atom type is: ', atom_type)
        # no_overlap.atoms.write('C:/Users/nasim/OneDrive/Desktop/cnvert/hola.pdb')
        combined.universe.trajectory.ts.dimensions[0] = SM.box_lx
        combined.universe.trajectory.ts.dimensions[1] = SM.box_ly
        combined.universe.trajectory.ts.dimensions[2] = SM.box_lz
        SM = active_window.universe_manager
        active_window.scene.rm(SM.sphere_actor)
        active_window.scene.rm(SM.bond_actor)
        active_window.load_universe(combined)
        active_window.render()
        return
        # import time
# start = time.time()
# load_data
# end = time.time() -start
# print(end)
# #