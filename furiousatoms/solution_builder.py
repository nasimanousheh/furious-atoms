import numpy as np
from numpy.linalg import norm
from math import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
import sys
from furiousatoms import io
import numpy as np
from fury import window
from PySide2 import QtWidgets
import MDAnalysis

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
        water_diameter = 3.16555789
        spacing_dia = 0.98
        spacing = spacing_dia * water_diameter
        total_water_inside = int((int(SM.box_lx/spacing) * int(SM.box_ly/spacing) * int(SM.box_lz/spacing)))
        water_amount = (SM.box_lx/spacing * SM.box_ly/spacing * SM.box_lz/spacing)
        water_concentration = water_amount / (0.602214076 * (SM.box_lx) * (SM.box_ly) * (SM.box_lz-water_diameter) * 0.001)
        print("water_concentration is: ", water_concentration)

        h2o = np.array([[ 0,        0,       0      ],  # oxygen
                        [ 0.95908, -0.02691, 0.03231],  # hydrogen
                        [-0.28004, -0.58767, 0.70556]]) # hydrogen
        coordinates = []
        coor=[]
        num_molecul_in_box_lx = int(SM.box_lx / spacing)
        num_molecul_in_box_ly = int(SM.box_ly / spacing)
        num_molecul_in_box_lz = int(SM.box_lz / spacing)
        mol = 0
        for i in range(num_molecul_in_box_lx):
            for j in range(num_molecul_in_box_ly):
                for k in range(num_molecul_in_box_lz):
                    if (mol < (total_water_inside*3)):
                        x = (-SM.box_lx/2 + 0.5 * (spacing)) + (i * spacing)
                        y = (-SM.box_ly/2 + 0.5 * (spacing)) + (j * spacing)
                        z = (-SM.box_lz/2 + 0.5 * (spacing)) + (k * spacing)
                        if (x > (SM.box_lx/2 - (0.5 * spacing)) or y > (SM.box_ly/2 - (0.5 * spacing)) or z > (SM.box_lz/2 - (0.5 * spacing))):
                           continue
                        xyz = np.array([x, y, z])
                        close = False
                        for p_n in range(len(SM.pos)):
                            dist = np.linalg.norm(xyz - SM.pos[p_n])
                            if dist < 2:
                                close = True
                                print("close")
                                break
                        if not close:
                            coordinates.extend(h2o + xyz.T)
                    mol = mol + 1

        coord_array = np.array(coordinates)
        n_atoms = len(coord_array)
        assert coord_array.shape == (n_atoms, 3)
        # n_residues = total_water_inside
        n_residues = int(n_atoms/3)
        # n_atoms = n_residues * 3
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

        sol.atoms.positions = coord_array
        assert not hasattr(sol, 'bonds')
        bonds = []
        for o in range(0, n_atoms, 3):
            bonds.extend([(o, o+1), (o, o+2)])
        sol.add_TopologyAttr('bonds', bonds)
        combined = MDAnalysis.Merge(SM.universe.atoms, sol.atoms)
        import tempfile
        dir_name = tempfile.mkdtemp(prefix='Furious_Atoms_')
        file_name = tempfile.mkstemp(suffix='.pdb', prefix='Solution', dir=dir_name)[1]
        combined.atoms.write(file_name)
        combined.universe.trajectory.ts.dimensions[0] = SM.box_lx
        combined.universe.trajectory.ts.dimensions[1] = SM.box_ly
        combined.universe.trajectory.ts.dimensions[2] = SM.box_lz
        SM = active_window.universe_manager
        active_window.scene.rm(SM.sphere_actor)
        active_window.scene.rm(SM.bond_actor)
        active_window.load_universe(combined)
        active_window.render()
        return file_name