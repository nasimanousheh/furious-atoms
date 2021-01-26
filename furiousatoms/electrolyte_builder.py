from MDAnalysis import *
import MDAnalysis as mda
import MDAnalysis.analysis.align
import numpy as np
import sys
import io
from furiousatoms import io
import vtk
import numpy as np
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from furiousatoms.sharedmem import SharedMemory
from furiousatoms.periodic_table import Ui_periodic

import sys

class Ui_electrolyte(QtWidgets.QMainWindow): #QWidget
    """ Ui_electrolyte class creates a widget for building electrolyte
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_electrolyte, self).__init__(parent)
        self.electrolyte = io.load_ui_widget("electrolyte.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.electrolyte)
        self.setCentralWidget(self.electrolyte)
        self.setLayout(self.v_layout)
        self.resize(620, 699)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.electrolyte.pushButton_cation_1.clicked.connect(self.get_periodic_table)

    def get_periodic_table(self):
        Ui_periodic.swnt = Ui_periodic()
        Ui_periodic.swnt.win = self
        Ui_periodic.swnt.show()
        Ui_periodic.swnt.showNormal()



# def Create_water(box_lx, box_ly, box_lz, water_diameter, spacing_dia, charge_hydrogen, charge_oxygen, charge_density, counter_type, counter_charge, type_cat_1, type_cat_2, type_cat_3, type_cat_4, type_cat_5,
#                 charge_cat_1, charge_cat_2, charge_cat_3, charge_cat_4, charge_cat_5, con_cat_1, con_cat_2, con_cat_3, con_cat_4, con_cat_5, type_an_1, type_an_2, type_an_3, type_an_4,
#                 type_an_5, charge_an_1, charge_an_2, charge_an_3, charge_an_4, charge_an_5, con_an_1, con_an_2, con_an_3, con_an_4, con_an_5):

    # box_lx = 50
    # box_ly = 50
    # box_lz = 30
    # water_diameter = 3.16555789
    # spacing_dia = 1
    # # charge of oxygen and hydrogen in SPC/E model:
    # charge_hydrogen = -0.8476
    # charge_oxygen = 0.4238
    # charge_density = -0.01
    # counter_type = 'Na'
    # counter_charge = 1

    # charge_cat_1 = 1
    # charge_cat_2 = 1
    # charge_cat_3 = 1
    # charge_cat_4 = 1
    # charge_cat_5 = 1
    # type_cat_1 = 'Na'
    # type_cat_2 = 'K'
    # type_cat_3 = 'Cs'
    # type_cat_4 = 'Rb'
    # type_cat_5 = 'H'
    # con_cat_1 = 0.1
    # con_cat_2 = 0.1
    # con_cat_3 = 0.1
    # con_cat_4 = 0.1
    # con_cat_5 = 0.1

    # charge_an_1 = -1
    # charge_an_2 = -1
    # charge_an_3 = -1
    # charge_an_4 = -1
    # charge_an_5 = -1
    # type_an_1 = 'Cl'
    # type_an_2 = 'B'
    # type_an_3 = 'I'
    # type_an_4 = 'F'
    # type_an_5 = 'At'
    # con_an_1 = 0.1
    # con_an_2 = 0.1
    # con_an_3 = 0.1
    # con_an_4 = 0.1
    # con_an_5 = 0.1

    # spacing = spacing_dia * water_diameter
    # volume_box = box_lx*box_ly*box_lz*0.001
    # total_water_inside = int(int(box_lx/spacing) * int(box_ly/spacing) * int(box_lz/spacing))
    # water_amount = float(box_lx/spacing * box_ly/spacing * box_lz/spacing)
    # water_concentration = water_amount / (0.602214076 * (box_lx) * (box_ly) * (box_lz-water_diameter) * 0.001)

    # total_pions_1 = int((con_cat_1 * 0.6022) * (volume_box))
    # total_pions_2 = int((con_cat_2 * 0.6022) * (volume_box))
    # total_pions_3 = int((con_cat_3 * 0.6022) * (volume_box))
    # total_pions_4 = int((con_cat_4 * 0.6022) * (volume_box))
    # total_pions_5 = int((con_cat_5 * 0.6022) * (volume_box))
    # total_pions_inside = int(total_pions_1 + total_pions_2 + total_pions_3 + total_pions_4 + total_pions_5)

    # total_nions_1 = int((con_an_1 * 0.6022) * (volume_box))
    # total_nions_2 = int((con_an_2 * 0.6022) * (volume_box))
    # total_nions_3 = int((con_an_3 * 0.6022) * (volume_box))
    # total_nions_4 = int((con_an_4 * 0.6022) * (volume_box))
    # total_nions_5 = int((con_an_5 * 0.6022) * (volume_box))
    # total_nions_inside = int(total_nions_1 + total_nions_2 + total_nions_3 + total_nions_4 + total_nions_5)
    # total_charge_pions = int((con_cat_1 * total_pions_1) + (con_cat_2 * total_pions_2) + (con_cat_3 * total_pions_3) + (con_cat_4 * total_pions_4) + (con_cat_5 * total_pions_5))
    # total_charge_nions = int((con_an_1 * total_nions_1) + (con_an_2 * total_nions_2) + (con_an_3 * total_nions_3) + (con_an_4 * total_nions_4) + (con_an_5 * total_nions_5))
    # if (total_charge_pions != total_charge_nions):
    #     print('The electrolyte is not electroneutral; Abortion')
    #     return

    # total_saltions_inside = int(total_nions_inside + total_pions_inside)

    # n_residues = total_water_inside
    # n_atoms = n_residues * 3
    # resindices = np.repeat(range(n_residues), 3)
    # assert len(resindices) == n_atoms

    # # all water molecules belong to 1 segment
    # segindices = [0] * n_residues
    # sol = mda.Universe.empty(n_atoms,
    #                         n_residues=n_residues,
    #                         atom_resindex=resindices,
    #                         residue_segindex=segindices,
    #                         trajectory=True)
    # sol.dimensions = [lx, ly, lz, 0, 0, 0]
    # sol.add_TopologyAttr('name', ['O', 'H1', 'H2']*n_residues)
    # sol.add_TopologyAttr('type', ['O', 'H', 'H']*n_residues)
    # sol.add_TopologyAttr('resname', ['SOL']*n_residues)
    # sol.add_TopologyAttr('resid', list(range(1, n_residues+1)))
    # sol.add_TopologyAttr('segid', ['SOL'])

    # h2o = np.array([[ 0,        0,       0      ],  # oxygen
    #                 [ 0.95908, -0.02691, 0.03231],  # hydrogen
    #                 [-0.28004, -0.58767, 0.70556]]) # hydrogen

    # coordinates = []
    # num_molecul_in_lx = int(lx / spacing)
    # num_molecul_in_ly = int(ly / spacing)
    # num_molecul_in_lz = int(lz / spacing)
    # for i in range(num_molecul_in_lx):
    #     for j in range(num_molecul_in_ly):
    #         for k in range(num_molecul_in_lz):
    #             x = (-lx/2 + (0.5*spacing)) + (i * spacing)
    #             y = (-ly/2 + (0.5*spacing)) + (j * spacing)
    #             z = (-lz/2 + (0.5*spacing)) + (k * spacing)
    #             xyz = np.array([x, y, z])
    #             coordinates.extend(h2o + xyz.T)
    # coord_array = np.array(coordinates)

    # assert coord_array.shape == (n_atoms, 3)
    # sol.atoms.positions = coord_array
    # assert not hasattr(sol, 'bonds')
    # bonds = []
    # for o in range(0, n_atoms, 3):
    #     bonds.extend([(o, o+1), (o, o+2)])
    # sol.add_TopologyAttr('bonds', bonds)

    # # water_center = sol.center_of_mass(pbc=True)
    # # dim = sol.dimensions
    # # box_center = np.sum(dim, axis=0)
    # # sol.atoms.translate(-lx/2 + spacing/2)
    # sol = MDAnalysis.Universe('C:/Users/nasim/OneDrive/Desktop/onlywater.pdb')
    # salt = MDAnalysis.Universe('C:/Users/nasim/OneDrive/Desktop/Fullerene_C720.pdb')


    # cog = sol.atoms.center_of_geometry()
    # print('Original solvent center of geometry: ', cog)
    # sol.atoms.positions -= cog
    # cog2 = sol.atoms.center_of_geometry()
    # print('New solvent center of geometry: ', cog2)
    # cog = salt.atoms.center_of_geometry()
    # print('Original solvent center of geometry: ', cog)
    # salt.atoms.positions -= cog
    # cog2 = salt.atoms.center_of_geometry()
    # print('New solvent center of geometry: ', cog2)


    # combined = MDAnalysis.Merge(salt.atoms, sol.atoms)
    # # u.select_atoms('same resid as (not around 30 salt)')
    # # point 5.0 5.0 5.0 3.5
    # no_overlap = combined.select_atoms("same resid as (not around 5 name C B)")

    # # no_overlap = combined.select_atoms("same resid as (not around 10 protein)")
    # u = mda.Merge(no_overlap)
    # print(len(u.atoms))
    # u.atoms.write('C:/Users/nasim/OneDrive/Desktop/solution.pdb')


