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
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from furiousatoms.io import create_universe, merged_universe_with_H
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

thre = 1e-10
vacuum = 4
SM = SharedMemory()
"""
    Ui_Graphene class creates a widget for building graphenes
"""
class Ui_graphene(QtWidgets.QMainWindow): #QWidget

    def __init__(self, app_path=None, parent=None):
        super(Ui_graphene, self).__init__(parent)
        self.graphene = io.load_ui_widget("graphene.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.graphene)
        self.setCentralWidget(self.graphene)
        self.setLayout(self.v_layout)
        self.resize(244, 226)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        bond_length_graphene = 1.421 # default value of C-C bond length
        self.graphene.lineEdit_bond_length_graphene.insert(str(bond_length_graphene))
        self.graphene.pushButton_build_graphene.clicked.connect(self.graphene_builder_callback)
        self.graphene.pushButton_build_graphene.clicked.connect(lambda:self.close())

    def graphene_builder_callback(self):
        H_termination_graphene = self.graphene.comboBox_H_termination_graphene.currentText()
        bond_length_graphene = self.graphene.lineEdit_bond_length_graphene.text()
        bond_length_graphene = float(bond_length_graphene)
        value_n_graphene = int(self.graphene.spinBox_chirality_N_graphene.text())
        value_m_graphene = int(self.graphene.spinBox_chirality_M_graphene.text())
        repeat_units_graphene = int(self.graphene.spinBox_repeat_units_graphene.text())
        graphene_type_1 = self.graphene.comboBox_type1_graphene.currentText()
        graphene_type_2 = self.graphene.comboBox_type2_graphene.currentText()
        structure_info = graphene_builder(H_termination_graphene, value_n_graphene, value_m_graphene, repeat_units_graphene, length=None, a=bond_length_graphene, species=(graphene_type_1, graphene_type_2), centered=True)
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(structure_info)
        window.show()

"""
  The numbers (n,m) show that your tube is obtained from taking one atom of the sheet and rolling it onto
  that atom that is at located na1+ma2 away from your original atom. N is the Number of hexagons in a unit cell. a1 and a2 are lattice vectors.
  (n,m=n) gives an “armchair” tube,e.g. (5,5). (n,m=0) gives an “zig-zag” tube, e.g. (6,0). Other tubes are “chiral”, e.g. (6,2)
"""


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
    all_bonds_graphene = np.array(fragments['bonds'])
    from furiousatoms.io import create_universe
    univ_graphene = create_universe(coord_array_graphene, all_bonds_graphene, atom_types_graphene)
    if H_termination_graphene == 'None':
        return univ_graphene
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

    coord_array_H_indice = np.array(H_coordinaes)
    num_hydrogen = 0
    bonds_hydrogen = []
    for x in range(num_atoms_graphene):
        indices_of_a = np.where(all_bonds_graphene == x)
        if len(indices_of_a[0]) == 1:
            bonds_hydrogen.extend([(x, num_atoms_graphene + num_hydrogen)])
            bonds_hydrogen.extend([(x, num_atoms_graphene + num_hydrogen + 1)])
            num_hydrogen = num_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            bonds_hydrogen.extend([(x, num_atoms_graphene + num_hydrogen)])
            num_hydrogen = num_hydrogen + 1
    atom_types_Hydrogen = list(['H']*num_hydrogen)
    merged_graphene_hydrogen = merged_universe_with_H(coord_array_graphene,  all_bonds_graphene, atom_types_graphene, coord_array_H_indice, bonds_hydrogen, atom_types_Hydrogen)


    # If the user chooses "All", hydrogenated graphene structure will be returned:
    if H_termination_graphene == 'All':
        return merged_graphene_hydrogen
