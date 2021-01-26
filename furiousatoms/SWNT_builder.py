from MDAnalysis import *
import MDAnalysis.analysis.align
import numpy as np
from numpy.linalg import norm
from fractions import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
from furiousatoms import io
import vtk
import numpy as np
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import sys


thre = 1e-10
vacuum = 4
SM = SharedMemory()
"""
    Ui_SWNT class creates a widget for building single-wall nanotube (SWNT)

"""

class Ui_SWNT(QtWidgets.QMainWindow): #QWidget
    def __init__(self, app_path=None, parent=None):
        super(Ui_SWNT, self).__init__(parent)
        self.SWNT = io.load_ui_widget("SWNT.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.SWNT)
        self.setCentralWidget(self.SWNT)
        self.setLayout(self.v_layout)
        self.resize(247, 285)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        bond_length_SWNT = 1.421 # default value of C-C bond length
        self.SWNT.lineEdit_bond_length_SWNT.insert(str(bond_length_SWNT))
        self.SWNT.lineEdit_bond_length_SWNT.textChanged.connect(self.SWNT_diameter_changed)
        self.SWNT.spinBox_chirality_N_SWNT.valueChanged.connect(self.SWNT_diameter_changed)
        self.SWNT.spinBox_chirality_M_SWNT.valueChanged.connect(self.SWNT_diameter_changed)
        self.SWNT.spinBox_repeat_units_SWNT.valueChanged.connect(self.SWNT_diameter_changed)
        self.SWNT.pushButton_build_SWNT.clicked.connect(self.SWNT_builder_callback)
        self.SWNT.comboBox_H_termination_SWNT.currentTextChanged.connect(self.SWNT_diameter_changed)

    def SWNT_diameter_changed(self):
        SM.H_termination_SWNT = self.SWNT.comboBox_H_termination_SWNT.currentText()
        bond_length_SWNT = self.SWNT.lineEdit_bond_length_SWNT.text()
        SM.bond_length_SWNT = float(bond_length_SWNT)
        value_n_SWNT = int(self.SWNT.spinBox_chirality_N_SWNT.text())
        value_m_SWNT = int(self.SWNT.spinBox_chirality_M_SWNT.text())
        repeat_units_SWNT = int(self.SWNT.spinBox_repeat_units_SWNT.text())
        a1 = np.array((np.sqrt(3)*SM.bond_length_SWNT, 0))
        a2 = np.array((np.sqrt(3)/2*SM.bond_length_SWNT, -3*SM.bond_length_SWNT/2))
        Ch = value_n_SWNT*a1+value_m_SWNT*a2
        d = gcd(value_n_SWNT, value_m_SWNT)
        dR = 3*d if (value_n_SWNT-value_m_SWNT) % (3*d) == 0 else d
        t1 = (2*value_m_SWNT+value_n_SWNT)//dR
        t2 = -(2*value_n_SWNT+value_m_SWNT)//dR
        T = t1*a1+t2*a2
        diameter_SWNT = float(np.linalg.norm(Ch)/np.pi)
        diameter_SWNT = "{:.2f}".format(diameter_SWNT)
        length_SWNT = np.linalg.norm(T) * repeat_units_SWNT
        length_SWNT = "{:.2f}".format(float(length_SWNT))
        self.SWNT.lineEdit_diameter_SWNT.setText(str(diameter_SWNT))
        self.SWNT.lineEdit_length_SWNT.setText(str(length_SWNT))

    def SWNT_builder_callback(self):
        value_n_SWNT = int(self.SWNT.spinBox_chirality_N_SWNT.text())
        value_m_SWNT = int(self.SWNT.spinBox_chirality_M_SWNT.text())
        repeat_units_SWNT = int(self.SWNT.spinBox_repeat_units_SWNT.text())
        SWNT_type_1 = self.SWNT.comboBox_type1_SWNT.currentText()
        SWNT_type_2 = self.SWNT.comboBox_type2_SWNT.currentText()
        universe = SWNT_builder(SM.H_termination_SWNT, value_n_SWNT, value_m_SWNT, repeat_units_SWNT, length=None, a=SM.bond_length_SWNT, species=(SWNT_type_1, SWNT_type_2), centered=True)
        file_name = 'fname.pdb'
        universe.atoms.write(file_name)
        self.win.process_load_file(fname=file_name)

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
