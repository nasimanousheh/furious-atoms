import numpy as np
from numpy.linalg import norm
from math import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule
from furiousatoms import io
import numpy as np
from fury import window
from PySide2 import QtWidgets
from furiousatoms.io import create_universe, merged_universe_with_H

thre = 1e-10
vacuum = 4
"""
    Ui_SWNT class creates a widget for building single-wall nanotube (SWNT)

"""

class Ui_SWNT(QtWidgets.QMainWindow):
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
        self.SWNT.pushButton_build_SWNT.clicked.connect(lambda:self.close())
        self.SWNT.SpinBox_lx.valueChanged.connect(self.initial_box_dim)
        self.SWNT.SpinBox_ly.valueChanged.connect(self.initial_box_dim)
        self.SWNT.SpinBox_lz.valueChanged.connect(self.initial_box_dim)
        self.SWNT.radioButton_desired_bond_length.toggled.connect(self.get_atom_type)
        self.SWNT.radioButton_bond_length.toggled.connect(self.get_atom_type)
        self.SWNT.comboBox_type1_SWNT.activated.connect(self.get_atom_type)
        self.SWNT.comboBox_type2_SWNT.activated.connect(self.get_atom_type)
    def get_atom_type(self):
        global bond_length_SWNT, bendFactor
        if self.SWNT.radioButton_bond_length.isChecked() == True:
            self.SWNT.SpinBox_desired_bond_length.setEnabled(False)
            SWNT_type_1 = self.SWNT.comboBox_type1_SWNT.currentText()
            SWNT_type_2 = self.SWNT.comboBox_type2_SWNT.currentText()
            if SWNT_type_1=="C" and SWNT_type_2=="C":
                bond_length_SWNT = 1.421 # default value of C-C bond length
            if SWNT_type_1=="N" and SWNT_type_2=="B":
                bond_length_SWNT = 1.47 # default value of N-B bond length
            if SWNT_type_1=="N" and SWNT_type_2=="Ga":
                bond_length_SWNT = 1.95 # default value of N-Ga bond length
            if SWNT_type_1=="N" and SWNT_type_2=="Al":
                bond_length_SWNT = 1.83 # default value of N-Al bond length
            if SWNT_type_1=="P" and SWNT_type_2=="Al":
                bond_length_SWNT = 2.3 # default value of P-Al bond length
            if SWNT_type_1=="P" and SWNT_type_2=="Ga":
                bond_length_SWNT = 2.28 # default value of P-Ga bond length
            if SWNT_type_1=="P" and SWNT_type_2=="C":
                bond_length_SWNT = 1.87 # default value of P-C bond length
            if SWNT_type_1=="N" and SWNT_type_2=="C":
                bond_length_SWNT = 1.47 # default value of N-C bond length
            if SWNT_type_1=="C" and SWNT_type_2=="B":
                bond_length_SWNT = 1.56 # default value of C-B bond length
            if SWNT_type_1=="C" and SWNT_type_2=="Al":
                bond_length_SWNT = 2.0 # default value of C-Al bond length
            if SWNT_type_1=="C" and SWNT_type_2=="Ga":
                bond_length_SWNT = 2.46 # default value of P-B bond length
            if SWNT_type_1=="P" and SWNT_type_2=="B":
                bond_length_SWNT = 1.74 # default value of P-B bond length
            self.SWNT.lineEdit_bond_length_SWNT.setText(str(bond_length_SWNT))

        elif self.SWNT.radioButton_desired_bond_length.isChecked() == True:
            self.SWNT.lineEdit_bond_length_SWNT.setText(str(' '))
            self.SWNT.SpinBox_desired_bond_length.setEnabled(True)
            bond_length_SWNT = float(self.SWNT.SpinBox_desired_bond_length.text())
        else:
            bond_length_SWNT = 1.421

    def SWNT_diameter_changed(self):
        global bond_length_SWNT, bendFactor
        try:
            bond_length_SWNT
        except NameError:
            bond_length_SWNT = 1.421
        # bond_length_SWNT = self.SWNT.lineEdit_bond_length_SWNT.text()
        # bond_length_SWNT = float(bond_length_SWNT)
        value_n_SWNT = int(self.SWNT.spinBox_chirality_N_SWNT.text())
        value_m_SWNT = int(self.SWNT.spinBox_chirality_M_SWNT.text())
        repeat_units_SWNT = int(self.SWNT.spinBox_repeat_units_SWNT.text())
        a1 = np.array((np.sqrt(3)*bond_length_SWNT, 0))
        a2 = np.array((np.sqrt(3)/2*bond_length_SWNT, -3*bond_length_SWNT/2))
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
        global bond_length_SWNT, bendFactor
        try:
            bendFactor
        except NameError:
            bendFactor = 1.0
        try:
            bond_length_SWNT
        except NameError:
            bond_length_SWNT = 1.421
        H_termination_SWNT = self.SWNT.comboBox_H_termination_SWNT.currentText()
        value_n_SWNT = int(self.SWNT.spinBox_chirality_N_SWNT.text())
        value_m_SWNT = int(self.SWNT.spinBox_chirality_M_SWNT.text())
        repeat_units_SWNT = int(self.SWNT.spinBox_repeat_units_SWNT.text())
        SWNT_type_1 = self.SWNT.comboBox_type1_SWNT.currentText()
        SWNT_type_2 = self.SWNT.comboBox_type2_SWNT.currentText()
        bendFactor = float(self.SWNT.doubleSpinBox_bend_factor.text())
        structure_info = SWNT_builder(H_termination_SWNT, value_n_SWNT, value_m_SWNT, repeat_units_SWNT, length=None, bond_length=bond_length_SWNT, species=(SWNT_type_1, SWNT_type_2), centered=True, bend = bendFactor)
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(structure_info)
        window.show()
    def initial_box_dim(self):
        global box_lx, box_ly, box_lz
        box_lx = float(self.SWNT.SpinBox_lx.text())
        box_ly = float(self.SWNT.SpinBox_ly.text())
        box_lz = float(self.SWNT.SpinBox_lz.text())
        return

"""
  The numbers (n,m) show that your tube is obtained from taking one atom of the sheet and rolling it onto
  that atom that is at located na1+ma2 away from your original atom. N is the Number of hexagons in a unit cell. a1 and a2 are lattice vectors.
  (n,m=n) gives an “armchair” tube,e.g. (5,5). (n,m=0) gives an “zig-zag” tube, e.g. (6,0). Other tubes are “chiral”, e.g. (6,2)
"""

def SWNT_builder(H_termination_SWNT, n, m, N, length, bond_length, species=('C', 'C'), centered=False, bend = 1):
    global box_lx, box_ly, box_lz, bendFactor
    bond_length_hydrogen = 1.1
    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*bond_length, 0))
    a2 = np.array((np.sqrt(3)/2*bond_length, -3*bond_length/2))
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
    # diameter_SWNT = np.linalg.norm(Ch)/np.pi
    diameter_SWNT = np.linalg.norm(Ch)/(np.pi*bendFactor)
    # This function converts the SWNT to nanotube with given diameter:
    def gr2tube(v):
        phi = np.pi + bendFactor * 2*np.pi*v.dot(Ch_proj)

        return np.array((diameter_SWNT/2*np.cos(phi),
                         diameter_SWNT/2*np.sin(phi),
                         v.dot(T_proj)*np.linalg.norm(T)))
    xyz = [gr2tube(v) for _, v in pts]
    atom_types_swnt = [v for v, _ in pts]
    m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    fragments = m.to_json(scale=1)
    # Number of atoms in SWNT:
    num_atoms_swnt = len(xyz)
    n_residues = 1
    # Atom coordinates:
    coord_array_swnt = np.array(xyz)
    assert coord_array_swnt.shape == (num_atoms_swnt, 3)
    # Bonds information connected the atoms:
    all_bonds_swnt = np.array(fragments['bonds'])
    try:
        box_lx or box_ly or box_lz
    except NameError:
        box_lx = box_ly = box_lz = 0.0
    univ_swnt = create_universe(coord_array_swnt, all_bonds_swnt, atom_types_swnt, box_lx, box_ly, box_lz)

    # If the user chooses "None", only SWNT structure without hydrogens will be returned:
    if H_termination_SWNT == 'None':
        return univ_swnt
   ##############################################Create Hydrogens at only one end of tube##############################################

   #To hydrogenate one end of tube, we devide the number of atoms in nanotube into two:
    scale_factor_H = bond_length_hydrogen/bond_length
    num_hydrogen = 0
    H_coordinaes = []
    length_SWNT = np.linalg.norm(T) * 4

    for a in range(num_atoms_swnt):

        if ((np.linalg.norm(coord_array_swnt[a])) < (length_SWNT/2)):

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
                H_coord_1 = xyz[a] + (- xyz[right_connection] + xyz[f_connec_to_end_atom_index]) * scale_factor_H
                H_coord_2 = xyz[a] + (- xyz[left_connection] + xyz[f_connec_to_end_atom_index])  * scale_factor_H
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
                    H_coord = xyz[end_atom_index] + (first_vector + second_vector)  * scale_factor_H
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
                    H_coord = xyz[end_atom_index] + (first_vector + second_vector) * scale_factor_H
                    H_coordinaes.extend([H_coord])
                    num_hydrogen = num_hydrogen + 1

# Here we define the bond information between the atoms of SWNT and hydrogen, if the number of hydrogen is not zero:
    pos_one_end_H = np.array(H_coordinaes)
    assert pos_one_end_H.shape == (num_hydrogen, 3)
    num_hydrogen = 0
    n_residues = 1
    one_end_bonds_H = []
    for x in range(num_atoms_swnt):
        if ((np.linalg.norm(coord_array_swnt[x])) < (length_SWNT/2)):
            indices_of_a = np.where(all_bonds_swnt == x)
            if len(indices_of_a[0]) == 1:
                one_end_bonds_H.extend([(x, num_atoms_swnt + num_hydrogen)])
                one_end_bonds_H.extend([(x, num_atoms_swnt + num_hydrogen + 1)])
                num_hydrogen = num_hydrogen + 2
            if len(indices_of_a[0]) == 2:
                one_end_bonds_H.extend([(x, num_atoms_swnt + num_hydrogen)])
                num_hydrogen = num_hydrogen + 1
    one_end_atom_types_H = list(['H']*num_hydrogen)
    merged_swnt_one_end_H = merged_universe_with_H(coord_array_swnt, all_bonds_swnt, atom_types_swnt, pos_one_end_H, one_end_bonds_H, one_end_atom_types_H, box_lx, box_ly, box_lz)

    # If the user chooses "One end", only SWNT structure with one end hydrogenated will be returned:
    if H_termination_SWNT == 'One end':
        return merged_swnt_one_end_H

#########################################################
    num_hydrogen = 0
    H_coordinaes = []
    for a in range(num_atoms_swnt):
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
            H_coord_1 = xyz[a] + (- xyz[right_connection] + xyz[f_connec_to_end_atom_index]) * scale_factor_H
            H_coord_2 = xyz[a] + (- xyz[left_connection] + xyz[f_connec_to_end_atom_index]) * scale_factor_H
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
                H_coord = xyz[end_atom_index] + (first_vector + second_vector) * scale_factor_H
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
                H_coord = xyz[end_atom_index] + (first_vector + second_vector) * scale_factor_H
                H_coordinaes.extend([H_coord])
                num_hydrogen = num_hydrogen + 1

    num_hydrogen = 0
    two_end_bonds_H = []
    n_residues = 1
    pos_two_end_H = np.array(H_coordinaes)
    # assert pos_two_end_H.shape == (num_hydrogen, 3)#######################
    for x in range(num_atoms_swnt):
        indices_of_a = np.where(all_bonds_swnt == x)
        if len(indices_of_a[0]) == 1:
            two_end_bonds_H.extend([(x, num_atoms_swnt + num_hydrogen)])
            two_end_bonds_H.extend([(x, num_atoms_swnt + num_hydrogen + 1)])
            num_hydrogen = num_hydrogen + 2
        if len(indices_of_a[0]) == 2:
            two_end_bonds_H.extend([(x, num_atoms_swnt + num_hydrogen)])
            num_hydrogen = num_hydrogen + 1
    # If the user chooses "Both ends", SWNT structure with both ends hydrogenated will be returned:
    two_end_atom_types_H = list(['H']*num_hydrogen)
    merged_swnt_two_end_H = merged_universe_with_H(coord_array_swnt, all_bonds_swnt, atom_types_swnt, pos_two_end_H, two_end_bonds_H, two_end_atom_types_H, box_lx, box_ly, box_lz)

    if H_termination_SWNT == 'Both ends':
        return merged_swnt_two_end_H