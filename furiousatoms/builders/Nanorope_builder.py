import numpy as np
from math import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule
import io
from furiousatoms import io
import numpy as np
from fury import window
from PySide2 import QtWidgets
from fury import window
from PySide2.QtGui import QIcon
from furiousatoms.molecular import MolecularStructure


thre = 1e-10
vacuum = 4
class Ui_NanoRope(QtWidgets.QMainWindow):
    """ Ui_NanoRope class creates a widget for building multple-walls nanotube (NanoRope)
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_NanoRope, self).__init__(parent)
        self.NanoRope = io.load_ui_widget("NanoRope.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.NanoRope)
        self.setCentralWidget(self.NanoRope)
        self.setLayout(self.v_layout)
        self.resize(248, 313)
        self.scene = window.Scene()
        self.setWindowIcon(QIcon(io.get_resources_file("splash.png")))
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        bond_length_NanoRope = 1.421 # default value of C-C bond length
        self.NanoRope.lineEdit_bond_length_NanoRope.insert(str(bond_length_NanoRope))
        # self.NanoRope.lineEdit_bond_length_NanoRope.textChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.spinBox_chirality_N_NanoRope.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.spinBox_chirality_M_NanoRope.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.spinBox_repeat_units_NanoRope.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.spinBox_nanotube_separation.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.pushButton_build_NanoRope.clicked.connect(self.NanoRope_builder_callback)
        self.NanoRope.spinBox_num_nanotubes.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.pushButton_build_NanoRope.clicked.connect(lambda:self.close())
        self.NanoRope.comboBox_H_termination_SWNT.currentTextChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.SpinBox_lx.valueChanged.connect(self.initial_box_dim)
        self.NanoRope.SpinBox_ly.valueChanged.connect(self.initial_box_dim)
        self.NanoRope.SpinBox_lz.valueChanged.connect(self.initial_box_dim)
        self.NanoRope.radioButton_desired_bond_length.toggled.connect(self.get_atom_type)
        self.NanoRope.radioButton_bond_length.toggled.connect(self.get_atom_type)
        self.NanoRope.comboBox_type1_NanoRope.activated.connect(self.get_atom_type)
        self.NanoRope.comboBox_type2_NanoRope.activated.connect(self.get_atom_type)
        self.NanoRope.SpinBox_desired_bond_length.valueChanged.connect(self.get_atom_type)

    def get_atom_type(self):
        global bendFactor, diameter_SWNT
        if self.NanoRope.radioButton_bond_length.isChecked() == True or self.NanoRope.radioButton_desired_bond_length.isChecked() == True:
            self.NanoRope.SpinBox_desired_bond_length.setEnabled(False)
            NanoRope_type_1 = self.NanoRope.comboBox_type1_NanoRope.currentText()
            NanoRope_type_2 = self.NanoRope.comboBox_type2_NanoRope.currentText()
            if NanoRope_type_1=="C" and NanoRope_type_2=="C":
                bond_length_NanoRope = 1.421 # default value of C-C bond length
                min = 1.3
                max = 2.0
            if NanoRope_type_1=="N" and NanoRope_type_2=="B":
                bond_length_NanoRope = 1.47 # default value of N-B bond length
                min = 1.3
                max = 2.0
            if NanoRope_type_1=="N" and NanoRope_type_2=="Ga":
                bond_length_NanoRope = 1.95 # default value of N-Ga bond length
                min = 1.9
                max = 2.5
            if NanoRope_type_1=="N" and NanoRope_type_2=="Al":
                bond_length_NanoRope = 1.83 # default value of N-Al bond length
                min = 1.8
                max = 2.5
            if NanoRope_type_1=="P" and NanoRope_type_2=="Al":
                bond_length_NanoRope = 2.3 # default value of P-Al bond length
                min = 1.8
                max = 2.6
            if NanoRope_type_1=="P" and NanoRope_type_2=="Ga":
                bond_length_NanoRope = 2.28 # default value of P-Ga bond length
                min = 1.9
                max = 2.6
            if NanoRope_type_1=="P" and NanoRope_type_2=="C":
                bond_length_NanoRope = 1.87 # default value of P-C bond length
                min = 1.6
                max = 2.3
            if NanoRope_type_1=="N" and NanoRope_type_2=="C":
                bond_length_NanoRope = 1.47 # default value of N-C bond length
                min = 1.3
                max = 1.9
            if NanoRope_type_1=="C" and NanoRope_type_2=="B":
                bond_length_NanoRope = 1.56 # default value of C-B bond length
                min = 1.3
                max = 2.0
            if NanoRope_type_1=="C" and NanoRope_type_2=="Al":
                bond_length_NanoRope = 2.0 # default value of C-Al bond length
                min = 1.8
                max = 2.5
            if NanoRope_type_1=="C" and NanoRope_type_2=="Ga":
                bond_length_NanoRope = 2.46 # default value of P-B bond length
                min = 1.9
                max = 2.6
            if NanoRope_type_1=="P" and NanoRope_type_2=="B":
                bond_length_NanoRope = 1.74 # default value of P-B bond length
                min = 1.6
                max = 2.4

            self.NanoRope.SpinBox_desired_bond_length.setRange(min, max)
            self.NanoRope.lineEdit_bond_length_NanoRope.setText(str(bond_length_NanoRope))

        if self.NanoRope.radioButton_desired_bond_length.isChecked() == True:
            self.NanoRope.SpinBox_desired_bond_length.setRange(min, max)
            self.NanoRope.lineEdit_bond_length_NanoRope.setText(str(' '))
            self.NanoRope.SpinBox_desired_bond_length.setEnabled(True)
            bond_length_NanoRope = float(self.NanoRope.SpinBox_desired_bond_length.text())

        return bond_length_NanoRope


    def NanoRope_diameter_changed(self):
        global bendFactor, diameter_SWNT
        bond_length_NanoRope = self.get_atom_type()
        value_n_NanoRope = int(self.NanoRope.spinBox_chirality_N_NanoRope.text())
        value_m_NanoRope = int(self.NanoRope.spinBox_chirality_M_NanoRope.text())
        repeat_units_NanoRope = int(self.NanoRope.spinBox_repeat_units_NanoRope.text())
        nanotube_separation = float(self.NanoRope.spinBox_nanotube_separation.text())
        a1 = np.array((np.sqrt(3)*bond_length_NanoRope, 0))
        a2 = np.array((np.sqrt(3)/2*bond_length_NanoRope, -3*bond_length_NanoRope/2))
        Ch = value_n_NanoRope*a1+value_m_NanoRope*a2
        d = gcd(value_n_NanoRope, value_m_NanoRope)
        dR = 3*d if (value_n_NanoRope-value_m_NanoRope) % (3*d) == 0 else d
        t1 = (2*value_m_NanoRope+value_n_NanoRope)//dR
        t2 = -(2*value_n_NanoRope+value_m_NanoRope)//dR
        T = t1*a1+t2*a2
        diameter_NanoRope = float(np.linalg.norm(Ch)/np.pi)
        diameter_NanoRope = "{:.2f}".format(diameter_NanoRope)
        length_NanoRope = np.linalg.norm(T) * repeat_units_NanoRope
        length_NanoRope = "{:.2f}".format(float(length_NanoRope))
        self.NanoRope.lineEdit_diameter_NanoRope.setText(str(diameter_NanoRope))
        self.NanoRope.lineEdit_length_NanoRope.setText(str(length_NanoRope))
        number_of_walls = 2
        if number_of_walls > 1:
            wall_separation = 3.43
        else:
            wall_separation = 0.0


    def initial_box_dim(self):
        global box_lx, box_ly, box_lz, bendFactor
        box_lx = float(self.NanoRope.SpinBox_lx.text())
        box_ly = float(self.NanoRope.SpinBox_ly.text())
        box_lz = float(self.NanoRope.SpinBox_lz.text())
        return

    def NanoRope_builder_callback(self):
        global bond_length_NanoRope, box_lx, box_ly, box_lz, bendFactor
        bond_length_NanoRope = self.get_atom_type()

        try:
            box_lx or box_ly or box_lz
        except NameError:
            box_lx = box_ly = box_lz = 0.0
        try:
            bendFactor
        except NameError:
            bendFactor = 1.0
        H_termination_SWNT = self.NanoRope.comboBox_H_termination_SWNT.currentText()
        value_n_NanoRope = int(self.NanoRope.spinBox_chirality_N_NanoRope.text())
        value_m_NanoRope = int(self.NanoRope.spinBox_chirality_M_NanoRope.text())
        repeat_units_NanoRope = int(self.NanoRope.spinBox_repeat_units_NanoRope.text())
        nanotube_separation = float(self.NanoRope.spinBox_nanotube_separation.text())
        number_tubes = int(self.NanoRope.spinBox_num_nanotubes.text())
        NanoRope_type_1 = self.NanoRope.comboBox_type1_NanoRope.currentText()
        NanoRope_type_2 = self.NanoRope.comboBox_type2_NanoRope.currentText()
        # bendFactor = float(self.NanoRope.doubleSpinBox_bend_factor.text())
        original_structure = SWNT_builder(H_termination_SWNT, value_n_NanoRope, value_m_NanoRope, repeat_units_NanoRope, length=None, bond_length=bond_length_NanoRope, species=(NanoRope_type_1, NanoRope_type_2), centered=True)
        entire_structure = MolecularStructure.create_empty()
        distance_between_tubes = 0
        for i in range(number_tubes):
            if i == 0:
                pos = original_structure.pos
            else:
                if (i > 0 and i <=6):
                    distance_between_tubes = diameter_SWNT + nanotube_separation
                    shifted_tube_i = np.array([np.cos(np.radians(60*i)) * distance_between_tubes, np.sin(np.radians(60*i))* distance_between_tubes, 0.0])
                    pos = original_structure.pos + shifted_tube_i

                if (i > 6 and i <= 19):
                    distance_between_tubes = 2 * (diameter_SWNT + nanotube_separation)
                    shifted_tube_i = np.array([np.cos(np.radians(30*i)) * distance_between_tubes, np.sin(np.radians(30*i))* distance_between_tubes, 0.0])
                    pos = original_structure.pos + shifted_tube_i

                if (i > 19 and i <= 38):
                    distance_between_tubes = 3 * (diameter_SWNT + nanotube_separation)
                    shifted_tube_i = np.array([np.cos(np.radians(20*i)) * distance_between_tubes, np.sin(np.radians(20*i))* distance_between_tubes, 0.0])
                    pos = original_structure.pos + shifted_tube_i

            structure_1 = MolecularStructure.create_empty()
            structure_1.pos = pos
            structure_1.atom_types = np.copy(original_structure.atom_types)
            structure_1.bonds = np.copy(original_structure.bonds)

            entire_structure = entire_structure.merge(structure_1, offset_bonds=True)
        entire_structure.center()
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_structure(entire_structure)
        window.show()

"""
  The numbers (n,m) show that your tube is obtained from taking one atom of the sheet and rolling it onto
  that atom that is at located na1+ma2 away from your original atom. N is the Number of hexagons in a unit cell. a1 and a2 are lattice vectors.
  (n,m=n) gives an “armchair” tube,e.g. (5,5). (n,m=0) gives an “zig-zag” tube, e.g. (6,0). Other tubes are “chiral”, e.g. (6,2)
"""

def SWNT_builder(H_termination_SWNT, n, m, N, length, bond_length, species=('C', 'C'), centered=False):
    global box_lx, box_ly, box_lz, bendFactor, diameter_SWNT
    bond_length_hydrogen = 1.0
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
    fragments = m.to_json(scale=1.3)
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
    box_size = [box_lx, box_ly, box_lz]
    univ_swnt = MolecularStructure(box_size, coord_array_swnt, all_bonds_swnt, atom_types_swnt)

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
    merged_swnt_one_end_H = MolecularStructure(box_size, coord_array_swnt, all_bonds_swnt, atom_types_swnt) \
            .merge(MolecularStructure(box_size, pos_one_end_H, one_end_bonds_H, one_end_atom_types_H))

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
    merged_swnt_two_end_H = MolecularStructure(box_size, coord_array_swnt, all_bonds_swnt, atom_types_swnt) \
            .merge(MolecularStructure(box_size, pos_two_end_H, two_end_bonds_H, two_end_atom_types_H))

    if H_termination_SWNT == 'Both ends':
        return merged_swnt_two_end_H