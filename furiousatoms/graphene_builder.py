import numpy as np
from numpy.linalg import norm
from math import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule
from furiousatoms import io
import numpy as np
from fury import window
from PySide2 import QtWidgets
from furiousatoms.io import merged_universe_with_H, create_universe
import MDAnalysis as mda


"""
    Ui_Graphene class creates a widget for building graphenes
"""
class Ui_graphene(QtWidgets.QMainWindow):

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
        self.graphene.pushButton_build_graphene_2.clicked.connect(self.graphene_shape_regtang)
        self.graphene.pushButton_build_graphene_1.clicked.connect(lambda:self.close())
        self.graphene.pushButton_build_graphene_1.clicked.connect(self.graphene_shape_diamond)
        self.graphene.pushButton_build_graphene_2.clicked.connect(lambda:self.close())
        self.graphene.SpinBox_lx.valueChanged.connect(self.initial_box_dim)
        self.graphene.SpinBox_ly.valueChanged.connect(self.initial_box_dim)
        self.graphene.SpinBox_lz.valueChanged.connect(self.initial_box_dim)
        self.graphene.radioButton_desired_bond_length.toggled.connect(self.get_atom_type)
        self.graphene.SpinBox_desired_bond_length.valueChanged.connect(self.get_atom_type)
        self.graphene.radioButton_bond_length.toggled.connect(self.get_atom_type)
        self.graphene.comboBox_type1_graphene.activated.connect(self.get_atom_type)
        self.graphene.comboBox_type2_graphene.activated.connect(self.get_atom_type)
        self.graphene.doubleSpinBox_unitcell_along_x.valueChanged.connect(self.initial_box_dim)
        self.graphene.doubleSpinBox_unitcell_along_y.valueChanged.connect(self.initial_box_dim)
        self.graphene.comboBox_type1_graphene.activated.connect(self.initial_box_dim)
        self.graphene.comboBox_type2_graphene.activated.connect(self.initial_box_dim)
        self.graphene.doubleSpinBox_unitcell_along_x.valueChanged.connect(self.get_atom_type)
        self.graphene.doubleSpinBox_unitcell_along_y.valueChanged.connect(self.get_atom_type)

    def get_atom_type(self):
        min = 1.3
        max = 2.6
        if self.graphene.radioButton_bond_length.isChecked() == True or self.graphene.radioButton_desired_bond_length.isChecked() == True:
            self.graphene.SpinBox_desired_bond_length.setEnabled(False)
            graphene_type_1 = self.graphene.comboBox_type1_graphene.currentText()
            graphene_type_2 = self.graphene.comboBox_type2_graphene.currentText()
            if graphene_type_1=="C" and graphene_type_2=="C":
                bond_length_graphene = 1.421 # default value of C-C bond length
                min = 1.3
                max = 2.0
            if graphene_type_1=="N" and graphene_type_2=="B":
                bond_length_graphene = 1.47 # default value of N-B bond length
                min = 1.3
                max = 2.0
            if graphene_type_1=="N" and graphene_type_2=="Ga":
                bond_length_graphene = 1.95 # default value of N-Ga bond length
                min = 1.9
                max = 2.5
            if graphene_type_1=="N" and graphene_type_2=="Al":
                bond_length_graphene = 1.83 # default value of N-Al bond length
                min = 1.8
                max = 2.5
            if graphene_type_1=="P" and graphene_type_2=="Al":
                bond_length_graphene = 2.3 # default value of P-Al bond length
                min = 1.8
                max = 2.6
            if graphene_type_1=="P" and graphene_type_2=="Ga":
                bond_length_graphene = 2.28 # default value of P-Ga bond length
                min = 1.9
                max = 2.6
            if graphene_type_1=="P" and graphene_type_2=="C":
                bond_length_graphene = 1.87 # default value of P-C bond length
                min = 1.6
                max = 2.3
            if graphene_type_1=="N" and graphene_type_2=="C":
                bond_length_graphene = 1.47 # default value of N-C bond length
                min = 1.3
                max = 1.9
            if graphene_type_1=="C" and graphene_type_2=="B":
                bond_length_graphene = 1.56 # default value of C-B bond length
                min = 1.3
                max = 2.0
            if graphene_type_1=="C" and graphene_type_2=="Al":
                bond_length_graphene = 2.0 # default value of C-Al bond length
                min = 1.8
                max = 2.5
            if graphene_type_1=="C" and graphene_type_2=="Ga":
                bond_length_graphene = 2.46 # default value of P-B bond length
                min = 1.9
                max = 2.6
            if graphene_type_1=="P" and graphene_type_2=="B":
                bond_length_graphene = 1.74 # default value of P-B bond length
                min = 1.6
                max = 2.4

            self.graphene.SpinBox_desired_bond_length.setRange(min, max)
            self.graphene.lineEdit_bond_length_graphene.setText(str(bond_length_graphene))

        if self.graphene.radioButton_desired_bond_length.isChecked() == True:
            self.graphene.SpinBox_desired_bond_length.setRange(min, max)
            self.graphene.lineEdit_bond_length_graphene.setText(str(' '))
            self.graphene.SpinBox_desired_bond_length.setEnabled(True)
            bond_length_graphene = float(self.graphene.SpinBox_desired_bond_length.text())

        return bond_length_graphene


    def graphene_shape_diamond(self):
        graphene_shape=True
        self.graphene_builder_callback(graphene_shape)
    def graphene_shape_regtang(self):
        graphene_shape=False
        self.graphene_builder_callback(graphene_shape)

    def graphene_builder_callback(self, graphene_shape):
        bond_length_graphene = self.get_atom_type()
        H_termination_graphene = self.graphene.comboBox_H_termination_graphene.currentText()
        value_n_graphene = int(self.graphene.spinBox_chirality_N_graphene.text())
        value_m_graphene = int(self.graphene.spinBox_chirality_M_graphene.text())
        repeat_units_graphene = int(self.graphene.spinBox_repeat_units_graphene.text())
        graphene_type_1 = self.graphene.comboBox_type1_graphene.currentText()
        graphene_type_2 = self.graphene.comboBox_type2_graphene.currentText()
        try:
            num_sheets = int(self.graphene.SpinBox_num_sheets.text())
        except NameError:
            num_sheets = 1

        try:
            sheet_separation = float(self.graphene.SpinBox_sheet_sep.text())
        except NameError:
            sheet_separation = 3.347

        if self.graphene.radioButton_armchair.isChecked() == True:
            edge_shape = 'armchair'
        else:
            edge_shape = 'zigzag'
        edge_length_x = float(self.graphene.doubleSpinBox_unitcell_along_x.text())
        edge_length_y = float(self.graphene.doubleSpinBox_unitcell_along_y.text())
        structure_info = graphene_builder(H_termination_graphene, value_n_graphene, value_m_graphene, repeat_units_graphene, length=None, bond_length_graphene=bond_length_graphene,
                                            species=(graphene_type_1, graphene_type_2), diamond_sheet=graphene_shape, edge_length_x=edge_length_x, edge_length_y=edge_length_y, edge_shape=edge_shape)
        if num_sheets > 1:
            structure_info = extend_the_sheets(structure_info, num_sheets, sheet_separation)

        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(structure_info)
        window.show()

    def initial_box_dim(self):
        global box_lx, box_ly, box_lz
        box_lx = float(self.graphene.SpinBox_lx.text())
        box_ly = float(self.graphene.SpinBox_ly.text())
        box_lz = float(self.graphene.SpinBox_lz.text())
        return

def extend_the_sheets(structure_info, num_sheets, sheet_separation):
    copied = []
    structure_info.dimensions = [sheet_separation, sheet_separation, sheet_separation, 90, 90, 90]
    box = structure_info.dimensions[:3]
    for x in range(1):
        for y in range(1):
            for z in range(num_sheets):
                u_ = structure_info.copy()
                move_by = box*(x, y, z)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)

        extended_universe = mda.Merge(*copied)
        extended_universe.dimensions = [box_lx, box_ly, box_lz, 90, 90, 90]
        cog = extended_universe.atoms.center_of_geometry()
        extended_universe.atoms.positions -= cog
        return extended_universe



def graphen_armchair(primitive_unitcell, edge_length_x, edge_length_y, bond_length_graphene):
    primitive_unitcell.bonds.to_indices()
    pos = primitive_unitcell.atoms.positions
    pos = pos.astype('float64')
    unit_cell_ly = np.linalg.norm(pos[2]-pos[3]) + bond_length_graphene
    unit_cell_lx = np.linalg.norm(pos[2]+pos[3]) * 2
    primitive_unitcell.dimensions = [unit_cell_lx, unit_cell_ly, unit_cell_lx, 90, 90, 90]
    box = primitive_unitcell.dimensions[:3]
    copied = []
    b = 0
    i = 0
    num_unitcell_in_lx = int(np.ceil(edge_length_x/unit_cell_lx))
    num_unitcell_in_ly = int(np.ceil(edge_length_y/unit_cell_ly))
    for x in range(num_unitcell_in_lx):
        i = 0
        for y in range(num_unitcell_in_ly):
            u_ = primitive_unitcell.copy()
            move_by = box*(x, y, 1)
            u_.atoms.translate(move_by)
            copied.append(u_.atoms)

    new_universe = mda.Merge(*copied)
    b = 0
    c = 0
    num_atoms_in_y_direction = 4 * num_unitcell_in_ly
    num_bonds_connect = (num_unitcell_in_lx-1)
    for b in range(num_unitcell_in_ly):
        b = c *4
        for i in range(num_bonds_connect):
            added_bonds_1 = np.array([[num_atoms_in_y_direction*i+3+b, num_atoms_in_y_direction*i+num_atoms_in_y_direction+b]])
            added_bonds_2 = np.array([[num_atoms_in_y_direction*i+2+b, num_atoms_in_y_direction*i+num_atoms_in_y_direction+b+1]])
            new_universe.add_bonds(added_bonds_1)
            new_universe.add_bonds(added_bonds_2)
        for j in range(num_unitcell_in_lx):
            if c < num_unitcell_in_ly-1:
                added_bonds_3 = np.array([[(num_atoms_in_y_direction*j)+3+b, (num_atoms_in_y_direction*j)+b+6]])
                new_universe.add_bonds(added_bonds_3)
        c = c + 1
    cog = new_universe.atoms.center_of_geometry()
    new_universe.atoms.positions -= cog
    return new_universe


def graphen_zigzag(primitive_unitcell, edge_length_x, edge_length_y, bond_length_graphene):
    primitive_unitcell.bonds.to_indices()
    pos = primitive_unitcell.atoms.positions
    pos = pos.astype('float64')
    unit_cell_ly = np.linalg.norm(pos[0]+pos[3]) *2
    unit_cell_lx = np.linalg.norm(pos[0]-pos[3]) + bond_length_graphene
    primitive_unitcell.dimensions = [unit_cell_lx, unit_cell_ly, unit_cell_lx, 90, 90, 90]
    box = primitive_unitcell.dimensions[:3]
    copied = []
    b = 0
    i = 0
    num_unitcell_in_lx = int(np.ceil(edge_length_x/unit_cell_lx))
    num_unitcell_in_ly = int(np.ceil(edge_length_y/unit_cell_ly))
    for x in range(num_unitcell_in_lx):
        i = 0
        for y in range(num_unitcell_in_ly):
            u_ = primitive_unitcell.copy()
            move_by = box*(x, y, 1)
            u_.atoms.translate(move_by)
            copied.append(u_.atoms)

    new_universe = mda.Merge(*copied)
    b = 0
    c = 0
    num_atoms_in_y_direction = 4 * num_unitcell_in_ly
    num_bonds_connect = (num_unitcell_in_lx-1)
    for b in range(num_unitcell_in_ly):
        b = c *4
        for i in range(num_bonds_connect):
            added_bonds_1 = np.array([[num_atoms_in_y_direction*i+3+b, num_atoms_in_y_direction*i+num_atoms_in_y_direction+b]])
            new_universe.add_bonds(added_bonds_1)
        for j in range(num_unitcell_in_lx):
            if c < num_unitcell_in_ly-1:
                added_bonds_2 = np.array([[(num_atoms_in_y_direction*j)+b, (num_atoms_in_y_direction*j)+b+5]])
                added_bonds_3 = np.array([[(num_atoms_in_y_direction*j)+3+b, (num_atoms_in_y_direction*j)+b+6]])
                new_universe.add_bonds(added_bonds_2)
                new_universe.add_bonds(added_bonds_3)
        c = c + 1
    cog = new_universe.atoms.center_of_geometry()
    new_universe.atoms.positions -= cog
    return new_universe

def graphene_builder(H_termination_graphene, n, m, N, length, bond_length_graphene, species=('C', 'C'), diamond_sheet=True, edge_length_x=1, edge_length_y=1, edge_shape = 'armchair'):
    global box_lx, box_ly, box_lz
    try:
        box_lx or box_ly or box_lz
    except NameError:
        box_lx = box_ly = box_lz = 0.0
    bond_length_hydrogen = 1.0
    if diamond_sheet is True:
        d = gcd(n, m)
        dR = 3*d if (n-m) % (3*d) == 0 else d
        t1 = (2*m+n)//dR
        t2 = -(2*n+m)//dR
        a1 = np.array((np.sqrt(3)*bond_length_graphene,0,0))
        a2 = np.array((np.sqrt(3)/2*bond_length_graphene, -3*bond_length_graphene/2,0))
        T = t1*a1+t2*a2
        basis = [np.array((0,0,0)), (a1+a2)/3]
        pts = []
        for i1, i2 in product(range(0, n+t1), range(t2, m)):
            shift = i1*a1+i2*a2
            for sp, b in zip(species, basis):
                pt = b+shift
                pts.append((sp, pt))

        xyz = [v for _, v in pts]
        atom_types_graphene = [v for v, _ in pts]
        m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
        fragments = m.to_json(scale =1.3)
        coord_array_graphene = np.array(xyz)
        num_atoms_graphene = len(coord_array_graphene)
        assert coord_array_graphene.shape == (num_atoms_graphene, 3)
        all_bonds_graphene = np.array(fragments['bonds'])
        univ_graphene = create_universe(coord_array_graphene, all_bonds_graphene, atom_types_graphene, box_lx, box_ly, box_lz)
    else:
        m = n = 1
        d = gcd(n, m)
        dR = 3*d if (n-m) % (3*d) == 0 else d
        t1 = (2*m+n)//dR
        t2 = -(2*n+m)//dR
        if edge_shape == 'armchair':
            a1 = np.array((np.sqrt(3)*bond_length_graphene,0,0))
            a2 = np.array((np.sqrt(3)/2*bond_length_graphene, -3*bond_length_graphene/2,0))
            T = t1*a1+t2*a2
            basis = [np.array((0,0,0)), (a1+a2)/3]
            pts = []
            for i1, i2 in product(range(0, n+t1), range(t2, m)):
                shift = i1*a1+i2*a2
                for sp, b in zip(species, basis):
                    pt = b+shift
                    pts.append((sp, pt))
            pts = pts[1:5]
            xyz = [v for _, v in pts]
            primitive_atom_types_graphene = [v for v, _ in pts]
            m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
            fragments = m.to_json(scale =1)
            primitive_coord_graphene = np.array(xyz)
            primitive_num_atoms_graphene = len(primitive_coord_graphene)
            assert primitive_coord_graphene.shape == (primitive_num_atoms_graphene, 3)
            primitive_bonds_graphene = np.array([[1,0],[2,1],[3,0]])
            primitive_graphene = create_universe(primitive_coord_graphene, primitive_bonds_graphene, primitive_atom_types_graphene, box_lx, box_ly, box_lz)
            univ_graphene = graphen_armchair(primitive_graphene, edge_length_x, edge_length_y, bond_length_graphene)

        if edge_shape == 'zigzag':
            a1 = np.array((3/2*bond_length_graphene, 1*np.sqrt(3)/2 * bond_length_graphene, 0))
            a2 = np.array((3/2*bond_length_graphene, -1*np.sqrt(3)/2 * bond_length_graphene, 0))
            T = t1*a1+t2*a2
            basis = [np.array((0,0,0)), (a1+a2)/3]
            pts = []
            for i1, i2 in product(range(0, n+t1), range(t2, m)):
                shift = i1*a1+i2*a2
                for sp, b in zip(species, basis):
                    pt = b+shift
                    pts.append((sp, pt))
            pts = pts[1:7]
            del pts[3]
            del pts[3]
            xyz = [v for _, v in pts]
            primitive_atom_types_graphene = [v for v, _ in pts]
            m = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
            fragments = m.to_json(scale =1)
            primitive_coord_graphene = np.array(xyz)
            primitive_num_atoms_graphene = len(primitive_coord_graphene)
            assert primitive_coord_graphene.shape == (primitive_num_atoms_graphene, 3)
            primitive_bonds_graphene = np.array([[1,0],[2,1],[2,3]])
            primitive_graphene = create_universe(primitive_coord_graphene, primitive_bonds_graphene, primitive_atom_types_graphene, box_lx, box_ly, box_lz)
            univ_graphene = graphen_zigzag(primitive_graphene, edge_length_x, edge_length_y, bond_length_graphene)

        all_bonds_graphene = univ_graphene.bonds.to_indices()
        xyz = univ_graphene.atoms.positions
        coord_array_graphene = xyz.astype('float64')
        num_atoms_graphene = len(coord_array_graphene)
        atom_types_graphene = univ_graphene.atoms.types

    if H_termination_graphene == 'None':
        return univ_graphene
   ##############################################Create Hydrogens at the end of graphene##############################################
    scale_factor_H = bond_length_hydrogen/bond_length_graphene
    num_hydrogen = 0
    H_coordinaes = []
    if diamond_sheet is True:
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
                H_coord_1_new = xyz[a] + (H_coord_1 - xyz[a]) * scale_factor_H
                H_coord_2_new = xyz[a] + (H_coord_2 - xyz[a]) * scale_factor_H
                H_coordinaes.extend([H_coord_1_new])
                H_coordinaes.extend([H_coord_2_new])
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
                    H_coord = xyz[end_atom_index] + (first_vector + second_vector) * scale_factor_H
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
                    H_coord = xyz[end_atom_index] + (first_vector + second_vector) * scale_factor_H
                    H_coordinaes.extend([H_coord])
                    num_hydrogen = num_hydrogen + 1
    else:
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
                if end_atom_indices[0][1] == a:
                    end_atom_indices = (all_bonds_graphene[indices_of_a])
                    f_connec_to_end_atom_index = end_atom_indices[0][0]
                    core_connections = all_bonds_graphene[np.where(all_bonds_graphene == f_connec_to_end_atom_index)[0]]
                    Both_connected_atoms_with_a = np.setdiff1d(core_connections, [f_connec_to_end_atom_index])
                    Both_connected_atoms = np.setdiff1d(Both_connected_atoms_with_a, a)
                    right_connection = Both_connected_atoms[0]
                H_coord_1 = xyz[a] - xyz[right_connection] + xyz[f_connec_to_end_atom_index]
                H_coord_3 = H_coord_1 + xyz[f_connec_to_end_atom_index] - xyz[a]
                H_coord_2 = -H_coord_3 + (2*xyz[a])

                H_coord_1_new = xyz[a] + (H_coord_1 - xyz[a]) * scale_factor_H
                H_coord_2_new = xyz[a] + (H_coord_2 - xyz[a]) * scale_factor_H
                H_coordinaes.extend([H_coord_1_new])
                H_coordinaes.extend([H_coord_2_new])
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
                    H_coord = xyz[end_atom_index] + (first_vector + second_vector) * scale_factor_H
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
                    H_coord = xyz[end_atom_index] + (first_vector + second_vector) * scale_factor_H
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
    merged_graphene_hydrogen = merged_universe_with_H(coord_array_graphene, all_bonds_graphene, atom_types_graphene, coord_array_H_indice, bonds_hydrogen, atom_types_Hydrogen, box_lx, box_ly, box_lz)
    # If the user chooses "All", hydrogenated graphene structure will be returned:
    if H_termination_graphene == 'All':
        return merged_graphene_hydrogen
