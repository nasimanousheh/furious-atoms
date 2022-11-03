from furiousatoms.io import create_universe, merged_two_universes
import numpy as np
from numpy.linalg import norm
from math import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule
import io
from furiousatoms import io
import numpy as np
from fury import window
from PySide2 import QtWidgets
from fury import window, utils
from furiousatoms.structure import bbox

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
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        bond_length_NanoRope = 1.421 # default value of C-C bond length
        self.NanoRope.lineEdit_bond_length_NanoRope.insert(str(bond_length_NanoRope))
        self.NanoRope.lineEdit_bond_length_NanoRope.textChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.spinBox_chirality_N_NanoRope.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.spinBox_chirality_M_NanoRope.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.spinBox_repeat_units_NanoRope.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.pushButton_build_NanoRope.clicked.connect(self.NanoRope_builder_callback)
        self.NanoRope.pushButton_build_NanoRope.clicked.connect(lambda:self.close())
        # self.NanoRope.spinBox_num_walls_NanoRope.valueChanged.connect(self.NanoRope_diameter_changed)
        self.NanoRope.SpinBox_lx.valueChanged.connect(self.initial_box_dim)
        self.NanoRope.SpinBox_lz.valueChanged.connect(self.initial_box_dim)

        self.NanoRope.radioButton_desired_bond_length.toggled.connect(self.get_atom_type)
        self.NanoRope.radioButton_bond_length.toggled.connect(self.get_atom_type)
        self.NanoRope.comboBox_type1_NanoRope.activated.connect(self.get_atom_type)
        self.NanoRope.comboBox_type2_NanoRope.activated.connect(self.get_atom_type)

    def get_atom_type(self):
        global bond_length_NanoRope, bendFactor
        if self.NanoRope.radioButton_bond_length.isChecked() == True:
            self.NanoRope.SpinBox_desired_bond_length.setEnabled(False)
            NanoRope_type_1 = self.NanoRope.comboBox_type1_NanoRope.currentText()
            NanoRope_type_2 = self.NanoRope.comboBox_type2_NanoRope.currentText()
            if NanoRope_type_1=="C" and NanoRope_type_2=="C":
                bond_length_NanoRope = 1.421 # default value of C-C bond length
            if NanoRope_type_1=="N" and NanoRope_type_2=="B":
                bond_length_NanoRope = 1.47 # default value of N-B bond length
            if NanoRope_type_1=="N" and NanoRope_type_2=="Ga":
                bond_length_NanoRope = 1.95 # default value of N-Ga bond length
            if NanoRope_type_1=="N" and NanoRope_type_2=="Al":
                bond_length_NanoRope = 1.83 # default value of N-Al bond length
            if NanoRope_type_1=="P" and NanoRope_type_2=="Al":
                bond_length_NanoRope = 2.3 # default value of P-Al bond length
            if NanoRope_type_1=="P" and NanoRope_type_2=="Ga":
                bond_length_NanoRope = 2.28 # default value of P-Ga bond length
            if NanoRope_type_1=="P" and NanoRope_type_2=="C":
                bond_length_NanoRope = 1.87 # default value of P-C bond length
            if NanoRope_type_1=="N" and NanoRope_type_2=="C":
                bond_length_NanoRope = 1.47 # default value of N-C bond length
            if NanoRope_type_1=="C" and NanoRope_type_2=="B":
                bond_length_NanoRope = 1.56 # default value of C-B bond length
            if NanoRope_type_1=="C" and NanoRope_type_2=="Al":
                bond_length_NanoRope = 2.0 # default value of C-Al bond length
            if NanoRope_type_1=="C" and NanoRope_type_2=="Ga":
                bond_length_NanoRope = 2.46 # default value of P-B bond length
            if NanoRope_type_1=="P" and NanoRope_type_2=="B":
                bond_length_NanoRope = 1.74 # default value of P-B bond length
            self.NanoRope.lineEdit_bond_length_NanoRope.setText(str(bond_length_NanoRope))

        elif self.NanoRope.radioButton_desired_bond_length.isChecked() == True:
            self.NanoRope.lineEdit_bond_length_NanoRope.setText(str(' '))
            self.NanoRope.SpinBox_desired_bond_length.setEnabled(True)
            bond_length_NanoRope = float(self.NanoRope.SpinBox_desired_bond_length.text())
        else:
            bond_length_NanoRope = 1.421


    def NanoRope_diameter_changed(self):
        global bond_length_NanoRope, bendFactor
        try:
            bond_length_NanoRope
        except NameError:
            bond_length_NanoRope = 1.421
        # bond_length_NanoRope = self.NanoRope.lineEdit_bond_length_NanoRope.text()
        # bond_length_NanoRope = float(bond_length_NanoRope)
        value_n_NanoRope = int(self.NanoRope.spinBox_chirality_N_NanoRope.text())
        value_m_NanoRope = int(self.NanoRope.spinBox_chirality_M_NanoRope.text())
        repeat_units_NanoRope = int(self.NanoRope.spinBox_repeat_units_NanoRope.text())
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

        # self.NanoRope.lineEdit_wall_separation_NanoRope.setText(str(wall_separation))

    def initial_box_dim(self):
        global box_lx, box_ly, box_lz, bendFactor
        box_lx = float(self.NanoRope.SpinBox_lx.text())
        box_ly = float(self.NanoRope.SpinBox_lx.text())
        box_lz = float(self.NanoRope.SpinBox_lz.text())
        self.NanoRope.lineEdit_ly.setText(str(box_ly))
        return

    def NanoRope_builder_callback(self):
        global bond_length_NanoRope, box_lx, box_ly, box_lz, bendFactor
        try:
            bond_length_NanoRope
        except NameError:
            bond_length_NanoRope = 1.421

        try:
            box_lx or box_ly or box_lz
        except NameError:
            box_lx = box_ly = box_lz = 0.0
        try:
            bendFactor
        except NameError:
            bendFactor = 1.0
        number_of_walls = 2
        value_n_NanoRope = int(self.NanoRope.spinBox_chirality_N_NanoRope.text())
        value_m_NanoRope = int(self.NanoRope.spinBox_chirality_M_NanoRope.text())
        repeat_units_NanoRope = int(self.NanoRope.spinBox_repeat_units_NanoRope.text())
        NanoRope_type_1 = self.NanoRope.comboBox_type1_NanoRope.currentText()
        NanoRope_type_2 = self.NanoRope.comboBox_type2_NanoRope.currentText()
        bendFactor = float(self.NanoRope.doubleSpinBox_bend_factor.text())

        universe = NanoRope_builder(value_n_NanoRope, value_m_NanoRope, repeat_units_NanoRope, a=bond_length_NanoRope, species=(NanoRope_type_1, NanoRope_type_2), centered=True, bend=bendFactor)
        universe_all = universe.copy()
        for i in range(1, number_of_walls):
            xyz = []
            type_atoms = []
            next_universe = NanoRope_builder(value_n_NanoRope, value_m_NanoRope, repeat_units_NanoRope, a=bond_length_NanoRope, species=(NanoRope_type_1, NanoRope_type_2), centered=True, bend=bendFactor)
            xyz.extend(next_universe.universe.atoms.positions)
            type_atoms.extend(next_universe.atoms.types)

            universe_all = merged_two_universes(universe_all.atoms.positions, universe_all.bonds.indices, universe_all.atoms.types, next_universe.atoms.positions, next_universe.bonds.indices, next_universe.atoms.types, box_lx, box_ly, box_lz)

        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(universe_all)
        window.show()

"""
  The numbers (n,m) show that your tube is obtained from taking one atom of the sheet and rolling it onto
  that atom that is at located na1+ma2 away from your original atom. N is the Number of hexagons in a unit cell. a1 and a2 are lattice vectors.
  (n,m=n) gives an “armchair” tube,e.g. (5,5). (n,m=0) gives an “zig-zag” tube, e.g. (6,0). Other tubes are “chiral”, e.g. (6,2)
"""

def NanoRope_builder(n, m, N, a, species=('B', 'C'), centered=False, bend=1):
    global box_lx, box_ly, box_lz, bendFactor
    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2))
    Ch = (n*a1+m*a2)
    T = t1*a1+t2*a2

    Ch_proj, T_proj = [(v/norm(v)**2) for v in [Ch, T]]
    basis = [np.array((0, 0)), ((a1+a2)/3)]
    pts = []
    xyz = []
    atom_types_swnt = []
    for i1, i2 in product(range(0, n+t1+1), range(t2, m+1)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            if all(-thre < pt.dot(v) < 1-thre for v in [Ch_proj, T_proj]):
                for k in range(N):
                    pts.append((sp, pt+k*T))

    diameter = (norm(Ch)/(np.pi*bendFactor))
    print(diameter)
    def gr2tube(v):
        # phi = 2*np.pi*v.dot(Ch_proj)
        phi = np.pi + bendFactor * 2*np.pi*v.dot(Ch_proj)
        return np.array((diameter/2*np.cos(phi),
                         diameter/2*np.sin(phi),
                         v.dot(T_proj)*norm(T)))
    xyz = [gr2tube(v) for _, v in pts]
    atom_types_swnt = [v for v, _ in pts]
    mol_1 = Molecule([Atom(sp, r) for (sp, _), r in zip(pts, xyz)])
    fragments = mol_1.to_json(scale =1)
    n_atoms_swnt = len(xyz)
    coord_array_swnt = np.array(xyz)
    assert coord_array_swnt.shape == (n_atoms_swnt, 3)
    all_bonds_swnt = np.array(fragments['bonds'])

    try:
        box_lx or box_ly or box_lz
    except NameError:
        box_lx = box_ly = box_lz = 0.0

    univ_swnt = create_universe(coord_array_swnt, all_bonds_swnt, atom_types_swnt, box_lx, box_ly, box_lz)
    return univ_swnt
