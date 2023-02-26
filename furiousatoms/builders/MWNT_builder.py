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
import MDAnalysis as mda
from PySide2.QtGui import QIcon

thre = 1e-10
vacuum = 4
class Ui_MWNT(QtWidgets.QMainWindow):
    """ Ui_MWNT class creates a widget for building multple-walls nanotube (MWNT)
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_MWNT, self).__init__(parent)
        self.MWNT = io.load_ui_widget("MWNT.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.MWNT)
        self.setCentralWidget(self.MWNT)
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
        bond_length_MWNT = 1.421 # default value of C-C bond length
        self.MWNT.lineEdit_bond_length_MWNT.insert(str(bond_length_MWNT))
        self.MWNT.spinBox_chirality_N_MWNT.valueChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.SpinBox_desired_bond_length.valueChanged.connect(self.get_atom_type)
        self.MWNT.spinBox_chirality_M_MWNT.valueChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.spinBox_repeat_units_MWNT.valueChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.pushButton_build_MWNT.clicked.connect(self.MWNT_builder_callback)
        self.MWNT.pushButton_build_MWNT.clicked.connect(lambda:self.close())
        self.MWNT.spinBox_num_walls_MWNT.valueChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.SpinBox_lx.valueChanged.connect(self.initial_box_dim)
        self.MWNT.SpinBox_ly.valueChanged.connect(self.initial_box_dim)
        self.MWNT.SpinBox_lz.valueChanged.connect(self.initial_box_dim)

        self.MWNT.radioButton_desired_bond_length.toggled.connect(self.get_atom_type)
        self.MWNT.radioButton_bond_length.toggled.connect(self.get_atom_type)
        self.MWNT.comboBox_type1_MWNT.activated.connect(self.get_atom_type)
        self.MWNT.comboBox_type2_MWNT.activated.connect(self.get_atom_type)

    def get_atom_type(self):
        global bendFactor
        bendFactor = 1.0
        # try:
        #     bendFactor = float(self.MWNT.doubleSpinBox_bend_factor.text())
        # except NameError:
            # bendFactor = 1.0
        if self.MWNT.radioButton_bond_length.isChecked() == True or self.MWNT.radioButton_desired_bond_length.isChecked() == True:
            self.MWNT.SpinBox_desired_bond_length.setEnabled(False)
            MWNT_type_1 = self.MWNT.comboBox_type1_MWNT.currentText()
            MWNT_type_2 = self.MWNT.comboBox_type2_MWNT.currentText()
            if MWNT_type_1=="C" and MWNT_type_2=="C":
                bond_length_MWNT = 1.421 # default value of C-C bond length
                min = 1.3
                max = 2.0
            if MWNT_type_1=="N" and MWNT_type_2=="B":
                bond_length_MWNT = 1.47 # default value of N-B bond length
                min = 1.3
                max = 2.0
            if MWNT_type_1=="N" and MWNT_type_2=="Ga":
                bond_length_MWNT = 1.95 # default value of N-Ga bond length
                min = 1.9
                max = 2.5
            if MWNT_type_1=="N" and MWNT_type_2=="Al":
                bond_length_MWNT = 1.83 # default value of N-Al bond length
                min = 1.8
                max = 2.5
            if MWNT_type_1=="P" and MWNT_type_2=="Al":
                bond_length_MWNT = 2.3 # default value of P-Al bond length
                min = 1.8
                max = 2.6
            if MWNT_type_1=="P" and MWNT_type_2=="Ga":
                bond_length_MWNT = 2.28 # default value of P-Ga bond length
                min = 1.9
                max = 2.6
            if MWNT_type_1=="P" and MWNT_type_2=="C":
                bond_length_MWNT = 1.87 # default value of P-C bond length
                min = 1.6
                max = 2.3
            if MWNT_type_1=="N" and MWNT_type_2=="C":
                bond_length_MWNT = 1.47 # default value of N-C bond length
                min = 1.3
                max = 1.9
            if MWNT_type_1=="C" and MWNT_type_2=="B":
                bond_length_MWNT = 1.56 # default value of C-B bond length
                min = 1.3
                max = 2.0
            if MWNT_type_1=="C" and MWNT_type_2=="Al":
                bond_length_MWNT = 2.0 # default value of C-Al bond length
                min = 1.8
                max = 2.5
            if MWNT_type_1=="C" and MWNT_type_2=="Ga":
                bond_length_MWNT = 2.46 # default value of P-B bond length
                min = 1.9
                max = 2.6
            if MWNT_type_1=="P" and MWNT_type_2=="B":
                bond_length_MWNT = 1.74 # default value of P-B bond length
                min = 1.6
                max = 2.4

            self.MWNT.SpinBox_desired_bond_length.setRange(min, max)
            self.MWNT.lineEdit_bond_length_MWNT.setText(str(bond_length_MWNT))

        if self.MWNT.radioButton_desired_bond_length.isChecked() == True:
            self.MWNT.SpinBox_desired_bond_length.setRange(min, max)
            self.MWNT.lineEdit_bond_length_MWNT.setText(str(' '))
            self.MWNT.SpinBox_desired_bond_length.setEnabled(True)
            bond_length_MWNT = float(self.MWNT.SpinBox_desired_bond_length.text())

        return bond_length_MWNT


    def MWNT_diameter_changed(self):
        global bendFactor
        bond_length_MWNT = self.get_atom_type()
        value_n_MWNT = int(self.MWNT.spinBox_chirality_N_MWNT.text())
        value_m_MWNT = int(self.MWNT.spinBox_chirality_M_MWNT.text())
        repeat_units_MWNT = int(self.MWNT.spinBox_repeat_units_MWNT.text())
        a1 = np.array((np.sqrt(3)*bond_length_MWNT, 0))
        a2 = np.array((np.sqrt(3)/2*bond_length_MWNT, -3*bond_length_MWNT/2))
        Ch = value_n_MWNT*a1+value_m_MWNT*a2
        d = gcd(value_n_MWNT, value_m_MWNT)
        dR = 3*d if (value_n_MWNT-value_m_MWNT) % (3*d) == 0 else d
        t1 = (2*value_m_MWNT+value_n_MWNT)//dR
        t2 = -(2*value_n_MWNT+value_m_MWNT)//dR
        T = t1*a1+t2*a2
        diameter_MWNT = float(np.linalg.norm(Ch)/np.pi)
        diameter_MWNT = "{:.2f}".format(diameter_MWNT)
        length_MWNT = np.linalg.norm(T) * repeat_units_MWNT
        length_MWNT = "{:.2f}".format(float(length_MWNT))
        self.MWNT.lineEdit_diameter_MWNT.setText(str(diameter_MWNT))
        self.MWNT.lineEdit_length_MWNT.setText(str(length_MWNT))
        number_of_walls = int(self.MWNT.spinBox_num_walls_MWNT.text())
        if number_of_walls > 1:
            wall_separation = 3.43
        else:
            wall_separation = 0.0

        self.MWNT.lineEdit_wall_separation_MWNT.setText(str(wall_separation))

    def initial_box_dim(self):
        global box_lx, box_ly, box_lz, bendFactor
        box_lx = float(self.MWNT.SpinBox_lx.text())
        box_ly = float(self.MWNT.SpinBox_ly.text())
        box_lz = float(self.MWNT.SpinBox_lz.text())
        return

    def MWNT_builder_callback(self):
        global box_lx, box_ly, box_lz, bendFactor
        bond_length_MWNT = self.get_atom_type()

        try:
            box_lx or box_ly or box_lz
        except NameError:
            box_lx = box_ly = box_lz = 0.0
        try:
            bendFactor
        except NameError:
            bendFactor = 1.0
        number_of_walls = int(self.MWNT.spinBox_num_walls_MWNT.text())
        value_n_MWNT = int(self.MWNT.spinBox_chirality_N_MWNT.text())
        value_m_MWNT = int(self.MWNT.spinBox_chirality_M_MWNT.text())
        repeat_units_MWNT = int(self.MWNT.spinBox_repeat_units_MWNT.text())
        MWNT_type_1 = self.MWNT.comboBox_type1_MWNT.currentText()
        MWNT_type_2 = self.MWNT.comboBox_type2_MWNT.currentText()

        list_universe = []
        for i in range(1, number_of_walls+1):
            universe_all = MWNT_builder(value_n_MWNT, value_m_MWNT, repeat_units_MWNT, a=bond_length_MWNT, species=(MWNT_type_1, MWNT_type_2), centered=True, wan = i)
            list_universe.append(universe_all.atoms)
        universe_all = mda.Merge(*list_universe)
        universe_all.trajectory.ts.dimensions = [box_lx, box_ly, box_lx, 90, 90, 90]
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(universe_all)
        window.show()

"""
  The numbers (n,m) show that your tube is obtained from taking one atom of the sheet and rolling it onto
  that atom that is at located na1+ma2 away from your original atom. N is the Number of hexagons in a unit cell. a1 and a2 are lattice vectors.
  (n,m=n) gives an “armchair” tube,e.g. (5,5). (n,m=0) gives an “zig-zag” tube, e.g. (6,0). Other tubes are “chiral”, e.g. (6,2)
"""

def MWNT_builder(n, m, N, a, species=('B', 'C'), centered=False, wan = 1):
    global box_lx, box_ly, box_lz, bendFactor
    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2))
    Ch = (n*a1+m*a2)
    T = t1*a1+t2*a2

    # diameter = (norm(Ch)/np.pi)*wan
    # Ch = (diameter*np.pi)/wan
    Ch_proj, T_proj = [(v/norm(v)**2) for v in [Ch*wan, T]]
    basis = [np.array((0, 0)), ((a1+a2)/3)]
    pts = []
    xyz = []
    atom_types_swnt = []
    for i1, i2 in product(range(0, wan*n+t1), range(t2, wan*m)):
        shift = i1*a1+i2*a2
        for sp, b in zip(species, basis):
            pt = b+shift
            if all(-thre < pt.dot(v) < 1-thre for v in [Ch_proj, T_proj]):
                for k in range(N):
                    pts.append((sp, pt+k*T))
    # diameter = (norm(Ch)/np.pi)*wan

    diameter = (norm(Ch)/(np.pi*bendFactor))*wan
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
    fragments = mol_1.to_json(scale =1.3)
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
