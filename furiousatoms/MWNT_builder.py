from MDAnalysis import *
import MDAnalysis.analysis.align
import numpy as np
from numpy.linalg import norm
from fractions import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
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

thre = 1e-10
vacuum = 4
SM = SharedMemory()

class Ui_MWNT(QtWidgets.QMainWindow): #QWidget
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
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        bond_length_MWNT = 1.421 # default value of C-C bond length
        self.MWNT.lineEdit_bond_length_MWNT.insert(str(bond_length_MWNT))
        self.MWNT.lineEdit_bond_length_MWNT.textChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.spinBox_chirality_N_MWNT.valueChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.spinBox_chirality_M_MWNT.valueChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.spinBox_repeat_units_MWNT.valueChanged.connect(self.MWNT_diameter_changed)
        self.MWNT.pushButton_build_MWNT.clicked.connect(self.MWNT_builder_callback)

    def MWNT_diameter_changed(self):
        SM.bond_length_MWNT = self.MWNT.lineEdit_bond_length_MWNT.text()
        SM.bond_length_MWNT = float(SM.bond_length_MWNT)
        value_n_MWNT = int(self.MWNT.spinBox_chirality_N_MWNT.text())
        value_m_MWNT = int(self.MWNT.spinBox_chirality_M_MWNT.text())
        repeat_units_MWNT = int(self.MWNT.spinBox_repeat_units_MWNT.text())
        a1 = np.array((np.sqrt(3)*SM.bond_length_MWNT, 0))
        a2 = np.array((np.sqrt(3)/2*SM.bond_length_MWNT, -3*SM.bond_length_MWNT/2))
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

    def MWNT_builder_callback(self):
        SM.number_of_walls = int(self.MWNT.spinBox_num_walls_MWNT.text())
        value_n_MWNT = int(self.MWNT.spinBox_chirality_N_MWNT.text())
        value_m_MWNT = int(self.MWNT.spinBox_chirality_M_MWNT.text())
        repeat_units_MWNT = int(self.MWNT.spinBox_repeat_units_MWNT.text())
        MWNT_type_1 = self.MWNT.comboBox_type1_MWNT.currentText()
        MWNT_type_2 = self.MWNT.comboBox_type2_MWNT.currentText()
        i = 1
        universe = MWNT_builder(value_n_MWNT, value_m_MWNT, repeat_units_MWNT, length=None, a=SM.bond_length_MWNT, species=(MWNT_type_1, MWNT_type_2), centered=True)
        for i in range(SM.number_of_walls):
            next_universe = MWNT_builder(value_n_MWNT + (6*i), value_m_MWNT + (6*i), repeat_units_MWNT, length=None, a=SM.bond_length_MWNT, species=(MWNT_type_1, MWNT_type_2), centered=True)
            structure_info = MDAnalysis.Merge(universe.atoms, next_universe.atoms)
        # file_name = 'fname.pdb'
        # universe.atoms.write(file_name)
        # self.win.process_load_file(fname=file_name)
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(structure_info)
        window.show()

"""
  The numbers (n,m) show that your tube is obtained from taking one atom of the sheet and rolling it onto
  that atom that is at located na1+ma2 away from your original atom. N is the Number of hexagons in a unit cell. a1 and a2 are lattice vectors.
  (n,m=n) gives an “armchair” tube,e.g. (5,5). (n,m=0) gives an “zig-zag” tube, e.g. (6,0). Other tubes are “chiral”, e.g. (6,2)
"""



def MWNT_builder(n, m, N=1, length=True, a=1.421, species=('B', 'C'), centered=False):
    d = gcd(n, m)
    dR = 3*d if (n-m) % (3*d) == 0 else d
    t1 = (2*m+n)//dR
    t2 = -(2*n+m)//dR
    a1 = np.array((np.sqrt(3)*a, 0))
    a2 = np.array((np.sqrt(3)/2*a, -3*a/2))
    Ch = n*a1+m*a2
    T = t1*a1+t2*a2
    # if length:
    #     N = int(np.ceil(length/norm(T)))
    Ch_proj, T_proj = [v/norm(v)**2 for v in [Ch, T]]
    basis = [np.array((0, 0)), (a1+a2)/3]
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
    diameter = norm(Ch)/np.pi
    print(diameter)
    def gr2tube(v):
        phi = 2*np.pi*v.dot(Ch_proj)
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
    swnt = MDAnalysis.Universe.empty(n_atoms_swnt, trajectory=True, n_residues=1)
    n_residues = 1
    swnt.atoms.positions = coord_array_swnt
    all_bonds_swnt = np.array(fragments['bonds'])
    swnt.add_TopologyAttr('name', atom_types_swnt)
    swnt.add_TopologyAttr('type', atom_types_swnt)
    swnt.add_TopologyAttr('resname', ['MOL']*n_residues)
    swnt.add_bonds(all_bonds_swnt)
    return swnt
