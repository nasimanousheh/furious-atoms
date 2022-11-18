from furiousatoms.io import create_universe, merged_two_universes
import numpy as np
from numpy.linalg import norm
from math import gcd
from itertools import product
import io
from furiousatoms import io
import numpy as np
from fury import window
from PySide2 import QtWidgets
from furiousatoms.io import  load_files
from fury import window, utils
from furiousatoms.structure import bbox

thre = 1e-10
vacuum = 4
class Ui_polymer(QtWidgets.QMainWindow):
    """ Ui_polymer class creates a widget for building multple-walls nanotube (polymer)
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_polymer, self).__init__(parent)
        self.polymer = io.load_ui_widget("polymer.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.polymer)
        self.setCentralWidget(self.polymer)
        self.setLayout(self.v_layout)
        self.resize(248, 313)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.polymer.comboBox_library.activated.connect(self.get_info_polymer)
        self.polymer.comboBox_repeat_unit.activated.connect(self.get_info_polymer)
        self.polymer.comboBox_tacticity.activated.connect(self.get_info_polymer)
        self.polymer.comboBox_chiral.activated.connect(self.get_info_polymer)
        self.polymer.comboBox_chain_length.activated.connect(self.get_info_polymer)
        self.polymer.comboBox_chain_number.activated.connect(self.get_info_polymer)
        self.polymer.pushButton_build_polymer.clicked.connect(self.polymer_builder_callback)
        self.polymer.pushButton_build_polymer.clicked.connect(lambda:self.close())

    def get_info_polymer(self):
        library = self.polymer.comboBox_library.currentText()
        # repeat_unit = self.polymer.comboBox_repeat_unit.currentText()
        # tacticity = self.polymer.comboBox_tacticity.currentText()
        # chiral = self.polymer.comboBox_chiral.currentText()
        # chain_length = int(self.polymer.comboBox_chain_length.text())
        # chain_number = int(self.polymer.comboBox_chain_number.text())

    def polymer_builder_callback(self):
        # library = self.polymer.comboBox_library.currentText()
        # repeat_unit = self.polymer.comboBox_repeat_unit.currentText()
        # tacticity = self.polymer.comboBox_tacticity.currentText()
        # chiral = self.polymer.comboBox_chiral.currentText()
        # chain_length = int(self.polymer.comboBox_chain_length.text())
        # chain_number = int(self.polymer.comboBox_chain_number.text())

        # fname_monomer = "polymers\acrylates\acrylic_acid.pdb"
        fname_monomer = 'C:/Users/nasim/Downloads/NewNS/gammaGraphyne_unitcell.pdb'
        universe_all = self.polymer_builder(fname_monomer)
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(universe_all)
        window.show()


    def polymer_builder(self, fname_monomer):
        load_file,_ = load_files(fname_monomer)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        first_carbon = pos[np.where(load_file.atoms.types=='C')][0]
        last_carbon = pos[np.where(load_file.atoms.types=='C')][-1]
        distance = max(pos[:, 0]) - min(pos[:, 0]) + 1.4
        # [ 4.97200012  2.09100008 -9.26099968]
        nx = load_file.dimensions [3]
        ny= load_file.dimensions [4]
        nz= load_file.dimensions [5]
        load_file.dimensions = [distance, distance, distance, nx, ny, nz]
        box = load_file.dimensions[:3]
        copied = []
        for x in range(3):
            for y in range(3):
                for z in range(1):
                    u_ = load_file.copy()
                    move_by = box*(x, y, z)
                    u_.atoms.translate(move_by)
                    copied.append(u_.atoms) #u_.bonds.indices
                    copied.add_bonds(universe_all.bonds.indices)

        import MDAnalysis as mda
        new_universe = mda.Merge(*copied)
        # new_box = box*(n_x, n_y, n_z)
        # new_universe.dimensions = list(new_box) + [90]*3
        return new_universe

        # for x in range(n_x):
        #     for y in range(n_y):
        #         for z in range(n_z):
        #             u_ = load_file.copy()
        #             move_by = box*(x, y, z)
        #             u_.atoms.translate(move_by)
        #             copied.append(u_.atoms)



        # new_universe = mda.Merge(*copied)
        # new_box = box*(n_x, n_y, n_z)
        # new_universe.dimensions = list(new_box) + [90]*3
        # atom_types_swnt = [v for v, _ in pts]

        # n_atoms_swnt = len(xyz)
        # coord_array_swnt = np.array(xyz)
        # assert coord_array_swnt.shape == (n_atoms_swnt, 3)
        # all_bonds_swnt = np.array(fragments['bonds'])
        # try:
        #     box_lx or box_ly or box_lz
        # except NameError:
        #     box_lx = box_ly = box_lz = 0.0
        # univ_polymer = create_universe(coord_array_polymer, all_bonds_swnt, atom_types_swnt, box_lx, box_ly, box_lz)
        # return univ_polymer
