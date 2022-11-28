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
        fname_monomer = 'C:/Users/nasim/Downloads/NewNS/6-6-12-graphyne_unitcell.pdb'

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
        distance = max(pos[:, 0]) - min(pos[:, 0]) + 1.23
        distance_2 = np.linalg.norm(pos[10]-pos[17]) + 1.4
        nx = load_file.dimensions [3]
        ny= load_file.dimensions [4]
        nz= load_file.dimensions [5]
        load_file.dimensions = [distance, distance_2, distance, nx, ny, nz]
        box = load_file.dimensions[:3]
        copied = []
        extend_in_x = 5
        extend_in_y = 5
        b = 0
        for x in range(extend_in_x):
            for y in range(extend_in_y):
                u_ = load_file.copy()
                move_by = box*(x, y, 1)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)

        import MDAnalysis as mda
        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        bond_connect = 18 * extend_in_y

        num_bonds_connect = (extend_in_x-1)

        for b in range(extend_in_y):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(bond_connect*i)+b, (bond_connect*(i+1))+5+b]])
                added_bonds_1 = np.array([[(bond_connect*i)+1+b, (bond_connect*(i+1))+4+b]])
                new_universe.add_bonds(added_bonds_1)
                if c < extend_in_y-1:
                    added_bonds_2 = np.array([[(bond_connect*i)+10+b, (bond_connect*i)+b+35]])
                    new_universe.add_bonds(added_bonds_2)
                    added_bonds_3 = np.array([[(bond_connect*i)+12+b, (bond_connect*i)+b+33]])
                    new_universe.add_bonds(added_bonds_3)
                new_universe.add_bonds(added_bonds)
                #This loop connects the atoms of the last coloumn
            for j in range(extend_in_x):
                if c < extend_in_y-1:
                    added_bonds_4 = np.array([[(bond_connect*j)+10+b, (bond_connect*j)+b+35]])
                    new_universe.add_bonds(added_bonds_4)
                    added_bonds_5 = np.array([[(bond_connect*j)+12+b, (bond_connect*j)+b+33]])
                    new_universe.add_bonds(added_bonds_5)

            c = c + 1


        copied_new = []
        for z in range(1):
            b_ = new_universe.copy()
            move_by = box*(x, y, z)
            b_.atoms.translate(move_by)
            copied_new.append(b_.atoms)

        new = mda.Merge(*copied_new)
        new_box = box*(nx, ny, nz)
        new_universe.dimensions = list(new_box) + [90]*3

        # new_universe.atoms.positions = xyz


        # def getPosOnBentLine(lineStart, lineEnd, t, bendFactor, pivot):
            # lineDir = lineEnd - lineStart
        # PI = 3.14
        # t=new_universe.atoms.positions
        # pivot=1
        # lineEnd = new_universe.atoms.positions[296]#[ 2.41700006 30.69400024 21.93000031]
        # bendFactor = 1
        # lineLength = 33.24999976158142 #len(lineDir)
        # lineStart = new_universe.atoms.positions[55] #[[ 2.41700006 30.69400024 21.93000031]]
        # circleRad = lineLength / (bendFactor * 2 * PI)
        # circleCenter = lineStart +  (lineEnd - lineStart)  * pivot + perp(lineDir) * circleRad

        # angle = PI + bendFactor * (1.0 - (t+pivot)) * 2 * PI
        # posOnCircle = circleCenter + (np.cos(angle), np.sin(angle)) * circleRad

        return new

        # return sol


        # for i in range(num_bonds_connect):
        #     added_bonds = np.array([[(bond_connect*i)+8, (bond_connect*(i+1))+7]])
        #     added_bonds_2 = np.array([[(bond_connect*i)+8+12, (bond_connect*(i+1))+7+12]])
        #     added_bonds_3 = np.array([[(bond_connect*i)+8+24, (bond_connect*(i+1))+7+24]])
        #     added_bonds_4 = np.array([[(bond_connect*i)+8+24+12, (bond_connect*(i+1))+7+24+12]])
            # new_universe.add_bonds(added_bonds)
            # new_universe.add_bonds(added_bonds_2)
            # new_universe.add_bonds(added_bonds_3)
            # new_universe.add_bonds(added_bonds_4)


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
