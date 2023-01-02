import numpy as np
from furiousatoms import io
from fury import window
from PySide2 import QtWidgets
from furiousatoms.io import  load_files
import MDAnalysis as mda

thre = 1e-10
vacuum = 4
class Ui_graphyne(QtWidgets.QMainWindow):
    """ Ui_graphyne class creates a widget for building multple-walls nanotube (graphyne)
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_graphyne, self).__init__(parent)
        self.graphyne = io.load_ui_widget("graphyne.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.graphyne)
        self.setCentralWidget(self.graphyne)
        self.setLayout(self.v_layout)
        self.resize(248, 313)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.graphyne.comboBox_graphyne.activated.connect(self.get_info_graphyne)
        self.graphyne.doubleSpinBox_lx_extent.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.doubleSpinBox_ly_extent.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.SpinBox_lx.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.SpinBox_ly.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.SpinBox_ly.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.SpinBox_num_sheets.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.pushButton_build_graphyne.clicked.connect(self.graphyne_builder_callback)
        self.graphyne.pushButton_build_graphyne.clicked.connect(lambda:self.close())

    def get_info_graphyne(self):
        num_unitcell_in_lx = float(self.graphyne.doubleSpinBox_lx_extent.text())
        num_unitcell_in_ly = float(self.graphyne.doubleSpinBox_ly_extent.text())
        try:
            box_lx = float(self.graphyne.SpinBox_lx.text())
            box_ly = float(self.graphyne.SpinBox_ly.text())
            box_lz = float(self.graphyne.SpinBox_lz.text())
        except:
            box_lx = box_lx = box_lz = 0

        try:
            num_sheets = int(self.graphyne.SpinBox_num_sheets.text())
        except:
            num_sheets = 1
        try:
            sheet_separation = float(self.graphyne.SpinBox_sheet_sep.text())
        except NameError:
            sheet_separation = 3.347

        return num_unitcell_in_lx, num_unitcell_in_ly, num_sheets

    def graphyne_builder_callback(self):
        num_unitcell_in_lx = float(self.graphyne.doubleSpinBox_lx_extent.text())
        num_unitcell_in_ly = float(self.graphyne.doubleSpinBox_ly_extent.text())
        num_unitcell_in_lx = int(num_unitcell_in_lx)
        num_unitcell_in_ly = int(num_unitcell_in_ly)

        try:
            num_sheets = int(self.graphyne.SpinBox_num_sheets.text())
        except:
            num_sheets = 1
        try:
            sheet_separation = float(self.graphyne.SpinBox_sheet_sep.text())
        except NameError:
            sheet_separation = 3.347

        graphyne_type = self.graphyne.comboBox_graphyne.currentText()
        if graphyne_type =="β-graphyne":
            fname = 'furiousatoms/graphyne_dataset/betaGraphyne_unitcell.pdb'
            structure_info = self.beta_graphyne_builder(fname, num_unitcell_in_lx, num_unitcell_in_ly, num_sheets)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)
        if graphyne_type =="γ-graphyne":
            fname = 'furiousatoms/graphyne_dataset/gammaGraphyne_unitcell.pdb'
            structure_info = self.gamma_graphyne_builder(fname, num_unitcell_in_lx, num_unitcell_in_ly, num_sheets)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)
        if graphyne_type =="graphyne_6_6_12":
            fname = 'furiousatoms/graphyne_dataset/6-6-12-graphyne_unitcell.pdb'
            structure_info = self.graphyne_6_6_12_builder(fname, num_unitcell_in_lx, num_unitcell_in_ly, num_sheets)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)
        cog = universe_all.atoms.center_of_geometry()
        print('Original solvent center of geometry: ', cog)
        universe_all.atoms.positions -= cog
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(universe_all)
        window.show()


    def beta_graphyne_builder(self, fname, num_unitcell_in_lx, num_unitcell_in_ly, num_sheets):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        distance = np.linalg.norm(pos[4]-pos[11]) + 1.4
        box_lx = load_file.dimensions[3]
        box_ly= load_file.dimensions[4]
        box_lz= load_file.dimensions[5]
        load_file.dimensions = [distance  ,  9.485-1.285, distance, box_lx, box_ly, box_lz]
        box = load_file.dimensions[:3]
        copied = []
        b = 0
        i = 0
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                u_ = load_file.copy()
                move_by = box*(x-i, y, 1)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)
                i = i-(0.5)

        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        bond_connect = 18 * num_unitcell_in_ly

        num_bonds_connect = (num_unitcell_in_lx-1)

        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(bond_connect*i)+4+b, (bond_connect*(i+1))+11+b]])
                if c < num_unitcell_in_ly:
                    added_bonds_2 = np.array([[(bond_connect*i)+2+b, (bond_connect*i)+b+(9+bond_connect)]])
                    new_universe.add_bonds(added_bonds_2)
                new_universe.add_bonds(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_3 = np.array([[(bond_connect*j)+16+b, (bond_connect*j)+b+32]])
                    new_universe.add_bonds(added_bonds_3)
                    added_bonds_4 = np.array([[(bond_connect*j)+17+b, (bond_connect*j)+b+33]])
                    new_universe.add_bonds(added_bonds_4)
            c = c + 1

        return new_universe

    def gamma_graphyne_builder(self, fname, num_unitcell_in_lx, num_unitcell_in_ly, num_sheets):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        distance_2 = max(pos[:, 0]) - min(pos[:, 0]) + 0.458509564
        distance = max(pos[:, 0]) - min(pos[:, 0]) + 1.4    #1.4000015258789062  #1.4000014682858655
        print(max(pos[:, 0]) - min(pos[:, 0]))
        nx = load_file.dimensions [3]
        ny= load_file.dimensions [4]
        nz= load_file.dimensions [5]
        load_file.dimensions = [distance, distance_2, distance, nx, ny, nz]
        box = load_file.dimensions[:3]
        copied = []
        b = 0
        i = 0
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                u_ = load_file.copy()
                move_by = box*(x-i, y, 1) #- (x * 2*1.73205080757) - (1.4*2)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)
                i = i+(0.5)

        import MDAnalysis as mda
        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        bond_connect = 12 * num_unitcell_in_ly

        num_bonds_connect = (num_unitcell_in_lx-1)

        for b in range(num_unitcell_in_ly):
            b = c *12
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(bond_connect*i)+8+b, (bond_connect*(i+1))+7+b]])
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(bond_connect*i)+3+b, (bond_connect*i)+b+(13+bond_connect)]])
                    new_universe.add_bonds(added_bonds_2)
                new_universe.add_bonds(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_1 = np.array([[(bond_connect*j)+2+b, (bond_connect*j)+b+12]])
                    new_universe.add_bonds(added_bonds_1)

            c = c + 1


        return new_universe

    def graphyne_6_6_12_builder(self, fname, num_unitcell_in_lx, num_unitcell_in_ly, num_sheets):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        distance = max(pos[:, 0]) - min(pos[:, 0]) + 1.23
        distance_2 = np.linalg.norm(pos[10]-pos[17]) + 1.4
        box_lx = load_file.dimensions [3]
        box_ly= load_file.dimensions [4]
        box_lz= load_file.dimensions [5]
        load_file.dimensions = [distance, distance_2, distance, box_lx, box_ly, box_lz]
        box = load_file.dimensions[:3]
        copied = []
        b = 0
        for x in range(num_unitcell_in_lx):
            for y in range(num_unitcell_in_ly):
                u_ = load_file.copy()
                move_by = box*(x, y, 1)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)

        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        bond_connect = 18 * num_unitcell_in_ly

        num_bonds_connect = (num_unitcell_in_lx-1)

        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(bond_connect*i)+b, (bond_connect*(i+1))+5+b]])
                added_bonds_1 = np.array([[(bond_connect*i)+1+b, (bond_connect*(i+1))+4+b]])
                new_universe.add_bonds(added_bonds_1)
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(bond_connect*i)+10+b, (bond_connect*i)+b+35]])
                    new_universe.add_bonds(added_bonds_2)
                    added_bonds_3 = np.array([[(bond_connect*i)+12+b, (bond_connect*i)+b+33]])
                    new_universe.add_bonds(added_bonds_3)
                new_universe.add_bonds(added_bonds)
                #This loop connects the atoms of the last coloumn
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_4 = np.array([[(bond_connect*j)+10+b, (bond_connect*j)+b+35]])
                    new_universe.add_bonds(added_bonds_4)
                    added_bonds_5 = np.array([[(bond_connect*j)+12+b, (bond_connect*j)+b+33]])
                    new_universe.add_bonds(added_bonds_5)
            c = c + 1

        return new_universe

    def extend_the_sheets(self, structure_info, num_sheets, sheet_separation):
        try:
            box_lx = float(self.graphyne.SpinBox_lx.text())
            box_ly = float(self.graphyne.SpinBox_ly.text())
            box_lz = float(self.graphyne.SpinBox_lz.text())
        except:
            box_lx = box_ly = box_lz = 0

        copied = []
        structure_info.dimensions = [sheet_separation, sheet_separation, sheet_separation, 90, 90, 90]
        box_sep = structure_info.dimensions[:3]
        for x in range(1):
            for y in range(1):
                for z in range(num_sheets):
                    u_ = structure_info.copy()
                    move_by = box_sep*(x, y, z)
                    u_.atoms.translate(move_by)
                    copied.append(u_.atoms)

            extended_universe = mda.Merge(*copied)
            extended_universe.dimensions = [box_lx, box_ly, box_lz, 90, 90, 90]
            cog = extended_universe.atoms.center_of_geometry()
            extended_universe.atoms.positions -= cog
            return extended_universe



