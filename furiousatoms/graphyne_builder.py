import numpy as np
from furiousatoms import io
from fury import window
from PySide2 import QtWidgets
from furiousatoms.io import  load_files
import MDAnalysis as mda
from PySide2.QtGui import QIcon
import os


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
        self.setWindowIcon(QIcon(io.get_resources_file("splash.png")))
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.graphyne.doubleSpinBox_lx_extent.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.doubleSpinBox_ly_extent.valueChanged.connect(self.get_info_graphyne)
        self.graphyne.pushButton_build_graphyne.clicked.connect(self.graphyne_builder_callback)
        self.graphyne.pushButton_build_graphyne.clicked.connect(lambda:self.close())

    def get_info_graphyne(self):
        edge_length_x = float(self.graphyne.doubleSpinBox_lx_extent.text())
        edge_length_y = float(self.graphyne.doubleSpinBox_ly_extent.text())
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

        return edge_length_x, edge_length_y, num_sheets

    def graphyne_builder_callback(self):
        edge_length_x = float(self.graphyne.doubleSpinBox_lx_extent.text())
        edge_length_y = float(self.graphyne.doubleSpinBox_ly_extent.text())

        try:
            num_sheets = int(self.graphyne.SpinBox_num_sheets.text())
        except:
            num_sheets = 1
        try:
            sheet_separation = float(self.graphyne.SpinBox_sheet_sep.text())
        except NameError:
            sheet_separation = 3.347

        graphyne_type = self.graphyne.comboBox_graphyne.currentText()
        dir_Graphyne_folder = io.get_frozen_path() if io.is_frozen() else io.get_application_path()
        Graphyne_folder = os.path.join(dir_Graphyne_folder, 'graphyne_dataset')
        if graphyne_type =="graphyne_12_12_12":
            fname= os.path.join(Graphyne_folder, 'betaGraphyne_unitcell.pdb')
            structure_info = self.beta_graphyne_builder(fname, edge_length_x, edge_length_y)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)
        if graphyne_type =="graphyne-1":
            fname=os.path.join(Graphyne_folder, 'gammaGraphyne_unitcell.pdb')
            structure_info = self.gamma_graphyne_builder(fname, edge_length_x, edge_length_y)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)
        if graphyne_type =="graphyne-2":
            fname=os.path.join(Graphyne_folder, 'graphdiyne_unitcell.pdb')
            structure_info = self.graphyne_2_builder(fname, edge_length_x, edge_length_y)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)
        if graphyne_type =="graphyne_6_6_12":
            fname=os.path.join(Graphyne_folder, '6-6-12-graphyne_unitcell.pdb')
            structure_info = self.graphyne_6_6_12_builder(fname, edge_length_x, edge_length_y)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)
        if graphyne_type =="twin Graphene":
            fname=os.path.join(Graphyne_folder, 'twinGraphene_unitcell.pdb')
            sheet_separation = sheet_separation + 2.034
            structure_info = self.twin_graphene_builder(fname, edge_length_x, edge_length_y)
            universe_all = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)

        cog = universe_all.atoms.center_of_geometry()
        universe_all.atoms.positions -= cog
        window = self.win.create_mdi_child()
        window.make_title()
        window.load_universe(universe_all)
        window.show()

    def beta_graphyne_builder(self, fname, edge_length_x, edge_length_y):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        unit_cell_lx = np.linalg.norm(pos[4]-pos[11]) + 1.4
        box_lx = load_file.dimensions[3]
        box_ly= load_file.dimensions[4]
        box_lz= load_file.dimensions[5]
        unit_cell_ly = 9.485-1.285
        load_file.dimensions = [unit_cell_lx  , unit_cell_ly, unit_cell_lx, box_lx, box_ly, box_lz]
        box = load_file.dimensions[:3]
        copied = []
        i = 0
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
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
        num_atoms_in_y_direction = 18 * num_unitcell_in_ly

        num_bonds_connect = (num_unitcell_in_lx-1)

        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+4+b, (num_atoms_in_y_direction*(i+1))+11+b]])
                if c < num_unitcell_in_ly:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+2+b, (num_atoms_in_y_direction*i)+b+(9+num_atoms_in_y_direction)]])
                    new_universe.add_bonds(added_bonds_2)
                new_universe.add_bonds(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_3 = np.array([[(num_atoms_in_y_direction*j)+16+b, (num_atoms_in_y_direction*j)+b+32]])
                    new_universe.add_bonds(added_bonds_3)
                    added_bonds_4 = np.array([[(num_atoms_in_y_direction*j)+17+b, (num_atoms_in_y_direction*j)+b+33]])
                    new_universe.add_bonds(added_bonds_4)
            c = c + 1

        return new_universe

    def gamma_graphyne_builder(self, fname, edge_length_x, edge_length_y):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        unit_cell_ly = max(pos[:, 0]) - min(pos[:, 0]) + 0.458509564
        unit_cell_lx = max(pos[:, 0]) - min(pos[:, 0]) + 1.4
        print(max(pos[:, 0]) - min(pos[:, 0]))
        nx = load_file.dimensions [3]
        ny= load_file.dimensions [4]
        nz= load_file.dimensions [5]
        load_file.dimensions = [unit_cell_lx, unit_cell_ly, unit_cell_lx, nx, ny, nz]
        box = load_file.dimensions[:3]
        copied = []
        i = 0
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                u_ = load_file.copy()
                move_by = box*(x-i, y, 1)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)
                i = i+(0.5)


        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        num_atoms_in_y_direction = 12 * num_unitcell_in_ly

        num_bonds_connect = (num_unitcell_in_lx-1)

        for b in range(num_unitcell_in_ly):
            b = c *12
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+8+b, (num_atoms_in_y_direction*(i+1))+7+b]])
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+3+b, (num_atoms_in_y_direction*i)+b+(13+num_atoms_in_y_direction)]])
                    new_universe.add_bonds(added_bonds_2)
                new_universe.add_bonds(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_1 = np.array([[(num_atoms_in_y_direction*j)+2+b, (num_atoms_in_y_direction*j)+b+12]])
                    new_universe.add_bonds(added_bonds_1)
            c = c + 1
        return new_universe

    def graphyne_2_builder(self, fname, edge_length_x, edge_length_y):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        unit_cell_ly = max(pos[:, 0]) - min(pos[:, 0]) -0.2
        unit_cell_lx = max(pos[:, 0]) - min(pos[:, 0]) + 1.42
        print(max(pos[:, 0]) - min(pos[:, 0]))
        nx = load_file.dimensions [3]
        ny= load_file.dimensions [4]
        nz= load_file.dimensions [5]
        load_file.dimensions = [unit_cell_lx, unit_cell_ly, unit_cell_lx, nx, ny, nz]
        box = load_file.dimensions[:3]
        copied = []
        i = 0
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                u_ = load_file.copy()
                move_by = box*(x-i, y, 1)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)
                i = i+(0.5)

        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        num_atoms_in_y_direction = 18 * num_unitcell_in_ly
        num_bonds_connect = (num_unitcell_in_lx-1)
        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+12+b, (num_atoms_in_y_direction*(i+1))+13+b]])
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+1+b, (num_atoms_in_y_direction*i)+b+(21+num_atoms_in_y_direction)]])
                    new_universe.add_bonds(added_bonds_2)
                new_universe.add_bonds(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_1 = np.array([[(num_atoms_in_y_direction*j)+b, (num_atoms_in_y_direction*j)+b+20]])
                    new_universe.add_bonds(added_bonds_1)
            c = c + 1
        return new_universe

    def graphyne_6_6_12_builder(self, fname, edge_length_x, edge_length_y):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        unit_cell_lx = max(pos[:, 0]) - min(pos[:, 0]) + 1.23
        unit_cell_ly = np.linalg.norm(pos[10]-pos[17]) + 1.4
        box_lx = load_file.dimensions [3]
        box_ly= load_file.dimensions [4]
        box_lz= load_file.dimensions [5]
        load_file.dimensions = [unit_cell_lx, unit_cell_ly, unit_cell_lx, box_lx, box_ly, box_lz]
        box = load_file.dimensions[:3]
        copied = []
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
        for x in range(num_unitcell_in_lx):
            for y in range(num_unitcell_in_ly):
                u_ = load_file.copy()
                move_by = box*(x, y, 1)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)

        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        num_atoms_in_y_direction = 18 * num_unitcell_in_ly
        num_bonds_connect = (num_unitcell_in_lx-1)
        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+b, (num_atoms_in_y_direction*(i+1))+5+b]])
                added_bonds_1 = np.array([[(num_atoms_in_y_direction*i)+1+b, (num_atoms_in_y_direction*(i+1))+4+b]])
                new_universe.add_bonds(added_bonds_1)
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+10+b, (num_atoms_in_y_direction*i)+b+35]])
                    new_universe.add_bonds(added_bonds_2)
                    added_bonds_3 = np.array([[(num_atoms_in_y_direction*i)+12+b, (num_atoms_in_y_direction*i)+b+33]])
                    new_universe.add_bonds(added_bonds_3)
                new_universe.add_bonds(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_4 = np.array([[(num_atoms_in_y_direction*j)+10+b, (num_atoms_in_y_direction*j)+b+35]])
                    new_universe.add_bonds(added_bonds_4)
                    added_bonds_5 = np.array([[(num_atoms_in_y_direction*j)+12+b, (num_atoms_in_y_direction*j)+b+33]])
                    new_universe.add_bonds(added_bonds_5)
            c = c + 1
        return new_universe

    def twin_graphene_builder(self, fname, edge_length_x, edge_length_y):
        load_file,_ = load_files(fname)
        load_file.bonds.to_indices()
        pos = load_file.atoms.positions
        pos = pos.astype('float64')
        unit_cell_lx = max(pos[:, 0]) - min(pos[:, 0]) + 1.421
        unit_cell_ly = max(pos[:, 0]) - min(pos[:, 0]) + 0.53
        box_lx = load_file.dimensions [3]
        box_ly= load_file.dimensions [4]
        box_lz= load_file.dimensions [5]
        load_file.dimensions = [unit_cell_lx, unit_cell_ly, unit_cell_lx, box_lx, box_ly, box_lz]
        box = load_file.dimensions[:3]
        copied = []
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                u_ = load_file.copy()
                move_by = box*(x-i, y, 1)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)
                i = i+(0.5)

        new_universe = mda.Merge(*copied)
        b = 0
        c = 0
        num_atoms_in_y_direction = 18 * num_unitcell_in_ly
        num_bonds_connect = (num_unitcell_in_lx-1)
        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+6+b, (num_atoms_in_y_direction*(i+1))+b]])
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+8+b, (num_atoms_in_y_direction*i)+b+(20+num_atoms_in_y_direction)]])
                    new_universe.add_bonds(added_bonds_2)
                new_universe.add_bonds(added_bonds)

            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_1 = np.array([[(num_atoms_in_y_direction*j)+10+b, (num_atoms_in_y_direction*j)+b+22]])
                    new_universe.add_bonds(added_bonds_1)
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



