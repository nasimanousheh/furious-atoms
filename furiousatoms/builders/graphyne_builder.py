import numpy as np
from furiousatoms import io
from fury import window
from PySide2 import QtWidgets
from furiousatoms.io import  load_files
from PySide2.QtGui import QIcon
import os
from furiousatoms.molecular import POS_DIM, BOND_DIM, MolecularStructure
from furiousatoms.builders.builder_util import copy_bonds


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
        graphyne_type = self.graphyne.comboBox_graphyne.currentText()

        graphyne = self.build_graphyne(edge_length_x, edge_length_y, graphyne_type)
        
        try:
            box_lx = float(self.graphyne.SpinBox_lx.text())
            box_ly = float(self.graphyne.SpinBox_ly.text())
            box_lz = float(self.graphyne.SpinBox_lz.text())
        except:
            box_lx = box_lx = box_lz = 0
        graphyne.box_size = [box_lx, box_ly, box_lz]

        window = self.win.create_mdi_child()
        window.make_title()
        window.load_structure(graphyne)
        window.show()

        

        
    def build_graphyne(self, edge_length_x, edge_length_y, graphyne_type):
        try:
            num_sheets = int(self.graphyne.SpinBox_num_sheets.text())
        except:
            num_sheets = 1
        try:
            sheet_separation = float(self.graphyne.SpinBox_sheet_sep.text())
        except NameError:
            sheet_separation = 3.347

        dir_Graphyne_folder = io.get_frozen_path() if io.is_frozen() else io.get_application_path()
        Graphyne_folder = os.path.join(dir_Graphyne_folder, 'graphyne_dataset')
        universe_all = None
        if graphyne_type =="graphyne_12_12_12":
            fname= os.path.join(Graphyne_folder, 'betaGraphyne_unitcell.pdb')
            builder = self.beta_graphyne_builder
        elif graphyne_type =="graphyne-1":
            fname=os.path.join(Graphyne_folder, 'gammaGraphyne_unitcell.pdb')
            builder = self.gamma_graphyne_builder
        elif graphyne_type =="graphyne-2":
            fname=os.path.join(Graphyne_folder, 'graphdiyne_unitcell.pdb')
            builder = self.graphyne_2_builder
        elif graphyne_type =="graphyne_6_6_12":
            fname=os.path.join(Graphyne_folder, '6-6-12-graphyne_unitcell.pdb')
            builder = self.graphyne_6_6_12_builder
        elif graphyne_type =="twin Graphene":
            fname=os.path.join(Graphyne_folder, 'twinGraphene_unitcell.pdb')
            sheet_separation = sheet_separation + 2.034
            builder = self.twin_graphene_builder
        else:
            raise ValueError("Illegal graphyne_type `" + graphyne_type + "`")

        structure_info = load_files(fname)
        structure_info = builder(structure_info, edge_length_x, edge_length_y)
        structure_info = self.extend_the_sheets(structure_info, num_sheets, sheet_separation)

        return structure_info




    def beta_graphyne_builder(self, s: MolecularStructure, edge_length_x, edge_length_y):
        unit_cell_lx = np.linalg.norm(s.pos[4] - s.pos[11]) + 1.4
        unit_cell_ly = 9.485-1.285
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
        
        UNIT_ATOM_COUNT = len(s.pos)
        NEW_ATOM_COUNT = UNIT_ATOM_COUNT * num_unitcell_in_lx * num_unitcell_in_ly
        copied_pos = np.zeros(shape=(NEW_ATOM_COUNT, 3))
        copied_atom_types = np.zeros(shape=(NEW_ATOM_COUNT), dtype=str)

        box = np.array([unit_cell_lx, unit_cell_ly, 0])
        i = 0
        atomId = 0
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                move_by = box*(x-i, y, 1)
                for j, atom in enumerate(s.pos): #j is an index in pos whereas atomId is an index in copied_pos
                    for k in range(0, 3):
                        copied_pos[atomId][k] = atom[k] + move_by[k]
                    copied_atom_types[atomId] = s.atom_types[j]
                    atomId += 1
                i = i-(0.5)

        s.pos = copied_pos
        s.atom_types = copied_atom_types
        
        s.bonds = copy_bonds(s.bonds, num_unitcell_in_lx, num_unitcell_in_ly, UNIT_ATOM_COUNT)

        b = 0
        c = 0
        new_bonds = []
        num_atoms_in_y_direction = UNIT_ATOM_COUNT * num_unitcell_in_ly

        num_bonds_connect = (num_unitcell_in_lx-1)

        for b in range(num_unitcell_in_ly):
            b = c * UNIT_ATOM_COUNT
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+4+b, (num_atoms_in_y_direction*(i+1))+11+b]])
                if c < num_unitcell_in_ly:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+2+b, (num_atoms_in_y_direction*i)+b+(9+num_atoms_in_y_direction)]])
                    new_bonds.append(added_bonds_2)
                new_bonds.append(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_3 = np.array([[(num_atoms_in_y_direction*j)+16+b, (num_atoms_in_y_direction*j)+b+32]])
                    new_bonds.append(added_bonds_3)
                    added_bonds_4 = np.array([[(num_atoms_in_y_direction*j)+17+b, (num_atoms_in_y_direction*j)+b+33]])
                    new_bonds.append(added_bonds_4)
            c = c + 1
        for bond in new_bonds:
            s.bonds = np.vstack((s.bonds, np.reshape(bond, (-1, 2)) ))

        return s
    
    def gamma_graphyne_builder(self, s: MolecularStructure, edge_length_x, edge_length_y):
        unit_cell_lx = max(s.pos[:, 0]) - min(s.pos[:, 0]) + 1.4
        unit_cell_ly = max(s.pos[:, 0]) - min(s.pos[:, 0]) + 0.458509564
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))

        UNIT_ATOM_COUNT = len(s.pos)
        NEW_ATOM_COUNT = UNIT_ATOM_COUNT * num_unitcell_in_lx * num_unitcell_in_ly
        copied_pos = np.zeros(shape=(NEW_ATOM_COUNT, 3))
        copied_atom_types = np.zeros(shape=(NEW_ATOM_COUNT), dtype=str)

        box = np.array([unit_cell_lx, unit_cell_ly, 0])
        i = 0
        atomId = 0
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                move_by = box*(x-i, y, 1)
                for j, atom in enumerate(s.pos): #j is an index in pos whereas atomId is an index in copied_pos
                    for k in range(0, 3):
                        copied_pos[atomId][k] = atom[k] + move_by[k]
                    copied_atom_types[atomId] = s.atom_types[j]
                    atomId += 1
                i = i+(0.5)

        s.pos = copied_pos
        s.atom_types = copied_atom_types
        
        s.bonds = copy_bonds(s.bonds, num_unitcell_in_lx, num_unitcell_in_ly, UNIT_ATOM_COUNT)

        b = 0
        c = 0
        num_atoms_in_y_direction = UNIT_ATOM_COUNT * num_unitcell_in_ly
        num_bonds_connect = (num_unitcell_in_lx - 1)
        new_bonds = []
        for b in range(num_unitcell_in_ly):
            b = c * UNIT_ATOM_COUNT
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+8+b, (num_atoms_in_y_direction*(i+1))+7+b]])
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+3+b, (num_atoms_in_y_direction*i)+b+(13+num_atoms_in_y_direction)]])
                    new_bonds.append(added_bonds_2)
                new_bonds.append(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_1 = np.array([[(num_atoms_in_y_direction*j)+2+b, (num_atoms_in_y_direction*j)+b+12]])
                    new_bonds.append(added_bonds_1)
            c = c + 1
        for bond in new_bonds:
            s.bonds = np.vstack((s.bonds, np.reshape(bond, (-1, 2)) ))

        return s

    def graphyne_2_builder(self, s: MolecularStructure, edge_length_x, edge_length_y):
        unit_cell_ly = max(s.pos[:, 0]) - min(s.pos[:, 0]) -0.2
        unit_cell_lx = max(s.pos[:, 0]) - min(s.pos[:, 0]) + 1.42
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
        
        UNIT_ATOM_COUNT = len(s.pos)
        NEW_ATOM_COUNT = UNIT_ATOM_COUNT * num_unitcell_in_lx * num_unitcell_in_ly
        copied_pos = np.zeros(shape=(NEW_ATOM_COUNT, 3))
        copied_atom_types = np.zeros(shape=(NEW_ATOM_COUNT), dtype=str)

        box = np.array([unit_cell_lx, unit_cell_ly, 0])
        i = 0
        atomId = 0
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                move_by = box*(x-i, y, 1)
                for j, atom in enumerate(s.pos): #j is an index in pos whereas atomId is an index in copied_pos
                    for k in range(0, 3):
                        copied_pos[atomId][k] = atom[k] + move_by[k]
                    copied_atom_types[atomId] = s.atom_types[j]
                    atomId += 1
                i = i+(0.5)

        s.pos = copied_pos
        s.atom_types = copied_atom_types
        
        s.bonds = copy_bonds(s.bonds, num_unitcell_in_lx, num_unitcell_in_ly, UNIT_ATOM_COUNT)

        b = 0
        c = 0
        num_atoms_in_y_direction = 18 * num_unitcell_in_ly
        num_bonds_connect = (num_unitcell_in_lx-1)
        new_bonds = []
        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+12+b, (num_atoms_in_y_direction*(i+1))+13+b]])
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+1+b, (num_atoms_in_y_direction*i)+b+(21+num_atoms_in_y_direction)]])
                    new_bonds.append(added_bonds_2)
                new_bonds.append(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_1 = np.array([[(num_atoms_in_y_direction*j)+b, (num_atoms_in_y_direction*j)+b+20]])
                    new_bonds.append(added_bonds_1)
            c = c + 1
        for bond in new_bonds:
            s.bonds = np.vstack((s.bonds, np.reshape(bond, (-1, 2)) ))

        return s

    def graphyne_6_6_12_builder(self, s: MolecularStructure, edge_length_x, edge_length_y):
        #TODO make unit_cell_lx/y more precise--use exact units 
        unit_cell_lx = max(s.pos[:, 0]) - min(s.pos[:, 0]) + 1.23
        unit_cell_ly = np.linalg.norm(s.pos[10]-s.pos[17]) + 1.4
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))

        UNIT_ATOM_COUNT = len(s.pos)
        NEW_ATOM_COUNT = UNIT_ATOM_COUNT * num_unitcell_in_lx * num_unitcell_in_ly
        copied_pos = np.zeros(shape=(NEW_ATOM_COUNT, 3))
        copied_atom_types = np.zeros(shape=(NEW_ATOM_COUNT), dtype=str)

        box = np.array([unit_cell_lx, unit_cell_ly, 0])
        i = 0
        atomId = 0
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                move_by = box*(x, y, 1)
                for j, atom in enumerate(s.pos): #j is an index in pos whereas atomId is an index in copied_pos
                    for k in range(0, 3):
                        copied_pos[atomId][k] = atom[k] + move_by[k]
                    copied_atom_types[atomId] = s.atom_types[j]
                    atomId += 1
                i = i+(0.5)

        s.pos = copied_pos
        s.atom_types = copied_atom_types

        s.bonds = copy_bonds(s.bonds, num_unitcell_in_lx, num_unitcell_in_ly, UNIT_ATOM_COUNT)
        
        b = 0
        c = 0
        num_atoms_in_y_direction = 18 * num_unitcell_in_ly
        num_bonds_connect = (num_unitcell_in_lx-1)
        new_bonds = []
        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+b, (num_atoms_in_y_direction*(i+1))+5+b]])
                added_bonds_1 = np.array([[(num_atoms_in_y_direction*i)+1+b, (num_atoms_in_y_direction*(i+1))+4+b]])
                new_bonds.append(added_bonds_1)
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+10+b, (num_atoms_in_y_direction*i)+b+35]])
                    new_bonds.append(added_bonds_2)
                    added_bonds_3 = np.array([[(num_atoms_in_y_direction*i)+12+b, (num_atoms_in_y_direction*i)+b+33]])
                    new_bonds.append(added_bonds_3)
                new_bonds.append(added_bonds)
            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_4 = np.array([[(num_atoms_in_y_direction*j)+10+b, (num_atoms_in_y_direction*j)+b+35]])
                    new_bonds.append(added_bonds_4)
                    added_bonds_5 = np.array([[(num_atoms_in_y_direction*j)+12+b, (num_atoms_in_y_direction*j)+b+33]])
                    new_bonds.append(added_bonds_5)
            c = c + 1
        for bond in new_bonds:
            s.bonds = np.vstack((s.bonds, np.reshape(bond, (-1, 2)) ))

        return s

    def twin_graphene_builder(self, s: MolecularStructure, edge_length_x, edge_length_y):
        unit_cell_lx = max(s.pos[:, 0]) - min(s.pos[:, 0]) + 1.421
        unit_cell_ly = max(s.pos[:, 0]) - min(s.pos[:, 0]) + 0.53
        num_unitcell_in_lx = int(np.floor(edge_length_x/unit_cell_lx))
        num_unitcell_in_ly = int(np.floor(edge_length_y/unit_cell_ly))
        
        UNIT_ATOM_COUNT = len(s.pos)
        NEW_ATOM_COUNT = UNIT_ATOM_COUNT * num_unitcell_in_lx * num_unitcell_in_ly
        copied_pos = np.zeros(shape=(NEW_ATOM_COUNT, 3))
        copied_atom_types = np.zeros(shape=(NEW_ATOM_COUNT), dtype=str)

        box = np.array([unit_cell_lx, unit_cell_ly, 0])
        i = 0
        atomId = 0
        for x in range(num_unitcell_in_lx):
            i = 0
            for y in range(num_unitcell_in_ly):
                move_by = box*(x-i, y, 1)
                for j, atom in enumerate(s.pos): #j is an index in pos whereas atomId is an index in copied_pos
                    for k in range(0, 3):
                        copied_pos[atomId][k] = atom[k] + move_by[k]
                    copied_atom_types[atomId] = s.atom_types[j]
                    atomId += 1
                i = i+(0.5)
                
        s.pos = copied_pos
        s.atom_types = copied_atom_types
        
        s.bonds = copy_bonds(s.bonds, num_unitcell_in_lx, num_unitcell_in_ly, UNIT_ATOM_COUNT)
        
        b = 0
        c = 0
        num_atoms_in_y_direction = 18 * num_unitcell_in_ly
        num_bonds_connect = (num_unitcell_in_lx-1)
        new_bonds = []
        for b in range(num_unitcell_in_ly):
            b = c *18
            for i in range(num_bonds_connect):
                added_bonds = np.array([[(num_atoms_in_y_direction*i)+6+b, (num_atoms_in_y_direction*(i+1))+b]])
                if c < num_unitcell_in_ly-1:
                    added_bonds_2 = np.array([[(num_atoms_in_y_direction*i)+8+b, (num_atoms_in_y_direction*i)+b+(20+num_atoms_in_y_direction)]])
                    new_bonds.append(added_bonds_2)
                new_bonds.append(added_bonds)

            for j in range(num_unitcell_in_lx):
                if c < num_unitcell_in_ly-1:
                    added_bonds_1 = np.array([[(num_atoms_in_y_direction*j)+10+b, (num_atoms_in_y_direction*j)+b+22]])
                    new_bonds.append(added_bonds_1)
            c = c + 1
        for bond in new_bonds:
            s.bonds = np.vstack((s.bonds, np.reshape(bond, (-1, 2)) ))
            
        return s

    def extend_the_sheets(self, structure_info, num_sheets, sheet_separation):
        s = structure_info
        if num_sheets > 1:
            ATOM_COUNT = len(s.pos)
            new_positions = np.zeros(shape=(ATOM_COUNT * num_sheets, POS_DIM))
            new_positions[:ATOM_COUNT] = s.pos #copy old data
            s.pos = new_positions
            for i in range(1, num_sheets):
                for j in range(ATOM_COUNT):
                    atom = s.pos[j]
                    for k in range(0, POS_DIM):
                        s.pos[i*ATOM_COUNT + j][k] = atom[k]
                    s.pos[i*ATOM_COUNT + j][2] -= sheet_separation * i

            new_atom_types = np.zeros(shape=(ATOM_COUNT * num_sheets), dtype=type(s.atom_types[0]))
            new_atom_types[:ATOM_COUNT] = s.atom_types
            s.atom_types = new_atom_types
            for i in range(1, num_sheets):
                for j in range(ATOM_COUNT):
                    s.atom_types[i*ATOM_COUNT + j] = s.atom_types[j]

            BOND_COUNT = len(s.bonds)
            new_bonds = np.zeros(shape=(BOND_COUNT * num_sheets, BOND_DIM), dtype='int')
            new_bonds[:BOND_COUNT] = s.bonds
            s.bonds = new_bonds
            for i in range(1, num_sheets):
                for j in range(BOND_COUNT):
                    for k in range(0, BOND_DIM):
                        #Make the bonds point to the new sheet's atoms
                        s.bonds[i*BOND_COUNT + j][k] = s.bonds[j][k] + (ATOM_COUNT * i)

        s.center()        
        return s