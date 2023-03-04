import numpy as np
from furiousatoms import io
import numpy as np
from fury import window, utils
from PySide2 import QtWidgets
from furiousatoms.structure import bbox
from furiousatoms.molecular import MolecularStructure
from PySide2.QtGui import QIcon

"""
    Ui_solution class creates a widget for building box and solution (e.g. water)
"""
class Ui_solution(QtWidgets.QMainWindow): #QWidget

    def __init__(self, app_path=None, parent=None):
        super(Ui_solution, self).__init__(parent)
        self.solution = io.load_ui_widget("solution.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.solution)
        self.setCentralWidget(self.solution)
        self.setLayout(self.v_layout)
        self.resize(280, 202)
        self.scene = window.Scene()
        self.setWindowIcon(QIcon(io.get_resources_file("splash.png")))
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.solution.SpinBox_lx.valueChanged.connect(self.box_builder_callback)
        self.solution.SpinBox_ly.valueChanged.connect(self.box_builder_callback)
        self.solution.SpinBox_lz.valueChanged.connect(self.box_builder_callback)
        self.solution.pushButton_build_solution.clicked.connect(self.solution_builder_callback)
        self.solution.pushButton_build_solution.clicked.connect(lambda:self.close())
        water_diameter = 3.1655
        self.solution.lineEdit_water_dia.setText(str(water_diameter))

    def initial_box_dim(self, box_lx, box_ly, box_lz):
        self.solution.SpinBox_lx.setValue(box_lx)
        self.solution.SpinBox_lz.setValue(box_lz)
        self.solution.SpinBox_ly.setValue(box_ly)

    def initial_values(self):
        water_diameter = 3.1655
        avoid_overlap = 2
        spacing_dia = 0.98
        spacing = spacing_dia * water_diameter
        space_format = "{:.4f}".format(spacing)
        self.solution.lineEdit_water_dia.setText(str(water_diameter))
        self.solution.SpinBox_avoid_overlap.setValue((avoid_overlap))

    def box_builder_callback(self):
        active_window = self.win.active_mdi_child()
        SM = active_window.universe_manager
        box_lx = float(self.solution.SpinBox_lx.text())
        box_ly = float(self.solution.SpinBox_ly.text())
        box_lz = float(self.solution.SpinBox_lz.text())
        SM.box_size = [box_lx, box_ly, box_lz]

    def solution_builder_callback(self):
        active_window = self.win.active_mdi_child()
        box_lx = float(self.solution.SpinBox_lx.text())
        box_ly = float(self.solution.SpinBox_ly.text())
        box_lz = float(self.solution.SpinBox_lz.text())
        if not (box_lx or box_ly or box_lz):
            #TODO maybe display a popup error for no water area?
            return
        SM = active_window.universe_manager
        box_size = [box_lx, box_ly, box_lz]
        water_diameter = 3.1655
        spacing_dia = 0.98
        try:
            avoid_overlap = float(self.solution.SpinBox_avoid_overlap.text())
        except:
            avoid_overlap = 2

        spacing = spacing_dia * water_diameter
        total_water_inside = int((int(box_lx/spacing) * int(box_ly/spacing) * int(box_lz/spacing)))
        water_amount = (box_lx/spacing * box_ly/spacing * box_lz/spacing)
        water_concentration = water_amount / (0.602214076 * (box_lx) * (box_ly) * (box_lz) * 0.001)
        print("water_concentration is: ", water_concentration)

        h2o = np.array([[ 0,        0,       0      ],  # oxygen
                        [ 0.95908, -0.02691, 0.03231],  # hydrogen
                        [-0.28004, -0.58767, 0.70556]]) # hydrogen
        coordinates = []
        num_molecul_in_box_lx = int(box_lx / spacing)
        num_molecul_in_box_ly = int(box_ly / spacing)
        num_molecul_in_box_lz = int(box_lz / spacing)
        mol = 0
        for i in range(num_molecul_in_box_lx):
            for j in range(num_molecul_in_box_ly):
                for k in range(num_molecul_in_box_lz):
                    if (mol < (total_water_inside*3)):
                        x = (-box_lx/2 + 0.5 * (spacing)) + (i * spacing)
                        y = (-box_ly/2 + 0.5 * (spacing)) + (j * spacing)
                        z = (-box_lz/2 + 0.5 * (spacing)) + (k * spacing)
                        if (x > (box_lx/2 - (0.5 * spacing)) or y > (box_ly/2 - (0.5 * spacing)) or z > (box_lz/2 - (0.5 * spacing))):
                           continue
                        xyz = np.array([x, y, z])
                        close = False
                        for p_n in range(len(SM.pos)):
                            dist = np.linalg.norm(xyz - SM.pos[p_n])
                            if dist < avoid_overlap:
                                close = True
                                break
                        if not close:
                            coordinates.extend(h2o + xyz.T)
                    mol = mol + 1

        coord_array = np.array(coordinates)
        n_atoms = len(coord_array)
        assert coord_array.shape == (n_atoms, 3)
        n_residues = int(n_atoms/3)
        resindices = np.repeat(range(n_residues), 3)
        assert len(resindices) == n_atoms
        box_size = [box_lx, box_ly, box_lz]
        atom_types = ['O', 'H', 'H']*n_residues
        atom_types = np.array(atom_types)

        bonds = []
        for o in range(0, n_atoms, 3):
            bonds.extend([(o, o+1), (o, o+2)])

        sol = MolecularStructure(box_size, coord_array, np.array(bonds, dtype=int), atom_types)
        sol.center()

        original = MolecularStructure(SM.box_size, SM.pos, SM.bonds, SM.atom_types)
        combined = original.merge(sol)
        combined.box_size = box_size

        import tempfile
        dir_name = tempfile.mkdtemp(prefix='Furious_Atoms_')
        file_name = tempfile.mkstemp(suffix='.pdb', prefix='Solution', dir=dir_name)[1]
        # combined.atoms.write(file_name)
        #TODO save universe after rwater added

        SM = active_window.universe_manager
        active_window.scene.rm(SM.bbox_actor)
        SM.bbox_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.scene.rm(SM.sphere_actor)
        active_window.scene.rm(SM.bond_actor)
        active_window.load_structure(combined)
        active_window.render()
        return combined