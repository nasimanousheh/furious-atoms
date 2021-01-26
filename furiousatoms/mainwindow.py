# Standard package
from json import load
import os
import fnmatch

# Local package
from furiousatoms import io
from furiousatoms.objects import mobius, box_edges
from fury import disable_warnings

disable_warnings()

# 3rd Party package
import vtk
import numpy as np
from fury import window, actor, utils, pick, ui
from fury.utils import numpy_to_vtk_points
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import MDAnalysis
from furiousatoms.sharedmem import SharedMemory
from numpy.linalg import norm
from fractions import gcd
import sys
import furiousatoms.forms.icons
from furiousatoms.periodic_table import Ui_periodic
from furiousatoms.SWNT_builder import  Ui_SWNT
from furiousatoms.graphene_builder import  Ui_graphene
from furiousatoms.MWNT_builder import  Ui_MWNT
from furiousatoms.electrolyte_builder import Ui_electrolyte

SM = SharedMemory()

class FuriousAtomsApp(QtWidgets.QMainWindow):
    """ Main Furious Atoms Class
    """
    def __init__(self, app_path=None, parent=None):
        super(FuriousAtomsApp, self).__init__(parent)
        self.ui = io.load_ui_widget("view.ui")
        self.setCentralWidget(self.ui)
        self.app_path = app_path or io.application_path()
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.qvtkwidget = QVTKRenderWindowInteractor(parent=self.ui.view_frame,
                                                     rw=self.showm.window,
                                                     iren=self.showm.iren)
        self.init_settings()
        self.init_interface()
        self.create_language_menu()
        self.create_connections()

    def init_settings(self):
        pass

    def init_interface(self):
        self.ui.widget_Animation.setVisible(self.ui.button_animation.isChecked())
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.qvtkwidget, stretch=1)
        self.ui.view_frame.setLayout(self.v_layout)
        self.scene.background((1, 1, 1))
        self.showm.initialize()

    def create_connections(self):
        self.ui.actionLoad_file.triggered.connect(self.open)
        self.ui.actionSave_file.triggered.connect(self.save)
        self.ui.actionExit.triggered.connect(self.quit_fired)
        self.ui.button_animation.toggled.connect(self.ui.widget_Animation.setVisible)
        self.ui.actionParticle.triggered.connect(self.delete_particles)
        self.ui.actionBond.triggered.connect(self.delete_bonds)
        self.ui.Button_bondcolor.clicked.connect(self.openColorDialog_bond)
        self.ui.Button_particlecolor.clicked.connect(self.openColorDialog_particle)
        self.ui.SpinBox_atom_radius.valueChanged.connect(self.update_particle_size)
        self.ui.horizontalSlider_metallicity_particle.valueChanged[int].connect(self.metallicity_particle)
        self.ui.horizontalSlider_metallicity_bond.valueChanged[int].connect(self.metallicity_bond)
        self.ui.Button_play.clicked.connect(self.play_movie)
        self.ui.Button_pause.clicked.connect(self.pause_movie)
        self.ui.Button_forward.clicked.connect(self.forward_movie)
        self.ui.Button_backward.clicked.connect(self.backward_movie)
        self.ui.horizontalSlider_animation.sliderMoved.connect(self.slider_changing)
        self.ui.comboBox_particleshape.currentTextChanged.connect(self.change_particle_shape)
        self.ui.comboBox_bondshape.currentTextChanged.connect(self.change_bond_shape)
        self.ui.actionSingle_Wall_Nanotube.triggered.connect(self.single_wall)
        self.ui.actionElectrolyte.triggered.connect(self.electrolyte)
        self.ui.actionGraphene_sheet.triggered.connect(self.graphene)
        self.ui.actionMulti_Wall_nanotube.triggered.connect(self.multiple_walls)
        self.ui.Button_cal_distance.clicked.connect(self.calculate_distance)

    def change_particle_shape(self):
        comboBox_particleshape = self.ui.comboBox_particleshape.currentText()
        if comboBox_particleshape == 'Point':
            SM.sphere_actor.VisibilityOff()
        #     all_vertices_radii = 1/np.repeat(SM.radii_spheres, SM.no_vertices_per_particle, axis=0)
        #     # all_vertices_radii = all_vertices_radii[:, None]
        #     SM.all_vertices_particles[:] = 0.1 * SM.all_vertices_particles[:] #- np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)) + np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)
        if comboBox_particleshape == 'Sphere':
            SM.sphere_actor.VisibilityOn()
        #     all_vertices_radii = 1/np.repeat(SM.radii_spheres, SM.no_vertices_per_particle, axis=0)
        #     # all_vertices_radii = all_vertices_radii[:, None]
        #     SM.all_vertices_particles[:] = 0.55 * SM.all_vertices_particles[:] #- np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)) + np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def change_bond_shape(self):
        comboBox_bondshape = self.ui.comboBox_bondshape.currentText()
        if comboBox_bondshape == 'Line':
            SM.bond_actor.VisibilityOff()
            SM.bond_actor = actor.line(SM.bonds, SM.bond_colors, linewidth=2)
            self.scene.add(SM.bond_actor)
        if comboBox_bondshape == 'Cylinder':
            SM.bond_actor = actor.streamtube(SM.bonds, SM.bond_colors, linewidth=SM.line_thickness)
            self.scene.add(SM.bond_actor)
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        SM.colors_backup_bond = SM.vcolors_bond
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def electrolyte(self):
        Ui_electrolyte.electrolyte = Ui_electrolyte()
        Ui_electrolyte.electrolyte.win = self
        Ui_electrolyte.electrolyte.show()
        Ui_electrolyte.electrolyte.showNormal()

    def single_wall(self):
        Ui_SWNT.swnt = Ui_SWNT()
        Ui_SWNT.swnt.win = self
        Ui_SWNT.swnt.show()
        Ui_SWNT.swnt.showNormal()

    def graphene(self):
        Ui_graphene.gr = Ui_graphene()
        Ui_graphene.gr.win = self
        Ui_graphene.gr.show()
        Ui_graphene.gr.showNormal()

    def multiple_walls(self):
        Ui_MWNT.smnt = Ui_MWNT()
        Ui_MWNT.smnt.win = self
        Ui_MWNT.smnt.show()
        Ui_MWNT.smnt.showNormal()

    def calculate_distance(self):
        object_indices_particles = np.where(SM.selected_particle == True)[0]
        SM.particle_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        if len(object_indices_particles)==2:
            distance_particle_particle = np.linalg.norm(SM.pos[object_indices_particles[0]] - SM.pos[object_indices_particles[1]])
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        distance = "{:.2f}".format(distance_particle_particle)
        avg = np.average(SM.pos[object_indices_particles], axis=0)
        label_actor = actor.label(text=(str(distance)), pos= avg, scale=(0.5, 0.5, 0.5), color=(0, 0, 0))
        first_pos_bond = SM.pos[object_indices_particles[0]]
        second_pos_bond = SM.pos[object_indices_particles[1]]
        bonds = np.hstack((first_pos_bond, second_pos_bond))
        bonds = bonds.reshape(1 , 2, 3)
        bonds += np.array([0, 0, 0.5])
        line_colors = (0, 0, 0, 1)
        line_actor = actor.line(bonds, line_colors, linewidth=5)
        self.scene.add(label_actor)
        line_actor.VisibilityOff()
        self.scene.add(line_actor)
        self.qvtkwidget.GetRenderWindow().Render()
        print('distance of particles', (distance_particle_particle))

    def slider_changing(self):
        SM.cnt = self.ui.horizontalSlider_animation.value()

    def play_movie(self):
        SM.play_factor = 1

    def pause_movie(self):
        SM.play_factor = 0

    def pause_movie(self):
        SM.play_factor = 0

    def forward_movie(self):
        SM.play_factor = -5

    def backward_movie(self):
        SM.play_factor = 5

    def update_particle_size(self, selected_value_radius):
        for i, atom_typ in enumerate(SM.unique_types):
            if self.h_box.itemAt(i).wid.isChecked():
                print(i, atom_typ, 'checked')
                SM.set_value_radius = SM.radii_spheres[SM.atom_type == atom_typ][0]
                # self.ui.SpinBox_atom_radius.setValue((SM.set_value_radius))
                self.ui.SpinBox_atom_radius.setValue(float((selected_value_radius)))
                all_vertices_radii = 1/np.repeat(SM.radii_spheres[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
                all_vertices_radii = all_vertices_radii[:, None]
                # all_vertices_radii = 1/np.repeat(SM.radii_spheres, SM.no_vertices_per_particle, axis=0)
                # all_vertices_radii = all_vertices_radii[:, None]
                # SM.all_vertices_particles[:] = all_vertices_radii * (SM.all_vertices_particles[:] - np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)) + np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
                SM.all_vertices_particles[:] = float(selected_value_radius) * all_vertices_radii * (SM.all_vertices_particles[:] - np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)) + np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def metallicity_particle(self, metallicity_degree_particle):
        SM.metallicCoefficient_particle = metallicity_degree_particle/100
        roughnessCoefficient_particle = 1.0 - metallicity_degree_particle/100
        SM.sphere_actor.GetProperty().SetMetallic(SM.metallicCoefficient_particle)
        SM.sphere_actor.GetProperty().SetRoughness(roughnessCoefficient_particle)
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def metallicity_bond(self, metallicity_degree_bond):
        metallicCoefficient_bond = metallicity_degree_bond/100
        roughnessCoefficient_bond = 1.0 - metallicity_degree_bond/100
        SM.bond_actor.GetProperty().SetMetallic(metallicCoefficient_bond)
        SM.bond_actor.GetProperty().SetRoughness(roughnessCoefficient_bond)
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def check_simulationcell(self, state):
        if (state == QtCore.Qt.Checked):
            SM.line_actor.VisibilityOn()
        else:
            SM.line_actor.VisibilityOff()
        utils.update_actor(SM.line_actor)
        SM.line_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def check_particles(self, state):
        if (state == QtCore.Qt.Checked):
            SM.sphere_actor.VisibilityOn()
        else:
            SM.sphere_actor.VisibilityOff()
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print('All particles are deleted')

    def check_bonds(self, state):
        if (state == QtCore.Qt.Checked):
            SM.bond_actor.VisibilityOn()
        else:
            SM.bond_actor.VisibilityOff()
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print('All bonds are deleted')

    def open(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.tr('Load'))
        if not fname:
            return
        self.process_load_file(fname)
        SM.enable_timer = True

    def save(self):
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, self.tr('Save'))
        # fname = open(fname, 'w')
        # MS.load_file = MDAnalysis.Universe(MS.load_file)
        # SM.universe = u.select_atoms("MS.load_file")
        SM.universe.atoms.write(fname)
        print('saving {}'.format(fname))

    def delete_particles(self):
        object_indices_particles = np.where(SM.selected_particle == True)[0]
        SM.particle_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        for object_index in object_indices_particles:
            SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.particle_color_add
        s = SM.universe.universe.atoms[:]
        t = SM.universe.universe.atoms[object_indices_particles]
        b = s.difference(t)
        SM.universe = MDAnalysis.core.groups.AtomGroup(b)
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print(object_indices_particles)
        # print(SM.pos[object_indices_particles])

    def delete_bonds(self):
        object_indices_bonds = np.where(SM.selected_bond == True)[0]
        SM.bond_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        for object_index_bond in object_indices_bonds:
            SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = SM.bond_color_add

        SM.universe.delete_bonds(SM.universe.bonds[object_indices_bonds])
        # SM.universe.delete_bonds(SM.universe.bonds.to_indices()) #delete if it is lammps
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print(len(SM.universe.bonds))


    def openColorDialog_particle(self):
        selected_color_particle = QtWidgets.QColorDialog.getColor()
        for i, atom_typ in enumerate(SM.unique_types):
            if self.h_box.itemAt(i).wid.isChecked():
                print(i, atom_typ, 'checked')
                object_indices_particles = np.where(SM.atom_type == atom_typ)[0]
                for object_index in object_indices_particles:
                    SM.colors_backup_particles[object_index] = selected_color_particle.getRgb()
                    SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.colors_backup_particles[object_index]
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def openColorDialog_bond(self):
        selected_color_bond = QtWidgets.QColorDialog.getColor()
        select_all_bonds = np.zeros(SM.no_bonds, dtype=np.bool)
        object_indices_bonds = np.where(select_all_bonds == False)[0]
        for object_index in object_indices_bonds:
            SM.colors_backup_bond[object_index] = selected_color_bond.getRgb()
            SM.vcolors_bond[object_index * SM.sec_bond: object_index * SM.sec_bond + SM.sec_bond] = SM.colors_backup_bond[object_index]
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def update_bonds_ui(self):
        if SM.no_atoms > 0:
            self.ui.Box_particles.setChecked(True)
            self.ui.Box_particles.stateChanged.connect(self.check_particles)
            self.ui.Edit_num_of_particles.insert(str(SM.no_atoms))
            self.ui.Edit_num_of_particle_types.insert(str(len(SM.unique_types)))
            self.ui.SpinBox_atom_radius.setValue((SM.set_value_radius))
        if ((SM.box_lx > 0) or (SM.box_ly > 0) or (SM.box_lz > 0)):
            self.ui.Box_simulationcell.setChecked(True)
            self.ui.Box_simulationcell.stateChanged.connect(self.check_simulationcell)
        if SM.no_bonds > 0:
            self.ui.Box_bonds.stateChanged.connect(self.check_bonds)
        if SM.n_frames > 1:
            self.ui.button_animation.setChecked(True)
        self.ui.Edit_num_of_bonds.insert(str(SM.no_bonds))
        self.ui.Edit_widthX.insert(str(SM.box_lx))
        self.ui.Edit_lengthY.insert(str(SM.box_ly))
        self.ui.Edit_hightZ.insert(str(SM.box_lz))
        self.ui.Edit_number_of_frames.insert(str(SM.n_frames))
        self.ui.Edit_directory.insert(str(SM.file_directory))
        self.ui.Edit_fileformat.insert((str(SM.extension)).upper())
        self.ui.Edit_currentfile.insert(str(SM.file_name))
        self.ui.horizontalSlider_animation.setMinimum(0)
        self.ui.horizontalSlider_animation.setMaximum(SM.n_frames)
        self.ui.horizontalSlider_animation.setSingleStep(1)
        self.ui.horizontalSlider_animation.setValue(SM.cnt)
        self.ui.Timer_animation.setMaximum(SM.n_frames)
        for i, atom_typ in enumerate(SM.unique_types):

            if self.h_box.itemAt(i).wid.isChecked():
                print(i, atom_typ, 'checked')


    def process_load_file(self, fname):
        SM.universe, SM.no_bonds = io.load_files(fname) #change it to universe
        SM.file_directory, SM.extension = os.path.splitext(fname)
        SM.file_name = os.path.basename(fname)
        self.process_universe(universe=SM.universe)

    def sky_box_effect(self, actor):
        self.scene.UseImageBasedLightingOn()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        cube_path = os.path.join(dir_path, 'skybox0')
        if not os.path.isdir(cube_path):
            print('This path does not exist:', cube_path)
            return
        cubemap = io.read_cubemap(cube_path, '/', '.jpg', 0)
        self.scene.SetEnvironmentTexture(cubemap)
        actor.GetProperty().SetInterpolationToPBR()
        SM.metallicCoefficient_particle = 0.5
        SM.roughnessCoefficient_particle = 0.1
        colors_sky = vtk.vtkNamedColors()
        actor.GetProperty().SetColor(colors_sky.GetColor3d('White'))
        actor.GetProperty().SetMetallic(SM.metallicCoefficient_particle)
        actor.GetProperty().SetRoughness(SM.roughnessCoefficient_particle)


    def process_universe(self, universe):
        SM.no_atoms = len(SM.universe.atoms)
        SM.n_frames = SM.universe.trajectory.n_frames
        SM.box = SM.universe.trajectory.ts.dimensions
        SM.box_lx = SM.universe.trajectory.ts.dimensions[0]
        SM.box_ly = SM.universe.trajectory.ts.dimensions[1]
        SM.box_lz = SM.universe.trajectory.ts.dimensions[2]
        SM.no_unique_types_particles = len(np.unique(SM.universe.atoms.types))
        box_centers = np.array([[0, 0, 25]])
        box_directions = np.array([[0, 1, 0]])
        SM.box_colors = np.array([[255, 255, 255, 255]])
        box_actor = actor.box(box_centers, box_directions, SM.box_colors,
                              scales=(SM.box_lx, SM.box_ly, SM.box_lz))
        box_actor.GetProperty().SetOpacity(0.)
        box_actor.GetProperty().SetRepresentationToWireframe()
        box_actor.GetProperty().SetLineWidth(10)

        SM.line_actor, _ = box_edges(SM.box_lx, SM.box_ly, SM.box_lz, colors=(0, 0, 0),
                                  linewidth=1, fake_tube=True)

        SM.atom_type = SM.universe.atoms.types
        SM.pos = SM.universe.trajectory[0].positions.copy().astype('f8')
        if SM.no_bonds > 0:
            self.ui.Box_bonds.setChecked(True)
            SM.pos = SM.universe.trajectory[0].positions.copy().astype('f8')
            bonds = SM.universe.bonds.to_indices()
            SM.no_bonds = len(SM.universe.bonds)
            first_pos_bond = SM.pos[(bonds[:, 0])]
            second_pos_bond = SM.pos[(bonds[:, 1])]
            bonds = np.hstack((first_pos_bond, second_pos_bond))
            SM.bonds = bonds.reshape((SM.no_bonds), 2, 3)
            SM.bond_colors = (0.8275, 0.8275, 0.8275, 1)
            SM.line_thickness = 0.2
            SM.bond_actor = actor.streamtube(SM.bonds, SM.bond_colors, linewidth=SM.line_thickness)
            SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
            SM.colors_backup_bond = SM.vcolors_bond.copy()
            SM.all_vertices_bonds = utils.vertices_from_actor(SM.bond_actor)
            SM.no_vertices_per_bond = len(SM.all_vertices_bonds) / SM.no_bonds
            SM.no_vertices_all_bonds = SM.all_vertices_bonds.shape[0]
            SM.sec_bond = np.int(SM.no_vertices_all_bonds / SM.no_bonds)
            SM.unique_types_bond = np.unique(SM.universe.bonds.types)
            SM.bond_actor.AddObserver("LeftButtonPressEvent", self.left_button_press_bond_callback)
        # if SM.no_bonds == 0:
        #     SM.pos_R = SM.universe.trajectory[0].positions.copy().astype('f8')
        #     SM.pos = MDAnalysis.lib.distances.transform_RtoS(SM.pos_R, SM.box,
        #                                                   backend='serial')
        #     open_initial_structur = True
        #     if open_initial_structur is True:
        #         SM.universe_initial_structure = MDAnalysis.Universe('C:/Users/nasim\Devel/furious-atoms/examples/Polymers/PVDF.data')
        #         SM.pos_initial_structure = SM.universe_initial_structure.trajectory[0].positions.copy().astype('f8')
        #         SM.bonds_initial_structure = SM.universe_initial_structure.bonds.to_indices()###
        #         SM.no_bonds = len(SM.universe_initial_structure.bonds)
        #         first_pos_bond = SM.pos_initial_structure[(SM.bonds_initial_structure[:, 0])]###
        #         second_pos_bond = SM.pos_initial_structure[(SM.bonds_initial_structure[:, 1])]###
        #         SM.bonds_initial_structure = np.hstack((first_pos_bond, second_pos_bond))##
        #         SM.bonds = SM.bonds_initial_structure.reshape((SM.no_bonds), 2, 3)
        #         SM.bond_colors = (0.8275, 0.8275, 0.8275, 1)
        #         SM.line_thickness = 0.2
        #         SM.bond_actor = actor.streamtube(SM.bonds, SM.bond_colors, linewidth=SM.line_thickness)
        #         SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        #         SM.colors_backup_bond = SM.vcolors_bond.copy()
        #         SM.all_vertices_bonds = utils.vertices_from_actor(SM.bond_actor)
        #         SM.no_vertices_per_bond = len(SM.all_vertices_bonds) / SM.no_bonds
        #         SM.no_vertices_all_bonds = SM.all_vertices_bonds.shape[0]
        #         SM.sec_bond = np.int(SM.no_vertices_all_bonds / SM.no_bonds)
        #         SM.unique_types_bond = np.unique(SM.universe_initial_structure.bonds.types)

        avg = np.average(SM.pos, axis=0)
        colors = np.ones((SM.no_atoms, 4))
        SM.unique_types = np.unique(SM.universe.atoms.types)
        SM.colors_unique_types = np.random.rand(len(SM.unique_types), 4)
        SM.colors_unique_types[:, 3] = 1
        for i, typ in enumerate(SM.unique_types):
            colors[SM.atom_type == typ] = SM.colors_unique_types[i]

        SM.radii_spheres = np.ones((SM.no_atoms))
        SM.radii_unique_types = 0.55 + np.zeros(len(SM.unique_types)) #np.random.rand(len(SM.unique_types))
        ##############Should be moved to up
        self.toggles = []
        self.lay = QtWidgets.QVBoxLayout()
        self.h_box = QtWidgets.QGridLayout()

        for i, typ in enumerate(SM.unique_types):
            SM.radii_spheres[SM.atom_type == typ] = SM.radii_unique_types[i]
            SM.set_value_radius = SM.radii_spheres[SM.atom_type == typ][0]
            self.ui.btn = QtWidgets.QRadioButton(str(typ) , self)
            self.ui.scrollArea_all_types_of_prticles.setLayout(self.lay)
            self.h_box.addWidget(self.ui.btn,i,1,1,1)
        self.lay.addLayout(self.h_box)

        ########################

        SM.selected_particle = np.zeros(SM.no_atoms, dtype=np.bool)
        SM.selected_bond = np.zeros(SM.no_bonds, dtype=np.bool)
        SM.sphere_actor = actor.sphere(centers=SM.pos,
                                    colors=colors,
                                    radii=SM.radii_spheres, theta=32, phi=32)
        SM.all_vertices_particles = utils.vertices_from_actor(SM.sphere_actor)
        SM.no_vertices_per_particle = len(SM.all_vertices_particles) / SM.no_atoms
        SM.initial_vertices_particles = SM.all_vertices_particles.copy() - np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)
        vertices_particle = utils.vertices_from_actor(SM.sphere_actor)
        SM.no_vertices_all_particles = vertices_particle.shape[0]
        SM.sec_particle = np.int(SM.no_vertices_all_particles / SM.no_atoms)
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        SM.colors_backup_particles = SM.vcolors_particle.copy()

        self.sky_box_effect(SM.sphere_actor)

        axes_actor = actor.axes(scale=(1, 1, 1), colorx=(1, 0, 0),
                                colory=(0, 1, 0), colorz=(0, 0, 1), opacity=1)

        if (SM.no_bonds > 0):
            self.scene.add(SM.bond_actor)

        self.scene.add(axes_actor)
        self.scene.add(SM.sphere_actor)
        self.scene.add(box_actor)
        self.scene.add(SM.line_actor)
        self.scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                              view_up=(0, 1, 0))

        SM.sphere_actor.AddObserver("LeftButtonPressEvent", self.left_button_press_particle_callback)
        self.update_bonds_ui()
        self.timer = QtCore.QTimer()
        duration = 200
        self.timer.start(duration)
        self.timer.timeout.connect(self.timer_callback)
        self.pickm = pick.PickingManager()
        SM.enable_timer = False

    def timer_callback(self):
        if SM.enable_timer is False:
            return
        if SM.cnt == SM.n_frames:
            return
        if SM.n_frames > 1:
            SM.pos_R = SM.universe.trajectory[SM.cnt].positions.copy().astype('f8')
            SM.pos = MDAnalysis.lib.distances.transform_RtoS(SM.pos_R, SM.box, backend='serial')
            SM.all_vertices_particles[:] = SM.initial_vertices_particles + \
                np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)
            utils.update_actor(SM.sphere_actor)
            # open_initial_structure = True
            # if open_initial_structure is True:
            #     SM.bonds_initial_structure = SM.universe_initial_structure.bonds.to_indices()
            #     first_pos_bond = SM.pos[(SM.bonds_initial_structure[:, 0])]
            #     second_pos_bond = SM.pos[(SM.bonds_initial_structure[:, 1])]
            #     SM.bonds_initial_structure = np.hstack((first_pos_bond, second_pos_bond))
            #     SM.bonds = SM.bonds_initial_structure.reshape((SM.no_bonds), 2, 3)
            #     points_array_bonds = np.vstack(SM.bonds)

            #     # Set Points to vtk array format
            #     vtk_points2 = numpy_to_vtk_points(points_array_bonds)
            #     SM.bond_actor.poly_data.SetPoints(vtk_points2)

            self.ui.Timer_animation.setValue(SM.cnt)
            self.ui.Edit_framenumber.setText(str(SM.cnt))
            self.ui.horizontalSlider_animation.setValue(SM.cnt)

        self.qvtkwidget.GetRenderWindow().Render()
        SM.cnt = SM.cnt + 1 * SM.play_factor


    def left_button_press_particle_callback(self, obj, event):

        event_pos = self.pickm.event_position(iren=self.showm.iren)
        picked_info = self.pickm.pick(event_pos, self.showm.scene)

        vertex_index_bond = picked_info['vertex']
        vertices = utils.vertices_from_actor(obj)
        SM.no_vertices_all_particles = vertices.shape[0]
        object_index = np.int(np.floor((vertex_index_bond / SM.no_vertices_all_particles) * SM.no_atoms))

        if not SM.selected_particle[object_index]:
            SM.particle_color_add = np.array([255, 0, 0, 255], dtype='uint8')
            SM.selected_particle[object_index] = True
        else:
            SM.particle_color_add = SM.colors_backup_particles[object_index]
            SM.selected_particle[object_index] = False

        SM.vcolors_particle = utils.colors_from_actor(obj, 'colors')
        SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.particle_color_add
        utils.update_actor(obj)
        obj.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()

    def left_button_press_bond_callback(self, obj, event):
        event_pos = self.pickm.event_position(iren=self.showm.iren)
        picked_info = self.pickm.pick(event_pos, self.showm.scene)
        vertex_index_bond = picked_info['vertex']
        vertices_bonds = utils.vertices_from_actor(obj)
        SM.no_vertices_all_bonds = vertices_bonds.shape[0]
        object_index_bond = np.int(np.floor((vertex_index_bond / SM.no_vertices_all_bonds) * SM.no_bonds))
        if not SM.selected_bond[object_index_bond]:
            SM.bond_color_add = np.array([255, 0, 0, 255], dtype='uint8')
            SM.selected_bond[object_index_bond] = True
        else:
            SM.bond_color_add = SM.colors_backup_bond[object_index_bond]
            SM.selected_bond[object_index_bond] = False

        SM.vcolors_bond = utils.colors_from_actor(obj, 'colors')
        SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = SM.bond_color_add
        utils.update_actor(obj)
        obj.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()

    def quit_fired(self):
        reply = QtWidgets.QMessageBox.warning(
            self,
            self.tr("Quit"),
            self.tr("Do you want to quit the software?"),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            QtCore.QCoreApplication.instance().quit()

    def retranslateUi(self):
        self.protocole_widget.retranslate_ui()
        self.data_selection_widget.retranslate_ui()
        self.processing_widget.retranslate_ui()
        self.hcs_widget.retranslate_ui()

    def closeEvent(self, event):
        reply = QtWidgets.QMessageBox.warning(
            self,
            self.tr("Quit"),
            self.tr("Do you want to quit the software?"),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def changeEvent(self, event):
        if event:
            if event.type() == QtCore.QEvent.LanguageChange:
                # self event is send if a translator is loaded
                # Todo: translation. see http://stackoverflow.com/questions/39734775/python-switch-translations-dynamically-in-a-pyqt4-application
                # self.ui.retranslateUi(self)
                self.retranslateUi()
            if event.type() == QtCore.QEvent.LocaleChange:
                # Todo: Change this part
                locale = QtCore.QLocale.system().name()
                locale.truncate(locale.lastIndexOf('_'))
                self.loadLanguage(locale)

        QtWidgets.QWidget.changeEvent(self, event)

    def create_language_menu(self):
        self.m_curr_lang = ""
        self.m_translator = QtCore.QTranslator(self)

        #  format systems language
        default_locale = QtCore.QLocale.system().name()  # e.g. "de_DE"
        default_locale = default_locale.split('_')[0]  # e.g. "de"

        # Get languages
        self.m_lang_path = io.get_languages_path()
        f_names = fnmatch.filter(os.listdir(self.m_lang_path), "*.qm")

        for f in f_names:
            #  get locale extracted by filename.  fr or fr_FR
            locale = os.path.splitext(f)[0]

            lang = QtCore.QLocale.languageToString(QtCore.QLocale(locale).language())
            ico = QtGui.QIcon("{0}/{1}.png".format(self.m_lang_path, locale))

            action = QtWidgets.QAction(ico, lang, self)
            action.setCheckable(True)
            action.setData(locale)

            self.langGroup.addAction(action)
            self.languageMenu.addAction(action)

            #  set default translators and language checked
            if default_locale == locale:
                action.setChecked(True)

            # Called every time, when a menu entry of the language menu
            # is called

    def slot_language_changed(self, action):
        if not action:
            return

        # load the language dependant on the action content
        self.loadLanguage(action.data())
        self.setWindowIcon(action.icon())

    def switch_translator(self, translator, f_name):
        """ Update the tranlator object

        Parameters
        -----------
        translator : Translator instance
            object which permit to change text
        f_name : str
            language file name
        """
        # Remove the old translator
        QtCore.QCoreApplication.instance().removeTranslator(translator)
        # Load the new translator
        if translator.load(f_name):
            QtCore.QCoreApplication.instance().installTranslator(translator)

    def load_language(self, r_language):
        """ Load the selected language

        Parameters
        -----------
        r_language: str
            language acronym
        """
        if self.m_curr_lang == r_language:
            return

        self.m_curr_lang = r_language
        locale = QtCore.QLocale(self.m_curr_lang)  # de-fr-en
        QtCore.QLocale.setDefault(locale)
        language_name = QtCore.QLocale.languageToString(locale.language())
        self.switch_translator(
            self.m_translator,
            QtCore.QDir(self.m_lang_path).filePath("{0}.qm").format(self.m_curr_lang))
        print(self.tr("Current Language changed to {0}").format(language_name))