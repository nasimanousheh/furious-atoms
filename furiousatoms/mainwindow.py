# Standard package
from json import load
import os
import fnmatch

# Local package
from furiousatoms import io
from fury import disable_warnings

disable_warnings()

# 3rd Party package
import vtk
import numpy as np
from fury import window, actor, utils, pick, ui, primitive
from fury.utils import numpy_to_vtk_points
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
import MDAnalysis
from numpy.linalg import norm
from fractions import gcd
import sys
import furiousatoms.forms.icons
from furiousatoms.viewer3d import Viewer3D
from furiousatoms.periodic_table import Ui_periodic
from furiousatoms.SWNT_builder import  Ui_SWNT
from furiousatoms.graphene_builder import  Ui_graphene
from furiousatoms.MWNT_builder import  Ui_MWNT
from furiousatoms.electrolyte_builder import Ui_electrolyte
from furiousatoms.fullerenes_database import load_CC1_file
from furiousatoms.viewer3d import Viewer3D


class FuriousAtomsApp(QtWidgets.QMainWindow):
    """ Main Furious Atoms Class
    """
    def __init__(self, app_path=None, parent=None):
        super(FuriousAtomsApp, self).__init__(parent)
        self.ui = io.load_ui_widget("view.ui")
        self.setCentralWidget(self.ui)
        self.app_path = app_path or io.application_path()

        self.init_settings()
        self.init_interface()
        self.create_language_menu()
        self.create_connections()

    def init_settings(self):
        pass

    def init_interface(self):
        # Hide by default animation
        self.ui.widget_Animation.setVisible(self.ui.button_animation.isChecked())
        # self.ui.mdiArea.subWindowActivated.connect(self.updateMenus)
        # self.windowMapper = QtCore.QSignalMapper(self)
        # self.windowMapper.mapped.connect(self.setActiveSubWindow)

    def create_connections(self):
        # File menu actions
        # self.ui.actionNew_file.triggered.connect(self.new_structure)
        self.ui.actionLoad_file.triggered.connect(self.open)
        self.ui.actionSave_file.triggered.connect(self.save)
        self.ui.actionExit.triggered.connect(self.quit_fired)

        self.ui.actionZoom_in.triggered.connect(self.slotZoomIn)
        self.ui.actionZoom_out.triggered.connect(self.slotZoomOut)
        self.ui.actionFit_to_Window.triggered.connect(self.slotZoomFit)
        # self.ui.button_animation.toggled.connect(self.ui.widget_Animation.setVisible)

        # View menu actions
        self.ui.actioncascade.triggered.connect(self.ui.mdiArea.cascadeSubWindows)
        self.ui.actiontiled.triggered.connect(self.ui.mdiArea.tileSubWindows)

        # Help menu actions

        # Modify menu actions
        self.ui.actionParticles.triggered.connect(self.delete_particles)
        self.ui.actionBonds.triggered.connect(self.delete_bonds)

        # Compute menu actions

        # Build menu actions
        self.ui.actionGraphene_sheet.triggered.connect(self.graphene)
        self.ui.actionSingle_Wall_Nanotube.triggered.connect(self.single_wall)
        self.ui.actionMulti_Wall_nanotube.triggered.connect(self.multiple_walls)
        self.ui.actionElectrolyte.triggered.connect(self.electrolyte)
        self.ui.actionFullerenes.triggered.connect(self.open_dataset_fullerene)
        self.ui.button_animation.toggled.connect(self.ui.widget_Animation.setVisible)
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
        self.ui.Button_cal_distance.clicked.connect(self.calculate_distance)
        self.ui.comboBox_particle_resolution.currentTextChanged.connect(self.change_particle_resolution)

        # General connections
        self.ui.mdiArea.subWindowActivated.connect(self.update_bonds_ui)

    def slotZoomIn(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        active_window.scene.zoom(1.25)
        active_window.render()

    def slotZoomOut(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        active_window.scene.zoom(0.75)
        active_window.render()
    def slotZoomFit(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        active_window.scene.ResetCamera()
        active_window.render()

    def change_particle_resolution(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        num = int(SM.no_vertices_per_particle)
        print('SM.no_vertices_per_particle', num)
        comboBox_particle_resolution = self.ui.comboBox_particle_resolution.currentText()

        colors = SM.colors_backup_particles[0::num].astype('f8').copy()/255
        active_window.scene.rm(SM.sphere_actor)
        vertices, faces = primitive.prim_sphere(name=comboBox_particle_resolution, gen_faces=False)
        res = primitive.repeat_primitive(vertices, faces, centers=SM.pos, colors=colors)#, dtype='uint8')
        big_verts, big_faces, big_colors, _ = res
        SM.sphere_actor = utils.get_actor_from_primitive(big_verts, big_faces, big_colors)
        SM.no_vertices_per_particle = len(vertices)
        SM.all_vertices_particles = SM.no_vertices_per_particle * SM.no_atoms

        print('colors.shape:' , colors.shape)

        # for i, atom_typ in enumerate(SM.unique_types):
        #     if self.ui.scrollArea_all_types_of_prticles.layout().itemAt(i).wid.isChecked():
        #         SM.set_value_radius = SM.radii_spheres[SM.atom_type == atom_typ][0]
        #         selected_value_radius = SM.radii_spheres[SM.atom_type == atom_typ]
        #         all_vertices_radii = 1/np.repeat(SM.radii_spheres[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
        #         all_vertices_radii = all_vertices_radii[:, None]
        #         selected_atom_mask = SM.atom_type == atom_typ
        #         all_vertices_mask = np.repeat(selected_atom_mask, SM.no_vertices_per_particle)
        #         SM.all_vertices_particles[all_vertices_mask] = float(selected_value_radius) * all_vertices_radii * (SM.all_vertices_particles[all_vertices_mask] - np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)) + np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
        # # SM.all_vertices_particles[:] = SM.initial_vertices_particles + \
        #         np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)




        # self.update_particle_size(selected_value_radius=0.2)
        active_window.scene.add(SM.sphere_actor)
        # SM.all_vertices_particles = utils.vertices_from_actor(SM.sphere_actor)
        # SM.no_vertices_per_particle = len(SM.all_vertices_particles) / SM.no_atoms
        # SM.initial_vertices_particles = SM.all_vertices_particles - np.repeat(SM.pos, SM.no_vertices_per_particle, axis=0)
        # vertices_particle = utils.vertices_from_actor(SM.sphere_actor)
        # SM.no_vertices_all_particles = vertices_particle.shape[0]
        # SM.sec_particle = np.int(SM.no_vertices_all_particles / SM.no_atoms)
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        SM.colors_backup_particles = SM.vcolors_particle.copy()



        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        utils.update_actor(SM.sphere_actor)
        # Viewer3D.display_universe(SM.sphere_actor)
        active_window.render()

    def change_particle_shape(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        comboBox_particleshape = self.ui.comboBox_particleshape.currentText()
        if comboBox_particleshape == 'Point':
            SM.sphere_actor.VisibilityOff()
        if comboBox_particleshape == 'Sphere':
            SM.sphere_actor.VisibilityOn()
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        active_window.render()

    def change_bond_shape(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        comboBox_bondshape = self.ui.comboBox_bondshape.currentText()
        if comboBox_bondshape == 'Line':
            SM.bond_actor.VisibilityOff()
            SM.bond_actor = actor.line(SM.bonds, SM.bond_colors, linewidth=2)
            active_window.scene.add(SM.bond_actor)
        if comboBox_bondshape == 'Cylinder':
            SM.bond_actor = actor.streamtube(SM.bonds, SM.bond_colors,
                                             linewidth=SM.line_thickness)
            active_window.scene.add(SM.bond_actor)
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        # SM.colors_backup_bond = SM.vcolors_bond
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        active_window.render()

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
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        object_indices_particles = np.where(SM.selected_particle == True)[0]
        SM.particle_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        SM.colors_backup_particles = SM.vcolors_particle.copy()
        if len(object_indices_particles)==2:
            distance_particle_particle = np.linalg.norm(SM.pos[object_indices_particles[0]] - SM.pos[object_indices_particles[1]])
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        distance = "{:.2f}".format(distance_particle_particle)
        avg = np.average(SM.pos[object_indices_particles], axis=0)
        label_actor = actor.label(text=(str(distance)), pos=avg,
                                  scale=(0.5, 0.5, 0.5), color=(0, 0, 0))
        first_pos_bond = SM.pos[object_indices_particles[0]]
        second_pos_bond = SM.pos[object_indices_particles[1]]
        bonds = np.hstack((first_pos_bond, second_pos_bond))
        bonds = bonds.reshape(1, 2, 3)
        bonds += np.array([0, 0, 0.5])
        line_colors = (0, 0, 0, 1)
        line_actor = actor.line(bonds, line_colors, linewidth=5)
        active_window.scene.add(label_actor)
        line_actor.VisibilityOff()
        active_window.scene.add(line_actor)
        active_window.render()
        print('distance of particles', (distance_particle_particle))

    def slider_changing(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.cnt = self.ui.horizontalSlider_animation.value()

    def play_movie(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.play_factor = 1

    def pause_movie(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.play_factor = 0

    def pause_movie(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.play_factor = 0

    def forward_movie(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.play_factor = -5

    def backward_movie(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.play_factor = 5

    def update_particle_size(self, selected_value_radius):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        for i, atom_typ in enumerate(SM.unique_types):
            if self.ui.scrollArea_all_types_of_prticles.layout().itemAt(i).wid.isChecked():
                print(i, atom_typ, 'checked')
                SM.set_value_radius = SM.radii_spheres[SM.atom_type == atom_typ][0]
                self.ui.SpinBox_atom_radius.setValue(float((selected_value_radius)))
                all_vertices_radii = 1/np.repeat(SM.radii_spheres[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
                all_vertices_radii = all_vertices_radii[:, None]
                selected_atom_mask = SM.atom_type == atom_typ
                all_vertices_mask = np.repeat(selected_atom_mask, SM.no_vertices_per_particle)
                SM.all_vertices_particles[all_vertices_mask] = float(selected_value_radius) * all_vertices_radii * (SM.all_vertices_particles[all_vertices_mask] - np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)) + np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
                SM.radii_spheres[SM.atom_type == atom_typ] = float(selected_value_radius)
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        active_window.render()

    def metallicity_particle(self, metallicity_degree_particle):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.metallicCoefficient_particle = metallicity_degree_particle/100
        roughnessCoefficient_particle = 1.0 - metallicity_degree_particle/100
        SM.sphere_actor.GetProperty().SetMetallic(SM.metallicCoefficient_particle)
        SM.sphere_actor.GetProperty().SetRoughness(roughnessCoefficient_particle)
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        active_window.render()

    def metallicity_bond(self, metallicity_degree_bond):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        metallicCoefficient_bond = metallicity_degree_bond/100
        roughnessCoefficient_bond = 1.0 - metallicity_degree_bond/100
        SM.bond_actor.GetProperty().SetMetallic(metallicCoefficient_bond)
        SM.bond_actor.GetProperty().SetRoughness(roughnessCoefficient_bond)
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        active_window.render()

    def check_simulationcell(self, state):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        if (state == QtCore.Qt.Checked):
            SM.bbox_actor.VisibilityOn()
        else:
            SM.bbox_actor.VisibilityOff()
        utils.update_actor(SM.bbox_actor)
        SM.bbox_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()

    def check_particles(self, state):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        if (state == QtCore.Qt.Checked):
            SM.sphere_actor.VisibilityOn()
        else:
            SM.sphere_actor.VisibilityOff()
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()
        print('All particles are deleted')

    def check_bonds(self, state):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        if SM.bond_actor is None:
            active_window.render()
            return

        if (state == QtCore.Qt.Checked):
            SM.bond_actor.VisibilityOn()
        else:
            SM.bond_actor.VisibilityOff()
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()
        print('All bonds are deleted')

    def new_window(self):
        child = self.create_mdi_child()
        child.make_title()
        child.show()

    def open(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.tr('Load'))
        if not fname:
            return

        existing = self.find_mdi_child(fname)
        if existing:
            self.ui.mdiArea.setActiveSubWindow(existing)
            return

        child = self.create_mdi_child()
        if child.load_file(fname):
            # self.statusBar().showMessage("File loaded", 2000)
            child.show()
        else:
            child.close()

    def open_dataset_fullerene(self):
        dir_fullerene_folder = os.path.dirname(os.path.realpath(__file__))
        fullerene_folder = os.path.join(dir_fullerene_folder,
                                        'fullerene_dataset')
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.tr('Load'),
                                                         fullerene_folder,
                                                         filter="*.cc1*")
        existing = self.find_mdi_child(fname)
        if existing:
            self.ui.mdiArea.setActiveSubWindow(existing)
            return

        child = self.create_mdi_child()

        # if child.load_CC1_file(fname):

        if child.load_fullerene_cc1_file(fname):
            child.show()
        else:
            child.close()

        # SM.enable_timer = True
    def save(self):
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, self.tr('Save'))
        # fname = open(fname, 'w')
        # MS.load_file = MDAnalysis.Universe(MS.load_file)
        # SM.universe = u.select_atoms("MS.load_file")
        if not fname:
            return

        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.universe.atoms.write(fname)
        # SM.atoms.write(fname)
        print('saving {}'.format(fname))

    def delete_particles(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        active_window.delete_particles()

    def delete_bonds(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        active_window.delete_bonds()

    def openColorDialog_particle(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        selected_color_particle = QtWidgets.QColorDialog.getColor()
        # if color = QtWidgets.QColorDialog.getColor(self.selectedColor, self)
        if selected_color_particle.isValid():
            for i, atom_typ in enumerate(SM.unique_types):
                if self.ui.scrollArea_all_types_of_prticles.layout().itemAt(i).wid.isChecked():
                    print(i, atom_typ, 'checked')
                    object_indices_particles = np.where(SM.atom_type == atom_typ)[0]
                    for object_index in object_indices_particles:
                        # SM.colors_backup_particles[object_index] = selected_color_particle.getRgb()
                        # SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.colors_backup_particles[object_index]
                        SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = selected_color_particle.getRgb()
                        test = 1
                        # [ 56,  45, 232, 255]
            SM.colors_backup_particles = SM.vcolors_particle.copy()
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()

    def openColorDialog_bond(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        selected_color_bond = QtWidgets.QColorDialog.getColor()
        select_all_bonds = np.zeros(SM.no_bonds, dtype=np.bool)
        object_indices_bonds = np.where(select_all_bonds == False)[0]
        for object_index in object_indices_bonds:
            SM.colors_backup_bond[object_index] = selected_color_bond.getRgb()
            SM.vcolors_bond[object_index * SM.sec_bond: object_index * SM.sec_bond + SM.sec_bond] = SM.colors_backup_bond[object_index]
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()

    def update_bonds_ui(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager

        scroll_layout = QtWidgets.QGridLayout()

        # reset the layout
        for i, typ in enumerate(SM.unique_types):

            SM.radii_spheres[SM.atom_type == typ] = SM.radii_unique_types[i]
            SM.set_value_radius = SM.radii_spheres[SM.atom_type == typ][0]
            btn = QtWidgets.QRadioButton(str(typ), self)
            scroll_layout.addWidget(btn, i, 1, 1, 1)

        self.ui.scrollArea_all_types_of_prticles.setLayout(scroll_layout)
        # Disconnect signal
        try:
            if SM.no_atoms > 0:
                self.ui.Box_particles.stateChanged.disconnect(self.check_particles)
        except RuntimeError:
            pass
        try:
            if SM.box_lx > 0 or SM.box_ly > 0 or SM.box_lz > 0:
                self.ui.Box_simulationcell.stateChanged.disconnect(self.check_simulationcell)
        except RuntimeError:
            pass
        try:
            if SM.no_bonds > 0:
                self.ui.Box_bonds.stateChanged.disconnect(self.check_bonds)
        except RuntimeError:
            pass

        # Setup checkbox
        self.ui.Box_particles.setChecked(SM.no_atoms > 0)
        self.ui.Box_simulationcell.setChecked(SM.box_lx > 0 or
                                              SM.box_ly > 0 or
                                              SM.box_lz > 0)
        self.ui.Box_bonds.setChecked(SM.no_bonds > 0)
        self.ui.button_animation.setChecked(SM.n_frames > 1)

        if SM.no_atoms > 0:
            self.ui.Box_particles.stateChanged.connect(self.check_particles)
            self.ui.Edit_num_of_particles.setText(str(SM.no_atoms))
            self.ui.Edit_num_of_particle_types.setText(str(len(SM.unique_types)))
            self.ui.SpinBox_atom_radius.setValue((SM.set_value_radius))
        if SM.box_lx > 0 or SM.box_ly > 0 or SM.box_lz > 0:
            self.ui.Box_simulationcell.stateChanged.connect(self.check_simulationcell)
        if SM.no_bonds > 0:
            self.ui.Box_bonds.stateChanged.connect(self.check_bonds)

        self.ui.Edit_num_of_bonds.setText(str(SM.no_bonds))
        self.ui.Edit_widthX.setText(str(SM.box_lx))
        self.ui.Edit_lengthY.setText(str(SM.box_ly))
        self.ui.Edit_hightZ.setText(str(SM.box_lz))
        self.ui.Edit_number_of_frames.setText(str(SM.n_frames))
        self.ui.Edit_directory.setText(str(active_window.current_filedir))
        self.ui.Edit_fileformat.setText((str(active_window.current_extension)).upper())
        self.ui.Edit_currentfile.setText(str(active_window.current_file))
        self.ui.horizontalSlider_animation.setMinimum(0)
        self.ui.horizontalSlider_animation.setMaximum(SM.n_frames)
        self.ui.horizontalSlider_animation.setSingleStep(1)
        self.ui.horizontalSlider_animation.setValue(SM.cnt)
        self.ui.Timer_animation.setMaximum(SM.n_frames)
        self.ui.tabWidget.update()

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

    def active_mdi_child(self):
        active_sub_window = self.ui.mdiArea.activeSubWindow()
        if active_sub_window:
            return active_sub_window.widget()
        return None

    def create_mdi_child(self, mdi_type="viewer"):
        if mdi_type == "viewer":
            child = Viewer3D()
        else:
            raise ValueError("Unknown MDI TYPE")

        self.ui.mdiArea.addSubWindow(child)
        # Todo: Make some connections
        return child

    def find_mdi_child(self, fname):
        canonical_filepath = os.path.abspath(fname)

        for window in self.ui.mdiArea.subWindowList():
            if window.widget().current_filepath == canonical_filepath:
                return window
        return None

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
            self.ui.mdiArea.closeAllSubWindows()
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