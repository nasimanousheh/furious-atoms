# Standard package
from json import load
import os
import fnmatch

# Local package
from furiousatoms import io
from fury import disable_warnings

disable_warnings()

# 3rd Party package
import numpy as np
from fury import window, actor, utils, pick, ui, primitive, material
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QActionEvent, QIcon
from PySide2 import QtWidgets
from PySide2.QtCore import QTimer
import MDAnalysis
from numpy.linalg import norm
import sys
import furiousatoms.forms.icons
from fury.io import save_image
from fury.lib import (RenderLargeImage, numpy_support)
from furiousatoms.viewer3d import Viewer3D, sky_box_effect_atom, sky_box_effect_bond
from furiousatoms.SWNT_builder import  Ui_SWNT
from furiousatoms.graphene_builder import  Ui_graphene
from furiousatoms.box_builder import  Ui_box
from furiousatoms.solution_builder import  Ui_solution
from furiousatoms.MWNT_builder import  Ui_MWNT
from furiousatoms.electrolyte_builder import Ui_electrolyte
from furiousatoms.fullerenes_builder import load_CC1_file
from fury.utils import (get_actor_from_primitive, normals_from_actor,
                        tangents_to_actor, update_polydata_normals,
                        tangents_from_direction_of_anisotropy)


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
        self.ui.actionNew_file.triggered.connect(self.new_window)
        self.ui.actionLoad_file.triggered.connect(self.open)
        self.ui.actionSave_file.triggered.connect(self.save)
        self.ui.actionSave_Image_File.triggered.connect(self.save_image)
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
        self.ui.actionAdd_Box.triggered.connect(self.box_builder)
        self.ui.actionAdd_solution.triggered.connect(self.solution_builder)
        self.ui.button_animation.toggled.connect(self.ui.widget_Animation.setVisible)
        self.ui.Button_bondcolor.clicked.connect(self.openColorDialog_bond)
        self.ui.Button_particlecolor.clicked.connect(self.openColorDialog_particle)
        self.ui.SpinBox_atom_radius.valueChanged.connect(self.update_particle_size)
        self.ui.Button_play.clicked.connect(self.play_movie)
        self.ui.Button_pause.clicked.connect(self.pause_movie)
        self.ui.Button_forward.clicked.connect(self.forward_movie)
        self.ui.Button_backward.clicked.connect(self.backward_movie)
        self.ui.horizontalSlider_animation.sliderMoved.connect(self.slider_changing)
        self.ui.comboBox_particleshape.currentTextChanged.connect(self.change_particle_shape)
        self.ui.comboBox_bondshape.currentTextChanged.connect(self.change_bond_shape)
        self.ui.Button_cal_distance.clicked.connect(self.calculate_distance)
        self.ui.comboBox_particle_resolution.currentTextChanged.connect(self.change_particle_resolution)
        self.ui.horizontalSlider_Opacity.valueChanged[int].connect(self.change_slice_opacity)
        self.ui.horizontalSlider_Metallic.valueChanged[int].connect(self.change_slice_metallic)
        self.ui.horizontalSlider_Roughness.valueChanged[int].connect(self.change_slice_roughness)
        self.ui.horizontalSlider_Anisotropic.valueChanged[int].connect(self.change_slice_anisotropic)
        self.ui.horizontalSlider_Anisotropic_rot.valueChanged[int].connect(self.change_slice_anisotropic_rot)
        self.ui.horizontalSlider_Anisotropic_X.valueChanged[int].connect(self.change_slice_anisotropic_X)
        self.ui.horizontalSlider_Anisotropic_Y.valueChanged[int].connect(self.change_slice_anisotropic_Y)
        self.ui.horizontalSlider_Anisotropic_Z.valueChanged[int].connect(self.change_slice_anisotropic_Z)
        self.ui.horizontalSlider_Coat_strength.valueChanged[int].connect(self.change_slice_Coat_strength)
        self.ui.horizontalSlider_Coat_roughness.valueChanged[int].connect(self.change_slice_Coat_rough)

        self.ui.treeWidget.setHeaderLabels(['color', 'Particle'])
        self.ui.treeWidget.itemClicked.connect(self.show_radius_value)
        # General connections
        self.ui.mdiArea.subWindowActivated.connect(self.update_information_ui)

    def show_radius_value(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        selected_item = self.ui.treeWidget.selectionModel()
        for i, atom_typ in enumerate(SM.unique_types):
            if selected_item.rowIntersectsSelection(i):
                SM.radii_spheres[SM.atom_type == atom_typ] = SM.radii_unique_types[i]
                set_value_radius = SM.radii_spheres[SM.atom_type == atom_typ][0]
                self.ui.SpinBox_atom_radius.setValue((set_value_radius))


    def change_slice_opacity(self, opacity_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.opacity = opacity_degree/100
        SM.sphere_actor.GetProperty().SetInterpolationToPBR()
        SM.sphere_actor.GetProperty().SetOpacity(SM.opacity)
        utils.update_actor(SM.sphere_actor)
        if SM.no_bonds > 0:
            SM.bond_actor.GetProperty().SetInterpolationToPBR()
            SM.bond_actor.GetProperty().SetOpacity(SM.opacity)
            utils.update_actor(SM.bond_actor)
        active_window.render()

    def change_slice_metallic(self, metallic_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.metallic = metallic_degree/100
        SM.pbr_params_atom.metallic = SM.metallic
        if SM.no_bonds > 0:
            SM.pbr_params_bond.metallic = SM.metallic
        active_window.render()

    def change_slice_roughness(self, roughness_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.roughness = roughness_degree/100
        SM.pbr_params_atom.roughness = SM.roughness
        if SM.no_bonds > 0:
            SM.pbr_params_bond.roughness = SM.roughness
        active_window.render()

    def change_slice_anisotropic(self, anisotropic_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.anisotropic = anisotropic_degree/100
        SM.pbr_params_atom.anisotropy = SM.anisotropic
        if SM.no_bonds > 0:
            SM.pbr_params_bond.anisotropy = SM.anisotropic
        active_window.render()

    def change_slice_anisotropic_rot(self, anisotropic_rot_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.anisotropic_rot = anisotropic_rot_degree/100
        SM.pbr_params_atom.anisotropy_rotation = SM.anisotropic_rot
        if SM.no_bonds > 0:
            SM.pbr_params_bond.anisotropy_rotation = SM.anisotropic_rot
        active_window.render()

    def change_slice_anisotropic_X(self, anisotropic_X_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.anisotropic_X = anisotropic_X_degree/100
        doa = [0, 1, .5]
        doa[0] = SM.anisotropic_X
        polydata = SM.sphere_actor.GetMapper().GetInput()
        update_polydata_normals(polydata)
        normals = normals_from_actor(SM.sphere_actor)
        tangents = tangents_from_direction_of_anisotropy(normals, doa)
        tangents_to_actor(SM.sphere_actor, tangents)
        if SM.no_bonds > 0:
            normals = normals_from_actor(SM.bond_actor)
            tangents = tangents_from_direction_of_anisotropy(normals, doa)
            tangents_to_actor(SM.bond_actor, tangents)
        active_window.render()

    def change_slice_anisotropic_Y(self, anisotropic_Y_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.anisotropic_Y = anisotropic_Y_degree/100
        doa = [0, 1, .5]
        doa[1] = SM.anisotropic_Y
        polydata = SM.sphere_actor.GetMapper().GetInput()
        update_polydata_normals(polydata)
        normals = normals_from_actor(SM.sphere_actor)
        tangents = tangents_from_direction_of_anisotropy(normals, doa)
        tangents_to_actor(SM.sphere_actor, tangents)
        if SM.no_bonds > 0:
            normals = normals_from_actor(SM.bond_actor)
            tangents = tangents_from_direction_of_anisotropy(normals, doa)
            tangents_to_actor(SM.bond_actor, tangents)
        active_window.render()

    def change_slice_anisotropic_Z(self, anisotropic_Z_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.anisotropic_Z = anisotropic_Z_degree/100
        doa = [0, 1, .5]
        doa[2] = SM.anisotropic_Z
        polydata = SM.sphere_actor.GetMapper().GetInput()
        update_polydata_normals(polydata)
        normals = normals_from_actor(SM.sphere_actor)
        tangents = tangents_from_direction_of_anisotropy(normals, doa)
        tangents_to_actor(SM.sphere_actor, tangents)
        if SM.no_bonds > 0:
            normals = normals_from_actor(SM.bond_actor)
            tangents = tangents_from_direction_of_anisotropy(normals, doa)
            tangents_to_actor(SM.bond_actor, tangents)
        active_window.render()

    def change_slice_clearcoat(self, clearcoat_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.clearcoat = clearcoat_degree/100
        utils.update_actor(SM.sphere_actor)
        utils.update_actor(SM.bond_actor)
        active_window.render()

    def change_slice_clearcoat_gloss(self, clearcoat_gloss_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.clearcoat_gloss = clearcoat_gloss_degree/100
        utils.update_actor(SM.sphere_actor)
        utils.update_actor(SM.bond_actor)
        active_window.render()
    def change_slice_sheen(self, sheen_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.sheen = sheen_degree/100
        utils.update_actor(SM.sphere_actor)
        utils.update_actor(SM.bond_actor)
        active_window.render()
    def change_slice_sheen_tint(self, sheen_tint_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.sheen_tint = sheen_tint_degree/100
        utils.update_actor(SM.sphere_actor)
        utils.update_actor(SM.bond_actor)
        active_window.render()

    def change_slice_Coat_rough(self, coat_rough_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.coat_rough = coat_rough_degree/100
        SM.pbr_params_atom.coat_roughness = SM.coat_rough
        if SM.no_bonds > 0:
            SM.pbr_params_bond.coat_roughness = SM.coat_rough
        active_window.render()
    def change_slice_Coat_strength(self, coat_strength_degree):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.coat_strength = coat_strength_degree/100
        SM.pbr_params_atom.coat_strength = SM.coat_strength
        if SM.no_bonds > 0:
            SM.pbr_params_bond.coat_strength = SM.coat_strength
        active_window.render()

    def slotZoomIn(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        active_window.scene.zoom(1.25)
        active_window.render()

    def slotZoomOut(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        active_window.scene.zoom(0.75)
        active_window.render()
    def slotZoomFit(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        active_window.scene.ResetCamera()
        active_window.render()

    def change_particle_resolution(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        num = int(SM.no_vertices_per_particle)
        comboBox_particle_resolution = self.ui.comboBox_particle_resolution.currentText()
        SM.colors = SM.colors_backup_particles[0::num].astype('f8').copy()/255
        active_window.scene.rm(SM.sphere_actor)
        vertices, faces = primitive.prim_sphere(name=comboBox_particle_resolution, gen_faces=False)
        res = primitive.repeat_primitive(vertices, faces, centers=SM.pos, colors=SM.colors, scales= SM.radii_spheres)
        big_verts, big_faces, big_colors, _ = res
        SM.sphere_actor = utils.get_actor_from_primitive(big_verts, big_faces, big_colors)
        SM.all_vertices_particles = utils.vertices_from_actor(SM.sphere_actor)
        SM.no_vertices_per_particle = len(SM.all_vertices_particles) / SM.no_atoms
        active_window.scene.add(SM.sphere_actor)
        vertices_particle = utils.vertices_from_actor(SM.sphere_actor)
        SM.no_vertices_all_particles = vertices_particle.shape[0]
        SM.sec_particle = np.int(SM.no_vertices_all_particles / SM.no_atoms)
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        SM.colors_backup_particles = SM.vcolors_particle.copy()

        for atom_typ in SM.unique_types:
            selected_value_radius = SM.radii_spheres[SM.atom_type == atom_typ][0]
            all_vertices_radii = 1/np.repeat(SM.radii_spheres[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
            all_vertices_radii = all_vertices_radii[:, None]
            selected_atom_mask = SM.atom_type == atom_typ
            all_vertices_mask = np.repeat(selected_atom_mask, SM.no_vertices_per_particle)
            SM.all_vertices_particles[all_vertices_mask] = float(selected_value_radius) * all_vertices_radii * (SM.all_vertices_particles[all_vertices_mask] - np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)) + np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
            SM.radii_spheres[SM.atom_type == atom_typ] = float(selected_value_radius)

        utils.update_actor(SM.sphere_actor)
        print('current value of radius: ',selected_value_radius)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        SM.pbr_params_atom = sky_box_effect_atom(active_window.scene, SM.sphere_actor, SM)
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
            SM.bond_colors = (0.8275, 0.8275, 0.8275, 1)
            SM.bond_actor = actor.line(SM.bonds, SM.bond_colors, linewidth=1)
            active_window.scene.add(SM.bond_actor)
        if comboBox_bondshape == 'Cylinder':
            SM.bond_colors = (0.8275, 0.8275, 0.8275, 1)
            SM.bond_actor = actor.streamtube(SM.bonds, SM.bond_colors, linewidth=0.5)
            active_window.scene.add(SM.bond_actor)
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        SM.colors_backup_bond = SM.vcolors_bond.copy()
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
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

    def box_builder(self):
        Ui_box.box = Ui_box()
        Ui_box.box.win = self
        Ui_box.box.show()
        Ui_box.box.showNormal()
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        Ui_box.box.initial_box_dim(SM.box_lx, SM.box_ly, SM.box_lz)

    def solution_builder(self):
        Ui_solution.sol = Ui_solution()
        Ui_solution.sol.win = self
        Ui_solution.sol.show()
        Ui_solution.sol.showNormal()
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        Ui_solution.sol.initial_box_dim(SM.box_lx, SM.box_ly, SM.box_lz)
        Ui_solution.sol.show()
        Ui_solution.sol.showNormal()

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

    def forward_movie(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.play_factor = 5

    def backward_movie(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.play_factor = -5

    def update_particle_size(self, selected_value_radius):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        selected_item = self.ui.treeWidget.selectionModel()
        for i, atom_typ in enumerate(SM.unique_types):
            if selected_item.rowIntersectsSelection(i):
                print(i, atom_typ, 'checked')
                # radii for each vertex of each atom with selected atom_type
                all_vertices_radii = 1/np.repeat(SM.radii_spheres[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
                all_vertices_radii = all_vertices_radii[:, None]
                selected_atom_mask = SM.atom_type == atom_typ
                all_vertices_mask = np.repeat(selected_atom_mask, SM.no_vertices_per_particle)
                SM.all_vertices_particles[all_vertices_mask] = float(selected_value_radius) * all_vertices_radii * \
                    (SM.all_vertices_particles[all_vertices_mask] - np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)) + \
                        np.repeat(SM.pos[SM.atom_type == atom_typ], SM.no_vertices_per_particle, axis=0)
                SM.radii_unique_types[i] = float(selected_value_radius)

        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPoints().GetData().Modified()
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

        if child.load_fullerene_cc1_file(fname):
            child.show()
        else:
            child.close()

    def save_image(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, self.tr('Save'), filter= 'PNG (*.png)')
        magnification = 4
        renderLarge = RenderLargeImage()
        renderLarge.SetInput(active_window.scene)
        renderLarge.SetMagnification(magnification)
        renderLarge.Update()
        renderLarge = RenderLargeImage()
        renderLarge.SetInput(active_window.scene)
        renderLarge.SetMagnification(magnification)
        renderLarge.Update()
        arr = numpy_support.vtk_to_numpy(renderLarge.GetOutput().GetPointData()
                                        .GetScalars())
        w, h, _ = renderLarge.GetOutput().GetDimensions()
        components = renderLarge.GetOutput().GetNumberOfScalarComponents()
        arr = arr.reshape((h, w, components))
        save_image(arr, filename)

    def save(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, self.tr('Save'), filter= 'LAMMPS (*.data);;PDB (*.pdb);;GROMACS (*.gro);;XYZ (*.xyz)')
        if not fname:
            return
        file_name = os.path.basename(fname)
        ext_file = os.path.splitext(file_name)[1]
        from furiousatoms.pdb2lmp import save_PDB2LMP
        if ext_file =='.data':
            if SM.universe_save is None:
                save_PDB2LMP(SM, fname, SM.universe)
            else:
                save_PDB2LMP(SM, fname, SM.universe_save)
        else:
            if SM.universe_save is None:
                num_atoms = SM.pos.shape[0]
                universe = MDAnalysis.Universe.empty(num_atoms, trajectory=True, n_residues=1)
                universe.atoms.positions = SM.pos
                n_residues = 1
                atom_types_list = list(SM.atom_type)
                if ((SM.box_lx==0) or (SM.box_ly==0) or (SM.box_lz==0)):
                    universe.dimensions[0] = 90
                    universe.dimensions[1] = 90
                    universe.dimensions[2] = 90
                else:
                    universe.dimensions[0] = SM.box_lx
                    universe.dimensions[1] = SM.box_ly
                    universe.dimensions[2] = SM.box_lz
                universe.add_TopologyAttr('name', atom_types_list)
                universe.add_TopologyAttr('type', atom_types_list)
                universe.add_TopologyAttr('resname', ['MOL']*n_residues)
                universe.add_TopologyAttr('masses')
                try:
                    universe.add_bonds(SM.universe.bonds.to_indices())
                except:
                    pass
                universe.atoms.write(fname)
            else:
                SM.universe_save.atoms.write(fname)

            with open(fname, 'r+') as fp:
                lines = fp.readlines()
                fp.seek(0)
                fp.truncate()
                if ((SM.box_lx==0) or (SM.box_ly==0) or (SM.box_lz==0)):
                    lines[0] = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom \n"
                else:
                    lines[0] = "Created by FURIOUS ATOMS.\n"
                fp.writelines(lines[:])

    def delete_particles(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.universe_save = active_window.delete_particles()
        self.ui.Edit_num_of_particles.setText(str(SM.no_atoms))

    def delete_bonds(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        active_window.delete_bonds()
        SM = active_window.universe_manager

    def openColorDialog_particle(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        selected_color_particle = QtWidgets.QColorDialog.getColor()
        selected_item = self.ui.treeWidget.selectionModel()

        if selected_color_particle.isValid():
            for i, atom_typ in enumerate(SM.unique_types):
                if selected_item.rowIntersectsSelection(i):
                    print(i, 'checked')
                    object_indices_particles = np.where(SM.atom_type == atom_typ)[0]
                    for object_index in object_indices_particles:
                        SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = selected_color_particle.getRgb()
                    r = int(selected_color_particle.getRgb()[0])
                    g = int(selected_color_particle.getRgb()[1])
                    b = int(selected_color_particle.getRgb()[2])
                    a = int(selected_color_particle.getRgb()[3])
                    items = self.ui.treeWidget.selectedItems()
                    for item in items:
                        item.setBackground(0,(QtGui.QBrush(QtGui.QColor(r, g, b, a))))
                    SM.colors_unique_types[i] = np.array([r/255, g/255, b/255, a/255], dtype='f8')
                    SM.colors[SM.atom_type == atom_typ] = SM.colors_unique_types[i]
            if SM.no_bonds > 0:
                bonds_indices = SM.universe.bonds.to_indices()
                bond_indices_1d = bonds_indices.ravel()
                half_bond = SM.sec_bond

                for k in object_indices_particles.tolist():
                    mem_block = np.where(bond_indices_1d==k)[0]
                    for j in mem_block:
                        SM.vcolors_bond[j * half_bond: j * half_bond + half_bond] = selected_color_particle.getRgb()
                SM.colors_backup_bond = SM.vcolors_bond.copy()
                utils.update_actor(SM.bond_actor)
                SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
            if SM.universe_save:
                SM.universe_save = active_window.delete_particles()


        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()


    # TODO: decide what you want to do for bonds that have two colors
    def openColorDialog_bond(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        selected_color_bond = QtWidgets.QColorDialog.getColor()
        if selected_color_bond.isValid():
            select_all_bonds = np.zeros(SM.no_bonds*2, dtype=np.bool)
            object_indices_bonds = np.where(select_all_bonds == False)[0]
            for object_index in object_indices_bonds:
                SM.colors_backup_bond[object_index] = selected_color_bond.getRgb()
                SM.vcolors_bond[object_index * SM.sec_bond: object_index * SM.sec_bond + SM.sec_bond] = SM.colors_backup_bond[object_index]
        SM.colors_backup_bond = SM.vcolors_bond.copy()
        if SM.universe_save:
            SM.universe_save = active_window.delete_bonds()
        utils.update_actor(SM.bond_actor)
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()

    def update_information_ui(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        SM.pbr_params_atom = sky_box_effect_atom(active_window.scene, SM.sphere_actor, SM)
        if SM.no_bonds > 0:
            SM.pbr_params_bond = sky_box_effect_bond(active_window.scene, SM.bond_actor, SM)
        self.ui.horizontalSlider_Opacity.setValue(SM.opacity*100)
        self.ui.horizontalSlider_Metallic.setValue(SM.metallic*100)
        self.ui.horizontalSlider_Roughness.setValue(SM.roughness*100)
        self.ui.horizontalSlider_Anisotropic.setValue(SM.anisotropic*100)
        self.ui.horizontalSlider_Anisotropic_rot.setValue(SM.anisotropic_rot*100)
        self.ui.horizontalSlider_Anisotropic_X.setValue(SM.anisotropic_X*100)
        self.ui.horizontalSlider_Anisotropic_Y.setValue(SM.anisotropic_Y*100)
        self.ui.horizontalSlider_Anisotropic_Z.setValue(SM.anisotropic_Z*100)
        self.ui.horizontalSlider_Coat_strength.setValue(SM.coat_strength*100)
        self.ui.horizontalSlider_Coat_roughness.setValue(SM.coat_rough*100)
        self.ui.treeWidget.clear()
        for i, typ in enumerate(SM.unique_types):
            SM.radii_spheres[SM.atom_type == typ] = SM.radii_unique_types[i]
            SM.colors[SM.atom_type == typ] = SM.colors_unique_types[i]
            set_value_radius = SM.radii_spheres[SM.atom_type == typ][0]
            self.ui.SpinBox_atom_radius.setValue((set_value_radius))
            r = (SM.colors[SM.atom_type == typ][0][0])*255
            g = (SM.colors[SM.atom_type == typ][0][1])*255
            b = (SM.colors[SM.atom_type == typ][0][2])*255
            a = (SM.colors[SM.atom_type == typ][0][3])*255
            cg = QtWidgets.QTreeWidgetItem(self.ui.treeWidget, [0, str(typ)])
            cg.setBackground(0,(QtGui.QBrush(QtGui.QColor(r, g, b, a))))

        self.ui.treeWidget.setCurrentItem(self.ui.treeWidget.topLevelItem(0))
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
        if SM.box_lx > 0 or SM.box_ly > 0 or SM.box_lz > 0:
            self.ui.Box_simulationcell.stateChanged.connect(self.check_simulationcell)
        if SM.no_bonds > 0:
            self.ui.Box_bonds.stateChanged.connect(self.check_bonds)

        self.ui.Edit_num_of_bonds.setText(str(SM.no_bonds))
        self.ui.Edit_widthX.setText(str(SM.box_lx))
        self.ui.Edit_widthY.setText(str(SM.box_ly))
        self.ui.Edit_widthZ.setText(str(SM.box_lz))
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
        self.timer = QTimer()
        self.timer.timeout.connect(self.timer_callback)
        duration = 200
        self.timer.start(duration)
        # MainWindow.iren.Initialize()

    def timer_callback(self):
        active_window = self.active_mdi_child()
        if not active_window:
            return
        SM = active_window.universe_manager
        # if SM.enable_timer is False:
        #     return
        if SM.cnt == SM.n_frames:
            return
        if SM.n_frames > 1:
            pos_R = SM.universe.trajectory[SM.cnt].positions.copy().astype('f8')
            pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, SM.box, backend='serial')
            SM.all_vertices_particles[:] = SM.initial_vertices_particles + \
                np.repeat(pos, SM.no_vertices_per_particle, axis=0)
            utils.update_actor(SM.sphere_actor)

            self.ui.Timer_animation.setValue(SM.cnt)
            self.ui.Edit_framenumber.setText(str(SM.cnt))
            self.ui.horizontalSlider_animation.setValue(SM.cnt)

        active_window.render()
        SM.cnt = SM.cnt + 1 * SM.play_factor

    def active_mdi_child(self):
        active_sub_window = self.ui.mdiArea.activeSubWindow()
        if active_sub_window:
            return active_sub_window.widget()
        else:
            self.ui.treeWidget.clear()
            self.ui.Edit_num_of_particles.setText(str())
            self.ui.Edit_num_of_particle_types.setText(str())
            self.ui.Edit_num_of_bonds.setText(str())
            self.ui.Edit_widthX.setText(str())
            self.ui.Edit_widthY.setText(str())
            self.ui.Edit_widthZ.setText(str())
            self.ui.Edit_number_of_frames.setText(str())
            self.ui.Edit_directory.setText(str())
            self.ui.Edit_fileformat.setText(str())
            self.ui.Edit_currentfile.setText(str())
            self.ui.SpinBox_atom_radius.setValue((0))

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
            app = QtGui.QGuiApplication.instance()
            app.closeAllWindows()
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