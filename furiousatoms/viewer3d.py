import os
import numpy as np
from furiousatoms.io import create_universe
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from furiousatoms.molecular import UniverseManager
from furiousatoms.fullerenes_builder import load_CC1_file
from furiousatoms import io
from fury import window, actor, utils, pick, ui, primitive, material
from fury.data import fetch_viz_cubemaps, read_viz_cubemap
from fury.io import load_cubemap_texture
from fury.utils import (normals_from_actor, tangents_to_actor,
                        tangents_from_direction_of_anisotropy, update_polydata_normals)
from fury.lib import numpy_support


def sky_box_effect_atom(scene, actor, universem):
    polydata = actor.GetMapper().GetInput()
    update_polydata_normals(polydata)
    normals = normals_from_actor(actor)
    doa = [0, 1, .5]
    tangents = tangents_from_direction_of_anisotropy(normals, doa)
    tangents_to_actor(actor, tangents)
    material_params = material.manifest_pbr(actor)
    return material_params

def sky_box_effect_bond(scene, actor, universem):
    doa = [0, 1, .5]
    normals = normals_from_actor(actor)
    tangents = tangents_from_direction_of_anisotropy(normals, doa)
    tangents_to_actor(actor, tangents)
    material_params = material.manifest_pbr(actor)
    return material_params

class Viewer3D(QtWidgets.QWidget):
    """ Basic 3D viewer widget
    """
    sequence_number = 1

    def __init__(self, app_path=None, parent=None):
        super(Viewer3D, self).__init__(parent)
        self.ui = io.load_ui_widget("viewer3d.ui")
        self.main_layout = QtWidgets.QVBoxLayout()
        self.main_layout.addWidget(self.ui)
        self.setLayout(self.main_layout)

        self.current_file = ""
        self.current_filepath = ""
        self.current_filedir = ""
        self.current_extension = ""

        self.init_settings()
        self.init_variables()
        self.init_interface()
        self.create_connections()

    def init_settings(self):
        pass

    def init_interface(self):
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self._qvtkwidget, stretch=1)
        self.ui.viewer_frame.setLayout(self.v_layout)
        self._scene.background((1, 1, 1))
        self._showm.initialize()

    def init_variables(self):
        self.is_untitled = True
        fetch_viz_cubemaps()
        textures = read_viz_cubemap('skybox')
        cubemap = load_cubemap_texture(textures)
        self._scene = window.Scene(skybox=cubemap)
        self._scene.skybox(visible=False)

        # self._scene = window.Scene()
        self._showm = window.ShowManager(scene=self._scene,
                                         order_transparent=True)
        self._qvtkwidget = QVTKRenderWindowInteractor(
            parent=self.ui.viewer_frame,
            rw=self._showm.window,
            iren=self._showm.iren)

        self.universe_manager = None
        self.pickm = pick.SelectionManager(select='faces')

    def create_connections(self):
        pass

    def make_title(self):
        self.is_untitled = True
        self.current_file = "viewer-%d" % Viewer3D.sequence_number
        Viewer3D.sequence_number += 1
        self.setWindowTitle(self.current_file)

    @property
    def scene(self):
        """The foo property."""
        return self._scene

    @property
    def showm(self):
        """The foo property."""
        return self._showm

    def render(self):
        self._qvtkwidget.GetRenderWindow().Render()

    def load_file(self, fname):
        success = False
        self.current_file = os.path.basename(fname)
        self.current_filepath = os.path.abspath(fname)
        self.current_filedir = os.path.dirname(self.current_filepath)
        self.current_extension = os.path.splitext(self.current_filepath)[1]
        self.is_untitled = False

        universe, no_bonds = io.load_files(fname)
        if not universe:
            return success
        self.load_universe(universe, no_bonds)
        success = True
        return success

    def load_universe(self, universe, no_bonds=0):
        self.universe_manager = UniverseManager(universe, no_bonds)
        self.create_universe_connections()
        self.display_universe()
        self.setWindowTitle(self.current_file)

    def create_universe_connections(self):
        if self.universe_manager.have_bonds:
            self.universe_manager.bond_actor.AddObserver(
                "LeftButtonPressEvent", self.left_button_press_bond_callback)
        self.universe_manager.sphere_actor.AddObserver(
            "LeftButtonPressEvent", self.left_button_press_particle_callback)

    def load_fullerene_cc1_file(self, fname):
        success = False
        self.current_file = os.path.basename(fname)
        self.current_filepath = os.path.abspath(fname)
        self.current_filedir = os.path.dirname(self.current_filepath)
        self.current_extension = os.path.splitext(self.current_filepath)[1]
        self.is_untitled = False
        un = load_CC1_file(fname)
        universe, no_bonds = io.load_files(un)
        if not universe:
            return success
        self.load_universe(universe, no_bonds)
        success = True
        return success

    def display_universe(self):
        axes_actor = actor.axes(scale=(1, 1, 1), colorx=(1, 0, 0),
                                colory=(0, 1, 0), colorz=(0, 0, 1), opacity=1)
        self.scene.add(axes_actor)
        for act in self.universe_manager.actors():
            self.scene.add(act)

        try:
            self.pbr_params_atom = sky_box_effect_atom(self.scene, self.universe_manager.sphere_actor, self.universe_manager)
        except:
            pass

        try:
            self.pbr_params_bond = sky_box_effect_bond(self.scene, self.universe_manager.bond_actor, self.universe_manager)
        except:
            pass

        try:
            self.scene.add(self.universe_manager.bbox_actor)
        except:
            pass

        self.scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                              view_up=(0, 1, 0))

    def delete_particles(self):
        SM = self.universe_manager
        SM.object_indices_particles = np.where(SM.selected_particle)[0].tolist()
        SM.object_indices_particles += np.where(SM.deleted_particles == True)[0].tolist()
        SM.object_indices_particles = np.asarray(SM.object_indices_particles)
        print('object_indices_particles: ', SM.object_indices_particles)
        print(SM.pos[SM.object_indices_particles])
        SM.particle_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        for object_index in SM.object_indices_particles:
            SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.particle_color_add
        utils.update_actor(SM.sphere_actor)
        # SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        final_pos = SM.pos.copy()
        final_pos_index = np.arange(final_pos.shape[0])
        final_pos = np.delete(final_pos, SM.object_indices_particles, axis=0)
        final_pos_index = np.delete(final_pos_index,
                                    SM.object_indices_particles)
        final_atom_types = SM.atom_type.copy()
        final_atom_types = np.delete(final_atom_types,
                                    SM.object_indices_particles)
        try:
            bonds_indices = SM.universe.bonds.to_indices()
            object_indices_bonds = np.where(SM.deleted_bonds == True)[0].tolist()
            for object_index in SM.object_indices_particles:
                object_indices_bonds += np.where(bonds_indices[:, 1] == object_index)[0].tolist()
                object_indices_bonds += np.where(bonds_indices[:, 0] == object_index)[0].tolist()
            bond_color_add = np.array([255, 0, 0, 0], dtype='uint8')
            SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
            for object_index_bond in object_indices_bonds:
                object_index_bond_n = object_index_bond * 2
                object_index_bond_n2 = object_index_bond * 2 + 1
                SM.vcolors_bond[object_index_bond_n * SM.sec_bond: object_index_bond_n * SM.sec_bond + SM.sec_bond] = bond_color_add
                SM.vcolors_bond[object_index_bond_n2 * SM.sec_bond: object_index_bond_n2 * SM.sec_bond + SM.sec_bond] = bond_color_add
            utils.update_actor(SM.bond_actor)
            SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
            self.render()
            final_bonds = bonds_indices.copy()
            final_bonds = np.delete(final_bonds, object_indices_bonds, axis=0)
            fb_shape = final_bonds.shape
            map_old_to_new = {}
            for i in range(final_pos.shape[0]):
                map_old_to_new[final_pos_index[i]] = i
            fb = final_bonds.ravel()
            for i in range(fb.shape[0]):
                fb[i] = map_old_to_new[fb[i]]

            final_bonds = fb.reshape(fb_shape)
            SM.deleted_bonds[object_indices_bonds] = True
        except:
            final_bonds = None

        SM.selected_particle[SM.object_indices_particles] = False
        SM.deleted_particles[SM.object_indices_particles] = True
        SM.universe_save = create_universe(final_pos, final_bonds, final_atom_types, SM.box_lx, SM.box_ly, SM.box_lz)
        return SM.universe_save
    def delete_bonds(self):
        SM = self.universe_manager
        SM.object_indices_particles = np.where(SM.selected_particle)[0].tolist()
        SM.object_indices_particles += np.where(SM.deleted_particles == True)[0].tolist()
        SM.object_indices_particles = np.asarray(SM.object_indices_particles)
        final_pos = SM.pos.copy()
        final_pos_index = np.arange(final_pos.shape[0])
        final_atom_types = SM.atom_type.copy()

        object_indices_bonds = np.where(SM.selected_bond)[0].tolist()
        object_indices_bonds += np.where(SM.deleted_bonds == True)[0].tolist()
        object_indices_bonds = np.asarray(object_indices_bonds)

        bonds_indices = SM.universe.bonds.to_indices()
        bond_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        for object_index_bond in object_indices_bonds * 2:
            if object_index_bond % 2 == 0:
                object_index_bond2 = object_index_bond + 1
            else:
                object_index_bond2 = object_index_bond - 1
            SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = bond_color_add
            SM.vcolors_bond[object_index_bond2 * SM.sec_bond: object_index_bond2 * SM.sec_bond + SM.sec_bond] = bond_color_add
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.render()
        final_bonds = bonds_indices.copy()
        final_bonds = np.delete(final_bonds, object_indices_bonds, axis=0)

        if len(SM.object_indices_particles)> 0:
            final_pos = np.delete(final_pos, SM.object_indices_particles, axis=0)
            final_pos_index = np.delete(final_pos_index,
                                        SM.object_indices_particles)
            final_atom_types = np.delete(final_atom_types,
                                        SM.object_indices_particles)
            fb_shape = final_bonds.shape
            map_old_to_new = {}
            for i in range(final_pos.shape[0]):
                map_old_to_new[final_pos_index[i]] = i
            fb = final_bonds.ravel()
            for i in range(fb.shape[0]):
                fb[i] = map_old_to_new[fb[i]]
            final_bonds = fb.reshape(fb_shape)


        SM.selected_bond[object_indices_bonds] = False
        SM.deleted_bonds[object_indices_bonds] = True
        SM.universe_save = create_universe(final_pos, final_bonds, final_atom_types, SM.box_lx, SM.box_ly, SM.box_lz)
        return SM.universe_save

    def left_button_press_particle_callback(self, obj, event):
        # print(obj)
        SM = self.universe_manager

        event_pos = self.pickm.event_position(iren=self.showm.iren)
        picked_info = self.pickm.pick(event_pos, self.showm.scene)
        polydata = obj.GetMapper().GetInput()
        faces = utils.get_polydata_triangles(polydata)

        face_indices = picked_info['face']
        vertex_index_particle = faces[face_indices[0]][0]
        vertices = utils.vertices_from_actor(obj)
        SM.no_vertices_all_particles = vertices.shape[0]
        object_index = np.int(np.floor((vertex_index_particle / SM.no_vertices_all_particles) * SM.no_atoms))
        print('left_button_press_particle_callback number of atoms is: ',SM.no_atoms)

        if (not SM.selected_particle[object_index]):
            SM.particle_color_add = np.array([255, 0, 0, 255], dtype='uint8')
            SM.selected_particle[object_index] = True
        else:
            # SM.particle_color_add = SM.colors_backup_particles[object_index]
            SM.particle_color_add = (SM.colors[object_index] * 255).astype('uint8')
            SM.selected_particle[object_index] = False

        SM.vcolors_particle = utils.colors_from_actor(obj, 'colors')
        SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.particle_color_add
        utils.update_actor(obj)
        obj.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()


    def left_button_press_bond_callback(self, obj, event):
        SM = self.universe_manager
        event_pos = self.pickm.event_position(iren=self.showm.iren)
        picked_info = self.pickm.pick(event_pos, self.showm.scene)
        polydata = obj.GetMapper().GetInput()
        faces = utils.get_polydata_triangles(polydata)

        face_indices = picked_info['face']
        vertex_index_bond = faces[face_indices[0]][0]
        vertices = utils.vertices_from_actor(obj)
        SM.no_vertices_all_bonds = vertices.shape[0]

        object_index_bond = np.int(np.floor((vertex_index_bond / SM.no_vertices_all_bonds) * 2 * SM.no_bonds))
        if object_index_bond % 2 == 0:
            object_index_bond2 = object_index_bond + 1
        else:
            object_index_bond2 = object_index_bond - 1

        if not SM.selected_bond[object_index_bond // 2 ]:
            bond_color_add = np.array([255, 0, 0, 255], dtype='uint8')
            bond_color_add2 = np.array([255, 0, 0, 255], dtype='uint8')
            SM.selected_bond[object_index_bond//2] = True
        else:
            bond_color_add = SM.colors_backup_bond[object_index_bond * SM.sec_bond]
            bond_color_add2 = SM.colors_backup_bond[object_index_bond2 * SM.sec_bond]
            SM.selected_bond[object_index_bond//2] = False

        SM.vcolors_bond = utils.colors_from_actor(obj, 'colors')

        SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = bond_color_add
        SM.vcolors_bond[object_index_bond2 * SM.sec_bond: object_index_bond2 * SM.sec_bond + SM.sec_bond] = bond_color_add2
        utils.update_actor(obj)
        obj.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()

    def process_universe(self, universe):
        self.timer = QtCore.QTimer()
        duration = 200
        self.timer.start(duration)
        self.timer.timeout.connect(self.timer_callback)
        self.enable_timer = False


if __name__ == "__main__":
    import sys

    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QtWidgets.QApplication(sys.argv)
    viewer = Viewer3D()
    viewer.scene.add(actor.axes())
    viewer.show()
    sys.exit(app.exec_())
