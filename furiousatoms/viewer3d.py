import os
import numpy as np
from furiousatoms.io import create_universe
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from furiousatoms.molecular import UniverseManager
from furiousatoms.fullerenes_builder import load_CC1_file
from furiousatoms import io

from fury.shaders import add_shader_callback, load, shader_to_actor
from fury.utils import (get_actor_from_polydata, get_polydata_triangles,
                        get_polydata_vertices, normals_from_v_f,
                        set_polydata_normals)
from fury import window, actor, utils, pick, ui, primitive


def sky_box_effect(scene, actor, universem):
    scene.UseImageBasedLightingOn()
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cube_path = os.path.join(dir_path, 'skybox0')
    if not os.path.isdir(cube_path):
        print('This path does not exist:', cube_path)
        return
    cubemap = io.read_cubemap(cube_path, '/', '.jpg', 0)
    scene.SetEnvironmentTexture(cubemap)
    actor.GetProperty().SetInterpolationToPBR()
    fs_dec_code = load('bxdf_dec.frag')
    fs_impl_code = load('bxdf_impl.frag')
    polydata = actor.GetMapper().GetInput()
    verts = get_polydata_vertices(polydata)
    faces = get_polydata_triangles(polydata)
    normals = normals_from_v_f(verts, faces)
    set_polydata_normals(polydata, normals)
    shader_to_actor(actor, 'fragment', decl_code=fs_dec_code)
    shader_to_actor(actor, 'fragment', impl_code=fs_impl_code,
                    block='light', debug=False)

    colors_sky = window.vtk.vtkNamedColors()
    actor.GetProperty().SetColor(colors_sky.GetColor3d('White'))
    actor.GetProperty().SetMetallic(universem.metallic)
    actor.GetProperty().SetRoughness(universem.roughness)
    actor.GetProperty().SetOpacity(universem.opacity)

    def uniforms_callback(_caller, _event, calldata=None):
        if calldata is not None:
            calldata.SetUniformf('anisotropic', universem.anisotropic)
            calldata.SetUniformf('clearcoat', universem.clearcoat)

    add_shader_callback(actor, uniforms_callback)

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
        self._scene = window.Scene()
        self._showm = window.ShowManager(scene=self._scene,
                                         order_transparent=True)
        self._qvtkwidget = QVTKRenderWindowInteractor(
            parent=self.ui.viewer_frame,
            rw=self._showm.window,
            iren=self._showm.iren)

        self.universe_manager = None
        self.pickm = pick.PickingManager()

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

        sky_box_effect(self.scene, self.universe_manager.sphere_actor, self.universe_manager)
        self.scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                              view_up=(0, 1, 0))

    def delete_particles(self):
        SM = self.universe_manager
        object_indices_particles = np.where(SM.selected_particle)[0]
        print('object_indices_particles: ', object_indices_particles)
        print(SM.pos[object_indices_particles])
        SM.particle_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        for object_index in object_indices_particles:
            SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.particle_color_add
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        bonds_indices = SM.universe.bonds.to_indices()
        object_indices_bonds = []
        for object_index in object_indices_particles:
            object_indices_bonds += np.where(bonds_indices[:, 1] == object_index)[0].tolist()
            object_indices_bonds += np.where(bonds_indices[:, 0] == object_index)[0].tolist()
        SM.bond_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        for object_index_bond in object_indices_bonds:
            SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = SM.bond_color_add
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.render()
        final_pos = SM.pos.copy()
        final_pos_index = np.arange(final_pos.shape[0])
        final_pos = np.delete(final_pos, object_indices_particles, axis=0)
        final_bonds = bonds_indices.copy()
        final_bonds = np.delete(final_bonds, object_indices_bonds, axis=0)
        final_atom_types = SM.atom_type
        final_atom_types = np.delete(final_atom_types, object_indices_particles)
        final_pos_index = np.delete(final_pos_index, object_indices_particles)
        fb_shape = final_bonds.shape
        map_old_to_new = {}
        for i in range(final_pos.shape[0]-1):
            map_old_to_new[final_pos_index[i]] = i

        fb = final_bonds.ravel()
        for i in range(fb.shape[0]):
            fb[i] = map_old_to_new[fb[i]]

        final_bonds = fb.reshape(fb_shape)
        SM.selected_particle[object_indices_particles] = False
        SM.universe = create_universe(final_pos, final_bonds, final_atom_types)
        return SM.universe

    def delete_bonds(self):
        SM = self.universe_manager
        bonds_indices = SM.universe.bonds.to_indices()
        object_indices_bonds = np.where(SM.selected_bond == True)[0]
        SM.bond_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        for object_index_bond in object_indices_bonds:
            SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = SM.bond_color_add
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.render()
        final_pos = SM.pos.copy()
        final_pos_index = np.arange(final_pos.shape[0])
        # final_pos = np.delete(final_pos, object_indices_particles, axis=0)
        final_bonds = bonds_indices.copy()
        final_bonds = np.delete(final_bonds, object_indices_bonds, axis=0)
        final_atom_types = SM.atom_type
        # final_atom_types = np.delete(final_atom_types, object_indices_particles)
        # final_pos_index = np.delete(final_pos_index, object_indices_particles)
        fb_shape = final_bonds.shape
        map_old_to_new = {}
        for i in range(final_pos.shape[0]-1):
            map_old_to_new[final_pos_index[i]] = i

        fb = final_bonds.ravel()
        for i in range(fb.shape[0]):
            fb[i] = map_old_to_new[fb[i]]

        final_bonds = fb.reshape(fb_shape)
        SM.selected_particle[object_indices_bonds] = False
        SM.universe = create_universe(final_pos, final_bonds, final_atom_types)
        return SM.universe


    def left_button_press_particle_callback(self, obj, event):
        SM = self.universe_manager
        event_pos = self.pickm.event_position(iren=self.showm.iren)
        picked_info = self.pickm.pick(event_pos, self.showm.scene)
        vertex_index_particle = picked_info['vertex']
        vertices = utils.vertices_from_actor(obj)
        SM.no_vertices_all_particles = vertices.shape[0]
        object_index = np.int(np.floor((vertex_index_particle / SM.no_vertices_all_particles) * SM.no_atoms))
        print('left_button_press_particle_callback number of atoms is: ',SM.no_atoms)
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
        SM = self.universe_manager
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
