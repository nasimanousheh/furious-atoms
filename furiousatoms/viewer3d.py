
import os
import numpy as np
import MDAnalysis

from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from furiousatoms.molecular import UniverseManager
from furiousatoms.fullerenes_database import load_CC1_file


from furiousatoms import io


class Viewer3D(QtWidgets.QWidget):
    """ Basic 3D viewer widget
    """
    def __init__(self, app_path=None, parent=None):
        super(Viewer3D, self).__init__(parent)
        self.ui = io.load_ui_widget("viewer3d.ui")
        self.main_layout = QtWidgets.QVBoxLayout()
        self.main_layout.addWidget(self.ui)
        self.setLayout(self.main_layout)

        self.current_file = ""
        self.current_filepath = ""

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

        universe, no_bonds = io.load_files(fname)
        if not universe:
            return success
        self.load_universe(universe, no_bonds)
        success = True
        return success

    def load_fullerene_cc1_file(self, fname):
        success = False
        self.current_file = os.path.basename(fname)
        self.current_filepath = os.path.abspath(fname)

        universe = load_CC1_file(fname)
        if not universe:
            return success
        self.load_universe(universe)
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

    def display_universe(self):
        axes_actor = actor.axes(scale=(1, 1, 1), colorx=(1, 0, 0),
                                colory=(0, 1, 0), colorz=(0, 0, 1), opacity=1)
        self.scene.add(axes_actor)
        for act in self.universe_manager.actors():
            self.scene.add(act)
        self.scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                              view_up=(0, 1, 0))

    def delete_particles(self):
        SM = self.universe_manager
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
        self.render()
        print(object_indices_particles)
        # print(SM.pos[object_indices_particles])

    def delete_bonds(self):
        SM = self.universe_manager
        object_indices_bonds = np.where(SM.selected_bond == True)[0]
        SM.bond_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        for object_index_bond in object_indices_bonds:
            SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = SM.bond_color_add

        SM.universe.delete_bonds(SM.universe.bonds[object_indices_bonds])
        # SM.universe.delete_bonds(SM.universe.bonds.to_indices()) #delete if it is lammps
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.render()
        print(len(SM.universe.bonds))

    def update_particle_size(self, selected_value_radius):
        SM = self.universe_manager
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
        self.render()

    def left_button_press_particle_callback(self, obj, event):
        SM = self.universe_manager
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
        self.sky_box_effect(SM.sphere_actor)

        self.update_bonds_ui()
        self.timer = QtCore.QTimer()
        duration = 200
        self.timer.start(duration)
        self.timer.timeout.connect(self.timer_callback)

        SM.enable_timer = False


if __name__ == "__main__":
    import sys

    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QtWidgets.QApplication(sys.argv)
    viewer = Viewer3D()
    viewer.scene.add(actor.axes())
    viewer.show()
    sys.exit(app.exec_())
