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

class ViewerVTK(QtWidgets.QWidget):
    """ Basic 3D viewer widget
    """
    sequence_number = 1

    def __init__(self, app_path=None, parent=None):
        super(ViewerVTK, self).__init__(parent)
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
        self._scene.background((0, 0, 0))
        self._showm.initialize()

    def init_variables(self):
        self.is_untitled = True
        # fetch_viz_cubemaps()
        # textures = read_viz_cubemap('skybox')
        # cubemap = load_cubemap_texture(textures)
        self._scene = window.Scene()#skybox=cubemap)
        # self._scene.skybox(visible=False)

        # self._scene = window.Scene()
        self._showm = window.ShowManager(scene=self._scene,
                                         order_transparent=True)
        self._qvtkwidget = QVTKRenderWindowInteractor(
            parent=self.ui.viewer_frame,
            rw=self._showm.window,
            iren=self._showm.iren)

        # self.universe_manager = None
        # self.pickm = pick.SelectionManager(select='faces')

    def create_connections(self):
        pass

    def make_title(self):
        self.is_untitled = True
        self.current_file = "viewer-%d" % ViewerVTK.sequence_number
        ViewerVTK.sequence_number += 1
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
    viewer = ViewerVTK()
    viewer.scene.add(actor.axes())
    viewer.show()
    sys.exit(app.exec_())
