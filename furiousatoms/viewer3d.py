
import numpy as np
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


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

        self._scene = window.Scene()
        self._showm = window.ShowManager(scene=self._scene,
                                         order_transparent=True)
        self._qvtkwidget = QVTKRenderWindowInteractor(
            parent=self.ui.viewer_frame,
            rw=self._showm.window,
            iren=self._showm.iren)

        self.init_settings()
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

    def create_connections(self):
        pass

    @property
    def scene(self):
        """The foo property."""
        return self._scene

    def add_to_scene(self, actors):
        self.scene.add(actors)


if __name__ == "__main__":
    import sys

    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QtWidgets.QApplication(sys.argv)
    viewer = Viewer3D()
    viewer.add_to_scene(actor.axes())
    viewer.show()
    sys.exit(app.exec_())
