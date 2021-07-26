import numpy as np
from numpy.linalg import norm
from fractions import gcd
from itertools import product
from furiousatoms.geomlib import Atom, Molecule, Crystal, getfragments
from furiousatoms.sharedmem import SharedMemory
import sys
from furiousatoms import io
import vtk
import numpy as np
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from furiousatoms.io import create_universe, merged_universe_with_H
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from furiousatoms.structure import bbox

SM = SharedMemory()
"""
    Ui_box_water class creates a widget for building box and water
"""

class Ui_box_water(QtWidgets.QMainWindow): #QWidget

    def __init__(self, app_path=None, parent=None):
        super(Ui_box_water, self).__init__(parent)
        self.box_water = io.load_ui_widget("box_water.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.box_water)
        self.setCentralWidget(self.box_water)
        self.setLayout(self.v_layout)
        self.resize(244, 226)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.box_water.pushButton_build_box.clicked.connect(self.box_builder_callback)
        self.box_water.pushButton_build_box.clicked.connect(lambda:self.close())
        self.box_water.pushButton_build_Solution.clicked.connect(lambda:self.close())


    def box_builder_callback(self):
        active_window = self.win.active_mdi_child()
        SM = active_window.universe_manager
        box_lx = float(self.box_water.SpinBox_lx.text())
        box_ly = float(self.box_water.SpinBox_ly.text())
        box_lz = float(self.box_water.SpinBox_lz.text())
        SM.universe.trajectory.ts.dimensions[0] = box_lx
        SM.universe.trajectory.ts.dimensions[1] = box_ly
        SM.universe.trajectory.ts.dimensions[2] = box_lz
        SM.bbox_actor, _ = bbox(box_lx, box_ly, box_lz, colors=(0, 0, 0), linewidth=1, fake_tube=True)
        active_window.scene.add(SM.bbox_actor)
        utils.update_actor(SM.bbox_actor)
        SM.bbox_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()
        active_window.show()