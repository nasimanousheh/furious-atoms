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
from fury import window, actor, utils, pick, ui, primitive
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from furiousatoms.io import create_universe, merged_universe_with_H
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from furiousatoms.structure import bbox
import MDAnalysis
import os

SM = SharedMemory()
"""
    Ui_box class creates a widget for building box and water
"""

class Ui_box(QtWidgets.QMainWindow): #QWidget

    def __init__(self, app_path=None, parent=None):
        super(Ui_box, self).__init__(parent)
        self.box = io.load_ui_widget("box.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.box)
        self.setCentralWidget(self.box)
        self.setLayout(self.v_layout)
        self.resize(225, 202)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.box.pushButton_build_box.clicked.connect(self.box_builder_callback)
        self.box.pushButton_build_box.clicked.connect(lambda:self.close())


    def box_builder_callback(self):
        active_window = self.win.active_mdi_child()
        SM = active_window.universe_manager
        active_window.scene.rm(SM.bbox_actor)
        box_lx = float(self.box.SpinBox_lx.text())
        box_ly = float(self.box.SpinBox_ly.text())
        box_lz = float(self.box.SpinBox_lz.text())
        SM.universe.trajectory.ts.dimensions[0] = box_lx
        SM.universe.trajectory.ts.dimensions[1] = box_ly
        SM.universe.trajectory.ts.dimensions[2] = box_lz
        SM.bbox_actor, _ = bbox(box_lx, box_ly, box_lz, colors=(0, 0, 0), linewidth=1, fake_tube=True)
        active_window.scene.add(SM.bbox_actor)
        utils.update_actor(SM.bbox_actor)
        SM.bbox_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        active_window.render()
        active_window.show()