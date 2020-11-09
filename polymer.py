import sys
import time
import vtk
import threading
import numpy as np
from PyQt5 import QtCore, QtWidgets
from fury import actor
from fury.window import Scene
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5.QtCore import QDir, QTimer
from PyQt5.QtWidgets import QAction, QApplication, QColorDialog, QFileDialog, QMenu, QMessageBox
from PyQt5.QtGui import QIcon, QImage
from lammps import load_lammps
import os
import itertools
from vtk.util import numpy_support
import pdb

import PySimpleGUI as SG


global sphere_actor, initial_vertices, no_vertices_per_sphere, total_no_frames


def get_vertices(act):
    all_vertices = np.array(numpy_support.vtk_to_numpy(
        act.GetMapper().GetInput().GetPoints().GetData()))
    return all_vertices


def get_vertices2(act):

    all_vertices = numpy_support.vtk_to_numpy(
        act.GetMapper().GetInput().GetPoints().GetData())
    return all_vertices


def set_vertices(act, num_arr):

    vtk_num_array = numpy_support.numpy_to_vtk(num_arr)
    act.GetMapper().GetInput().GetPoints().SetData(vtk_num_array)


def modified(act):
    act.GetMapper().GetInput().GetPoints().GetData().Modified()
    act.GetMapper().GetInput().ComputeBounds()


def process_lammps_file(fname):
    global sphere_actor, initial_vertices, no_vertices_per_sphere, total_no_frames, lammps_dix
    frames_cnt = 0
    lammps_dix = load_lammps(fname)
    total_no_frames = len(lammps_dix['index'])
    print('Total number of frames', total_no_frames)

    no_atoms = int(lammps_dix['index'][frames_cnt]['no_atoms'])
    box_lx = float(lammps_dix['index'][frames_cnt]['box'][0])
    box_ly = float(lammps_dix['index'][frames_cnt]['box'][1])
    box_lz = float(lammps_dix['index'][frames_cnt]['box'][2])

    box_centers = np.array([[0, 0, 0.]])
    box_directions = np.array([[0, 1., 0]])
    box_colors = np.array([[1, 0, 0, 0.2]])

    box_actor = actor.box(box_centers, box_directions, box_colors,
                          scales=(box_lx * 2, box_ly * 2, box_lz * 2))

    box_actor.GetProperty().SetRepresentationToWireframe()
    box_actor.GetProperty().SetLineWidth(10)
    atom_types = lammps_dix['index'][frames_cnt]['coords'][:, 1]
    pos = lammps_dix['index'][frames_cnt]['coords'][:, 2:5]
    colors = np.ones((no_atoms, 3))



    # pdb.set_trace()
    colors[atom_types == 1] = np.array([1., 0, 0])
    colors[atom_types == -1] = np.array([0, 0, 1.])

    radii = 0.1 * np.ones(no_atoms)

    print(no_atoms)
    print(pos.shape)
    print(colors.shape)
    print(radii.shape)

    # pdb.set_trace()

    radii[atom_types == 1] = 0.5
    radii[atom_types == -1] = 0.66

    sphere_actor = actor.sphere(centers=pos,
                                colors=colors,
                                radii=radii)

    # global all_vertices
    all_vertices = get_vertices(sphere_actor)
    initial_vertices = all_vertices.copy()
    no_vertices_per_sphere = len(all_vertices) / no_atoms

    window.ren.add(sphere_actor)
    window.ren.add(box_actor)

    print('Processed done!')


global cnt, enable_timer
cnt = 0
enable_timer = False

def timer_callback():
    global enable_timer
    if enable_timer is False:
        return

    global sphere_actor, cnt, initial_vertices, lammps_dix
    global no_vertices_per_sphere, window, no_atoms

    if cnt == total_no_frames:
        return

    pos = lammps_dix['index'][cnt]['coords'][:, 2:5]

    all_vertices = np.array(numpy_support.vtk_to_numpy(
        sphere_actor.GetMapper().GetInput().GetPoints().GetData()))

    all_vertices = get_vertices2(sphere_actor)
    all_vertices[:] = initial_vertices + \
        np.repeat(pos, no_vertices_per_sphere, axis=0)

    # all_vertices[:] += np.random.rand(*all_vertices.shape)

    # set_vertices(sphere_actor, all_vertices)
    modified(sphere_actor)
    cnt += 1
    window.vtkWidget.GetRenderWindow().Render()


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self,parent = None):
        QtWidgets.QMainWindow.__init__(self, parent)

        self.frame = QtWidgets.QFrame()
        self.setWindowIcon(QIcon('Photo_Water.jpg'))
        self.setWindowTitle('Furious Atoms (Open Visualization Tool)')

        self.vl = QtWidgets.QVBoxLayout()

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)
        bar = self.menuBar()
        # color = QColorDialog.getColor()
        self.createActions()
        self.createMenus()
        tb = self.addToolBar("Load File")
        new = QAction(QIcon("Open_file.bmp"),"Load File",self)
        tb.addAction(new)

        edit_menu = bar.addMenu('Edit')

        # adding actions to edit menu
        undo_action = QtWidgets.QAction('Undo', self)
        redo_action = QtWidgets.QAction('Redo', self)
        edit_menu.addAction(undo_action)
        edit_menu.addAction(redo_action)
        # self.fileMenu = QAction(QIcon('Photo_StavrosGaryfallidis.jpg'), '&File', self)

        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)
        self.ren = Scene()
        self.ren.background((1, 1, 1))
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        # self.iren.AddObserver("TimerEvent", timer_callback)
        # timer_id = self.iren.CreateRepeatingTimer(5000)
        self.timer = QTimer()
        self.timer.timeout.connect(timer_callback)
        duration = 200
        self.timer.start(duration)

        self.resize(2000, 1000)
        self.iren.Initialize()
        self.show()

    # def newFile(self):
    #     self.infoLabel.setText("Invoked <b>File|New</b>")

    def openColorDialog(self):
        color = QColorDialog.getColor()

    def open(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File', filter = "*.lammp*")

        process_lammps_file(fname)
        global enable_timer
        enable_timer = True

    def save(self):
        fname, _ = QFileDialog.getSaveFileName(self, 'Save File', filter = "*.lammp*")

    def print_(self):
        dialog = QPrintDialog(self.printer, self)
        if dialog.exec_():
            painter = QPainter(self.printer)
            rect = painter.viewport()
            size = self.imageLabel.pixmap().size()
            size.scale(rect.size(), Qt.KeepAspectRatio)
            painter.setViewport(rect.x(), rect.y(), size.width(), size.height())
            painter.setWindow(self.imageLabel.pixmap().rect())
            painter.drawPixmap(0, 0, self.imageLabel.pixmap())

    def zoomIn(self):
        self.scaleImage(1.25)

    def zoomOut(self):
        self.scaleImage(0.8)

    def normalSize(self):
        self.imageLabel.adjustSize()
        self.scaleFactor = 1.0

    def fitToWindow(self):
        fitToWindow = self.fitToWindowAct.isChecked()
        self.scrollArea.setWidgetResizable(fitToWindow)
        if not fitToWindow:
            self.normalSize()

        self.updateActions()

    def about(self):
        QMessageBox.about(self, "About Furious Atoms",
                "<p>The <b>Image Viewer</b> example shows how to combine "
                "QLabel and QScrollArea to display an image. QLabel is "
                "typically used for displaying text, but it can also display "
                "an image. QScrollArea provides a scrolling view around "
                "another widget. If the child widget exceeds the size of the "
                "frame, QScrollArea automatically provides scroll bars.</p>"
                "<p>The example demonstrates how QLabel's ability to scale "
                "its contents (QLabel.scaledContents), and QScrollArea's "
                "ability to automatically resize its contents "
                "(QScrollArea.widgetResizable), can be used to implement "
                "zooming and scaling features.</p>"
                "<p>In addition the example shows how to use QPainter to "
                "print an image.</p>")

    def createActions(self):

        self.openColorDialog = QAction("&Color...", self, shortcut="Ctrl+O",
                triggered=self.open)
        self.openAct = QAction("&Load File", self, shortcut="Ctrl+O",
                               triggered=self.open)

        self.printAct = QAction("&Print...", self, shortcut="Ctrl+P",
                enabled=False, triggered=self.print_)

        self.exitAct = QAction("E&xit", self, shortcut="Ctrl+Q",
                triggered=self.close)

        self.zoomInAct = QAction("Zoom &In (25%)", self, shortcut="Ctrl++",
                enabled=False, triggered=self.zoomIn)

        self.saveAct = QAction("&Export File", self, shortcut="Ctrl+S",
                               triggered=self.save)
        self.newAct = QAction("&New", self, shortcut="Ctrl+S",
                               triggered=self.open)

        self.zoomOutAct = QAction("Zoom &Out (25%)", self, shortcut="Ctrl+-",
                enabled=False, triggered=self.zoomOut)

        self.normalSizeAct = QAction("&Normal Size", self, shortcut="Ctrl+S",
                enabled=False, triggered=self.normalSize)

        self.fitToWindowAct = QAction("&Fit to Window", self, enabled=False,
                checkable=True, shortcut="Ctrl+F", triggered=self.fitToWindow)

        self.aboutAct = QAction("&About", self, triggered=self.about)


    def createMenus(self):
        self.fileMenu = QMenu("&File", self)
        self.fileMenu.addAction(self.newAct)
        self.fileMenu.addAction(self.openAct)
        self.fileMenu.addAction(self.saveAct)
        self.fileMenu.addAction(self.printAct)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        self.viewMenu = QMenu("&View", self)
        self.viewMenu.addAction(self.zoomInAct)
        self.viewMenu.addAction(self.zoomOutAct)
        self.viewMenu.addAction(self.normalSizeAct)
        self.viewMenu.addSeparator()
        self.viewMenu.addAction(self.fitToWindowAct)

        self.helpMenu = QMenu("&Help", self)
        self.helpMenu.addAction(self.aboutAct)

        self.menuBar().addMenu(self.fileMenu)
        self.menuBar().addMenu(self.viewMenu)
        self.menuBar().addMenu(self.helpMenu)

    def updateActions(self):
        self.zoomInAct.setEnabled(not self.fitToWindowAct.isChecked())
        self.zoomOutAct.setEnabled(not self.fitToWindowAct.isChecked())
        self.normalSizeAct.setEnabled(not self.fitToWindowAct.isChecked())

    def scaleImage(self, factor):
        self.scaleFactor *= factor
        self.imageLabel.resize(self.scaleFactor * self.imageLabel.pixmap().size())

        self.adjustScrollBar(self.scrollArea.horizontalScrollBar(), factor)
        self.adjustScrollBar(self.scrollArea.verticalScrollBar(), factor)

        self.zoomInAct.setEnabled(self.scaleFactor < 3.0)
        self.zoomOutAct.setEnabled(self.scaleFactor > 0.333)

    def adjustScrollBar(self, scrollBar, factor):
        scrollBar.setValue(int(factor * scrollBar.value()
                                + ((factor - 1) * scrollBar.pageStep()/2)))


def calc_prod(arr):

    for i in range(10**6):
        np.prod(arr)


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)

    global window
    window = MainWindow()

    arr = np.random.rand(10)
    # t1 = threading.Thread(target=calc_prod, args=(arr, ))

    # t1.start()

    sys.exit(app.exec_())

    # t1.join()