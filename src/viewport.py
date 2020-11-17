# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'windows.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from loadfile import load_files
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
from PyQt5.QtGui import QIcon
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from fury.window import Scene


import sys
import vtk
import numpy as np
from PyQt5 import QtWidgets
from fury import window, actor, ui, utils, pick
from PyQt5.QtWidgets import QAction, QColorDialog, QFileDialog, QMenu, QMessageBox
from pylammpsmpi import LammpsLibrary
import os
from numpy.linalg import norm
import MDAnalysis as mda
import MDAnalysis.analysis.align
from MDAnalysis import *
from vtk.util import numpy_support
from fury import disable_warnings
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF,DCD
disable_warnings()




global sphere_actor, initial_vertices, no_vertices_per_sphere, n_frames

def GetMobius():
    uResolution = 51
    vResolution = 51
    surface = vtk.vtkParametricMobius()
    surface.SetMinimumV(-0.25)
    surface.SetMaximumV(0.25)

    source = vtk.vtkParametricFunctionSource()
    source.SetUResolution(uResolution)
    source.SetVResolution(vResolution)
    source.SetParametricFunction(surface)
    source.Update()


def ReadCubeMap(folderRoot, fileRoot, ext, key):
    """
    Read the cube map.
    :param folderRoot: The folder where the cube maps are stored.
    :param fileRoot: The root of the individual cube map file names.
    :param ext: The extension of the cube map files.
    :param key: The key to data used to build the full file name.
    :return: The cubemap texture.
    """
    # A map of cube map naming conventions and the corresponding file name
    # components.
    fileNames = {
        0: ['right', 'left', 'top', 'bottom', 'front', 'back'],
        1: ['posx', 'negx', 'posy', 'negy', 'posz', 'negz'],
        2: ['-px', '-nx', '-py', '-ny', '-pz', '-nz'],
        3: ['0', '1', '2', '3', '4', '5']}
    if key in fileNames:
        fns = fileNames[key]
    else:
        print('ReadCubeMap(): invalid key, unable to continue.')
        sys.exit()
    texture = vtk.vtkTexture()
    texture.CubeMapOn()
    # Build the file names.
    for i in range(0, len(fns)):
        fns[i] = folderRoot + fileRoot + fns[i] + ext
        if not os.path.isfile(fns[i]):
            print('Nonexistent texture file:', fns[i])
            return texture
    i = 0
    for fn in fns:
        # Read the images
        readerFactory = vtk.vtkImageReader2Factory()
        imgReader = readerFactory.CreateImageReader2(fn)
        imgReader.SetFileName(fn)

        flip = vtk.vtkImageFlip()
        flip.SetInputConnection(imgReader.GetOutputPort())
        flip.SetFilteredAxis(1)  # flip y axis
        texture.SetInputConnection(i, flip.GetOutputPort(0))
        i += 1
    return texture


def box_edges(box_lx, box_ly, box_lz):

    edge1 = 0.5 * np.array([[box_lx, box_ly, box_lz],
                            [box_lx, box_ly, -box_lz],
                            [-box_lx, box_ly, -box_lz],
                            [-box_lx, box_ly, box_lz],
                            [box_lx, box_ly, box_lz]])
    edge2 = 0.5 * np.array([[box_lx, box_ly, box_lz],
                            [box_lx, -box_ly, box_lz]])
    edge3 = 0.5 * np.array([[box_lx, box_ly, -box_lz],
                            [box_lx, -box_ly, -box_lz]])
    edge4 = 0.5 * np.array([[-box_lx, box_ly, -box_lz],
                            [-box_lx, -box_ly, -box_lz]])
    edge5 = 0.5 * np.array([[-box_lx, box_ly, box_lz],
                            [-box_lx, -box_ly, box_lz]])
    lines = [edge1, -edge1, edge2, edge3, edge4, edge5]
    return lines


def process_load_file(fname):
    global sphere_actor, initial_vertices, no_vertices_per_sphere, n_frames, no_bonds, bond_actor, no_vertices_per_bond,pos, load_file, all_vertices, box, MainWindow
    frames_cnt = 0
    load_file,no_bonds = load_files(fname)
    n_frames = load_file.trajectory.n_frames
    str_no_bonds = str(no_bonds)
    MainWindow.dcfdLineEdit.insert(str_no_bonds)
    str_n_frames = str(n_frames)
    MainWindow.dcfdLineEdit_2.insert(str_n_frames)

    no_atoms = len(load_file.atoms)
    box = load_file.trajectory.ts.dimensions
    box_lx = load_file.trajectory.ts.dimensions[0]
    box_ly = load_file.trajectory.ts.dimensions[1]
    box_lz = load_file.trajectory.ts.dimensions[2]
    str_no_atoms = str(no_atoms)
    MainWindow.particleLineEdit.insert(str_no_atoms)
    box_centers = np.array([[0, 0, 25]])
    box_directions = np.array([[0, 1, 0]])
    box_colors = np.array([[255, 255, 255]])
    box_actor = actor.box(box_centers, box_directions, box_colors,
                          scales=(box_lx, box_ly, box_lz))
    box_actor.GetProperty().SetOpacity(0.)
    str_box_lx = str(box_lx)
    str_box_ly = str(box_ly)
    str_box_lz = str(box_lz)
    MainWindow.atomTypeLineEdit.insert(str_box_lx)
    MainWindow.atomTypeLineEdit_2.insert(str_box_ly)
    MainWindow.atomTypeLineEdit_3.insert(str_box_lz)


    lines = box_edges(box_lx, box_ly, box_lz)
    line_actor = actor.line(lines, colors=(0, 0, 0), linewidth=1, fake_tube=True)

    box_actor.GetProperty().SetRepresentationToWireframe()
    box_actor.GetProperty().SetLineWidth(10)
    np.unique(load_file.atoms.types)
    atom_type = load_file.atoms.types

    colors = np.ones((no_atoms, 3))
    if no_bonds == 0:
        pos_R = load_file.trajectory[0].positions.copy().astype('f8')
        pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, box, backend='serial')
        MainWindow.particleLineEdit_3.insert('0')

    if no_bonds > 0:
        pos = load_file.trajectory[0].positions.copy().astype('f8')
        bonds = load_file.bonds.to_indices()
        first_pos_bond = pos[(bonds[:,0])]
        second_pos_bond = pos[(bonds[:,1])]
        bonds = np.hstack((first_pos_bond, second_pos_bond))
        bonds = bonds.reshape(no_bonds, 2, 3)
        bond_actor = actor.streamtube(bonds, colors=(0.8275, 0.8275, 0.8275), linewidth=0.2)
        if MainWindow.CheckBox.isChecked() == True:
            MainWindow.ren.add(bond_actor)
        unique_types_bond = np.unique(load_file.bonds.types)
        str_no_unique_types_bond = str(len(unique_types_bond))
        MainWindow.particleLineEdit_3.insert(str_no_unique_types_bond)

    avg = np.average(pos, axis=0)

    radii = 0.5* np.ones(no_atoms)
    unique_types = np.unique(load_file.atoms.types)
    colors_unique_types = np.random.rand(len(unique_types), 3)
    str_no_unique_types = str(len(unique_types))
    MainWindow.particleLineEdit_2.insert(str_no_unique_types)


    for i, typ in enumerate(unique_types):
        colors[atom_type == typ] = colors_unique_types[i]

    sphere_actor = actor.sphere(centers=pos,
                                colors=colors,
                                radii=radii,theta=32, phi=32)
    all_vertices = utils.vertices_from_actor(sphere_actor)
    num_vertices = all_vertices.shape[0]
    num_objects = no_atoms
    no_vertices_per_sphere = len(all_vertices) / no_atoms
    initial_vertices = all_vertices.copy() - np.repeat(pos, no_vertices_per_sphere, axis=0)

    vcolors = utils.colors_from_actor(sphere_actor, 'colors')

    sphere_actor.GetProperty().SetInterpolationToPBR()
    # Lets use a smooth metallic surface
    # Build the pipeline
    mapper = vtk.vtkPolyDataMapper()
    source = GetMobius()
    mapper.SetInputData(source)
    vtk.vtkActor().SetMapper(mapper)

    MainWindow.ren.UseImageBasedLightingOn()
    cube_path = 'c:/Users/nasim/Devel/furious-atoms/src/skybox0'
    surface = 'boy'
    if not os.path.isdir(cube_path):
        print('This path does not exist:', cube_path)
        return
    surface = surface.lower()
    cubemap = ReadCubeMap(cube_path, '/', '.jpg', 0)
    MainWindow.ren.SetEnvironmentTexture(cubemap)
    sphere_actor.GetProperty().SetInterpolationToPBR()

    # configure the basic properties
    metallicCoefficient = 0.5 #0
    roughnessCoefficient = 0.1 #1
    colors_sky = vtk.vtkNamedColors()
    sphere_actor.GetProperty().SetColor(colors_sky.GetColor3d('White'))
    sphere_actor.GetProperty().SetMetallic(metallicCoefficient)
    sphere_actor.GetProperty().SetRoughness(roughnessCoefficient)

    basic_passes = vtk.vtkRenderStepsPass()
    ssao = vtk.vtkSSAOPass()
    sceneSize = 2000
    ssao.SetRadius(0.1 * sceneSize)
    ssao.SetBias(0.001 * sceneSize)
    ssao.SetKernelSize(128)
    ssao.BlurOff()
    ssao.SetDelegatePass(basic_passes)
    glrenderer = vtk.vtkOpenGLRenderer.SafeDownCast(Scene())
    glrenderer.SetPass(ssao)
    axes_actor = actor.axes(scale=(1, 1, 1), colorx=(1, 0, 0), colory=(0, 1, 0), colorz=(0, 0, 1), opacity=1)
    MainWindow.ren.add(axes_actor)
    # if MainWindow.CheckBox.isChecked() == True:
    MainWindow.ren.add(sphere_actor)

    # if MainWindow.CheckBox_2.isChecked() == True:
    MainWindow.ren.add(box_actor)
    MainWindow.ren.add(line_actor)
    MainWindow.ren.set_camera(focal_point=avg)

    print('Processed done!')
global cnt, enable_timer
cnt = 0
enable_timer = False


def timer_callback():
    global enable_timer
    if enable_timer is False:
        return

    global sphere_actor, cnt, initial_vertices
    global no_vertices_per_sphere, window, no_atoms, no_bonds, load_file, all_vertices, box, n_frames

    if cnt == n_frames:
        return

    if no_bonds==0:
        if MainWindow.CheckBox_3.isChecked() == True:
            pos_R = load_file.trajectory[cnt].positions.copy().astype('f8')
            pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, box, backend='serial')
            all_vertices[:] = initial_vertices + \
                np.repeat(pos, no_vertices_per_sphere, axis=0)
            utils.update_actor(sphere_actor)
            MainWindow.Timer.setRange(cnt, n_frames)
            str_cnt = str(cnt)
            MainWindow.horizontalSlider.setValue(cnt)

    MainWindow.vtkWidget.GetRenderWindow().Render()
    cnt += 1

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(2990, 1302)
        MainWindow.setMinimumSize(QtCore.QSize(0, 0))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        MainWindow.setFont(font)
        MainWindow.setMouseTracking(False)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/main.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setIconSize(QtCore.QSize(40, 40))
        MainWindow.setAnimated(True)
        MainWindow.setDocumentMode(False)
        MainWindow.setTabShape(QtWidgets.QTabWidget.Triangular)
        MainWindow.setUnifiedTitleAndToolBarOnMac(True)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.widget = QtWidgets.QWidget(self.centralwidget)
        self.widget.setGeometry(QtCore.QRect(1, 0, 361, 261))
        self.widget.setObjectName("widget")
        MainWindow.CheckBox_3 = QtWidgets.QCheckBox(self.widget)
        MainWindow.CheckBox_3.setGeometry(QtCore.QRect(0, 180, 82, 17))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        MainWindow.CheckBox_3.setFont(font)
        MainWindow.CheckBox_3.setObjectName("CheckBox_3")
        self.particleLabel_2 = QtWidgets.QLabel(self.widget)
        self.particleLabel_2.setGeometry(QtCore.QRect(208, 55, 91, 22))
        self.particleLabel_2.setObjectName("particleLabel_2")
        MainWindow.dcfdLineEdit_2 = QtWidgets.QLineEdit(self.widget)
        MainWindow.dcfdLineEdit_2.setGeometry(QtCore.QRect(130, 201, 61, 22))
        MainWindow.dcfdLineEdit_2.setObjectName("dcfdLineEdit_2")
        self.dcfdLabel = QtWidgets.QLabel(self.widget)
        self.dcfdLabel.setGeometry(QtCore.QRect(10, 85, 111, 22))
        self.dcfdLabel.setObjectName("dcfdLabel")
        MainWindow.particleLineEdit_3 = QtWidgets.QLineEdit(self.widget)
        MainWindow.particleLineEdit_3.setGeometry(QtCore.QRect(300, 85, 61, 22))
        MainWindow.particleLineEdit_3.setObjectName("particleLineEdit_3")
        MainWindow.particleLineEdit = QtWidgets.QLineEdit(self.widget)
        MainWindow.particleLineEdit.setGeometry(QtCore.QRect(125, 55, 61, 22))
        MainWindow.particleLineEdit.setObjectName("particleLineEdit")
        MainWindow.particleLineEdit.setReadOnly(True)###########
        MainWindow.particleLineEdit_2 = QtWidgets.QLineEdit(self.widget)
        MainWindow.particleLineEdit_2.setGeometry(QtCore.QRect(300, 55, 61, 22))
        MainWindow.particleLineEdit_2.setObjectName("particleLineEdit_2")
        self.atomTypeLabel_2 = QtWidgets.QLabel(self.widget)
        self.atomTypeLabel_2.setGeometry(QtCore.QRect(110, 144, 21, 22))
        self.atomTypeLabel_2.setObjectName("atomTypeLabel_2")
        self.particleLabel_3 = QtWidgets.QLabel(self.widget)
        self.particleLabel_3.setGeometry(QtCore.QRect(208, 85, 91, 22))
        self.particleLabel_3.setObjectName("particleLabel_3")
        self.particleLabel = QtWidgets.QLabel(self.widget)
        self.particleLabel.setGeometry(QtCore.QRect(10, 55, 111, 22))
        self.particleLabel.setObjectName("particleLabel")
        self.atomTypeLabel = QtWidgets.QLabel(self.widget)
        self.atomTypeLabel.setGeometry(QtCore.QRect(10, 144, 21, 22))
        self.atomTypeLabel.setObjectName("atomTypeLabel")
        self.atomTypeLabel_3 = QtWidgets.QLabel(self.widget)
        self.atomTypeLabel_3.setGeometry(QtCore.QRect(206, 144, 21, 22))
        self.atomTypeLabel_3.setObjectName("atomTypeLabel_3")
        MainWindow.atomTypeLineEdit = QtWidgets.QLineEdit(self.widget)
        MainWindow.atomTypeLineEdit.setGeometry(QtCore.QRect(30, 144, 61, 22))
        MainWindow.atomTypeLineEdit.setObjectName("atomTypeLineEdit")
        self.dcfdLabel_2 = QtWidgets.QLabel(self.widget)
        self.dcfdLabel_2.setGeometry(QtCore.QRect(10, 201, 111, 22))
        self.dcfdLabel_2.setObjectName("dcfdLabel_2")
        MainWindow.dcfdLineEdit = QtWidgets.QLineEdit(self.widget)
        MainWindow.dcfdLineEdit.setGeometry(QtCore.QRect(125, 85, 61, 22))
        MainWindow.dcfdLineEdit.setObjectName("dcfdLineEdit")
        MainWindow.CheckBox = QtWidgets.QCheckBox(self.widget)
        MainWindow.CheckBox.setGeometry(QtCore.QRect(2, 35, 82, 17))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        MainWindow.CheckBox.setFont(font)
        MainWindow.CheckBox.setObjectName("CheckBox")
        MainWindow.CheckBox_2 = QtWidgets.QCheckBox(self.widget)
        MainWindow.CheckBox_2.setGeometry(QtCore.QRect(1, 124, 82, 17))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        MainWindow.CheckBox_2.setFont(font)
        MainWindow.CheckBox_2.setObjectName("CheckBox_2")
        MainWindow.atomTypeLineEdit_2 = QtWidgets.QLineEdit(self.widget)
        MainWindow.atomTypeLineEdit_2.setGeometry(QtCore.QRect(133, 144, 61, 22))
        MainWindow.atomTypeLineEdit_2.setObjectName("atomTypeLineEdit_2")
        MainWindow.atomTypeLineEdit_3 = QtWidgets.QLineEdit(self.widget)
        MainWindow.atomTypeLineEdit_3.setGeometry(QtCore.QRect(230, 144, 61, 22))
        MainWindow.atomTypeLineEdit_3.setObjectName("atomTypeLineEdit_3")
        self.label = QtWidgets.QLabel(self.widget)
        self.label.setGeometry(QtCore.QRect(0, 0, 361, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setStyleSheet("background-color: rgb(204, 204, 204);")
        self.label.setObjectName("label")
        self.widget_2 = QtWidgets.QWidget(self.centralwidget)
        self.widget_2.setGeometry(QtCore.QRect(430, 1275, 2591, 71))
        self.widget_2.setObjectName("widget_2")
        MainWindow.pushButton_8 = QtWidgets.QPushButton(self.widget_2)
        MainWindow.pushButton_8.setGeometry(QtCore.QRect(160, 10, 51, 51))
        MainWindow.pushButton_8.setText("")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/play/play.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.pushButton_8.setIcon(icon1)
        MainWindow.pushButton_8.setIconSize(QtCore.QSize(45, 45))
        MainWindow.pushButton_8.setObjectName("pushButton_8")
        MainWindow.Timer = QtWidgets.QSpinBox(self.widget_2)
        MainWindow.Timer.setGeometry(QtCore.QRect(0, 10, 61, 51))
        MainWindow.Timer.setObjectName("Timer")
        MainWindow.pushButton_10 = QtWidgets.QPushButton(self.widget_2)
        MainWindow.pushButton_10.setGeometry(QtCore.QRect(60, 10, 51, 51))
        MainWindow.pushButton_10.setText("")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/play/backward.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.pushButton_10.setIcon(icon2)
        MainWindow.pushButton_10.setIconSize(QtCore.QSize(40, 40))
        MainWindow.pushButton_10.setObjectName("pushButton_10")
        MainWindow.horizontalSlider = QtWidgets.QSlider(self.widget_2)
        MainWindow.horizontalSlider.setGeometry(QtCore.QRect(270, 10, 2291, 51))
        MainWindow.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        MainWindow.horizontalSlider.setObjectName("horizontalSlider")
        MainWindow.horizontalSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        MainWindow.pushButton_5 = QtWidgets.QPushButton(self.widget_2)
        MainWindow.pushButton_5.setGeometry(QtCore.QRect(110, 10, 51, 51))
        MainWindow.pushButton_5.setText("")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/play/pause.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.pushButton_5.setIcon(icon3)
        MainWindow.pushButton_5.setIconSize(QtCore.QSize(45, 45))
        MainWindow.pushButton_5.setObjectName("pushButton_5")
        MainWindow.pushButton_6 = QtWidgets.QPushButton(self.widget_2)
        MainWindow.pushButton_6.setGeometry(QtCore.QRect(210, 10, 51, 51))
        MainWindow.pushButton_6.setText("")
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/play/forward.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.pushButton_6.setIcon(icon4)
        MainWindow.pushButton_6.setIconSize(QtCore.QSize(40, 40))
        MainWindow.pushButton_6.setObjectName("pushButton_6")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 2990, 21))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuView = QtWidgets.QMenu(self.menubar)
        self.menuView.setObjectName("menuView")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        self.menuModify = QtWidgets.QMenu(self.menubar)
        self.menuModify.setObjectName("menuModify")
        self.menuCompute = QtWidgets.QMenu(self.menubar)
        self.menuCompute.setObjectName("menuCompute")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionLoad_file = QtWidgets.QAction(MainWindow)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/file/load.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionLoad_file.setIcon(icon5)
        self.actionLoad_file.setObjectName("actionLoad_file")
        self.actionSave_file = QtWidgets.QAction(MainWindow)
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/file/save.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSave_file.setIcon(icon6)
        self.actionSave_file.setObjectName("actionSave_file")
        self.actionExit = QtWidgets.QAction(MainWindow)
        icon7 = QtGui.QIcon()
        icon7.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/file/exit.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionExit.setIcon(icon7)
        self.actionExit.setObjectName("actionExit")
        self.actionZoom_in = QtWidgets.QAction(MainWindow)
        icon8 = QtGui.QIcon()
        icon8.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/view/zoomin.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionZoom_in.setIcon(icon8)
        self.actionZoom_in.setObjectName("actionZoom_in")
        self.actionZoom_Out = QtWidgets.QAction(MainWindow)
        icon9 = QtGui.QIcon()
        icon9.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/view/zoomout.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionZoom_Out.setIcon(icon9)
        self.actionZoom_Out.setObjectName("actionZoom_Out")
        self.actionFit_to_Window = QtWidgets.QAction(MainWindow)
        icon10 = QtGui.QIcon()
        icon10.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/view/fit-to-width.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionFit_to_Window.setIcon(icon10)
        self.actionFit_to_Window.setObjectName("actionFit_to_Window")
        self.actionAbout = QtWidgets.QAction(MainWindow)
        icon11 = QtGui.QIcon()
        icon11.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/about/about.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionAbout.setIcon(icon11)
        self.actionAbout.setObjectName("actionAbout")
        self.actionNew = QtWidgets.QAction(MainWindow)
        icon12 = QtGui.QIcon()
        icon12.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/file/new.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionNew.setIcon(icon12)
        self.actionNew.setObjectName("actionNew")
        self.actionSave_as = QtWidgets.QAction(MainWindow)
        icon13 = QtGui.QIcon()
        icon13.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/file/save as.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSave_as.setIcon(icon13)
        self.actionSave_as.setObjectName("actionSave_as")
        self.actionRemove_Particle = QtWidgets.QAction(MainWindow)
        self.actionRemove_Particle.setObjectName("actionRemove_Particle")
        self.actionReplace_Particle = QtWidgets.QAction(MainWindow)
        self.actionReplace_Particle.setObjectName("actionReplace_Particle")
        self.actionRemove = QtWidgets.QAction(MainWindow)
        self.actionRemove.setObjectName("actionRemove")
        self.actionReplace = QtWidgets.QAction(MainWindow)
        self.actionReplace.setObjectName("actionReplace")
        self.actionPCF = QtWidgets.QAction(MainWindow)
        self.actionPCF.setObjectName("actionPCF")
        self.actionPair_Correlation_Function = QtWidgets.QAction(MainWindow)
        self.actionPair_Correlation_Function.setObjectName("actionPair_Correlation_Function")
        self.actionMean_Square_Desplacement = QtWidgets.QAction(MainWindow)
        self.actionMean_Square_Desplacement.setObjectName("actionMean_Square_Desplacement")
        self.actionParticle = QtWidgets.QAction(MainWindow)
        icon14 = QtGui.QIcon()
        icon14.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/shape/sphere.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionParticle.setIcon(icon14)
        self.actionParticle.setObjectName("actionParticle")
        self.actionBond = QtWidgets.QAction(MainWindow)
        icon15 = QtGui.QIcon()
        icon15.addPixmap(QtGui.QPixmap("C:/Users/nasim/OneDrive/Desktop/image-icons/shape/bond.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionBond.setIcon(icon15)
        self.actionBond.setObjectName("actionBond")
        self.menuFile.addAction(self.actionNew)
        self.menuFile.addAction(self.actionLoad_file)
        self.menuFile.addAction(self.actionSave_file)
        self.menuFile.addAction(self.actionSave_as)
        self.menuFile.addAction(self.actionExit)
        self.menuView.addAction(self.actionZoom_in)
        self.menuView.addAction(self.actionZoom_Out)
        self.menuView.addAction(self.actionFit_to_Window)
        self.menuHelp.addAction(self.actionAbout)
        self.menuModify.addAction(self.actionParticle)
        self.menuModify.addAction(self.actionBond)
        self.menuCompute.addAction(self.actionPCF)
        self.menuCompute.addAction(self.actionPair_Correlation_Function)
        self.menuCompute.addAction(self.actionMean_Square_Desplacement)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menubar.addAction(self.menuModify.menuAction())
        self.menubar.addAction(self.menuCompute.menuAction())
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
##############################

        self.actionLoad_file.triggered.connect(self.open)
        self.actionSave_file.triggered.connect(self.save)
        # MainWindow.pushButton_8.

    # def click(self,text):
    #     self.particleLineEdit.setText(text)

    def open(self):
        fname, _ = QFileDialog.getOpenFileName(MainWindow, 'Load')#, filter = "*.lammp*")
        process_load_file(fname)
        global enable_timer
        enable_timer = True

    def save(self):
        fname, _ = QFileDialog.getSaveFileName(MainWindow, 'Save')#, filter = "*.lammp*")

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Furious Atoms (Open Visualization Tool)"))
        MainWindow.CheckBox_3.setText(_translate("MainWindow", "Animation"))
        self.particleLabel_2.setText(_translate("MainWindow", "Types of Atom:"))
        self.dcfdLabel.setText(_translate("MainWindow", "Number of Bonds:"))
        self.atomTypeLabel_2.setText(_translate("MainWindow", "Ly"))
        self.particleLabel_3.setText(_translate("MainWindow", "Types of Bond:"))
        self.particleLabel.setText(_translate("MainWindow", "Number of Atoms:"))
        self.atomTypeLabel.setText(_translate("MainWindow", "Lx"))
        self.atomTypeLabel_3.setText(_translate("MainWindow", "Lz"))
        self.dcfdLabel_2.setText(_translate("MainWindow", "Number of Frames"))
        MainWindow.CheckBox.setText(_translate("MainWindow", "Particle"))
        MainWindow.CheckBox_2.setText(_translate("MainWindow", "Simulation Cell"))
        self.label.setText(_translate("MainWindow", "Properties"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.menuModify.setTitle(_translate("MainWindow", "Modify"))
        self.menuCompute.setTitle(_translate("MainWindow", "Compute"))
        self.actionLoad_file.setText(_translate("MainWindow", "Load"))
        self.actionLoad_file.setShortcut(_translate("MainWindow", "Ctrl+O"))
        self.actionSave_file.setText(_translate("MainWindow", "Save"))
        self.actionSave_file.setShortcut(_translate("MainWindow", "Ctrl+S"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))
        self.actionExit.setShortcut(_translate("MainWindow", "Ctrl+Q"))
        self.actionZoom_in.setText(_translate("MainWindow", "Zoom In"))
        self.actionZoom_in.setShortcut(_translate("MainWindow", "Ctrl++"))
        self.actionZoom_Out.setText(_translate("MainWindow", "Zoom Out"))
        self.actionZoom_Out.setShortcut(_translate("MainWindow", "Ctrl+-"))
        self.actionFit_to_Window.setText(_translate("MainWindow", "Fit to Window"))
        self.actionFit_to_Window.setShortcut(_translate("MainWindow", "Ctrl+F"))
        self.actionAbout.setText(_translate("MainWindow", "About"))
        self.actionNew.setText(_translate("MainWindow", "New"))
        self.actionNew.setShortcut(_translate("MainWindow", "Ctrl+N"))
        self.actionSave_as.setText(_translate("MainWindow", "Save as"))
        self.actionRemove_Particle.setText(_translate("MainWindow", "Remove Particle"))
        self.actionReplace_Particle.setText(_translate("MainWindow", "Replace Particle"))
        self.actionRemove.setText(_translate("MainWindow", "Remove"))
        self.actionReplace.setText(_translate("MainWindow", "Replace"))
        self.actionPCF.setText(_translate("MainWindow", "Density"))
        self.actionPair_Correlation_Function.setText(_translate("MainWindow", "Pair Correlation Function"))
        self.actionMean_Square_Desplacement.setText(_translate("MainWindow", "Mean Square Displacement"))
        self.actionParticle.setText(_translate("MainWindow", "Particle"))
        self.actionBond.setText(_translate("MainWindow", "Bond"))
        # FURY

        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setGeometry(QtCore.QRect(355, -10, 3090, 1302))

        # self.frame = QtWidgets.QFrame()
        self.vl = QtWidgets.QVBoxLayout()
        self.frame.setLayout(self.vl)
        # MainWindow.setCentralWidget(self.frame)
        MainWindow.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(MainWindow.vtkWidget)
        MainWindow.ren = Scene()

        # self.ren.background((0, 0, 0))
        # self.ren.background((.3, .4, .5))
        MainWindow.ren.background((1.0, 1.0, 1.0))

        MainWindow.vtkWidget.GetRenderWindow().AddRenderer(MainWindow.ren)
        self.iren = MainWindow.vtkWidget.GetRenderWindow().GetInteractor()
        self.timer = QTimer()
        self.timer.timeout.connect(timer_callback)
        duration = 200
        self.timer.start(duration)

        MainWindow.resize(2000, 1000)
        self.iren.Initialize()
        # MainWindow.show()




if __name__ == "__main__":
    import sys
    global MainWindow
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

