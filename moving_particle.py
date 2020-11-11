import sys
import vtk
import numpy as np
from PyQt5 import QtWidgets
from fury import window, actor, ui, utils, pick
from fury.window import Scene
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QAction, QColorDialog, QFileDialog, QMenu, QMessageBox
from PyQt5.QtGui import QIcon
from lammps import load_lammps
from pylammpsmpi import LammpsLibrary
import os
from numpy.linalg import norm
import MDAnalysis as mda
import MDAnalysis.analysis.align
from MDAnalysis import *
from vtk.util import numpy_support
from fury import disable_warnings
import MDAnalysis as mda
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
    print(box_lx, box_ly, box_lz)
    return lines


def process_lammps_file(fname):
    global sphere_actor, initial_vertices, no_vertices_per_sphere, n_frames, no_bonds, bond_actor, no_vertices_per_bond,pos, lammps_file,all_vertices,box
    frames_cnt = 0
    lammps_file,no_bonds = load_lammps(fname)
    n_frames = lammps_file.trajectory.n_frames
    print('Total number of frames', n_frames)

    no_atoms = len(lammps_file.atoms)
    box = lammps_file.trajectory.ts.dimensions
    box_lx = lammps_file.trajectory.ts.dimensions[0]
    box_ly = lammps_file.trajectory.ts.dimensions[1]
    box_lz = lammps_file.trajectory.ts.dimensions[2]

    box_centers = np.array([[0, 0, 25]])
    box_directions = np.array([[0, 1, 0]])
    box_colors = np.array([[255, 255, 255]])
    box_actor = actor.box(box_centers, box_directions, box_colors,
                          scales=(box_lx, box_ly, box_lz))
    box_actor.GetProperty().SetOpacity(0.)

    lines = box_edges(box_lx, box_ly, box_lz)
    line_actor = actor.line(lines, colors=(0.5, 0, 0.5), linewidth=10, fake_tube=True)

    box_actor.GetProperty().SetRepresentationToWireframe()
    box_actor.GetProperty().SetLineWidth(10)
    np.unique(lammps_file.atoms.types)
    atom_type = lammps_file.atoms.types
    colors = np.ones((no_atoms, 3))
    if no_bonds == 0:
        pos_R = lammps_file.trajectory[0].positions.copy().astype('f8')
        pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, box, backend='serial')


    if no_bonds > 0:
        pos = lammps_file.trajectory[0].positions.copy().astype('f8')
        bonds = lammps_file.bonds.to_indices()
        first_pos_bond = pos[(bonds[:,0])]
        second_pos_bond = pos[(bonds[:,1])]

        bonds = np.hstack((first_pos_bond, second_pos_bond))
        bonds = bonds.reshape(no_bonds, 2, 3)
        bond_actor = actor.streamtube(bonds, colors=(0.8275, 0.8275, 0.8275), linewidth=0.2)
        window.ren.add(bond_actor)
    avg = np.average(pos, axis=0)

    radii = 0.5* np.ones(no_atoms)
    unique_types = np.unique(lammps_file.atoms.types)
    colors_unique_types = np.random.rand(len(unique_types), 3)

    for i, typ in enumerate(unique_types):
        # print(i, typ)
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

    window.ren.UseImageBasedLightingOn()
    cube_path = 'c:/Users/nasim/Devel/furious-atoms/skybox0'
    surface = 'boy'
    if not os.path.isdir(cube_path):
        print('This path does not exist:', cube_path)
        return
    surface = surface.lower()
    cubemap = ReadCubeMap(cube_path, '/', '.jpg', 0)
    window.ren.SetEnvironmentTexture(cubemap)
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
    window.ren.add(axes_actor)
    window.ren.add(sphere_actor)
    window.ren.add(box_actor)
    window.ren.add(line_actor)
    window.ren.set_camera(focal_point=avg)

    print('Processed done!')
global cnt, enable_timer
cnt = 0
enable_timer = False


def timer_callback():
    global enable_timer
    if enable_timer is False:
        return

    global sphere_actor, cnt, initial_vertices
    global no_vertices_per_sphere, window, no_atoms, no_bonds, lammps_file, all_vertices, box, n_frames

    if cnt == n_frames:
        return

    if no_bonds==0:
        pos_R = lammps_file.trajectory[cnt].positions.copy().astype('f8')
        pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, box, backend='serial')

        all_vertices[:] = initial_vertices + \
            np.repeat(pos, no_vertices_per_sphere, axis=0)

        utils.update_actor(sphere_actor)


    window.vtkWidget.GetRenderWindow().Render()
    cnt += 1


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

        # self.ren.background((0, 0, 0))
        # self.ren.background((.3, .4, .5))
        self.ren.background((1.0, 1.0, 1.0))

        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.timer = QTimer()
        self.timer.timeout.connect(timer_callback)
        duration = 200
        self.timer.start(duration)

        self.resize(2000, 1000)
        self.iren.Initialize()
        self.show()

    def openColorDialog(self):
        color = QColorDialog.getColor()

    def open(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')#, filter = "*.lammp*")

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

        self.zoomOutAct = QAction("Zoom &Out (25%)", self, shortcut="Ctrl+-",
                enabled=False, triggered=self.zoomOut)

        self.normalSizeAct = QAction("&Normal Size", self, shortcut="Ctrl+S",
                enabled=False, triggered=self.normalSize)

        self.fitToWindowAct = QAction("&Fit to Window", self, enabled=False,
                checkable=True, shortcut="Ctrl+F", triggered=self.fitToWindow)

        self.aboutAct = QAction("&About", self, triggered=self.about)


    def createMenus(self):
        self.fileMenu = QMenu("&File", self)
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
    sys.exit(app.exec_())
