# Standard package
from json import load
import os
import fnmatch

# Local package
from furiousatoms import io
from furiousatoms.objects import mobius, box_edges
from fury import disable_warnings

disable_warnings()

# 3rd Party package
import vtk
import numpy as np
from fury import window, actor, utils
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import MDAnalysis


class FuriousAtomsApp(QtWidgets.QMainWindow):
    """

    """

    def __init__(self, app_path=None, parent=None):
        super(FuriousAtomsApp, self).__init__(parent)

        self.ui = io.load_ui_widget("view.ui")
        self.setCentralWidget(self.ui)

        self.app_path = app_path or io.application_path()
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.qvtkwidget = QVTKRenderWindowInteractor(parent=self.ui.view_frame,
                                                     rw=self.showm.window,
                                                     iren=self.showm.iren)
        self.init_settings()
        self.init_interface()
        self.create_language_menu()
        self.create_connections()

    def init_settings(self):
        pass

    def init_interface(self):
        self.ui.widget_Animation.setVisible(self.ui.button_animation.isChecked())
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.qvtkwidget, stretch=1)
        self.ui.view_frame.setLayout(self.v_layout)
        self.scene.background((1, 1, 1))
        self.showm.initialize()

    def create_connections(self):
        self.ui.actionLoad_file.triggered.connect(self.open)
        self.ui.actionSave_file.triggered.connect(self.save)
        self.ui.actionExit.triggered.connect(self.quit_fired)
        self.ui.button_animation.toggled.connect(self.ui.widget_Animation.setVisible)
        # self.actionBond.triggered.connect(delete_bonds)
        # self.actionParticle.triggered.connect(delete_particles)

    def open(self):
        # , filter = "*.lammp*")
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.tr('Load'))

        if not fname:
            return

        self.process_load_file(fname)
        global enable_timer
        enable_timer = True

    def save(self):
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, self.tr('Save'))
        print('saving {}'.format(fname))

    def update_bonds_ui(self, load_file, no_bonds, box_shape,
                        no_unique_types_bond):
        if len(box_shape) != 3:
            msg = "Wrong Box Shape size. Should be 3"
            raise ValueError(msg)
        n_frames = load_file.trajectory.n_frames
        no_atoms = len(load_file.atoms)
        no_unique_types = len(np.unique(load_file.atoms.types))
        # TODO: Update widget name
        # self.ui.dcfdLineEdit.insert(str(no_bonds))
        # self.ui.dcfdLineEdit_2.insert(str(n_frames))
        # self.ui.particleLineEdit.insert(str(no_atoms))
        # self.ui.atomTypeLineEdit.insert(str(box_shape[0]))
        # self.ui.atomTypeLineEdit_2.insert(str(box_shape[1]))
        # self.ui.atomTypeLineEdit_3.insert(str(box_shape[2]))
        # self.ui.particleLineEdit_2.insert(str(no_unique_types))
        # self.ui.particleLineEdit_3.insert(str(no_unique_types_bond))

    def process_load_file(self, fname):
        load_file, no_bonds = io.load_files(fname)

        no_atoms = len(load_file.atoms)
        box = load_file.trajectory.ts.dimensions
        box_lx = load_file.trajectory.ts.dimensions[0]
        box_ly = load_file.trajectory.ts.dimensions[1]
        box_lz = load_file.trajectory.ts.dimensions[2]

        box_centers = np.array([[0, 0, 25]])
        box_directions = np.array([[0, 1, 0]])
        box_colors = np.array([[255, 255, 255]])
        box_actor = actor.box(box_centers, box_directions, box_colors,
                              scales=(box_lx, box_ly, box_lz))
        box_actor.GetProperty().SetOpacity(0.)
        box_actor.GetProperty().SetRepresentationToWireframe()
        box_actor.GetProperty().SetLineWidth(10)

        line_actor, _ = box_edges(box_lx, box_ly, box_lz, colors=(0, 0, 0),
                                  linewidth=1, fake_tube=True)

        # np.unique(load_file.atoms.types)
        atom_type = load_file.atoms.types

        colors = np.ones((no_atoms, 4))
        pos = load_file.trajectory[0].positions.copy().astype('f8')
        #######################
        # s = load_file.universe.atoms[:]
        # t = load_file.universe.atoms[0:1]
        # b = s.difference(t)
        # load_file = MDAnalysis.core.groups.AtomGroup(b)
        ########################

        unique_types_bond = 0
        if no_bonds == 0:
            pos_R = load_file.trajectory[0].positions.copy().astype('f8')
            pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, box,
                                                          backend='serial')

        if no_bonds > 0:
            pos = load_file.trajectory[0].positions.copy().astype('f8')
            # if MainWindow.CheckBox.isChecked() == True:
            # load_file.delete_bonds(load_file.bonds[first_index_bond:first_index_bond+1])
            # load_file.delete_bonds(load_file.bonds.to_indices())
            bonds = load_file.bonds.to_indices()
            no_bonds = len(load_file.bonds)
            first_pos_bond = pos[(bonds[:, 0])]
            second_pos_bond = pos[(bonds[:, 1])]
            bonds = np.hstack((first_pos_bond, second_pos_bond))
            bonds = bonds.reshape((no_bonds), 2, 3)
            bond_colors = (0.8275, 0.8275, 0.8275, 1)
            bond_actor = actor.streamtube(bonds, bond_colors, linewidth=0.2,
                                          opacity=0.995)
            vcolors_bond = utils.colors_from_actor(bond_actor, 'colors')
            colors_backup_bond = vcolors_bond.copy()
            all_vertices_bonds = utils.vertices_from_actor(bond_actor)
            no_vertices_per_bond = len(all_vertices_bonds) / no_bonds
            # initial_vertices_bonds = all_vertices_bonds.copy() - \
            #   np.repeat(bonds, no_vertices_per_bonds, axis=0)

            self.scene.add(bond_actor)
            unique_types_bond = np.unique(load_file.bonds.types)
            str_no_unique_types_bond = str(len(unique_types_bond))

        avg = np.average(pos, axis=0)

        radii = 0.5 * np.ones(no_atoms)
        unique_types = np.unique(load_file.atoms.types)
        colors_unique_types = np.random.rand(len(unique_types), 4)
        colors_unique_types[:, 3] = 1

        for i, typ in enumerate(unique_types):
            colors[atom_type == typ] = colors_unique_types[i]

        selected = np.zeros(no_atoms, dtype=np.bool)
        selected_bond = np.zeros(no_bonds, dtype=np.bool)
        sphere_actor = actor.sphere(centers=pos,
                                    colors=colors,
                                    radii=radii, theta=32, phi=32)
        all_vertices = utils.vertices_from_actor(sphere_actor)
        no_vertices_per_sphere = len(all_vertices) / no_atoms
        initial_vertices = all_vertices.copy() - np.repeat(pos, no_vertices_per_sphere, axis=0)

        self.update_bonds_ui(load_file, no_bonds,
                             box_shape=[box_lx, box_ly, box_lz],
                             no_unique_types_bond=unique_types_bond)

        vcolors = utils.colors_from_actor(sphere_actor, 'colors')
        colors_backup = vcolors.copy()

        sphere_actor.GetProperty().SetInterpolationToPBR()
        # Lets use a smooth metallic surface
        # Build the pipeline
        # mapper = vtk.vtkPolyDataMapper()
        # source = mobius()
        # utils.set_input(mapper, source)
        # # mapper.SetInputData(source)
        # vtk.vtkActor().SetMapper(mapper)

        self.scene.UseImageBasedLightingOn()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        cube_path = os.path.join(dir_path, 'skybox0')
        # surface = 'boy'
        if not os.path.isdir(cube_path):
            print('This path does not exist:', cube_path)
            return
        # surface = surface.lower()
        cubemap = io.read_cubemap(cube_path, '/', '.jpg', 0)
        self.scene.SetEnvironmentTexture(cubemap)
        sphere_actor.GetProperty().SetInterpolationToPBR()

        # configure the basic properties
        metallicCoefficient = 0.5  # 0
        roughnessCoefficient = 0.1  # 1
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
        # glrenderer = vtk.vtkOpenGLRenderer.SafeDownCast(self.scene)
        # glrenderer.SetPass(ssao)
        axes_actor = actor.axes(scale=(1, 1, 1), colorx=(1, 0, 0),
                                colory=(0, 1, 0), colorz=(0, 0, 1), opacity=1)

        self.scene.add(axes_actor)
        self.scene.add(sphere_actor)
        self.scene.add(box_actor)
        self.scene.add(line_actor)

        self.scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                              view_up=(0, 1, 0))

    def quit_fired(self):
        reply = QtWidgets.QMessageBox.warning(
            self,
            self.tr("Quit"),
            self.tr("Do you want to quit the software?"),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            QtCore.QCoreApplication.instance().quit()

    def retranslateUi(self):
        self.protocole_widget.retranslate_ui()
        self.data_selection_widget.retranslate_ui()
        self.processing_widget.retranslate_ui()
        self.hcs_widget.retranslate_ui()

    def closeEvent(self, event):
        reply = QtWidgets.QMessageBox.warning(
            self,
            self.tr("Quit"),
            self.tr("Do you want to quit the software?"),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def changeEvent(self, event):
        if event:
            if event.type() == QtCore.QEvent.LanguageChange:
                # self event is send if a translator is loaded
                # Todo: translation. see http://stackoverflow.com/questions/39734775/python-switch-translations-dynamically-in-a-pyqt4-application
                # self.ui.retranslateUi(self)
                self.retranslateUi()
            if event.type() == QtCore.QEvent.LocaleChange:
                # Todo: Change this part
                locale = QtCore.QLocale.system().name()
                locale.truncate(locale.lastIndexOf('_'))
                self.loadLanguage(locale)

        QtWidgets.QWidget.changeEvent(self, event)

    def create_language_menu(self):
        self.m_curr_lang = ""
        self.m_translator = QtCore.QTranslator(self)

        #  format systems language
        default_locale = QtCore.QLocale.system().name()  # e.g. "de_DE"
        default_locale = default_locale.split('_')[0]  # e.g. "de"

        # Get languages
        self.m_lang_path = io.get_languages_path()
        f_names = fnmatch.filter(os.listdir(self.m_lang_path), "*.qm")

        for f in f_names:
            #  get locale extracted by filename.  fr or fr_FR
            locale = os.path.splitext(f)[0]

            lang = QtCore.QLocale.languageToString(QtCore.QLocale(locale).language())
            ico = QtGui.QIcon("{0}/{1}.png".format(self.m_lang_path, locale))

            action = QtWidgets.QAction(ico, lang, self)
            action.setCheckable(True)
            action.setData(locale)

            self.langGroup.addAction(action)
            self.languageMenu.addAction(action)

            #  set default translators and language checked
            if default_locale == locale:
                action.setChecked(True)

            # Called every time, when a menu entry of the language menu
            # is called

    def slot_language_changed(self, action):
        if not action:
            return

        # load the language dependant on the action content
        self.loadLanguage(action.data())
        self.setWindowIcon(action.icon())

    def switch_translator(self, translator, f_name):
        """ Update the tranlator object

        Parameters
        -----------
        translator : Translator instance
            object which permit to change text
        f_name : str
            language file name
        """
        # Remove the old translator
        QtCore.QCoreApplication.instance().removeTranslator(translator)
        # Load the new translator
        if translator.load(f_name):
            QtCore.QCoreApplication.instance().installTranslator(translator)

    def load_language(self, r_language):
        """ Load the selected language

        Parameters
        -----------
        r_language: str
            language acronym
        """
        if self.m_curr_lang == r_language:
            return

        self.m_curr_lang = r_language
        locale = QtCore.QLocale(self.m_curr_lang)  # de-fr-en
        QtCore.QLocale.setDefault(locale)
        language_name = QtCore.QLocale.languageToString(locale.language())
        self.switch_translator(
            self.m_translator,
            QtCore.QDir(self.m_lang_path).filePath("{0}.qm").format(self.m_curr_lang))
        print(self.tr("Current Language changed to {0}").format(language_name))
