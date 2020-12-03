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
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import MDAnalysis
from furiousatoms.sharedmem import SharedMemory

SM = SharedMemory()

class FuriousAtomsApp(QtWidgets.QMainWindow):
    """ Main Furious Atoms Class
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
        self.ui.actionParticle.triggered.connect(self.delete_particles)
        self.ui.actionBond.triggered.connect(self.delete_bonds)
        self.ui.Button_bondcolor.clicked.connect(self.openColorDialog_bond)
        self.ui.Button_particlecolor.clicked.connect(self.openColorDialog_particle)
        # self.ui.LineEdit_atom_radius.textChanged.connect(self.change_particle_size)
        # selected_value = self.ui.LineEdit_atom_radius.text()


    # def change_particle_size(self, selected_value):

    #     # object_indices_particles = np.where(SM.selected_particle == True)[0]
    #     SM.radii_spheres = float(selected_value)
    #     SM.sphere_actor = actor.sphere(centers=SM.pos,
    #                                 colors=SM.colors_backup_particles,
    #                                 radii=SM.radii_spheres, theta=32, phi=32)
    #     SM.all_vertices_particles = utils.vertices_from_actor(SM.sphere_actor)

    #     # SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
    #     # for object_index in object_indices_particles:
    #     #     SM.all_vertices_particles[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = utils.vertices_from_actor(SM.sphere_actor)
    #     utils.update_actor(SM.sphere_actor)
    #     SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
    #     self.qvtkwidget.GetRenderWindow().Render()



    # SM.radii_spheres = np.ones((SM.no_atoms, 3))
    # SM.radii_unique_types = np.random.rand(len(unique_types), 3)
    # for i, typ in enumerate(unique_types):
    #     SM.radii_spheres[atom_type == typ] = SM.radii_unique_types[i]



    def check_simulationcell(self, state):
        if (state == QtCore.Qt.Checked):
            SM.line_actor.VisibilityOn()
        else:
            SM.line_actor.VisibilityOff()
        utils.update_actor(SM.line_actor)
        SM.line_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()


    def check_particles(self, state):
        if (state == QtCore.Qt.Checked):
            SM.sphere_actor.VisibilityOn()
        else:
            SM.sphere_actor.VisibilityOff()
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print('All particles are deleted')

    def check_bonds(self, state):
        if (state == QtCore.Qt.Checked):
            SM.bond_actor.VisibilityOn()
        else:
            SM.bond_actor.VisibilityOff()
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print('All bonds are deleted')

    def open(self):
        # , filter = "*.lammp*")
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.tr('Load'))

        if not fname:
            return

        self.process_load_file(fname)
        SM.enable_timer = True


    def save(self):
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, self.tr('Save'))
        print('saving {}'.format(fname))


    def delete_particles(self):
        object_indices_particles = np.where(SM.selected_particle == True)[0]
        print(object_indices_particles)
        SM.particle_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        for object_index in object_indices_particles:
            SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.particle_color_add
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print('The particle is deleted')


    def delete_bonds(self):
        object_indices_bonds = np.where(SM.selected_bond == True)[0]
        SM.bond_color_add = np.array([255, 0, 0, 0], dtype='uint8')
        SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
        for object_index_bond in object_indices_bonds:
            SM.vcolors_bond[object_index_bond * SM.sec_bond: object_index_bond * SM.sec_bond + SM.sec_bond] = SM.bond_color_add
        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()
        print('The bond is deleted')


    def openColorDialog_particle(self):
        selected_color_particle = QtWidgets.QColorDialog.getColor()
        select_all_particles = np.zeros(SM.no_atoms, dtype=np.bool)
        object_indices_particles = np.where(select_all_particles == False)[0]
        for object_index in object_indices_particles:
            SM.particle_color_add = selected_color_particle.getRgb()
            SM.vcolors_particle[object_index * SM.sec_particle: object_index * SM.sec_particle + SM.sec_particle] = SM.particle_color_add
        utils.update_actor(SM.sphere_actor)
        SM.sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()



    def openColorDialog_bond(self):
        selected_color_bond = QtWidgets.QColorDialog.getColor()
        select_all_bonds = np.zeros(SM.no_bonds, dtype=np.bool)
        object_indices_bonds = np.where(select_all_bonds == False)[0]
        for object_index in object_indices_bonds:
            SM.bond_color_add = selected_color_bond.getRgb()
            SM.vcolors_bond[object_index * SM.sec_bond: object_index * SM.sec_bond + SM.sec_bond] = SM.bond_color_add

        utils.update_actor(SM.bond_actor)
        SM.bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
        self.qvtkwidget.GetRenderWindow().Render()

    def update_bonds_ui(self, load_file, no_bonds, box_shape,
                        no_unique_types_bond):
        if len(box_shape) != 3:
            msg = "Wrong Box Shape size. Should be 3"
            raise ValueError(msg)
        SM.n_frames = SM.load_file.trajectory.n_frames
        SM.no_atoms = len(SM.load_file.atoms)
        SM.no_unique_types_particles = len(np.unique(SM.load_file.atoms.types))
        self.ui.Edit_num_of_particletypes.insert(str(SM.no_unique_types_particles))
        # TODO: Update widget name
        # self.ui.dcfdLineEdit.insert(str(no_bonds))
        # self.ui.dcfdLineEdit_2.insert(str(n_frames))
        # self.ui.particleLineEdit.insert(str(no_atoms))
        # self.ui.atomTypeLineEdit.insert(str(box_shape[0]))
        # self.ui.atomTypeLineEdit_2.insert(str(box_shape[1]))
        # self.ui.atomTypeLineEdit_3.insert(str(box_shape[2]))
        # self.ui.particleLineEdit_2.insert(str(no_unique_types_particles))
        # self.ui.particleLineEdit_3.insert(str(no_unique_types_bond))

    def process_load_file(self, fname):

        SM.load_file, SM.no_bonds = io.load_files(fname)
        file_directory, extension = os.path.splitext(fname)
        self.ui.Edit_directory.insert(str(file_directory))
        file_name = os.path.basename(fname)
        self.ui.Edit_currentfile.insert(str(file_name))
        self.ui.Edit_fileformat.insert((str(extension)).upper())
        SM.no_atoms = len(SM.load_file.atoms)
        SM.n_frames = SM.load_file.trajectory.n_frames
        self.ui.Edit_number_of_frames.insert(str(SM.n_frames))
        if SM.no_atoms > 0:
            self.ui.Box_particles.setChecked(True)
            self.ui.Box_particles.stateChanged.connect(self.check_particles)

        SM.box = SM.load_file.trajectory.ts.dimensions
        box_lx = SM.load_file.trajectory.ts.dimensions[0]
        box_ly = SM.load_file.trajectory.ts.dimensions[1]
        box_lz = SM.load_file.trajectory.ts.dimensions[2]
        self.ui.Edit_widthX.insert(str(box_lx))
        self.ui.Edit_lengthY.insert(str(box_ly))
        self.ui.Edit_hightZ.insert(str(box_lz))
        SM.no_unique_types_particles = len(np.unique(SM.load_file.atoms.types))
        print(SM.no_unique_types_particles)
        if box_lx > 0:
            self.ui.Box_simulationcell.setChecked(True)
            self.ui.Box_simulationcell.stateChanged.connect(self.check_simulationcell)

        box_centers = np.array([[0, 0, 25]])
        box_directions = np.array([[0, 1, 0]])
        SM.box_colors = np.array([[255, 255, 255, 255]])
        box_actor = actor.box(box_centers, box_directions, SM.box_colors,
                              scales=(box_lx, box_ly, box_lz))
        # SM.vcolors_box = utils.colors_from_actor(SM.box_actor, 'colors')
        SM.colors_backup_box = SM.vcolors_box
        box_actor.GetProperty().SetOpacity(0.)
        box_actor.GetProperty().SetRepresentationToWireframe()
        box_actor.GetProperty().SetLineWidth(10)

        SM.line_actor, _ = box_edges(box_lx, box_ly, box_lz, colors=(0, 0, 0),
                                  linewidth=1, fake_tube=True)

        # np.unique(load_file.atoms.types)
        atom_type = SM.load_file.atoms.types
        pos = SM.load_file.trajectory[0].positions.copy().astype('f8')
        #######################
        # s = load_file.universe.atoms[:]
        # t = load_file.universe.atoms[0:1]
        # b = s.difference(t)
        # load_file = MDAnalysis.core.groups.AtomGroup(b)
        ########################

        unique_types_bond = 0
        if SM.no_bonds == 0:
            pos_R = SM.load_file.trajectory[0].positions.copy().astype('f8')
            pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, SM.box,
                                                          backend='serial')

        if SM.no_bonds > 0:
            self.ui.Box_bonds.setChecked(True)
            pos = SM.load_file.trajectory[0].positions.copy().astype('f8')
            bonds = SM.load_file.bonds.to_indices()
            SM.no_bonds = len(SM.load_file.bonds)
            first_pos_bond = pos[(bonds[:, 0])]
            second_pos_bond = pos[(bonds[:, 1])]
            bonds = np.hstack((first_pos_bond, second_pos_bond))
            bonds = bonds.reshape((SM.no_bonds), 2, 3)
            SM.bond_colors = (0.8275, 0.8275, 0.8275, 1)
            SM.bond_actor = actor.streamtube(bonds, SM.bond_colors, linewidth=0.2)
            SM.vcolors_bond = utils.colors_from_actor(SM.bond_actor, 'colors')
            SM.colors_backup_bond = SM.vcolors_bond.copy()
            SM.all_vertices_bonds = utils.vertices_from_actor(SM.bond_actor)
            SM.no_vertices_per_bond = len(SM.all_vertices_bonds) / SM.no_bonds
            # SM.initial_vertices_bonds = SM.all_vertices_bonds.copy() - np.repeat(pos, SM.no_vertices_per_bond, axis=0)

            vertices_bonds = utils.vertices_from_actor(SM.bond_actor)
            SM.no_vertices_all_bonds = vertices_bonds.shape[0]
            SM.sec_bond = np.int(SM.no_vertices_all_bonds / SM.no_bonds)
            self.update_bonds_ui(SM.load_file, SM.no_bonds,
                                box_shape=[box_lx, box_ly, box_lz],
                                no_unique_types_bond=unique_types_bond)

            self.scene.add(SM.bond_actor)
            self.ui.Box_bonds.stateChanged.connect(self.check_bonds)
            unique_types_bond = np.unique(SM.load_file.bonds.types)
            str_no_unique_types_bond = str(len(unique_types_bond))
            SM.bond_actor.AddObserver("LeftButtonPressEvent", self.left_button_press_bond_callback)

        avg = np.average(pos, axis=0)

        colors = np.ones((SM.no_atoms, 4))
        unique_types = np.unique(SM.load_file.atoms.types)
        SM.colors_unique_types = np.random.rand(len(unique_types), 4)
        SM.colors_unique_types[:, 3] = 1

        for i, typ in enumerate(unique_types):
            colors[atom_type == typ] = SM.colors_unique_types[i]

        # radii = 0.5 * np.ones(SM.no_atoms)
        # SM.radii_spheres =  0.5 * np.ones(SM.no_atoms)
        SM.radii_spheres = np.ones((SM.no_atoms, 3))
        SM.radii_unique_types = np.random.rand(len(unique_types), 3)
        for i, typ in enumerate(unique_types):
            SM.radii_spheres[atom_type == typ] = SM.radii_unique_types[i]

        SM.selected_particle = np.zeros(SM.no_atoms, dtype=np.bool)
        SM.selected_bond = np.zeros(SM.no_bonds, dtype=np.bool)
        SM.sphere_actor = actor.sphere(centers=pos,
                                    colors=colors,
                                    radii=SM.radii_spheres, theta=32, phi=32)
        print("radius value is: ", SM.radii_spheres)
        SM.all_vertices_particles = utils.vertices_from_actor(SM.sphere_actor)
        SM.no_vertices_per_particle = len(SM.all_vertices_particles) / SM.no_atoms
        SM.initial_vertices_particles = SM.all_vertices_particles.copy() - np.repeat(pos, SM.no_vertices_per_particle, axis=0)
        vertices_particle = utils.vertices_from_actor(SM.sphere_actor)
        SM.no_vertices_all_particles = vertices_particle.shape[0]
        SM.sec_particle = np.int(SM.no_vertices_all_particles / SM.no_atoms)
        SM.vcolors_particle = utils.colors_from_actor(SM.sphere_actor, 'colors')
        SM.colors_backup_particles = SM.vcolors_particle.copy()

        SM.sphere_actor.GetProperty().SetInterpolationToPBR()
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
        SM.sphere_actor.GetProperty().SetInterpolationToPBR()

        # configure the basic properties
        metallicCoefficient = 0.5  # 0
        roughnessCoefficient = 0.1  # 1
        colors_sky = vtk.vtkNamedColors()
        SM.sphere_actor.GetProperty().SetColor(colors_sky.GetColor3d('White'))
        SM.sphere_actor.GetProperty().SetMetallic(metallicCoefficient)
        SM.sphere_actor.GetProperty().SetRoughness(roughnessCoefficient)
        axes_actor = actor.axes(scale=(1, 1, 1), colorx=(1, 0, 0),
                                colory=(0, 1, 0), colorz=(0, 0, 1), opacity=1)

        self.scene.add(axes_actor)
        self.scene.add(SM.sphere_actor)
        self.scene.add(box_actor)
        self.scene.add(SM.line_actor)
        self.scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                              view_up=(0, 1, 0))

        SM.sphere_actor.AddObserver("LeftButtonPressEvent", self.left_button_press_particle_callback)


        if self.ui.button_animation.isChecked():
            self.timer = QtCore.QTimer()
            duration = 200
            self.timer.start(duration)
            self.timer.timeout.connect(self.timer_callback)
            self.showm.add_timer_callback(True, 1, self.timer_callback)

        self.pickm = pick.PickingManager()

    SM.enable_timer = False

    def timer_callback(self):
        # print('number of frames: ' ,SM.n_frames)
        if SM.enable_timer is False:
            return
        if SM.cnt == SM.n_frames:
            return
        if SM.n_frames > 1:
            # self.ui.button_animation.setChecked(True)
            pos_R = SM.load_file.trajectory[SM.cnt].positions.copy().astype('f8')
            pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, SM.box, backend='serial')
            SM.all_vertices_particles[:] = SM.initial_vertices_particles + \
                np.repeat(pos, SM.no_vertices_per_particle, axis=0)
            utils.update_actor(SM.sphere_actor)

            # SM.all_vertices_bonds[:] = SM.initial_vertices_bonds + \
            #     np.repeat(pos, SM.no_vertices_per_bond, axis=0)
            # utils.update_actor(SM.bond_actor)
            # self.timer = QtCore.QTimer()
            # self.Timer.setRange(SM.cnt, SM.n_frames)
            # if (self.ui.Button_play.isClicked()):
        # self.ui.Edit_framenumber.insert(str(SM.cnt))
        self.ui.horizontalSlider_animation.setMinimum(0)
        self.ui.horizontalSlider_animation.setMaximum(SM.n_frames)
        self.ui.horizontalSlider_animation.setSingleStep(1)
        self.ui.horizontalSlider_animation.setValue(SM.cnt)
        self.ui.horizontalSlider_animation.setTickInterval(1)
# Timer_animation
        self.qvtkwidget.GetRenderWindow().Render()
        SM.cnt += 1


    def left_button_press_particle_callback(self, obj, event):

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
        event_pos = self.pickm.event_position(iren=self.showm.iren)
        picked_info = self.pickm.pick(event_pos, self.showm.scene)
        vertex_index_particle = picked_info['vertex']
        vertices_bonds = utils.vertices_from_actor(obj)
        SM.no_vertices_all_bonds = vertices_bonds.shape[0]
        object_index_bond = np.int(np.floor((vertex_index_particle / SM.no_vertices_all_bonds) * SM.no_bonds))
        # Find how many vertices correspond to each object
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
