from MDAnalysis import *
import MDAnalysis as mda
import MDAnalysis.analysis.align
import numpy as np
import sys
import io
from furiousatoms import io
import vtk
import numpy as np
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from furiousatoms.sharedmem import SharedMemory
from furiousatoms.periodic_table import Ui_periodic
from furiousatoms.warning_message import Ui_warning_charge_salt_1, Ui_warning_oppos_charge
import sys
import math
SM = SharedMemory()
class Ui_electrolyte(QtWidgets.QMainWindow): #QWidget
    """ Ui_electrolyte class creates a widget for building electrolyte
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_electrolyte, self).__init__(parent)
        self.elect = io.load_ui_widget("electrolyte.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.elect)
        self.setCentralWidget(self.elect)
        self.setLayout(self.v_layout)
        self.resize(631, 380)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):

        # connect counterion, ions and anions icons to periodic table to get information of ion (type, valency and mass):
        self.elect.pushButton_type_counterion.clicked.connect(self.info_counterion)
        self.elect.lineEdit_valency_counterion.textChanged[str].connect(self.check_neutrality)
        self.elect.lineEdit_mass_counterion.textChanged[str].connect(self.check_neutrality)
        self.elect.pushButton_cation_salt_1.clicked.connect(self.info_cation_salt_1)
        self.elect.pushButton_anion_salt_1.clicked.connect(self.info_anion_salt_1)
        self.elect.pushButton_cation_salt_2.clicked.connect(self.info_cation_salt_2)
        self.elect.pushButton_anion_salt_2.clicked.connect(self.info_anion_salt_2)
        self.elect.pushButton_cation_salt_3.clicked.connect(self.info_cation_salt_3)
        self.elect.pushButton_anion_salt_3.clicked.connect(self.info_anion_salt_3)
        self.elect.pushButton_cation_salt_4.clicked.connect(self.info_cation_salt_4)
        self.elect.pushButton_anion_salt_4.clicked.connect(self.info_anion_salt_4)
        self.elect.SpinBox_con_cation_salt_1.valueChanged.connect(self.add_salt)
        self.elect.SpinBox_con_cation_salt_2.valueChanged.connect(self.add_salt)
        self.elect.SpinBox_con_cation_salt_3.valueChanged.connect(self.add_salt)
        self.elect.SpinBox_con_cation_salt_4.valueChanged.connect(self.add_salt)
        self.elect.pushButton_build_electrolyte.clicked.connect(self.build_electrolyte)
        self.elect.pushButton_build_electrolyte.clicked.connect(lambda:self.close())
        self.elect.lineEdit_type_counterion.textChanged[str].connect(self.dummy_wall)
        self.elect.lineEdit_valency_counterion.textChanged[str].connect(self.dummy_wall)
        self.elect.lineEdit_mass_counterion.textChanged[str].connect(self.dummy_wall)
        self.elect.SpinBox_lx.valueChanged.connect(self.dummy_wall)
        self.elect.SpinBox_lz.valueChanged.connect(self.dummy_wall)
        self.elect.SpinBox_lx.valueChanged.connect(self.add_water)
        self.elect.SpinBox_lz.valueChanged.connect(self.add_water)
        self.elect.lineEdit_surface_charge.textChanged[str].connect(self.dummy_wall)
        self.elect.SpinBox_lx.valueChanged.connect(self.add_salt)
        self.elect.SpinBox_lz.valueChanged.connect(self.add_salt)
        self.elect.lineEdit_type_cation_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_cation_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_cation_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_type_cation_salt_2.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_cation_salt_2.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_cation_salt_2.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_type_cation_salt_3.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_cation_salt_3.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_cation_salt_3.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_type_cation_salt_4.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_cation_salt_4.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_cation_salt_4.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_type_anion_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_anion_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_anion_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_type_anion_salt_2.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_anion_salt_2.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_anion_salt_2.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_type_anion_salt_3.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_anion_salt_3.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_anion_salt_3.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_type_anion_salt_4.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_anion_salt_4.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_anion_salt_4.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_water_diameter.textChanged[str].connect(self.add_water)
        self.elect.lineEdit_space_diameter.textChanged[str].connect(self.add_water)
        self.elect.lineEdit_valency_anion_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.comboBox_water_model.currentTextChanged.connect(self.add_water)
        self.elect.checkBox_dummy_wall.stateChanged.connect(self.check_box_dummy)
        self.elect.checkBox_add_water.stateChanged.connect(self.check_box_add_water)
        self.elect.frame_dummy_wall.setEnabled(False)
        self.elect.frame_dummy_wall_charge.setEnabled(False)
        self.elect.frame_water_properties.setEnabled(False)
        self.elect.frame_water_model.setEnabled(False)
        self.elect.frame_dummy_wall_properties_final.setEnabled(False)
        self.elect.frame_cation_properties_final.setEnabled(False)
        self.elect.water_properties_final.setEnabled(False)
        SM.water_diameter = self.elect.lineEdit_water_diameter.setText(str(3.16555789))
        spacing_dia = self.elect.lineEdit_space_diameter.setText(str(1.0))
        SM.spacing = self.elect.lineEdit_space_bet_water.setText(str(3.16555789))
        SM.charge_density = self.elect.lineEdit_surface_charge.setText(str(0.0))
    def show_periodic_table(self):
        Ui_periodic.pt = Ui_periodic()
        Ui_periodic.pt.win = self
        Ui_periodic.pt.show()
    def info_counterion(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_counterion
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_counterion
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_counterion
    def info_cation_salt_1(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_cation_salt_1
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_cation_salt_1
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_cation_salt_1
    def info_anion_salt_1(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_anion_salt_1
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_anion_salt_1
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_anion_salt_1
    def info_cation_salt_2(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_cation_salt_2
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_cation_salt_2
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_cation_salt_2
    def info_anion_salt_2(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_anion_salt_2
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_anion_salt_2
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_anion_salt_2
    def info_cation_salt_3(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_cation_salt_3
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_cation_salt_3
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_cation_salt_3
    def info_anion_salt_3(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_anion_salt_3
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_anion_salt_3
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_anion_salt_3
    def info_cation_salt_4(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_cation_salt_4
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_cation_salt_4
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_cation_salt_4
    def info_anion_salt_4(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_anion_salt_4
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_anion_salt_4
        Ui_periodic.pt.current_edit_mass = self.elect.lineEdit_mass_anion_salt_4

    def show_warning_salt_1(self):
        Ui_warning_charge_salt_1.msg = Ui_warning_charge_salt_1()
        Ui_warning_charge_salt_1.msg.win = self
        Ui_warning_charge_salt_1.msg.show()
    def show_warning_opposite_charge(self):
        Ui_warning_oppos_charge.msg_coun = Ui_warning_oppos_charge()
        Ui_warning_oppos_charge.msg_coun.win = self
        Ui_warning_oppos_charge.msg_coun.show()

    def add_water(self):
        SM.water_diameter = self.elect.lineEdit_water_diameter.text()
        SM.spacing_dia = self.elect.lineEdit_space_diameter.text()
        try:
            SM.spacing_dia = float( SM.spacing_dia)
        except ValueError:
            SM.spacing_dia = 0
        SM.water_diameter = float( SM.water_diameter)
        SM.spacing =  SM.spacing_dia *  SM.water_diameter
        self.elect.lineEdit_space_bet_water.setText(str( SM.spacing))
        self.elect.lineEdit_oxygen_charge.text()
        self.elect.lineEdit_hydrogen_charge.text()
        try:
            water_amount = float( SM.box_lx/ SM.spacing *  SM.box_ly/ SM.spacing *  SM.box_lz/ SM.spacing)
            SM.total_water_inside = int(int( SM.box_lx/ SM.spacing) * int( SM.box_ly/ SM.spacing) * int( SM.box_lz/ SM.spacing))
            water_concentration = water_amount / (0.602214076 * ( SM.box_lx) * ( SM.box_ly) * ( SM.box_lz- SM.water_diameter) * 0.001)
        except ZeroDivisionError:
            SM.total_water_inside = water_amount = water_concentration = 0
        self.elect.lineEdit_total_num_water.setText(str(SM.total_water_inside))
        self.elect.lineEdit_total_num_water_final.setText(str(SM.total_water_inside))
        self.elect.lineEdit_water_con.setText(str(str("{:.5f}".format(water_concentration))))
        self.elect.lineEdit_water_con_final.setText(str("{:.5f}".format(water_concentration)))
# Choose the water model in the system:
        comboBox_water_model = self.elect.comboBox_water_model.currentText()
        if comboBox_water_model == 'SPC':
            charge_oxygen = self.elect.lineEdit_oxygen_charge.setText(str(-0.82))
            charge_hydrogen = self.elect.lineEdit_hydrogen_charge.setText(str(0.41))
        if comboBox_water_model == 'SPC/E':
            charge_oxygen = self.elect.lineEdit_oxygen_charge.setText(str(-0.8476))
            charge_hydrogen =self.elect.lineEdit_hydrogen_charge.setText(str(0.4238))
        if comboBox_water_model == 'TIPS':
            charge_oxygen = self.elect.lineEdit_oxygen_charge.setText(str(-0.80))
            charge_hydrogen =self.elect.lineEdit_hydrogen_charge.setText(str(0.40))
        if comboBox_water_model == 'TIP3P (Jorgensen)':
            charge_oxygen = self.elect.lineEdit_oxygen_charge.setText(str(-0.834))
            charge_hydrogen =self.elect.lineEdit_hydrogen_charge.setText(str(0.417))
        if comboBox_water_model == 'TIP3P (Price)':
            charge_oxygen = self.elect.lineEdit_oxygen_charge.setText(str(-0.830))
            charge_hydrogen = self.elect.lineEdit_hydrogen_charge.setText(str(0.415))

###################STEP 1###################
    def check_neutrality(self):
        try:
            SM.valency_counterion = float(self.elect.lineEdit_valency_counterion.text())
        except ValueError:
            SM.valency_counterion = 0
        if SM.charge_density != 0:
            SM.total_surface_charge = int( SM.charge_density *  SM.box_lx *  SM.box_ly / 16) # lx and ly are in A
        try:
            SM.counterions = int(2.0 * abs(SM.total_surface_charge/SM.valency_counterion))
        except ZeroDivisionError:
            SM.counterions = 0
        print("valency_counterion ", SM.valency_counterion)
        print("counterions ", SM.counterions)
        print("total surface charge ", SM.total_surface_charge)
        charge_system = SM.valency_counterion * SM.counterions + SM.total_surface_charge * 2.0
        print("total charge in the system ", charge_system)
        if (charge_system == 0):
            print("system is charge neutral")
        else:
            print("system is not electroneutral; aborting..." )
            self.show_warning_opposite_charge()
    def check_box_add_water(self, state):
        if (state == QtCore.Qt.Checked):
            self.elect.frame_water_properties.setEnabled(True)
            self.elect.frame_water_model.setEnabled(True)
            self.elect.water_properties_final.setEnabled(True)
        else:
            self.elect.frame_water_properties.setEnabled(False)
            self.elect.frame_water_model.setEnabled(False)
            self.elect.water_properties_final.setEnabled(False)
    def check_box_dummy(self, state):
        if (state == QtCore.Qt.Checked):
            self.elect.frame_dummy_wall.setEnabled(True)
            self.elect.frame_dummy_wall_charge.setEnabled(True)
            self.elect.frame_dummy_wall_properties_final.setEnabled(True)
            self.elect.frame_cation_properties_final.setEnabled(True)
        else:
            self.elect.frame_dummy_wall.setEnabled(False)
            self.elect.frame_dummy_wall_charge.setEnabled(False)
            self.elect.frame_dummy_wall_properties_final.setEnabled(False)
            self.elect.frame_cation_properties_final.setEnabled(False)

    def dummy_wall(self):
        SM.box_lx = SM.box_ly = self.elect.SpinBox_lx.value()
        SM.box_lz = self.elect.SpinBox_lz.value()
        self.elect.lineEdit_lx_final.setText(str(SM.box_lx))
        self.elect.lineEdit_ly_final.setText(str(SM.box_ly))
        self.elect.lineEdit_lz_final.setText(str(SM.box_lz))
        self.elect.lineEdit_ly.setText(str(SM.box_ly))
        SM.type_counter = str(self.elect.lineEdit_type_counterion.text())
        SM.mass_counter = str(self.elect.lineEdit_mass_counterion.text())
        try:
            SM.valency_counterion = float(self.elect.lineEdit_valency_counterion.text())
        except ValueError:
            SM.valency_counterion = 0.0
        width = 2  #A ; note CG width 3 A with lx = 15 nm
# Create two dummy walls on z direction:
        coordinates_wallR = []
        SM.box_lx = float(SM.box_lx)
        SM.box_ly = float(SM.box_ly)
        SM.box_lz = float(SM.box_lz)
        nx = int( SM.box_lx / width)
        ny = int( SM.box_ly / width)
        for i in range(ny):
            for j in range(nx):
                x = float(-0.5* SM.box_lx+(0.5)*width+i*width)
                y = float(-0.5* SM.box_ly+(0.5)*width+j*width)
                z = -0.5 *  SM.box_lz
                xyz_wallR = np.array([x, y, z])
                coordinates_wallR.append(xyz_wallR)
        SM.wallR = np.array(coordinates_wallR)

# Define surface charge density on the walls and number of counterions to keep the system electroneutral:
        SM.charge_density = self.elect.lineEdit_surface_charge.text()
        try:
            SM.charge_density = float(SM.charge_density)
        except ValueError:
            SM.charge_density = 0
            # self.elect.lineEdit_surface_charge_final.setText(str(SM.charge_density))
        if SM.charge_density == 0:
            SM.type_counter = None
            SM.mass_counterion = 0.0
            SM.valency_counterion = 0.0
            SM.counterions = 0
        # self.elect.lineEdit_type_counterion.setText(SM.type_counter)
        self.elect.lineEdit_type_counter_final.setText(SM.type_counter)
        # self.elect.lineEdit_valency_counterion.setText(str(SM.valency_counterion))
        # self.elect.lineEdit_mass_counterion.setText(str(SM.mass_counterion))
        unitcharge = 1.60217646 * pow(10.0 , -19)
        SM.total_surface_charge = int( SM.charge_density *  SM.box_lx *  SM.box_ly / 16) # lx and ly are in A
        surface_area =  SM.box_lx *  SM.box_ly * 0.01 * pow(10.0 , -18) # in unit of squared meter;
        try:
            SM.counterions = int(2.0 * abs(SM.total_surface_charge/SM.valency_counterion))
        except ZeroDivisionError:
            SM.counterions = 0
        if (len(SM.wallR) > 0):
            num_type_wall = 2
            SM.charge_meshpoint = (SM.total_surface_charge * 1.0) / len(SM.wallR)
        else:
            num_type_wall = 0
            SM.charge_meshpoint = 0

        self.elect.lineEdit_grid_size.setText(str(width))
        self.elect.lineEdit_mesh_charge.setText(str("{:.5f}".format(SM.charge_meshpoint)))
        # self.elect.lineEdit_mesh_charge.setText(str(SM.charge_meshpoint))

        self.elect.lineEdit_surface_charge_final.setText(str(SM.charge_density))
        self.elect.lineEdit_num_counterion.setText(str(SM.counterions))
        self.elect.lineEdit_num_mesh_points.setText(str(len(SM.wallR)))
        self.elect.lineEdit_num_mesh_points_final.setText(str(len(SM.wallR)))

        self.elect.lineEdit_total_counter_final.setText(str(SM.counterions))
        self.elect.lineEdit_type_counter_final.setText(str(SM.type_counter))
        self.elect.lineEdit_valency_counter_final.setText(str(SM.valency_counterion))
        # self.elect.lineEdit_surface_charge.setText(str(SM.charge_density))

###################STEP 2###################
    def add_salt(self):
        #Get the valiables from widget:
        SM.type_cation_salt_1 = str(self.elect.lineEdit_type_cation_salt_1.text())
        SM.charge_cation_salt_1 = self.elect.lineEdit_valency_cation_salt_1.text()
        SM.mass_cation_salt_1 = self.elect.lineEdit_mass_cation_salt_1.text()
        SM.type_anion_salt_1 = str(self.elect.lineEdit_type_anion_salt_1.text())
        SM.charge_anion_salt_1 = self.elect.lineEdit_valency_anion_salt_1.text()
        SM.mass_anion_salt_1 = self.elect.lineEdit_mass_anion_salt_1.text()
        SM.type_cation_salt_2 = str(self.elect.lineEdit_type_cation_salt_2.text())
        SM.charge_cation_salt_2 = self.elect.lineEdit_valency_cation_salt_2.text()
        SM.mass_cation_salt_2 = self.elect.lineEdit_mass_cation_salt_2.text()
        SM.type_anion_salt_2 = str(self.elect.lineEdit_type_anion_salt_2.text())
        SM.charge_anion_salt_2 = self.elect.lineEdit_valency_anion_salt_2.text()
        SM.mass_anion_salt_2 = self.elect.lineEdit_mass_anion_salt_2.text()
        SM.type_cation_salt_3 = str(self.elect.lineEdit_type_cation_salt_3.text())
        SM.charge_cation_salt_3 = self.elect.lineEdit_valency_cation_salt_3.text()
        SM.mass_cation_salt_3 = self.elect.lineEdit_mass_cation_salt_3.text()
        SM.type_anion_salt_3 = str(self.elect.lineEdit_type_anion_salt_3.text())
        SM.charge_anion_salt_3 = self.elect.lineEdit_valency_anion_salt_3.text()
        SM.mass_anion_salt_3 = self.elect.lineEdit_mass_anion_salt_3.text()
        SM.type_cation_salt_4 = str(self.elect.lineEdit_type_cation_salt_4.text())
        SM.charge_cation_salt_4 = self.elect.lineEdit_valency_cation_salt_4.text()
        SM.mass_cation_salt_4 = self.elect.lineEdit_mass_cation_salt_4.text()
        SM.type_anion_salt_4 = str(self.elect.lineEdit_type_anion_salt_4.text())
        SM.charge_anion_salt_4 = self.elect.lineEdit_valency_anion_salt_4.text()
        SM.mass_anion_salt_4 = self.elect.lineEdit_mass_anion_salt_4.text()

        # Here we add the final types of ions in page 3:
        self.elect.lineEdit_type_cation_salt_1_final.setText(str(SM.type_cation_salt_1))
        self.elect.lineEdit_type_anion_salt_1_final.setText(str(SM.type_anion_salt_1))
        self.elect.lineEdit_type_cation_salt_2_final.setText(str(SM.type_cation_salt_2))
        self.elect.lineEdit_type_anion_salt_2_final.setText(str(SM.type_anion_salt_2))
        self.elect.lineEdit_type_cation_salt_3_final.setText(str(SM.type_cation_salt_3))
        self.elect.lineEdit_type_anion_salt_3_final.setText(str(SM.type_anion_salt_3))
        self.elect.lineEdit_type_cation_salt_4_final.setText(str(SM.type_cation_salt_4))
        self.elect.lineEdit_type_anion_salt_4_final.setText(str(SM.type_anion_salt_4))

        num_type_cation_salt_1 = num_type_anion_salt_1 = num_type_cation_salt_2 = num_type_anion_salt_2 = num_type_cation_salt_3 = num_type_cation_salt_3 = num_type_anion_salt_3 = num_type_cation_salt_4 = num_type_anion_salt_4 = 0

        # We get the values of cation concentration from spin box. If the spin box does not have value, to prevent getting error,
        # we use Errors and Exceptions
        try:
            SM.con_cation_salt_1 = float(self.elect.SpinBox_con_cation_salt_1.value())
        except TypeError:
            SM.con_cation_salt_1 = 0
        try:
            SM.con_cation_salt_2 = float(self.elect.SpinBox_con_cation_salt_2.value())
        except TypeError:
            SM.con_cation_salt_2 = 0
        try:
            SM.con_cation_salt_3 = float(self.elect.SpinBox_con_cation_salt_3.value())
        except TypeError:
            SM.con_cation_salt_3 = 0
        try:
            SM.con_cation_salt_4 = float(self.elect.SpinBox_con_cation_salt_4.value())
        except TypeError:
            SM.con_cation_salt_4 = 0

    # We get the values of cation and anion charges from widget. If the there is no value, to prevent getting error,
    # we use Errors and Exceptions
        try:
            SM.charge_cation_salt_1 = float(SM.charge_cation_salt_1)
        except ValueError:
            SM.charge_cation_salt_1 = 0
        try:
            SM.charge_anion_salt_1 = float(SM.charge_anion_salt_1)
        except ValueError:
            SM.charge_anion_salt_1 = 0
        try:
            SM.charge_cation_salt_2 = float(SM.charge_cation_salt_2)
        except ValueError:
            SM.charge_cation_salt_2 = 0
        try:
            SM.charge_anion_salt_2 = float(SM.charge_anion_salt_2)
        except ValueError:
            SM.charge_anion_salt_2 = 0
        try:
            SM.charge_cation_salt_3 = float(SM.charge_cation_salt_3)
        except ValueError:
            SM.charge_cation_salt_3 = 0
        try:
            SM.charge_anion_salt_3 = float(SM.charge_anion_salt_3)
        except ValueError:
            SM.charge_anion_salt_3 = 0
        try:
            SM.charge_cation_salt_4 = float(SM.charge_cation_salt_4)
        except ValueError:
            SM.charge_cation_salt_4 = 0
        try:
            SM.charge_anion_salt_4 = float(SM.charge_anion_salt_4)
        except ValueError:
            SM.charge_anion_salt_4 = 0
    # We define the concentration of anions based on concentration of cations.
        volume_box =  SM.box_lx* SM.box_ly* SM.box_lz*0.001
        SM.total_cation_salt_1 =int((SM.con_cation_salt_1 * 0.6022) * (volume_box))
        try:
            if ((SM.total_cation_salt_1 % SM.charge_anion_salt_1) !=0):
                SM.total_cation_salt_1 = int(abs(SM.total_cation_salt_1 - (SM.total_cation_salt_1 % SM.charge_anion_salt_1) + SM.charge_anion_salt_1))
            SM.total_anion_salt_1 = int(abs(SM.charge_cation_salt_1 * SM.total_cation_salt_1 / SM.charge_anion_salt_1))
            SM.con_anion_salt_1  = SM.total_anion_salt_1 / (0.6022 * (volume_box))
        except ZeroDivisionError:
            SM.charge_cation_salt_1 = SM.total_cation_salt_1 = 0
        SM.con_anion_salt_1 = round(SM.con_anion_salt_1, 1)
        self.elect.lineEdit_con_anion_salt_1.setText(str("{:.1f}".format(SM.con_anion_salt_1)))
        self.elect.lineEdit_num_cation_salt_1.setText(str(SM.total_cation_salt_1))
        self.elect.lineEdit_num_anion_salt_1.setText(str(SM.total_anion_salt_1))
        SM.total_cation_salt_2 =int((SM.con_cation_salt_2 * 0.6022) * (volume_box))
        try:
            if ((SM.total_cation_salt_2 % SM.charge_anion_salt_2) !=0):
                SM.total_cation_salt_2 = SM.total_cation_salt_2 - (SM.total_cation_salt_2 % SM.charge_anion_salt_2) + SM.charge_anion_salt_2
            SM.total_anion_salt_2 = int(abs(SM.charge_cation_salt_2) * SM.total_cation_salt_2 / SM.charge_anion_salt_2)
            SM.con_anion_salt_2  = SM.total_anion_salt_2 / (0.6022 * (volume_box))
        except ZeroDivisionError:
            SM.charge_cation_salt_2 = SM.total_cation_salt_2 = 0
        SM.con_anion_salt_2 = round(SM.con_anion_salt_2, 1)
        self.elect.lineEdit_con_anion_salt_2.setText(str("{:.1f}".format(SM.con_anion_salt_2)))
        self.elect.lineEdit_num_cation_salt_2.setText(str(SM.total_cation_salt_2))
        self.elect.lineEdit_num_anion_salt_2.setText(str(SM.total_anion_salt_2))
        SM.total_cation_salt_3 =int((SM.con_cation_salt_3 * 0.6022) * (volume_box))
        try:
            if ((SM.total_cation_salt_3 % SM.charge_anion_salt_3) !=0):
                SM.total_cation_salt_3 = SM.total_cation_salt_3 - (SM.total_cation_salt_3 % SM.charge_anion_salt_3) + SM.charge_anion_salt_3
            SM.total_anion_salt_3 = int(abs(SM.charge_cation_salt_3) * SM.total_cation_salt_3 / SM.charge_anion_salt_3)
            SM.con_anion_salt_3  = SM.total_anion_salt_3 / (0.6022 * (volume_box))
        except ZeroDivisionError:
            SM.charge_cation_salt_3 = SM.total_cation_salt_3 = 0
        SM.con_anion_salt_3 = round(SM.con_anion_salt_3, 1)
        self.elect.lineEdit_con_anion_salt_3.setText(str("{:.1f}".format(SM.con_anion_salt_3)))
        self.elect.lineEdit_num_cation_salt_3.setText(str(SM.total_cation_salt_3))
        self.elect.lineEdit_num_anion_salt_3.setText(str(SM.total_anion_salt_3))
        SM.total_cation_salt_4 =int((SM.con_cation_salt_4 * 0.6022) * (volume_box))
        try:
            if ((SM.total_cation_salt_4 % SM.charge_anion_salt_4) !=0):
                SM.total_cation_salt_4 = SM.total_cation_salt_4 - (SM.total_cation_salt_4 % SM.charge_anion_salt_4) + SM.charge_anion_salt_4
            SM.total_anion_salt_4 = int(abs(SM.charge_cation_salt_4) * SM.total_cation_salt_4 / SM.charge_anion_salt_4)
            SM.con_anion_salt_4  = SM.total_anion_salt_4 / (0.6022 * (volume_box))
        except ZeroDivisionError:
            SM.charge_cation_salt_4 = SM.total_cation_salt_4 = 0
        SM.con_anion_salt_4 = round(SM.con_anion_salt_4, 1)
        self.elect.lineEdit_con_anion_salt_4.setText(str("{:.1f}".format(SM.con_anion_salt_4)))
        self.elect.lineEdit_num_cation_salt_4.setText(str(SM.total_cation_salt_4))
        self.elect.lineEdit_num_anion_salt_4.setText(str(SM.total_anion_salt_4))
        # total_saltions_inside = int(total_nions_inside + total_pions_inside + counterions)
        SM.total_cation_salt_2 =int((SM.con_anion_salt_2 * 0.6022) * (volume_box))
        SM.total_anion_salt_2 =int((SM.con_anion_salt_2 * 0.6022) * (volume_box))
        SM.total_cation_salt_3 =int((SM.con_anion_salt_3 * 0.6022) * (volume_box))
        SM.total_anion_salt_3 =int((SM.con_anion_salt_3 * 0.6022) * (volume_box))
        SM.total_cation_salt_4 =int((SM.con_anion_salt_4 * 0.6022) * (volume_box))
        SM.total_anion_salt_4 =int((SM.con_anion_salt_4 * 0.6022) * (volume_box))
        self.elect.lineEdit_num_cation_salt_2.setText(str(SM.total_cation_salt_2))
        self.elect.lineEdit_num_anion_salt_2.setText(str(SM.total_anion_salt_2))
        self.elect.lineEdit_num_cation_salt_3.setText(str(SM.total_cation_salt_3))
        self.elect.lineEdit_num_anion_salt_3.setText(str(SM.total_anion_salt_3))
        self.elect.lineEdit_num_cation_salt_4.setText(str(SM.total_cation_salt_4))
        self.elect.lineEdit_num_anion_salt_4.setText(str(SM.total_anion_salt_4))
        if SM.con_anion_salt_1 > 0:
            num_type_cation_salt_1 = 1
        if SM.con_anion_salt_1 > 0:
            num_type_anion_salt_1 = 1
        if SM.con_anion_salt_2 > 0:
            num_type_cation_salt_2 = 1
        if SM.con_anion_salt_2 > 0:
            num_type_anion_salt_2 = 1

        if SM.con_anion_salt_3 > 0:
            num_type_cation_salt_3 = 1
        if SM.con_anion_salt_3 > 0:
            num_type_anion_salt_3 = 1
        if SM.con_anion_salt_4 > 0:
            num_type_cation_salt_4 = 1
        if SM.con_anion_salt_4 > 0:
            num_type_anion_salt_4 = 1

        SM.total_num_salt_types = num_type_cation_salt_1 + num_type_anion_salt_1 + num_type_cation_salt_2 + num_type_anion_salt_2 + num_type_cation_salt_3 + num_type_cation_salt_3 + num_type_anion_salt_3 + num_type_cation_salt_4 + num_type_anion_salt_4
        SM.total_saltions_inside = int(SM.total_cation_salt_1 + SM.total_anion_salt_1 + SM.total_cation_salt_2 + SM.total_anion_salt_2 + SM.total_cation_salt_3 + SM.total_anion_salt_3 + SM.total_cation_salt_4 + SM.total_anion_salt_4)
        self.elect.lineEdit_total_ion_number_final.setText(str(SM.total_saltions_inside))
        try:
            SM.total_ions_concentration = ((SM.total_saltions_inside + SM.counterions) / (0.6022 * volume_box))
        except ZeroDivisionError:
            SM.total_ions_concentration = 0

        self.elect.lineEdit_total_ion_con_final.setText(str("{:.5f}".format(SM.total_ions_concentration)))

    def build_electrolyte(self):
        try:
            SM.charge_cation_salt_1 = float(SM.charge_cation_salt_1)
        except ValueError:
            SM.charge_cation_salt_1 = 0
        try:
            SM.charge_anion_salt_1 = float(SM.charge_anion_salt_1)
        except ValueError:
            SM.charge_anion_salt_1 = 0
        try:
            SM.charge_cation_salt_2 = float(SM.charge_cation_salt_2)
        except ValueError:
            SM.charge_cation_salt_2 = 0
        try:
            SM.charge_anion_salt_2 = float(SM.charge_anion_salt_2)
        except ValueError:
            SM.charge_anion_salt_2 = 0
        try:
            SM.charge_cation_salt_3 = float(SM.charge_cation_salt_3)
        except ValueError:
            SM.charge_cation_salt_3 = 0
        try:
            SM.charge_anion_salt_3 = float(SM.charge_anion_salt_3)
        except ValueError:
            SM.charge_anion_salt_3 = 0
        try:
            SM.charge_cation_salt_4 = float(SM.charge_cation_salt_4)
        except ValueError:
            SM.charge_cation_salt_4 = 0
        try:
            SM.charge_anion_salt_4 = float(SM.charge_anion_salt_4)
        except ValueError:
            SM.charge_anion_salt_4 = 0
        #test if the system is electroneutral:
        total_charge_ions = int((SM.charge_cation_salt_1 * SM.total_cation_salt_1) + (SM.charge_anion_salt_1 * SM.total_anion_salt_1) + (SM.charge_cation_salt_2 * SM.total_cation_salt_2) + (SM.charge_anion_salt_2 * SM.total_anion_salt_2) + (SM.charge_cation_salt_3 * SM.total_cation_salt_3) + (SM.charge_anion_salt_3 * SM.total_anion_salt_3) + (SM.charge_cation_salt_4 * SM.total_cation_salt_4) + (SM.charge_anion_salt_4 * SM.total_anion_salt_4))
        # total_charge_nions = int((SM.charge_cation_salt_3 * SM.total_cation_salt_3) + (SM.charge_anion_salt_3 * SM.total_anion_salt_3) + (SM.charge_cation_salt_4 * SM.total_cation_salt_4) + (SM.charge_anion_salt_4 * SM.total_anion_salt_4))
        if (total_charge_ions != 0):
            print('SM.charge_cation_salt_1', SM.charge_cation_salt_1, 'SM.total_cation_salt_1',SM.total_cation_salt_1)
            print('SM.charge_anion_salt_4', SM.charge_anion_salt_4, 'SM.total_anion_salt_4',SM.total_anion_salt_4)
            print('The electrolyte is not electroneutral; Abortion', total_charge_ions)
            self.show_warning_message_salt()
            return
            # return
        SM.charge_density = self.elect.lineEdit_surface_charge.text()
        try:
            SM.charge_density = float(SM.charge_density)
        except ValueError:
            SM.charge_density = 0
        SM.total_surface_charge = int( SM.charge_density *  SM.box_lx *  SM.box_ly / 16) # lx and ly are in A

        print("total surface charge ", SM.total_surface_charge)
        charge_system = SM.valency_counterion * SM.counterions + SM.total_surface_charge * 2.0

        print("total charge in the system ", charge_system)
        if (charge_system == 0):
            print("system is charge neutral")
        else:
            print("system is not electroneutral; aborting..., charge_system is", charge_system)
            self.show_warning_message_salt()
            return
        SM.total_saltions_inside = SM.total_saltions_inside + SM.counterions
        if SM.total_water_inside > 0 :
            SM.total_water_inside = SM.total_water_inside - SM.total_saltions_inside

        print("SM.total_water_inside:  ", SM.total_water_inside)
        print("total SM.counterions is: ", SM.counterions)
        print("total surface charge density is: ",  SM.charge_density, " C.m^-2")
        #   long double
        if ( SM.box_lx== 0 or  SM.box_ly== 0 or  SM.box_lz == 0):
            final_water_concentration = 0
        else:
            final_water_concentration = (SM.total_water_inside) /(0.6022 * ( SM.box_lx) * ( SM.box_ly) * ( SM.box_lz- SM.water_diameter) * 0.001)
        print("final water conc (in M): ", final_water_concentration)
        try:
            num_molecul_in_lx = int( SM.box_lx /  SM.spacing)
            num_molecul_in_ly = int( SM.box_ly /  SM.spacing)
            num_molecul_in_lz = int( SM.box_lz /  SM.spacing)
        except ZeroDivisionError:
            num_molecul_in_lx = num_molecul_in_ly = num_molecul_in_lz = 0


        if self.elect.frame_dummy_wall.isEnabled():
            pass
        else:
            SM.wallR = np.array([])
            SM.charge_density = SM.charge_meshpoint = SM.counterions = SM.valency_counterion = 0
        if self.elect.frame_water_properties.isEnabled():
            pass
        else:
            SM.total_water_inside = SM.num_type_atom_water = SM.final_water_concentration = 0
        coordinates_ions = []

        for i in range(num_molecul_in_lx):
            for j in range(num_molecul_in_ly):
                for k in range(num_molecul_in_lz):
                    if (len(coordinates_ions) < (SM.total_saltions_inside + SM.total_water_inside)):
                        x = float((- SM.box_lx/2 + (0.5* SM.spacing)) + i *  SM.spacing)
                        y = float((- SM.box_ly/2 + (0.5* SM.spacing)) + j *  SM.spacing)
                        z = float((- SM.box_lz/2 + (0.5* SM.spacing)) + k *  SM.spacing)
                        if ((x > ( SM.box_lx/2 - 0.5 *  SM.spacing)) or (y > ( SM.box_ly/2 - 0.5 *  SM.spacing)) or (z > ( SM.box_lz/2 - 0.5 *  SM.spacing))):
                            continue
                        position = np.array([x, y, z])
                        coordinates_ions.append(position)
        ions = np.array(coordinates_ions)

        total_atoms = (SM.total_water_inside * 3) + SM.total_saltions_inside + (2 * len(SM.wallR))
        if SM.total_water_inside > 0:
            num_type_atom_water = 2
        else:
            num_type_atom_water = 0
        if SM.counterions > 0:
            num_type_counterion = 1
        else:
            num_type_counterion = 0
        if (len(SM.wallR) > 0):
            num_type_wall = 2
        else:
            num_type_wall = 0

        total_type_atom_in_system = SM.total_num_salt_types + num_type_atom_water + num_type_counterion + num_type_wall
        mass_hydrogen ='1.008'
        type_hydrogen = 'H'
        mass_oxygen ='15.9994'
        type_oxygen = 'O'
        type_Rwall = 'R'
        type_Lwall = 'L'
        mass_Rwall = '0.00001'
        mass_Lwall = '0.00001'
        charge_oxygen = float(self.elect.lineEdit_oxygen_charge.text())
        charge_hydrogen = float(self.elect.lineEdit_hydrogen_charge.text())
        print("total mesh points on each surface: ", len(SM.wallR))
        print("charge on each meshpoint ", SM.charge_meshpoint, " e")
        charge_system = SM.charge_meshpoint * 2 * len(SM.wallR) + SM.counterions * SM.valency_counterion
        # charge_system = SM.valency_counterion * SM.counterions + SM.total_surface_charge * 2.0
        print("recheck: total charge in the system ", SM.counterions)
        if charge_system == 0:
            print("system is charge neutral again")
        else:
            print("system is not electroneutral; aborting...")
            self.show_warning_message_salt()
            return
        file_name = 'many_particle.DATA'
        outdump = open(file_name, "w")
        outdump.write("LAMMPS data file\n\n")
        outdump.write("{}\t{} \n".format(total_atoms, ' atoms'))
        outdump.write("{}\t{} \n".format(SM.total_water_inside * 2, ' bonds'))
        outdump.write("{}\t{} \n".format(SM.total_water_inside, ' angles'))
        outdump.write("{}\t{} \n".format('0', ' dihedrals'))
        outdump.write("{}\t{} \n\n".format('0', ' impropers'))
        outdump.write("{}\t{} \n".format(total_type_atom_in_system, 'atom types'))
        if SM.total_water_inside > 0:
            outdump.write("{}\t{} \n".format('1', 'bond types'))
            outdump.write("{}\t{} \n\n".format('1', 'angle types'))
        else:
            outdump.write("{}\t{} \n".format('0', 'bond types'))
            outdump.write("{}\t{} \n\n".format('0', 'angle types'))

        outdump.write("{}\t{}\t{} \n".format(-0.5* SM.box_lx, 0.5* SM.box_lx, ' xlo xhi'))
        outdump.write("{}\t{}\t{} \n".format(-0.5* SM.box_ly, 0.5* SM.box_ly, ' ylo yhi'))
        outdump.write("{}\t{}\t{} \n\n".format((-0.5* SM.box_lz) - 0.0005, (0.5* SM.box_lz) + 0.0005, ' zlo zhi'))
        outdump.write("Masses\n\n")

        if SM.total_water_inside > 0:
            outdump.write("{}\t{} \n".format(type_hydrogen, mass_hydrogen))
            outdump.write("{}\t{} \n".format(type_oxygen, mass_oxygen))
        if SM.counterions > 0:
            outdump.write("{}\t{} \n".format(SM.type_counter, SM.mass_counter))
        if SM.total_cation_salt_1 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cation_salt_1, SM.mass_cation_salt_1))
        if SM.total_anion_salt_1 > 0:
            outdump.write("{}\t{} \n".format(SM.type_anion_salt_1, SM.mass_anion_salt_1))
        if SM.total_cation_salt_2 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cation_salt_2, SM.mass_cation_salt_2))
        if SM.total_anion_salt_2 > 0:
            outdump.write("{}\t{} \n".format(SM.type_anion_salt_2, SM.mass_anion_salt_2))


        if SM.total_cation_salt_3 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cation_salt_3, SM.mass_cation_salt_3))
        if SM.total_anion_salt_3 > 0:
            outdump.write("{}\t{} \n".format(SM.type_anion_salt_3, SM.mass_anion_salt_3))
        if SM.total_cation_salt_4 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cation_salt_4, SM.mass_cation_salt_4))
        if SM.total_anion_salt_4 > 0:
            outdump.write("{}\t{} \n".format(SM.type_anion_salt_4, SM.mass_anion_salt_4))
        if len(SM.wallR) > 0:
            outdump.write("{}\t{} \n".format(type_Rwall, mass_Rwall))
            outdump.write("{}\t{} \n".format(type_Lwall, mass_Lwall))

        outdump.write("\n")
        outdump.write("Atoms          # full\n\n")
        if SM.total_water_inside > 0:
            num_molecules = 0
            for i in range(SM.total_water_inside):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 1 + num_molecules + num_molecules), (num_molecules + 1), type_oxygen, charge_oxygen, ions[i][0], ions[i][1], (ions[i][2]), '0   0   0'))
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 2 + num_molecules+ num_molecules), (num_molecules + 1), type_hydrogen, charge_hydrogen, ions[i][0] +  0.95908, (ions[i][1] + (-0.02691)), (ions[i][2]) + 0.03231, ' 0   0   0 '))
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 3 + num_molecules+ num_molecules), (num_molecules + 1), type_hydrogen, charge_hydrogen, ions[i][0]+(-0.28004), (ions[i][1] + (-0.58767)), (ions[i][2]) + 0.70556, '0   0   0 '))
                num_molecules = num_molecules + 1

        if SM.counterions > 0:
            num_molecules = 0
            for j in range(SM.total_water_inside, (SM.counterions+SM.total_water_inside)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_counter, SM.valency_counterion, ions[j][0], ions[j][1], ions[j][2],  '0   0   0 '))

        # crystal pack of counter ions (in this system, it is the same as positive ions):
        if SM.total_cation_salt_1 > 0:
            for j in range((SM.counterions+SM.total_water_inside), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cation_salt_1, SM.charge_cation_salt_1, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_anion_salt_1 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_cation_salt_1), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_anion_salt_1, SM.charge_anion_salt_1, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_cation_salt_2 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1 +SM.total_cation_salt_2)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cation_salt_2, SM.charge_cation_salt_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_anion_salt_2 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1+SM.total_cation_salt_2), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1 +SM.total_cation_salt_2+SM.total_anion_salt_2)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_anion_salt_2, SM.charge_anion_salt_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_cation_salt_3 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1+SM.total_cation_salt_2+SM.total_anion_salt_2), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1 +SM.total_cation_salt_2+SM.total_anion_salt_2+SM.total_cation_salt_3)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cation_salt_3, SM.charge_cation_salt_3, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_anion_salt_3 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1+SM.total_cation_salt_2+SM.total_anion_salt_2+SM.total_cation_salt_3), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1 +SM.total_cation_salt_2+SM.total_anion_salt_2+SM.total_cation_salt_3+SM.total_anion_salt_3)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_anion_salt_3, SM.charge_anion_salt_3, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_cation_salt_4 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1+SM.total_cation_salt_2+SM.total_anion_salt_2+SM.total_cation_salt_3+SM.total_anion_salt_3), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1 +SM.total_cation_salt_2+SM.total_anion_salt_2+SM.total_cation_salt_3+SM.total_anion_salt_3+SM.total_cation_salt_4)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cation_salt_4, SM.charge_cation_salt_4, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_anion_salt_4 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1+SM.total_cation_salt_2+SM.total_anion_salt_2+SM.total_cation_salt_3+SM.total_anion_salt_3+SM.total_cation_salt_4), (SM.counterions+SM.total_water_inside+SM.total_cation_salt_1+SM.total_anion_salt_1 +SM.total_cation_salt_2+SM.total_anion_salt_2+SM.total_cation_salt_3+SM.total_anion_salt_3+SM.total_cation_salt_4+SM.total_anion_salt_4)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_anion_salt_4, SM.charge_anion_salt_4, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))

        if (len(SM.wallR) > 0):
            for r  in range(len(SM.wallR)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((r + 1 + (3 * SM.total_water_inside) + SM.total_saltions_inside),  "0" , type_Rwall, SM.charge_meshpoint, SM.wallR[r][0], SM.wallR[r][1], SM.wallR[r][2], '0   0   0 '))
            # mesh points on left wall; the same as right wall with oppositee sign in z direction
            for l  in range(len(SM.wallR)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((l + 1 + (3 * SM.total_water_inside) + SM.total_saltions_inside + len(SM.wallR)),  "0" , type_Lwall, SM.charge_meshpoint, SM.wallR[l][0], SM.wallR[l][1], SM.wallR[l][2] * -1, '0   0   0 '))
        if SM.total_water_inside > 0:
            outdump.write("\n")
            outdump.write("Bonds\n\n")
            num_molecules = 0
            for g in range(SM.total_water_inside):
                outdump.write("{}\t{}\t{}\t{}\n".format((g + 1 + num_molecules), '1', (g + 1 + num_molecules + num_molecules), (g + 2 + num_molecules + num_molecules)))
                outdump.write("{}\t{}\t{}\t{}\n".format((g + 2 + num_molecules), '1', (g + 1 + num_molecules + num_molecules), (g + 3 + num_molecules + num_molecules)))
                num_molecules = num_molecules + 1
            outdump.write("\n")
            outdump.write("Angles\n\n")
            num_molecules = 0
            for n in range(SM.total_water_inside):
                outdump.write("{}\t{}\t{}\t{}\t{}\n".format( n + 1, '1', (n + 2 + num_molecules + num_molecules), (n + 1 + num_molecules + num_molecules), (n + 3 + num_molecules + num_molecules)))
                num_molecules = num_molecules + 1


        window = self.win.create_mdi_child()
        window.load_file(fname=file_name)
        window.show()
        return file_name










































        #  SM.charge_density = -0.01
        # counter_type = 'Na'
        # SM.valency_counterion = 1

        # SM.charge_cation_salt_1 = 1
        # SM.charge_anion_salt_1 = 1
        # SM.charge_cation_salt_2 = 1
        # SM.charge_anion_salt_2 = 1
        # SM.charge_cation_salt_3 = 1
        # SM.type_cation_salt_1 = 'Na'
        # SM.type_anion_salt_1 = 'K'
        # SM.type_cation_salt_2 = 'Cs'
        # SM.type_anion_salt_2 = 'Rb'
        # SM.type_cation_salt_3 = 'H'
        # SM.con_anion_salt_1 = 0.1
        # SM.con_anion_salt_1 = 0.1
        # SM.con_anion_salt_2 = 0.1
        # SM.con_anion_salt_2 = 0.1
        # SM.con_anion_salt_3 = 0.1

        # SM.charge_cation_salt_3 = -1
        # SM.charge_anion_salt_3 = -1
        # SM.charge_cation_salt_4 = -1
        # SM.charge_anion_salt_4 = -1
        # SM.charge_an_5 = -1
        # SM.type_cation_salt_3 = 'Cl'
        # SM.type_anion_salt_3 = 'B'
        # SM.type_cation_salt_4 = 'I'
        # SM.type_anion_salt_4 = 'F'
        # SM.type_an_5 = 'At'
        # SM.con_anion_salt_3 = 0.1
        # SM.con_anion_salt_3 = 0.1
        # SM.con_anion_salt_4 = 0.1
        # SM.con_anion_salt_4 = 0.1
        # con_an_5 = 0.1



        # SM.total_cation_salt_1 = int((SM.con_anion_salt_1 * 0.6022) * (volume_box))
        # SM.total_anion_salt_1 = int((SM.con_anion_salt_1 * 0.6022) * (volume_box))
        # SM.total_cation_salt_2 = int((SM.con_anion_salt_2 * 0.6022) * (volume_box))
        # SM.total_anion_salt_2 = int((SM.con_anion_salt_2 * 0.6022) * (volume_box))
        # SM.total_cation_salt_3 = int((SM.con_anion_salt_3 * 0.6022) * (volume_box))
        # SM.total_saltions_inside = int(SM.total_cation_salt_1 + SM.total_anion_salt_1 + SM.total_cation_salt_2 + SM.total_anion_salt_2 + SM.total_cation_salt_3)

        # SM.total_cation_salt_3 = int((SM.con_anion_salt_3 * 0.6022) * (volume_box))
        # SM.total_anion_salt_3 = int((SM.con_anion_salt_3 * 0.6022) * (volume_box))
        # SM.total_cation_salt_4 = int((SM.con_anion_salt_4 * 0.6022) * (volume_box))
        # SM.total_anion_salt_4 = int((SM.con_anion_salt_4 * 0.6022) * (volume_box))
        # SM.total_ncation_salt_3 = int((con_an_5 * 0.6022) * (volume_box))
        # SM.total_nions_inside = int(SM.total_cation_salt_3 + SM.total_anion_salt_3 + SM.total_cation_salt_4 + SM.total_anion_salt_4 + SM.total_ncation_salt_3)
        # total_charge_pions = int((SM.con_anion_salt_1 * SM.total_cation_salt_1) + (SM.con_anion_salt_1 * SM.total_anion_salt_1) + (SM.con_anion_salt_2 * SM.total_cation_salt_2) + (SM.con_anion_salt_2 * SM.total_anion_salt_2) + (SM.con_anion_salt_3 * SM.total_cation_salt_3))
        # total_charge_nions = int((SM.con_anion_salt_3 * SM.total_cation_salt_3) + (SM.con_anion_salt_3 * SM.total_anion_salt_3) + (SM.con_anion_salt_4 * SM.total_cation_salt_4) + (SM.con_anion_salt_4 * SM.total_anion_salt_4) + (con_an_5 * SM.total_ncation_salt_3))
        # if (total_charge_pions != total_charge_nions):
        #     print('The electrolyte is not electroneutral; Abortion')
        #     return

        # SM.total_saltions_inside = int(SM.total_nions_inside + SM.total_saltions_inside)

        # n_residues = SM.total_water_inside
        # n_atoms = n_residues * 3
        # resindices = np.repeat(range(n_residues), 3)
        # assert len(resindices) == n_atoms

        # # all water molecules belong to 1 segment
        # segindices = [0] * n_residues
        # sol = mda.Universe.empty(n_atoms,
        #                         n_residues=n_residues,
        #                         atom_resindex=resindices,
        #                         residue_segindex=segindices,
        #                         trajectory=True)
        # sol.dimensions = [lx, ly, lz, 0, 0, 0]
        # sol.add_TopologyAttr('name', ['O', 'H1', 'H2']*n_residues)
        # sol.add_TopologyAttr('type', ['O', 'H', 'H']*n_residues)
        # sol.add_TopologyAttr('resname', ['SOL']*n_residues)
        # sol.add_TopologyAttr('resid', list(range(1, n_residues+1)))
        # sol.add_TopologyAttr('segid', ['SOL'])

        # h2o = np.array([[ 0,        0,       0      ],  # oxygen
        #                 [ 0.95908, -0.02691, 0.03231],  # hydrogen
        #                 [-0.28004, -0.58767, 0.70556]]) # hydrogen

        # coordinates = []
        # num_molecul_in_lx = int(lx /  SM.spacing)
        # num_molecul_in_ly = int(ly /  SM.spacing)
        # num_molecul_in_lz = int(lz /  SM.spacing)
        # for i in range(num_molecul_in_lx):
        #     for j in range(num_molecul_in_ly):
        #         for k in range(num_molecul_in_lz):
        #             x = (-lx/2 + (0.5* SM.spacing)) + (i *  SM.spacing)
        #             y = (-ly/2 + (0.5* SM.spacing)) + (j *  SM.spacing)
        #             z = (-lz/2 + (0.5* SM.spacing)) + (k *  SM.spacing)
        #             xyz = np.array([x, y, z])
        #             coordinates.extend(h2o + xyz.T)
        # coord_array = np.array(coordinates)

        # assert coord_array.shape == (n_atoms, 3)
        # sol.atoms.positions = coord_array
        # assert not hasattr(sol, 'bonds')
        # bonds = []
        # for o in range(0, n_atoms, 3):
        #     bonds.extend([(o, o+1), (o, o+2)])
        # sol.add_TopologyAttr('bonds', bonds)

        # # water_center = sol.center_of_mass(pbc=True)
        # # dim = sol.dimensions
        # # box_center = np.sum(dim, axis=0)
        # # sol.atoms.translate(-lx/2 +  SM.spacing/2)
        # sol = MDAnalysis.Universe('C:/Users/nasim/OneDrive/Desktop/onlywater.pdb')
        # salt = MDAnalysis.Universe('C:/Users/nasim/OneDrive/Desktop/Fullerene_C720.pdb')


        # cog = sol.atoms.center_of_geometry()
        # print('Original solvent center of geometry: ', cog)
        # sol.atoms.positions -= cog
        # cog2 = sol.atoms.center_of_geometry()
        # print('New solvent center of geometry: ', cog2)
        # cog = salt.atoms.center_of_geometry()
        # print('Original solvent center of geometry: ', cog)
        # salt.atoms.positions -= cog
        # cog2 = salt.atoms.center_of_geometry()
        # print('New solvent center of geometry: ', cog2)


        # combined = MDAnalysis.Merge(salt.atoms, sol.atoms)
        # # u.select_atoms('same resid as (not around 30 salt)')
        # # point 5.0 5.0 5.0 3.5
        # no_overlap = combined.select_atoms("same resid as (not around 5 name C B)")

        # # no_overlap = combined.select_atoms("same resid as (not around 10 protein)")
        # u = mda.Merge(no_overlap)
        # print(len(u.atoms))
        # u.atoms.write('C:/Users/nasim/OneDrive/Desktop/solution.pdb')









    # def create_electrolyte( SM.box_lx,  SM.box_ly,  SM.box_lz,  SM.water_diameter,  spacing_dia, charge_hydrogen, charge_oxygen,  SM.charge_density, counter_type, SM.valency_counterion, SM.type_cation_salt_1, SM.type_anion_salt_1, SM.type_cation_salt_2, SM.type_anion_salt_2, SM.type_cation_salt_3,
    #                 SM.charge_cation_salt_1, SM.charge_anion_salt_1, SM.charge_cation_salt_2, SM.charge_anion_salt_2, SM.charge_cation_salt_3, SM.con_anion_salt_1, SM.con_anion_salt_1, SM.con_anion_salt_2, SM.con_anion_salt_2, SM.con_anion_salt_3, SM.type_cation_salt_3, SM.type_anion_salt_3, SM.type_cation_salt_4, SM.type_anion_salt_4,
    #                 SM.type_an_5, SM.charge_cation_salt_3, SM.charge_anion_salt_3, SM.charge_cation_salt_4, SM.charge_anion_salt_4, SM.charge_an_5, SM.con_anion_salt_3, SM.con_anion_salt_3, SM.con_anion_salt_4, SM.con_anion_salt_4, con_an_5):





