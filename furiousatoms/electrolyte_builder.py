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
import sys
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
        # self.elect.pushButton_type_counterion.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_type_counterion.clicked.connect(self.info_counterion)


        self.elect.pushButton_cation_1.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_cation_2.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_cation_3.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_cation_4.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_cation_5.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_anion_1.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_anion_2.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_anion_3.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_anion_4.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_anion_5.clicked.connect(self.show_periodic_table)
        self.elect.pushButton_cation_1.clicked.connect(self.info_cation_1)
        self.elect.pushButton_cation_2.clicked.connect(self.info_cation_2)
        self.elect.pushButton_cation_3.clicked.connect(self.info_cation_3)
        self.elect.pushButton_cation_4.clicked.connect(self.info_cation_4)
        self.elect.pushButton_cation_5.clicked.connect(self.info_cation_5)
        self.elect.pushButton_anion_1.clicked.connect(self.info_anion_1)
        self.elect.pushButton_anion_2.clicked.connect(self.info_anion_2)
        self.elect.pushButton_anion_3.clicked.connect(self.info_anion_3)
        self.elect.pushButton_anion_4.clicked.connect(self.info_anion_4)
        self.elect.pushButton_anion_5.clicked.connect(self.info_anion_5)

        self.elect.SpinBox_con_cation_1.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_cation_2.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_cation_3.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_cation_4.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_cation_5.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_anion_1.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_anion_2.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_anion_3.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_anion_4.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_con_anion_5.valueChanged.connect(self.calculate_parameters)
        self.elect.pushButton_build_electrolyte.clicked.connect(self.build_electrolyte)
        self.elect.lineEdit_type_counterion.textChanged[str].connect(self.calculate_parameters)###############
        # self.elect.pushButton_type_counterion.clicked.connect(self.calculate_parameters)


        SM.water_diameter = self.elect.lineEdit_water_diameter.setText(str(3.16555789))
        spacing_dia = self.elect.lineEdit_space_diameter.setText(str(1.0))
        SM.spacing = self.elect.lineEdit_space_bet_water.setText(str(3.16555789))
        SM.charge_density = self.elect.lineEdit_surface_charge.setText(str(0.0))

        self.elect.lineEdit_water_diameter.textChanged[str].connect(self.calculate_parameters)
        self.elect.lineEdit_space_diameter.textChanged[str].connect(self.calculate_parameters)
        self.elect.lineEdit_surface_charge.textChanged[str].connect(self.calculate_parameters)
        self.elect.SpinBox_lx.valueChanged.connect(self.calculate_parameters)
        self.elect.SpinBox_lz.valueChanged.connect(self.calculate_parameters)
        self.elect.comboBox_water_model.currentTextChanged.connect(self.calculate_parameters)


    def info_counterion(self):
        self.show_periodic_table()
        Ui_periodic.pt.current_edit_symbol = self.elect.lineEdit_type_counterion
        Ui_periodic.pt.current_edit_valency = self.elect.lineEdit_valency_counterion
        # TEST = self.elect.lineEdit_type_counterion.textChanged()
        # print(TEST, 'HOOOOOOOOOOOOO')




        # SM.mass_counter = float(SM.all_information_element['mass'])
        # SM.charge_counter = int(SM.all_information_element['valency'])
        # SM.type_counter =  str(SM.all_information_element['symbol'])

        # SM.charge_density = float(self.elect.lineEdit_surface_charge.text())
        # SM.box_lx = SM.box_ly = self.elect.SpinBox_lx.value()
        # SM.box_lz = self.elect.SpinBox_lz.value()
        # SM.total_surface_charge = int( SM.charge_density *  SM.box_lx *  SM.box_ly / 16) # lx and ly are in A
        # try:
        #     SM.counterions = int(2.0 * abs(SM.total_surface_charge)/SM.charge_counter)
        # except ZeroDivisionError:
        #     SM.counterions = 0
        # self.elect.lineEdit_type_counterion.setText(str(SM.type_counter))
        # self.elect.lineEdit_valency_counterion.setText(str(SM.charge_counter))
        # self.elect.lineEdit_num_counterion.setText(str(SM.counterions))
        # self.elect.lineEdit_total_counter_final.setText(str(SM.counterions))
        # self.elect.lineEdit_type_counter_final.setText(str(SM.type_counter))
        # self.elect.lineEdit_valency_counter_final.setText(str(SM.charge_counter))

    def info_cation_1(self):
        try:
            SM.all_information_element = Ui_periodic.pt.show_atom_type()
            SM.mass_cat_1 = SM.all_information_element['mass']
            SM.charge_cat_1 = SM.all_information_element['valency']
            SM.type_cat_1 =  str(SM.all_information_element['symbol'])
        except TypeError:
            SM.mass_cat_1 = 0
            SM.charge_cat_1 = 0
            SM.type_cat_1 = None
        self.elect.lineEdit_type_cat_1.setText(str(SM.type_cat_1))#, ',' ,str(SM.charge_cat_1))


    def info_cation_2(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_cat_2 = SM.all_information_element['mass']
        SM.charge_cat_2 = SM.all_information_element['valency']
        SM.type_cat_2 =  str(SM.all_information_element['symbol'])
    def info_cation_3(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_cat_3 = SM.all_information_element['mass']
        SM.charge_cat_3 = SM.all_information_element['valency']
        SM.type_cat_3 =  str(SM.all_information_element['symbol'])
    def info_cation_4(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_cat_4 = SM.all_information_element['mass']
        SM.charge_cat_4 = SM.all_information_element['valency']
        SM.type_cat_4 =  str(SM.all_information_element['symbol'])
    def info_cation_5(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_cat_5 = SM.all_information_element['mass']
        SM.charge_cat_5 = SM.all_information_element['valency']
        SM.type_cat_5 =  str(SM.all_information_element['symbol'])


    def info_anion_1(self):
        try:
            SM.all_information_element = Ui_periodic.pt.show_atom_type()
            SM.mass_an_1 = SM.all_information_element['mass']
            SM.charge_an_1 = SM.all_information_element['valency']
            SM.type_an_1 =  str(SM.all_information_element['symbol'])
        except TypeError:
            SM.mass_an_1 = 0
            SM.charge_an_1 = 0
            SM.type_an_1 = None
        # self.elect.lineEdit_type_an_1.setText(str(SM.type_an_1))

    def info_anion_2(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_an_2 = SM.all_information_element['mass']
        SM.charge_an_2 = SM.all_information_element['valency']
        SM.type_an_2 =  str(SM.all_information_element['symbol'])
    def info_anion_3(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_an_3 = SM.all_information_element['mass']
        SM.charge_an_3 = SM.all_information_element['valency']
        SM.type_an_3 =  str(SM.all_information_element['symbol'])
    def info_anion_4(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_an_4 = SM.all_information_element['mass']
        SM.charge_an_4 = SM.all_information_element['valency']
        SM.type_an_4 =  str(SM.all_information_element['symbol'])
    def info_anion_5(self):
        SM.all_information_element = Ui_periodic.pt.show_atom_type()
        SM.mass_an_5 = SM.all_information_element['mass']
        SM.charge_an_5 = SM.all_information_element['valency']
        SM.type_an_5 =  str(SM.all_information_element['symbol'])
    def show_periodic_table(self):
        Ui_periodic.pt = Ui_periodic()
        Ui_periodic.pt.win = self
        Ui_periodic.pt.show()


    def calculate_parameters(self):
        SM.type_counter =  str(self.elect.lineEdit_type_counterion.text())

        try:
            SM.charge_counter = float(self.elect.lineEdit_valency_counterion.text())
        except ValueError:
            SM.charge_counter  = 0
        try:
            SM.charge_density = float(self.elect.lineEdit_surface_charge.text())
        except ValueError:
            SM.charge_density = 0

        SM.box_lx = SM.box_ly = self.elect.SpinBox_lx.value()
        SM.box_lz = self.elect.SpinBox_lz.value()
        SM.total_surface_charge = int( SM.charge_density *  SM.box_lx *  SM.box_ly / 16) # lx and ly are in A
        print( ' SM.total_surface_charge is:  ', SM.total_surface_charge)
        print( ' SM.charge_counter is:  ', SM.charge_counter)
        try:
            SM.counterions = int(2.0 * abs(SM.total_surface_charge)/SM.charge_counter)
        except ZeroDivisionError:
            SM.counterions = 0

        self.elect.lineEdit_num_counterion.setText(str(SM.counterions))
        self.elect.lineEdit_total_counter_final.setText(str(SM.counterions))
        self.elect.lineEdit_type_counter_final.setText(str(SM.type_counter))
        self.elect.lineEdit_valency_counter_final.setText(str(SM.charge_counter))
        self.elect.lineEdit_type_cat_1.setText(str(SM.type_cat_1)) #, ',' ,str(SM.charge_cat_1))
        self.elect.lineEdit_type_an_1.setText(str(SM.type_an_1))
        self.elect.lineEdit_oxygen_charge.text()
        self.elect.lineEdit_hydrogen_charge.text()
        SM.water_diameter = self.elect.lineEdit_water_diameter.text()
        SM.spacing_dia = self.elect.lineEdit_space_diameter.text()
        SM.spacing_dia = float( SM.spacing_dia)
        SM.water_diameter = float( SM.water_diameter)
        SM.spacing =  SM.spacing_dia *  SM.water_diameter
        self.elect.lineEdit_space_bet_water.setText(str( SM.spacing))
###################SIMULATION BOX###################
        self.elect.lineEdit_lx_final.setText(str(SM.box_lx))
        self.elect.lineEdit_ly_final.setText(str(SM.box_ly))
        self.elect.lineEdit_lz_final.setText(str(SM.box_lz))
        self.elect.lineEdit_ly.setText(str(SM.box_ly))
        SM.charge_density = self.elect.lineEdit_surface_charge.text()
        self.elect.lineEdit_surface_charge_final.setText(str(SM.charge_density))
################### Number of mesh points###################
        width = 2  #A ; note CG width 3 A with lx = 15 nm
        self.elect.lineEdit_grid_size.setText(str(width))
################### Number of mesh points###################
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
        self.elect.lineEdit_num_mesh_points.setText(str(len(SM.wallR)))
################### Number ofcounterions###################
        try:
            SM.charge_density = float(SM.charge_density)
        except ValueError:
            SM.charge_density = 0
            self.elect.lineEdit_surface_charge_final.setText(str(SM.charge_density))
        SM.total_surface_charge = int( SM.charge_density *  SM.box_lx *  SM.box_ly / 16) # lx and ly are in A
        if SM.charge_density==0 or SM.charge_counter==0:
            SM.counterions = 0
        else:
            unitcharge = 1.60217646 * pow(10.0 , -19)
            surface_area =  SM.box_lx *  SM.box_ly * 0.01 * pow(10.0 , -18) # in unit of squared meter;
            SM.counterions = int(2.0 * abs(SM.total_surface_charge)/SM.charge_counter)

################### Charge of mesh points###################
        if (len(SM.wallR) > 0):
            num_type_wall = 2
            SM.charge_meshpoint = (SM.total_surface_charge * 1.0) / len(SM.wallR)
        else:
            num_type_wall = 0
            SM.charge_meshpoint = 0
        self.elect.lineEdit_mesh_charge.setText(str(SM.charge_meshpoint))
        volume_box =  SM.box_lx* SM.box_ly* SM.box_lz*0.001
        try:
            water_amount = float( SM.box_lx/ SM.spacing *  SM.box_ly/ SM.spacing *  SM.box_lz/ SM.spacing)
            SM.total_water_inside = int(int( SM.box_lx/ SM.spacing) * int( SM.box_ly/ SM.spacing) * int( SM.box_lz/ SM.spacing))
            water_concentration = water_amount / (0.602214076 * ( SM.box_lx) * ( SM.box_ly) * ( SM.box_lz- SM.water_diameter) * 0.001)
        except ZeroDivisionError:
            SM.total_water_inside = water_amount = water_concentration = 0

        self.elect.lineEdit_total_num_water.setText(str(SM.total_water_inside))
        self.elect.lineEdit_water_con.setText(str(water_concentration))
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

        SM.charge_density = self.elect.lineEdit_surface_charge.text()
        self.elect.lineEdit_surface_charge_final.setText(str(SM.charge_density))

        ###############create electrolyte

        num_type_cat_1=num_type_cat_2=num_type_cat_3=num_type_cat_4=num_type_cat_5=num_type_an_1=num_type_an_2=num_type_an_3=num_type_an_4=num_type_an_5=0
        con_cat_1 = float(self.elect.SpinBox_con_cation_1.text())
        con_cat_2 = float(self.elect.SpinBox_con_cation_2.text())
        con_cat_3 = float(self.elect.SpinBox_con_cation_3.text())
        con_cat_4 = float(self.elect.SpinBox_con_cation_4.text())
        con_cat_5 = float(self.elect.SpinBox_con_cation_5.text())
        con_an_1 = float(self.elect.SpinBox_con_anion_1.text())
        con_an_2 = float(self.elect.SpinBox_con_anion_2.text())
        con_an_3 = float(self.elect.SpinBox_con_anion_3.text())
        con_an_4 = float(self.elect.SpinBox_con_anion_4.text())
        con_an_5 = float(self.elect.SpinBox_con_anion_5.text())
        if con_cat_1 > 0:
            num_type_cat_1 = 1
        if con_cat_2 > 0:
            num_type_cat_2 = 1
        if con_cat_3 > 0:
            num_type_cat_3 = 1
        if con_cat_4 > 0:
            num_type_cat_4 = 1
        if con_cat_5 > 0:
            num_type_cat_5 = 1

        if con_an_1 > 0:
            num_type_an_1 = 1
        if con_an_2 > 0:
            num_type_an_2 = 1
        if con_an_3 > 0:
            num_type_an_3 = 1
        if con_an_4 > 0:
            num_type_an_4 = 1
        if con_an_5 > 0:
            num_type_an_5 = 1
        SM.total_num_salt_types = num_type_cat_1 + num_type_cat_2 + num_type_cat_3 + num_type_cat_4 + num_type_cat_5 + num_type_an_1 + num_type_an_2 + num_type_an_3 + num_type_an_4 + num_type_an_5
        SM.total_pions_1 = int((con_cat_1 * 0.6022) * (volume_box))
        SM.total_pions_2 = int((con_cat_2 * 0.6022) * (volume_box))
        SM.total_pions_3 = int((con_cat_3 * 0.6022) * (volume_box))
        SM.total_pions_4 = int((con_cat_4 * 0.6022) * (volume_box))
        SM.total_pions_5 = int((con_cat_5 * 0.6022) * (volume_box))
        SM.total_pions_inside = int(SM.total_pions_1 + SM.total_pions_2 + SM.total_pions_3 + SM.total_pions_4 + SM.total_pions_5)
        self.elect.lineEdit_num_cation_1.setText(str(SM.total_pions_1))
        self.elect.lineEdit_num_cation_2.setText(str(SM.total_pions_2))
        self.elect.lineEdit_num_cation_3.setText(str(SM.total_pions_3))
        self.elect.lineEdit_num_cation_4.setText(str(SM.total_pions_4))
        self.elect.lineEdit_num_cation_5.setText(str(SM.total_pions_5))

        SM.total_nions_1 = int((con_an_1 * 0.6022) * (volume_box))
        SM.total_nions_2 = int((con_an_2 * 0.6022) * (volume_box))
        SM.total_nions_3 = int((con_an_3 * 0.6022) * (volume_box))
        SM.total_nions_4 = int((con_an_4 * 0.6022) * (volume_box))
        SM.total_nions_5 = int((con_an_5 * 0.6022) * (volume_box))
        SM.total_nions_inside = int(SM.total_nions_1 + SM.total_nions_2 + SM.total_nions_3 + SM.total_nions_4 + SM.total_nions_5)

        self.elect.lineEdit_num_anion_1.setText(str(SM.total_nions_1))
        self.elect.lineEdit_num_anion_2.setText(str(SM.total_nions_2))
        self.elect.lineEdit_num_anion_3.setText(str(SM.total_nions_3))
        self.elect.lineEdit_num_anion_4.setText(str(SM.total_nions_4))
        self.elect.lineEdit_num_anion_5.setText(str(SM.total_nions_5))


    def build_electrolyte(self):
        #test if the system is electroneutral:
        total_charge_pions = int((SM.charge_cat_1 * SM.total_pions_1) + (SM.charge_cat_2 * SM.total_pions_2) + (SM.charge_cat_3 * SM.total_pions_3) + (SM.charge_cat_4 * SM.total_pions_4) + (SM.charge_cat_5 * SM.total_pions_5))
        total_charge_nions = int((SM.charge_an_1 * SM.total_nions_1) + (SM.charge_an_2 * SM.total_nions_2) + (SM.charge_an_3 * SM.total_nions_3) + (SM.charge_an_4 * SM.total_nions_4) + (SM.charge_an_5 * SM.total_nions_5))
        if (total_charge_pions + total_charge_nions != 0):
            print('The electrolyte is not electroneutral; Abortion')
            # return
        total_saltions_inside = int(SM.total_nions_inside + SM.total_pions_inside)

        print("total surface charge ", SM.total_surface_charge)
        charge_system = SM.charge_counter * SM.counterions + SM.total_surface_charge * 2.0

        print("total charge in the system ", charge_system)
        if (charge_system == 0):
            print("system is charge neutral")
        else:
            print("system is not electroneutral; aborting..." )
            # return

        total_saltions_inside = total_saltions_inside + SM.counterions
        SM.total_water_inside = SM.total_water_inside - total_saltions_inside

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
        coordinates_ions = []
        for i in range(num_molecul_in_lx):
            for j in range(num_molecul_in_ly):
                for k in range(num_molecul_in_lz):
                    if (len(coordinates_ions) < (total_saltions_inside + SM.total_water_inside)):
                        x = float((- SM.box_lx/2 + (0.5* SM.spacing)) + i *  SM.spacing)
                        y = float((- SM.box_ly/2 + (0.5* SM.spacing)) + j *  SM.spacing)
                        z = float((- SM.box_lz/2 + (0.5* SM.spacing)) + k *  SM.spacing)
                        if ((x > ( SM.box_lx/2 - 0.5 *  SM.spacing)) or (y > ( SM.box_ly/2 - 0.5 *  SM.spacing)) or (z > ( SM.box_lz/2 - 0.5 *  SM.spacing))):
                            continue
                        position = np.array([x, y, z])
                        coordinates_ions.append(position)#.T)
        ions = np.array(coordinates_ions)

        total_atoms = (SM.total_water_inside * 3) + total_saltions_inside + (2 * len(SM.wallR))
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
        self.elect.lineEdit_mesh_charge.setText(str(SM.charge_meshpoint))
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
        charge_system = SM.charge_meshpoint * 2 * len(SM.wallR) + SM.counterions
        print("recheck: total charge in the system ", charge_system)
        if charge_system == 0:
            print("system is charge neutral again")
        else:
            print("system is not electroneutral; aborting...")
        file_name = 'many_particle.data'
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
        if SM.total_pions_1 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cat_1, SM.mass_cat_1))
        if SM.total_pions_2 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cat_2, SM.charge_cat_2))
        if SM.total_pions_3 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cat_3, SM.mass_cat_3))
        if SM.total_pions_4 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cat_4, SM.mass_cat_4))
        if SM.total_pions_5 > 0:
            outdump.write("{}\t{} \n".format(SM.type_cat_5, SM.mass_cat_5))

        if SM.total_nions_1 > 0:
            outdump.write("{}\t{} \n".format(SM.type_an_1, SM.mass_an_1))
        if SM.total_nions_2 > 0:
            outdump.write("{}\t{} \n".format(SM.type_an_2, SM.mass_an_2))
        if SM.total_nions_3 > 0:
            outdump.write("{}\t{} \n".format(SM.type_an_3, SM.mass_an_3))
        if SM.total_nions_4 > 0:
            outdump.write("{}\t{} \n".format(SM.type_an_4, SM.mass_an_4))
        if SM.total_nions_5 > 0:
            outdump.write("{}\t{} \n".format(SM.type_an_5, SM.mass_an_5))

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
            for j in range(SM.total_water_inside, (SM.total_pions_inside+SM.total_water_inside)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_counter, SM.charge_counter, ions[j][0], ions[j][1], ions[j][2],  '0   0   0 '))

        # crystal pack of counter ions (in this system, it is the same as positive ions):
        if SM.total_pions_1 > 0:
            for j in range((SM.counterions+SM.total_water_inside), (SM.counterions+SM.total_water_inside+SM.total_pions_1)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cat_1, SM.charge_cat_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_pions_2 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_1), (SM.counterions+SM.total_water_inside+SM.total_pions_1+SM.total_pions_2)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cat_2, SM.charge_cat_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_pions_3 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_1+SM.total_pions_2), (SM.counterions+SM.total_water_inside+SM.total_pions_1+SM.total_pions_2 +SM.total_pions_3)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cat_3, SM.charge_cat_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_pions_4 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_1+SM.total_pions_2+SM.total_pions_3), (SM.counterions+SM.total_water_inside+SM.total_pions_1+SM.total_pions_2 +SM.total_pions_3+SM.total_pions_4)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cat_4, SM.charge_cat_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_pions_5 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_1+SM.total_pions_2+SM.total_pions_3+SM.total_pions_4), (SM.counterions+SM.total_water_inside+SM.total_pions_1+SM.total_pions_2 +SM.total_pions_3+SM.total_pions_4+SM.total_pions_5)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_cat_5, SM.charge_cat_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))


        # # crystal pack of negative ions: //the diameter of Cl is 5 A, we should pack them with more space;

        if SM.total_nions_1 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_inside), (SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_an_1, SM.charge_an_1, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_nions_2 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1), (SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1+SM.total_nions_2)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_an_2, SM.charge_an_2, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_nions_3 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1+SM.total_nions_2), (SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1+SM.total_nions_2 +SM.total_nions_3)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_an_3, SM.charge_an_3, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_nions_4 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1+SM.total_nions_2+SM.total_nions_3), (SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_pions_1+SM.total_pions_2 +SM.total_pions_3+SM.total_pions_4)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_an_4, SM.charge_an_4, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if SM.total_nions_5 > 0:
            for j in range((SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1+SM.total_nions_2+SM.total_nions_3+SM.total_nions_4), (SM.counterions+SM.total_water_inside+SM.total_pions_inside+SM.total_nions_1+SM.total_nions_2 +SM.total_nions_3+SM.total_nions_4+SM.total_nions_5)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + SM.total_water_inside*2), num_molecules, SM.type_an_5, SM.charge_an_5, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))

        if (len(SM.wallR) > 0):
            for r  in range(len(SM.wallR)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((r + 1 + (3 * SM.total_water_inside) + total_saltions_inside),  "0" , type_Rwall, SM.charge_meshpoint, SM.wallR[r][0], SM.wallR[r][1], SM.wallR[r][2], '0   0   0 '))
            # mesh points on left wall; the same as right wall with opposite sign in z direction
            for l  in range(len(SM.wallR)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((l + 1 + (3 * SM.total_water_inside) + total_saltions_inside + len(SM.wallR)),  "0" , type_Lwall, SM.charge_meshpoint, SM.wallR[l][0], SM.wallR[l][1], SM.wallR[l][2] * -1, '0   0   0 '))
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
        self.win.process_load_file(fname=file_name)











































        #  SM.charge_density = -0.01
        # counter_type = 'Na'
        # SM.charge_counter = 1

        # SM.charge_cat_1 = 1
        # SM.charge_cat_2 = 1
        # SM.charge_cat_3 = 1
        # SM.charge_cat_4 = 1
        # SM.charge_cat_5 = 1
        # SM.type_cat_1 = 'Na'
        # SM.type_cat_2 = 'K'
        # SM.type_cat_3 = 'Cs'
        # SM.type_cat_4 = 'Rb'
        # SM.type_cat_5 = 'H'
        # con_cat_1 = 0.1
        # con_cat_2 = 0.1
        # con_cat_3 = 0.1
        # con_cat_4 = 0.1
        # con_cat_5 = 0.1

        # SM.charge_an_1 = -1
        # SM.charge_an_2 = -1
        # SM.charge_an_3 = -1
        # SM.charge_an_4 = -1
        # SM.charge_an_5 = -1
        # SM.type_an_1 = 'Cl'
        # SM.type_an_2 = 'B'
        # SM.type_an_3 = 'I'
        # SM.type_an_4 = 'F'
        # SM.type_an_5 = 'At'
        # con_an_1 = 0.1
        # con_an_2 = 0.1
        # con_an_3 = 0.1
        # con_an_4 = 0.1
        # con_an_5 = 0.1



        # SM.total_pions_1 = int((con_cat_1 * 0.6022) * (volume_box))
        # SM.total_pions_2 = int((con_cat_2 * 0.6022) * (volume_box))
        # SM.total_pions_3 = int((con_cat_3 * 0.6022) * (volume_box))
        # SM.total_pions_4 = int((con_cat_4 * 0.6022) * (volume_box))
        # SM.total_pions_5 = int((con_cat_5 * 0.6022) * (volume_box))
        # SM.total_pions_inside = int(SM.total_pions_1 + SM.total_pions_2 + SM.total_pions_3 + SM.total_pions_4 + SM.total_pions_5)

        # SM.total_nions_1 = int((con_an_1 * 0.6022) * (volume_box))
        # SM.total_nions_2 = int((con_an_2 * 0.6022) * (volume_box))
        # SM.total_nions_3 = int((con_an_3 * 0.6022) * (volume_box))
        # SM.total_nions_4 = int((con_an_4 * 0.6022) * (volume_box))
        # SM.total_nions_5 = int((con_an_5 * 0.6022) * (volume_box))
        # SM.total_nions_inside = int(SM.total_nions_1 + SM.total_nions_2 + SM.total_nions_3 + SM.total_nions_4 + SM.total_nions_5)
        # total_charge_pions = int((con_cat_1 * SM.total_pions_1) + (con_cat_2 * SM.total_pions_2) + (con_cat_3 * SM.total_pions_3) + (con_cat_4 * SM.total_pions_4) + (con_cat_5 * SM.total_pions_5))
        # total_charge_nions = int((con_an_1 * SM.total_nions_1) + (con_an_2 * SM.total_nions_2) + (con_an_3 * SM.total_nions_3) + (con_an_4 * SM.total_nions_4) + (con_an_5 * SM.total_nions_5))
        # if (total_charge_pions != total_charge_nions):
        #     print('The electrolyte is not electroneutral; Abortion')
        #     return

        # total_saltions_inside = int(SM.total_nions_inside + SM.total_pions_inside)

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









    # def create_electrolyte( SM.box_lx,  SM.box_ly,  SM.box_lz,  SM.water_diameter,  spacing_dia, charge_hydrogen, charge_oxygen,  SM.charge_density, counter_type, SM.charge_counter, SM.type_cat_1, SM.type_cat_2, SM.type_cat_3, SM.type_cat_4, SM.type_cat_5,
    #                 SM.charge_cat_1, SM.charge_cat_2, SM.charge_cat_3, SM.charge_cat_4, SM.charge_cat_5, con_cat_1, con_cat_2, con_cat_3, con_cat_4, con_cat_5, SM.type_an_1, SM.type_an_2, SM.type_an_3, SM.type_an_4,
    #                 SM.type_an_5, SM.charge_an_1, SM.charge_an_2, SM.charge_an_3, SM.charge_an_4, SM.charge_an_5, con_an_1, con_an_2, con_an_3, con_an_4, con_an_5):





