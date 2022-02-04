import numpy as np
import io
from furiousatoms import io
import numpy as np
from fury import window
from PySide2 import QtCore
from PySide2 import QtWidgets
from furiousatoms.sharedmem import SharedMemory
from furiousatoms.periodic_table import Ui_periodic_cation
from furiousatoms.periodic_table import Ui_periodic_anion
ShM= SharedMemory()
class Ui_electrolyte(QtWidgets.QMainWindow):
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
        self.elect.SpinBox_con_cation_salt_1.valueChanged.connect(self.add_salt)
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
        self.elect.lineEdit_type_anion_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_valency_anion_salt_1.textChanged[str].connect(self.add_salt)
        self.elect.lineEdit_mass_anion_salt_1.textChanged[str].connect(self.add_salt)
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
        self.elect.lineEdit_grid_size.textChanged[str].connect(self.dummy_wall)
        water_diameter = self.elect.lineEdit_water_diameter.setText(str(3.16555789))
        spacing_dia = self.elect.lineEdit_space_diameter.setText(str(1.0))
        spacing = self.elect.lineEdit_space_bet_water.setText(str(3.16555789))
        ShM.charge_density = self.elect.lineEdit_surface_charge.setText(str(0.0))
    def show_periodic_cation_table(self):
        Ui_periodic_cation.pt = Ui_periodic_cation()
        Ui_periodic_cation.pt.win = self
        Ui_periodic_cation.pt.show()

    def show_periodic_anion_table(self):
        Ui_periodic_anion.pt = Ui_periodic_anion()
        Ui_periodic_anion.pt.win = self
        Ui_periodic_anion.pt.show()

    def info_counterion(self):
        ShM.charge_density = float(self.elect.lineEdit_surface_charge.text())
        if ShM.charge_density > 0:
            self.show_periodic_anion_table()
            Ui_periodic_anion.pt.current_edit_symbol = self.elect.lineEdit_type_counterion
            Ui_periodic_anion.pt.current_edit_valency = self.elect.lineEdit_valency_counterion
            Ui_periodic_anion.pt.current_edit_mass = self.elect.lineEdit_mass_counterion
        elif ShM.charge_density < 0:
            self.show_periodic_cation_table()
            Ui_periodic_cation.pt.current_edit_symbol = self.elect.lineEdit_type_counterion
            Ui_periodic_cation.pt.current_edit_valency = self.elect.lineEdit_valency_counterion
            Ui_periodic_cation.pt.current_edit_mass = self.elect.lineEdit_mass_counterion

    def info_cation_salt_1(self):
        self.show_periodic_cation_table()
        Ui_periodic_cation.pt.current_edit_symbol = self.elect.lineEdit_type_cation_salt_1
        Ui_periodic_cation.pt.current_edit_valency = self.elect.lineEdit_valency_cation_salt_1
        Ui_periodic_cation.pt.current_edit_mass = self.elect.lineEdit_mass_cation_salt_1
    def info_anion_salt_1(self):
        self.show_periodic_anion_table()
        Ui_periodic_anion.pt.current_edit_symbol = self.elect.lineEdit_type_anion_salt_1
        Ui_periodic_anion.pt.current_edit_valency = self.elect.lineEdit_valency_anion_salt_1
        Ui_periodic_anion.pt.current_edit_mass = self.elect.lineEdit_mass_anion_salt_1

    def add_water(self):
        water_diameter = self.elect.lineEdit_water_diameter.text()
        spacing_dia = self.elect.lineEdit_space_diameter.text()
        try:
            spacing_dia = float(spacing_dia)
        except ValueError:
            spacing_dia = 0
        water_diameter = float(water_diameter)
        spacing =  spacing_dia *  water_diameter
        self.elect.lineEdit_space_bet_water.setText(str(spacing))
        self.elect.lineEdit_oxygen_charge.text()
        self.elect.lineEdit_hydrogen_charge.text()
        try:
            ShM.total_water_inside = (int( ShM.box_lx/ spacing) * int( ShM.box_ly/ spacing) * int( ShM.box_lz/ spacing))
            water_concentration = ShM.total_water_inside  / (0.602214076 * ( ShM.box_lx) * ( ShM.box_ly) * ( ShM.box_lz - water_diameter) * 0.001)
        except ZeroDivisionError:
            ShM.total_water_inside = water_amount = water_concentration = 0
        self.elect.lineEdit_total_num_water.setText(str(ShM.total_water_inside))
        self.elect.lineEdit_total_num_water_final.setText(str(ShM.total_water_inside))
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
            ShM.valency_counterion = float(self.elect.lineEdit_valency_counterion.text())
        except ValueError:
            ShM.valency_counterion = 0
        if ShM.charge_density != 0:
            ShM.total_surface_charge = int( ShM.charge_density *  ShM.box_lx *  ShM.box_ly / 16) # lx and ly are in A
        else:
            ShM.total_surface_charge = 0
        try:
            ShM.counterions = int(2.0 * abs(ShM.total_surface_charge/ShM.valency_counterion))
        except ZeroDivisionError:
            ShM.counterions = 0
        print("valency_counterion ", ShM.valency_counterion)
        print("counterions ", ShM.counterions)
        print("total surface charge ", ShM.total_surface_charge)
        charge_system = ShM.valency_counterion * ShM.counterions + ShM.total_surface_charge * 2.0
        print("total charge in the system ", charge_system)
        if (charge_system == 0):
            print("system is charge neutral")
        else:
            print("system is not electroneutral; aborting..." )
            # self.show_warning_opposite_charge()
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
        ShM.box_lx = ShM.box_ly = self.elect.SpinBox_lx.value()
        ShM.box_lz = self.elect.SpinBox_lz.value()
        self.elect.lineEdit_lx_final.setText(str(ShM.box_lx))
        self.elect.lineEdit_ly_final.setText(str(ShM.box_ly))
        self.elect.lineEdit_lz_final.setText(str(ShM.box_lz))
        self.elect.lineEdit_ly.setText(str(ShM.box_ly))
        ShM.type_counter = str(self.elect.lineEdit_type_counterion.text())
        ShM.mass_counter = str(self.elect.lineEdit_mass_counterion.text())
        try:
            ShM.valency_counterion = float(self.elect.lineEdit_valency_counterion.text())
        except ValueError:
            ShM.valency_counterion = 0.0
        try:
            width = int(self.elect.lineEdit_grid_size.text())
        except:
            width = 1
            self.elect.lineEdit_grid_size.setText(str(width))

        coordinates_wallR = []
        ShM.box_lx = float(ShM.box_lx)
        ShM.box_ly = float(ShM.box_ly)
        ShM.box_lz = float(ShM.box_lz)
        nx = int( ShM.box_lx / width)
        ny = int( ShM.box_ly / width)
        for i in range(ny):
            for j in range(nx):
                x = float(-0.5* ShM.box_lx+(0.5)*width+i*width)
                y = float(-0.5* ShM.box_ly+(0.5)*width+j*width)
                z = -0.5 *  ShM.box_lz
                xyz_wallR = np.array([x, y, z])
                coordinates_wallR.append(xyz_wallR)
        ShM.wallR = np.array(coordinates_wallR)

# Define surface charge density on the walls and number of counterions to keep the system electroneutral:
        ShM.charge_density = self.elect.lineEdit_surface_charge.text()
        try:
            ShM.charge_density = float(ShM.charge_density)
        except ValueError:
            ShM.charge_density = 0
            # self.elect.lineEdit_surface_charge_final.setText(str(ShM.charge_density))
        if ShM.charge_density == 0:
            ShM.type_counter = None
            ShM.mass_counterion = 0.0
            ShM.valency_counterion = 0.0
            ShM.counterions = 0
        # self.elect.lineEdit_type_counterion.setText(ShM.type_counter)
        self.elect.lineEdit_type_counter_final.setText(ShM.type_counter)
        unitcharge = 1.60217646 * pow(10.0 , -19)
        ShM.total_surface_charge = int( ShM.charge_density *  ShM.box_lx *  ShM.box_ly / 16) # lx and ly are in A
        surface_area =  ShM.box_lx *  ShM.box_ly * 0.01 * pow(10.0 , -18) # in unit of squared meter;
        try:
            ShM.counterions = int(2.0 * abs(ShM.total_surface_charge/ShM.valency_counterion))
        except ZeroDivisionError:
            ShM.counterions = 0
        if (len(ShM.wallR) > 0):
            num_type_wall = 2
            ShM.charge_meshpoint = (ShM.total_surface_charge * 1.0) / len(ShM.wallR)
        else:
            num_type_wall = 0
            ShM.charge_meshpoint = 0

        self.elect.lineEdit_mesh_charge.setText(str("{:.5f}".format(ShM.charge_meshpoint)))
        self.elect.lineEdit_surface_charge_final.setText(str(ShM.charge_density))
        self.elect.lineEdit_num_counterion.setText(str(ShM.counterions))
        self.elect.lineEdit_num_mesh_points.setText(str(len(ShM.wallR)))
        self.elect.lineEdit_num_mesh_points_final.setText(str(len(ShM.wallR)))

        self.elect.lineEdit_total_counter_final.setText(str(ShM.counterions))
        self.elect.lineEdit_type_counter_final.setText(str(ShM.type_counter))
        self.elect.lineEdit_valency_counter_final.setText(str(ShM.valency_counterion))

###################STEP 2###################
    def add_salt(self):
        #Get the valiables from widget:
        ShM.type_cation_salt_1 = str(self.elect.lineEdit_type_cation_salt_1.text())
        charge_cation_salt_1 = self.elect.lineEdit_valency_cation_salt_1.text()
        ShM.mass_cation_salt_1 = self.elect.lineEdit_mass_cation_salt_1.text()
        ShM.type_anion_salt_1 = str(self.elect.lineEdit_type_anion_salt_1.text())
        charge_anion_salt_1 = self.elect.lineEdit_valency_anion_salt_1.text()
        ShM.mass_anion_salt_1 = self.elect.lineEdit_mass_anion_salt_1.text()
        self.elect.lineEdit_type_cation_salt_1_final.setText(str(ShM.type_cation_salt_1))
        self.elect.lineEdit_type_anion_salt_1_final.setText(str(ShM.type_anion_salt_1))
        num_type_cation_salt_1 = num_type_anion_salt_1 = 0

        # We get the values of cation concentration from spin box. If the spin box does not have value, to prevent getting error,
        # we use Errors and Exceptions
        try:
            ShM.con_cation_salt_1 = float(self.elect.SpinBox_con_cation_salt_1.value())
        except TypeError:
            ShM.con_cation_salt_1 = 0

    # We get the values of cation and anion charges from widget. If the there is no value, to prevent getting error,
    # we use Errors and Exceptions
        try:
            charge_cation_salt_1 = float(charge_cation_salt_1)
        except ValueError:
            charge_cation_salt_1 = 0
        try:
            charge_anion_salt_1 = float(charge_anion_salt_1)
        except ValueError:
            charge_anion_salt_1 = 0

    # We define the concentration of anions based on concentration of cations.
        volume_box =  ShM.box_lx* ShM.box_ly* ShM.box_lz*0.001
        ShM.total_cation_salt_1 =int((ShM.con_cation_salt_1 * 0.6022) * (volume_box))
        try:
            if ((ShM.total_cation_salt_1 % charge_anion_salt_1) !=0):
                ShM.total_cation_salt_1 = int(abs(ShM.total_cation_salt_1 - (ShM.total_cation_salt_1 % charge_anion_salt_1) + charge_anion_salt_1))
            ShM.total_anion_salt_1 = int(abs(charge_cation_salt_1 * ShM.total_cation_salt_1 / charge_anion_salt_1))
            # ShM.con_anion_salt_1  = ShM.total_anion_salt_1 / (0.6022 * (volume_box))
            ShM.con_anion_salt_1 = ShM.con_cation_salt_1
        except ZeroDivisionError:
            charge_cation_salt_1 = ShM.total_cation_salt_1 = 0
        self.elect.lineEdit_con_anion_salt_1.setText(str("{:.1f}".format(ShM.con_anion_salt_1)))
        self.elect.lineEdit_num_cation_salt_1.setText(str(ShM.total_cation_salt_1))
        self.elect.lineEdit_num_anion_salt_1.setText(str(ShM.total_anion_salt_1))

        if ShM.con_anion_salt_1 > 0:
            num_type_cation_salt_1 = 1
        if ShM.con_anion_salt_1 > 0:
            num_type_anion_salt_1 = 1

        ShM.total_num_salt_types = num_type_cation_salt_1 + num_type_anion_salt_1
        ShM.total_saltions_inside = int(ShM.total_cation_salt_1 + ShM.total_anion_salt_1)
        self.elect.lineEdit_total_ion_number_final.setText(str(ShM.total_saltions_inside))
        try:
            ShM.total_ions_concentration = ((ShM.total_saltions_inside + ShM.counterions) / (0.6022 * volume_box))
        except ZeroDivisionError:
            ShM.total_ions_concentration = 0

        self.elect.lineEdit_total_ion_con_final.setText(str("{:.5f}".format(ShM.total_ions_concentration)))

    def build_electrolyte(self):
        type_cation_salt_1 = str(self.elect.lineEdit_type_cation_salt_1.text())
        charge_cation_salt_1 = self.elect.lineEdit_valency_cation_salt_1.text()
        mass_cation_salt_1 = self.elect.lineEdit_mass_cation_salt_1.text()
        type_anion_salt_1 = str(self.elect.lineEdit_type_anion_salt_1.text())
        charge_anion_salt_1 = self.elect.lineEdit_valency_anion_salt_1.text()
        mass_anion_salt_1 = self.elect.lineEdit_mass_anion_salt_1.text()
        type_counter = str(self.elect.lineEdit_type_counterion.text())
        mass_counter = str(self.elect.lineEdit_mass_counterion.text())
        water_diameter = self.elect.lineEdit_water_diameter.text()
        spacing_dia = self.elect.lineEdit_space_diameter.text()
        charge_cation_salt_1 = self.elect.lineEdit_valency_cation_salt_1.text()
        charge_anion_salt_1 = self.elect.lineEdit_valency_anion_salt_1.text()
        box_lx = box_ly = self.elect.SpinBox_lx.value()
        box_lz = self.elect.SpinBox_lz.value()
        box_lx = float(box_lx)
        box_ly = float(box_ly)
        box_lz = float(box_lz)

        try:
            valency_counterion = float(self.elect.lineEdit_valency_counterion.text())
        except ValueError:
            valency_counterion = 0.0
        try:
            spacing_dia = float(spacing_dia)
        except ValueError:
            spacing_dia = 0
        water_diameter = float(water_diameter)
        spacing =  spacing_dia *  water_diameter

        try:
            charge_cation_salt_1 = int(charge_cation_salt_1)
        except ValueError:
            charge_cation_salt_1 = 0
        try:
            charge_anion_salt_1 = int(charge_anion_salt_1)
        except ValueError:
            charge_anion_salt_1 = 0
        try:
            con_cation_salt_1 = float(self.elect.SpinBox_con_cation_salt_1.value())
        except TypeError:
            con_cation_salt_1 = 0
        con_anion_salt_1 = con_cation_salt_1
        volume_box =  box_lx* box_ly* box_lz*0.001
        total_cation_salt_1 =int((con_cation_salt_1 * 0.6022) * (volume_box))
        try:
            if ((total_cation_salt_1 % charge_anion_salt_1) !=0):
                total_cation_salt_1 = int(abs(total_cation_salt_1 - (total_cation_salt_1 % charge_anion_salt_1) + charge_anion_salt_1))
            total_anion_salt_1 = int(abs(charge_cation_salt_1 * total_cation_salt_1 / charge_anion_salt_1))
            # con_anion_salt_1  = total_anion_salt_1 / (0.6022 * (volume_box))
        except ZeroDivisionError:
            charge_cation_salt_1 = total_cation_salt_1 = total_anion_salt_1 = 0
        self.elect.lineEdit_con_anion_salt_1.setText(str("{:.1f}".format(con_anion_salt_1)))
        self.elect.lineEdit_num_cation_salt_1.setText(str(total_cation_salt_1))


        #test if the system is electroneutral:
        total_charge_ions = int((charge_cation_salt_1 * total_cation_salt_1) + (charge_anion_salt_1 * total_anion_salt_1))
        if (total_charge_ions != 0):
            print('The electrolyte is not electroneutral; Abortion', total_charge_ions)
            return

        charge_density = self.elect.lineEdit_surface_charge.text()
        try:
            charge_density = float(charge_density)
        except ValueError:
            charge_density = 0
            type_counter = None
            mass_counterion = 0.0
            valency_counterion = 0.0
            counterions = 0

        self.elect.lineEdit_type_counter_final.setText(type_counter)
        unitcharge = 1.60217646 * pow(10.0 , -19)
        total_surface_charge = int( charge_density *  box_lx *  box_ly / 16) # lx and ly are in A
        surface_area =  box_lx *  box_ly * 0.01 * pow(10.0 , -18) # in unit of squared meter;
        try:
            counterions = int(2.0 * abs(total_surface_charge/valency_counterion))
        except ZeroDivisionError:
            counterions = 0

        print("total surface charge ", total_surface_charge)
        charge_system = valency_counterion * counterions + total_surface_charge * 2.0

        print("total charge in the system ", charge_system)
        if (charge_system == 0):
            print("system is charge neutral")
        else:
            print("system is not electroneutral; aborting..., charge_system is", charge_system)
           # self.show_warning_message_salt()
            return
        total_saltions_inside = int(ShM.total_cation_salt_1 + ShM.total_anion_salt_1)
        total_saltions_inside = total_saltions_inside + counterions
        try:
            total_water_inside = (int(box_lx/ spacing) * int(box_ly/ spacing) * int(box_lz/ spacing))
            water_amount = float((box_lx/ spacing) *  (box_ly/ spacing) *  (box_lz/ spacing))
            water_concentration = water_amount / (0.602214076 * (box_lx) * (box_ly) * (box_lz- water_diameter) * 0.001)
        except ZeroDivisionError:
            total_water_inside = water_amount = water_concentration = 0
        if total_water_inside > 0 :
            total_water_inside = total_water_inside - total_saltions_inside

        print("total_water_inside:  ", total_water_inside)
        print("total counterions is: ", counterions)
        print("total surface charge density is: ",  charge_density, " C.m^-2")
        #   long double
        if ( box_lx== 0 or  box_ly== 0 or  box_lz == 0):
            final_water_concentration = 0
        else:
            final_water_concentration = (total_water_inside) /(0.6022 * (box_lx) * (box_ly) * (box_lz- water_diameter) * 0.001)
        print("final water conc (in M): ", final_water_concentration)
        try:
            num_molecul_in_lx = int( box_lx /  spacing)
            num_molecul_in_ly = int( box_ly /  spacing)
            num_molecul_in_lz = int( box_lz /  spacing)
        except ZeroDivisionError:
            num_molecul_in_lx = num_molecul_in_ly = num_molecul_in_lz = 0


        if self.elect.frame_dummy_wall.isEnabled():
            pass
        else:
            ShM.wallR = np.array([])
            charge_density = ShM.charge_meshpoint = counterions = valency_counterion = 0
        if self.elect.frame_water_properties.isEnabled():
            pass
        else:
            total_water_inside = num_type_atom_water = final_water_concentration = 0
        coordinates_ions = []
        for i in range(num_molecul_in_lx):
            for j in range(num_molecul_in_ly):
                for k in range(num_molecul_in_lz):
                    if (len(coordinates_ions) < (total_saltions_inside + total_water_inside)):
                        x = float((- box_lx/2 + (0.5* spacing)) + i *  spacing)
                        y = float((- box_ly/2 + (0.5* spacing)) + j *  spacing)
                        z = float((- box_lz/2 + (0.5* spacing)) + k *  spacing)
                        if (x > (box_lx/2 - (0.5 *  spacing)) or y > ( box_ly/2 - (0.5 *  spacing)) or z > ( box_lz/2 - (0.5 *  spacing))):
                           continue
                        position = np.array([x, y, z])
                        coordinates_ions.append(position)
        ions = np.array(coordinates_ions)

        total_atoms = (total_water_inside * 3) + total_saltions_inside + (2 * len(ShM.wallR))
        if total_water_inside > 0:
            num_type_atom_water = 2
        else:
            num_type_atom_water = 0
        if counterions > 0:
            num_type_counterion = 1
        else:
            num_type_counterion = 0
        if (len(ShM.wallR) > 0):
            num_type_wall = 2
        else:
            num_type_wall = 0
        if con_anion_salt_1 > 0:
            num_type_cation_salt_1 = num_type_anion_salt_1 = 1
        else:
            num_type_cation_salt_1 = num_type_anion_salt_1 = 0


        total_num_salt_types = num_type_cation_salt_1 + num_type_anion_salt_1
        total_type_atom_in_system = total_num_salt_types + num_type_atom_water + num_type_counterion + num_type_wall
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
        print("total mesh points on each surface: ", len(ShM.wallR))
        print("charge on each meshpoint ", ShM.charge_meshpoint, " e")
        charge_system = ShM.charge_meshpoint * 2 * len(ShM.wallR) + counterions * valency_counterion
        print("recheck: total charge in the system ", counterions)
        if charge_system == 0:
            print("system is charge neutral again")
        else:
            print("system is not electroneutral; aborting...")
            self.show_warning_message_salt()
            return
        file_name, _ = QtWidgets.QFileDialog.getSaveFileName(self, self.tr('Save'), filter= 'LAMMPS (*.data)') #;;PDB (*.pdb);;GROMACS (*.gro);;XYZ (*.xyz)')

        # import tempfile
        # dir_name = tempfile.mkdtemp(prefix='Furious_Atoms_')
        # file_name = tempfile.mkstemp(suffix='.DATA', prefix='Electrolyte', dir=dir_name)[1]
        outdump = open(file_name, "w")
        outdump.write("LAMMPS data file\n\n")
        outdump.write("{}\t{} \n".format(total_atoms, ' atoms'))
        outdump.write("{}\t{} \n".format(total_water_inside * 2, ' bonds'))
        outdump.write("{}\t{} \n".format(total_water_inside, ' angles'))
        outdump.write("{}\t{} \n".format('0', ' dihedrals'))
        outdump.write("{}\t{} \n\n".format('0', ' impropers'))
        outdump.write("{}\t{} \n".format(total_type_atom_in_system, 'atom types'))
        if total_water_inside > 0:
            outdump.write("{}\t{} \n".format('1', 'bond types'))
            outdump.write("{}\t{} \n\n".format('1', 'angle types'))
        else:
            outdump.write("{}\t{} \n".format('0', 'bond types'))
            outdump.write("{}\t{} \n\n".format('0', 'angle types'))

        outdump.write("{}\t{}\t{} \n".format(-0.5* box_lx, 0.5* box_lx, ' xlo xhi'))
        outdump.write("{}\t{}\t{} \n".format(-0.5* box_ly, 0.5* box_ly, ' ylo yhi'))
        outdump.write("{}\t{}\t{} \n\n".format((-0.5* box_lz) - 0.0005, (0.5* box_lz) + 0.0005, ' zlo zhi'))
        outdump.write("Masses\n\n")

        if total_water_inside > 0:
            outdump.write("{}\t{} \n".format(type_hydrogen, mass_hydrogen))
            outdump.write("{}\t{} \n".format(type_oxygen, mass_oxygen))
        if counterions > 0:
            outdump.write("{}\t{} \n".format(type_counter, mass_counter))
        if total_cation_salt_1 > 0:
            outdump.write("{}\t{} \n".format(type_cation_salt_1, mass_cation_salt_1))
        if total_anion_salt_1 > 0:
            outdump.write("{}\t{} \n".format(type_anion_salt_1, mass_anion_salt_1))
        if len(ShM.wallR) > 0:
            outdump.write("{}\t{} \n".format(type_Rwall, mass_Rwall))
            outdump.write("{}\t{} \n".format(type_Lwall, mass_Lwall))

        outdump.write("\n")
        outdump.write("Atoms          # full\n\n")
        if total_water_inside > 0:
            num_molecules = 0
            for i in range(total_water_inside):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 1 + num_molecules + num_molecules), (num_molecules + 1), type_oxygen, charge_oxygen, ions[i][0], ions[i][1], (ions[i][2]), '0   0   0'))
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 2 + num_molecules+ num_molecules), (num_molecules + 1), type_hydrogen, charge_hydrogen, ions[i][0] +  0.95908, (ions[i][1] + (-0.02691)), (ions[i][2]) + 0.03231, ' 0   0   0 '))
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 3 + num_molecules+ num_molecules), (num_molecules + 1), type_hydrogen, charge_hydrogen, ions[i][0]+(-0.28004), (ions[i][1] + (-0.58767)), (ions[i][2]) + 0.70556, '0   0   0 '))
                num_molecules = num_molecules + 1


        if counterions > 0:
            num_molecules = 0
            for j in range(total_water_inside, (counterions+total_water_inside)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + total_water_inside*2), num_molecules, type_counter, valency_counterion, ions[j][0], ions[j][1], ions[j][2],  '0   0   0 '))

        # crystal pack of counter ions (in this system, it is the same as positive ions):
        num_molecules = 0
        if total_cation_salt_1 > 0:
            for j in range((counterions+total_water_inside), (counterions+total_water_inside+total_cation_salt_1)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + total_water_inside*2), num_molecules, type_cation_salt_1, charge_cation_salt_1, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))
        if total_anion_salt_1 > 0:
            for j in range((counterions+total_water_inside+total_cation_salt_1), (counterions+total_water_inside+total_cation_salt_1+total_anion_salt_1)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((j + 1 + total_water_inside*2), num_molecules, type_anion_salt_1, charge_anion_salt_1, ions[j][0], ions[j][1], ions[j][2], '0   0   0 '))


        if (len(ShM.wallR) > 0):
            for r  in range(len(ShM.wallR)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((r + 1 + (3 * total_water_inside) + total_saltions_inside),  "0" , type_Rwall, ShM.charge_meshpoint, ShM.wallR[r][0], ShM.wallR[r][1], ShM.wallR[r][2], '0   0   0 '))
            # mesh points on left wall; the same as right wall with oppositee sign in z direction
            for l  in range(len(ShM.wallR)):
                outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((l + 1 + (3 * total_water_inside) + total_saltions_inside + len(ShM.wallR)),  "0" , type_Lwall, ShM.charge_meshpoint, ShM.wallR[l][0], ShM.wallR[l][1], ShM.wallR[l][2] * -1, '0   0   0 '))
        if total_water_inside > 0:
            outdump.write("\n")
            outdump.write("Bonds\n\n")
            num_molecules = 0
            for g in range(total_water_inside):
                outdump.write("{}\t{}\t{}\t{}\n".format((g + 1 + num_molecules), '1', (g + 1 + num_molecules + num_molecules), (g + 2 + num_molecules + num_molecules)))
                outdump.write("{}\t{}\t{}\t{}\n".format((g + 2 + num_molecules), '1', (g + 1 + num_molecules + num_molecules), (g + 3 + num_molecules + num_molecules)))
                num_molecules = num_molecules + 1
            outdump.write("\n")
            outdump.write("Angles\n\n")
            num_molecules = 0
            for n in range(total_water_inside):
                outdump.write("{}\t{}\t{}\t{}\t{}\n".format( n + 1, '1', (n + 2 + num_molecules + num_molecules), (n + 1 + num_molecules + num_molecules), (n + 3 + num_molecules + num_molecules)))
                num_molecules = num_molecules + 1

        window = self.win.create_mdi_child()
        window.load_file(fname=file_name)
        window.show()