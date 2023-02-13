from furiousatoms import io
from fury import window
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
from furiousatoms.sharedmem import SharedMemory
import csv
from io import StringIO
from PySide2.QtGui import QIcon
from furiousatoms.element_lookup import elems_csv

SM = SharedMemory()
class Ui_periodic_cation(QtWidgets.QMainWindow):
    """ Ui_periodic_cation class creates a widget for building periodic_cation table of elements
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_periodic_cation, self).__init__(parent)
        self.periodic_cation = io.load_ui_widget("periodic_cation_table.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.periodic_cation)
        self.setCentralWidget(self.periodic_cation)
        self.setLayout(self.v_layout)
        self.resize(1074, 488)
        self.setWindowIcon(QIcon(io.get_resources_file("splash.png")))
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.periodic_cation.buttonGroup_all_elements.buttonClicked.connect(self.get_element_info)
        self.periodic_cation.buttonGroup_all_elements.buttonClicked.connect(lambda:self.close())


    def get_element_info(self, button):
        SM.info_element = (button.text())
        spl_word = '\n'
        SM.info_element_number = SM.info_element.partition(spl_word)[0]
        SM.info_element_symbol = SM.info_element.partition(spl_word)[2]
        elems = [row for row in csv.DictReader(StringIO(elems_csv))]
        for row in elems:
            for key in row:
                try:
                    row[key] = int(row[key])
                    continue
                except ValueError:
                    pass
                try:
                    row[key] = float(row[key])
                except ValueError:
                    pass
        elems_symbol = {row['symbol']: row for row in elems}
        elems_number = {row['number']: row for row in elems}
        SM.all_information_element = elems_symbol[SM.info_element_symbol]
        self.current_edit_symbol.setText(str(SM.all_information_element['symbol']))
        self.current_edit_valency.setText(str(SM.all_information_element['valency']))
        self.current_edit_mass.setText(str(SM.all_information_element['mass']))
        return SM.all_information_element

class Ui_periodic_anion(QtWidgets.QMainWindow):
    """ Ui_periodic_anion class creates a widget for building periodic_anion table of elements
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_periodic_anion, self).__init__(parent)
        self.periodic_anion = io.load_ui_widget("periodic_anion_table.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.periodic_anion)
        self.setCentralWidget(self.periodic_anion)
        self.setLayout(self.v_layout)
        self.resize(1074, 488)
        self.scene = window.Scene()
        self.setWindowIcon(QIcon(io.get_resources_file("splash.png")))
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.periodic_anion.buttonGroup_all_elements.buttonClicked.connect(self.get_element_info)
        self.periodic_anion.buttonGroup_all_elements.buttonClicked.connect(lambda:self.close())


    def get_element_info(self, button):
        SM.info_element = (button.text())
        spl_word = '\n'
        SM.info_element_number = SM.info_element.partition(spl_word)[0]
        SM.info_element_symbol = SM.info_element.partition(spl_word)[2]
        elems = [row for row in csv.DictReader(StringIO(elems_csv))]
        for row in elems:
            for key in row:
                try:
                    row[key] = int(row[key])
                    continue
                except ValueError:
                    pass
                try:
                    row[key] = float(row[key])
                except ValueError:
                    pass
        elems_symbol = {row['symbol']: row for row in elems}
        elems_number = {row['number']: row for row in elems}
        SM.all_information_element = elems_symbol[SM.info_element_symbol]
        self.current_edit_symbol.setText(str(SM.all_information_element['symbol']))
        self.current_edit_valency.setText(str(SM.all_information_element['valency']))
        self.current_edit_mass.setText(str(SM.all_information_element['mass']))
        return SM.all_information_element