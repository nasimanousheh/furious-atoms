from furiousatoms import io
import numpy as np
from fury import window, actor, utils, pick, ui
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QIcon
from PySide2 import QtWidgets
import sys
from io import StringIO
from functools import cmp_to_key

class Ui_warning_build_elect(QtWidgets.QMainWindow):
    """ Ui_warning_message class
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_warning_build_elect, self).__init__(parent)
        self.msg = io.load_ui_widget("warning_message.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.msg)
        self.setCentralWidget(self.msg)
        self.setLayout(self.v_layout)
        self.resize(657, 49)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.msg.pushButton_warning.clicked.connect(lambda:self.close())


class Ui_warning_oppos_charge(QtWidgets.QMainWindow):
    """ Ui_warning_message class
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_warning_oppos_charge, self).__init__(parent)
        self.msg = io.load_ui_widget("warning_opposit_charge.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.msg)
        self.setCentralWidget(self.msg)
        self.setLayout(self.v_layout)
        self.resize(622, 72)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.msg.pushButton_Ok.clicked.connect(self.open_period_table)
        self.msg.pushButton_Ok.clicked.connect(lambda:self.close())

    def open_period_table(self):
        self.msg_coun.win.elect.pushButton_type_counterion.click()
        # self.msg.close()


class Ui_warning_charge_salt_1(QtWidgets.QMainWindow):
    """ Ui_warning_message class
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_warning_charge_salt_1, self).__init__(parent)
        self.msg = io.load_ui_widget("warning_opposit_charge.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.msg)
        self.setCentralWidget(self.msg)
        self.setLayout(self.v_layout)
        self.resize(622, 72)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.msg.pushButton_Ok.clicked.connect(self.open_period_table)
        self.msg.pushButton_Ok.clicked.connect(lambda:self.close())

    def open_period_table(self):
        self.msg_coun.win.elect.pushButton_type_cation_salt_1.click()




class Ui_warning_atom_delete(QtWidgets.QMainWindow):
    """ Ui_warning_message class
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_warning_atom_delete, self).__init__(parent)
        self.msg = io.load_ui_widget("warning_message_atom_selection.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.msg)
        self.setCentralWidget(self.msg)
        self.setLayout(self.v_layout)
        self.resize(460, 80)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.msg.pushButton_warning.clicked.connect(lambda:self.close())
        self.msg.pushButton_Ok.clicked.connect(lambda:self.close())

class Ui_warning_bond_delete(QtWidgets.QMainWindow):
    """ Ui_warning_message class
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_warning_bond_delete, self).__init__(parent)
        self.msg = io.load_ui_widget("warning_message_bond_selection.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.msg)
        self.setCentralWidget(self.msg)
        self.setLayout(self.v_layout)
        self.resize(460, 80)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        self.msg.pushButton_warning.clicked.connect(lambda:self.close())
        self.msg.pushButton_Ok.clicked.connect(lambda:self.close())
