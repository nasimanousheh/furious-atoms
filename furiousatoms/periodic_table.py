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
import sys
import csv
from io import StringIO
from functools import cmp_to_key
elems_csv = """\
number,symbol,name,vdw radius,covalent radius,mass,ionization energy
1,H,hydrogen,1.2,0.38,1.0079,13.5984
2,He,helium,1.4,0.32,4.0026,24.5874
3,Li,lithium,1.82,1.34,6.941,5.3917
4,Be,beryllium,1.53,0.9,9.0122,9.3227
5,B,boron,1.92,0.82,10.811,8.298
6,C,carbon,1.7,0.77,12.0107,11.2603
7,N,nitrogen,1.55,0.75,14.0067,14.5341
8,O,oxygen,1.52,0.73,15.9994,13.6181
9,F,fluorine,1.47,0.71,18.9984,17.4228
10,Ne,neon,1.54,0.69,20.1797,21.5645
11,Na,sodium,2.27,1.54,22.9897,5.1391
12,Mg,magnesium,1.73,1.3,24.305,7.6462
13,Al,aluminium,1.84,1.18,26.9815,5.9858
14,Si,silicon,2.1,1.11,28.0855,8.1517
15,P,phosphorus,1.8,1.06,30.9738,10.4867
16,S,sulfur,1.8,1.02,32.065,10.36
17,Cl,chlorine,1.75,0.99,35.453,12.9676
18,Ar,argon,1.88,0.97,39.948,15.7596
19,K,potassium,2.75,1.96,39.0983,4.3407
20,Ca,calcium,2.31,1.74,40.078,6.1132
21,Sc,scandium,2.11,1.44,44.9559,6.5615
22,Ti,titanium,,1.36,47.867,6.8281
23,V,vanadium,,1.25,50.9415,6.7462
24,Cr,chromium,,1.27,51.9961,6.7665
25,Mn,manganese,,1.39,54.938,7.434
26,Fe,iron,,1.25,55.845,7.9024
27,Co,cobalt,,1.26,58.9332,7.881
28,Ni,nickel,1.63,1.21,58.6934,7.6398
29,Cu,copper,1.4,1.38,63.546,7.7264
30,Zn,zinc,1.39,1.31,65.39,9.3942
31,Ga,gallium,1.87,1.26,69.723,5.9993
32,Ge,germanium,2.11,1.22,72.64,7.8994
33,As,arsenic,1.85,1.19,74.9216,9.7886
34,Se,selenium,1.9,1.16,78.96,9.7524
35,Br,bromine,1.85,1.14,79.904,11.8138
36,Kr,krypton,2.02,1.1,83.8,13.9996
37,Rb,rubidium,3.03,2.11,85.4678,4.1771
38,Sr,strontium,2.49,1.92,87.62,5.6949
39,Y,yttrium,,1.62,88.9059,6.2173
40,Zr,zirconium,,1.48,91.224,6.6339
41,Nb,niobium,,1.37,92.9064,6.7589
42,Mo,molybdenum,,1.45,95.94,7.0924
43,Tc,technetium,,1.56,98,7.28
44,Ru,ruthenium,,1.26,101.07,7.3605
45,Rh,rhodium,,1.35,102.9055,7.4589
46,Pd,palladium,1.63,1.31,106.42,8.3369
47,Ag,silver,1.72,1.53,107.8682,7.5762
48,Cd,cadmium,1.58,1.48,112.411,8.9938
49,In,indium,1.93,1.44,114.818,5.7864
50,Sn,tin,2.17,1.41,118.71,7.3439
51,Sb,antimony,2.06,1.38,121.76,8.6084
52,Te,tellurium,2.06,1.35,127.6,9.0096
53,I,iodine,1.98,1.33,126.9045,10.4513
54,Xe,xenon,2.16,1.3,131.293,12.1298
55,Cs,caesium,3.43,2.25,132.9055,3.8939
56,Ba,barium,2.68,1.98,137.327,5.2117
57,La,lanthanum,,1.69,138.9055,5.5769
58,Ce,cerium,,,140.116,5.5387
59,Pr,praseodymium,,,140.9077,5.473
60,Nd,neodymium,,,144.24,5.525
61,Pm,promethium,,,145,5.582
62,Sm,samarium,,,150.36,5.6437
63,Eu,europium,,,151.964,5.6704
64,Gd,gadolinium,,,157.25,6.1501
65,Tb,terbium,,,158.9253,5.8638
66,Dy,dysprosium,,,162.5,5.9389
67,Ho,holmium,,,164.9303,6.0215
68,Er,erbium,,,167.259,6.1077
69,Tm,thulium,,,168.9342,6.1843
70,Yb,ytterbium,,,173.04,6.2542
71,Lu,lutetium,,1.6,174.967,5.4259
72,Hf,hafnium,,1.5,178.49,6.8251
73,Ta,tantalum,,1.38,180.9479,7.5496
74,W,tungsten,,1.46,183.84,7.864
75,Re,rhenium,,1.59,186.207,7.8335
76,Os,osmium,,1.28,190.23,8.4382
77,Ir,iridium,,1.37,192.217,8.967
78,Pt,platinum,1.75,1.28,195.078,8.9587
79,Au,gold,1.66,1.44,196.9665,9.2255
80,Hg,mercury,1.55,1.49,200.59,10.4375
81,Tl,thallium,1.96,1.48,204.3833,6.1082
82,Pb,lead,2.02,1.47,207.2,7.4167
83,Bi,bismuth,2.07,1.46,208.9804,7.2856
84,Po,polonium,1.97,,209,8.417
85,At,astatine,2.02,,210,9.3
86,Rn,radon,2.2,1.45,222,10.7485
87,Fr,francium,3.48,,223,4.0727
88,Ra,radium,2.83,,226,5.2784
89,Ac,actinium,,,227,5.17
90,Th,thorium,,,232.0381,6.3067
91,Pa,protactinium,,,231.0359,5.89
92,U,uranium,1.86,,238.0289,6.1941
"""

class Ui_periodic(QtWidgets.QMainWindow): #QWidget
    """ Ui_periodic class creates a widget for building periodic table of elements
    """
    def __init__(self, app_path=None, parent=None):
        super(Ui_periodic, self).__init__(parent)
        self.periodic = io.load_ui_widget("periodic_table.ui")
        self.v_layout = QtWidgets.QVBoxLayout()
        self.v_layout.addWidget(self.periodic)
        self.setCentralWidget(self.periodic)
        self.setLayout(self.v_layout)
        self.resize(1074, 488)
        self.scene = window.Scene()
        self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        self.init_settings()
        self.create_connections()

    def init_settings(self):
        pass

    def create_connections(self):
        # pushButton_cation_1.clicked.connect(self.get_element)
        self.periodic.H.clicked.connect(self.get_element_info)
        self.periodic.buttonGroup_all_elements.checkedButton()

    def get_element_info(self):
        # Mass =
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
        print(elems_symbol['H'])
    #     self.periodic.lineEdit_bond_length_periodic.textChanged.connect(self.periodic_diameter_changed)
    #     self.periodic.spinBox_chirality_N_periodic.valueChanged.connect(self.periodic_diameter_changed)
    #     self.periodic.spinBox_chirality_M_periodic.valueChanged.connect(self.periodic_diameter_changed)
    #     self.periodic.spinBox_repeat_units_periodic.valueChanged.connect(self.periodic_diameter_changed)
    #     self.periodic.pushButton_build_periodic.clicked.connect(self.periodic_builder_callback)

    # def periodic_diameter_changed(self):
    #     # H_termination_periodic = self.periodic.comboBox_H_termination_periodic.currentText()
    #     SM.bond_length_periodic = self.periodic.lineEdit_bond_length_periodic.text()
    #     SM.bond_length_periodic = float(SM.bond_length_periodic)
    #     value_n_periodic = int(self.periodic.spinBox_chirality_N_periodic.text())
    #     value_m_periodic = int(self.periodic.spinBox_chirality_M_periodic.text())
    #     # SM.number_of_walls = int(self.periodic.spinBox_num_walls_periodic.text())
    #     repeat_units_periodic = int(self.periodic.spinBox_repeat_units_periodic.text())
    #     a1 = np.array((np.sqrt(3)*SM.bond_length_periodic, 0))
    #     a2 = np.array((np.sqrt(3)/2*SM.bond_length_periodic, -3*SM.bond_length_periodic/2))
    #     Ch = value_n_periodic*a1+value_m_periodic*a2
    #     d = gcd(value_n_periodic, value_m_periodic)
    #     dR = 3*d if (value_n_periodic-value_m_periodic) % (3*d) == 0 else d
    #     t1 = (2*value_m_periodic+value_n_periodic)//dR
    #     t2 = -(2*value_n_periodic+value_m_periodic)//dR
    #     T = t1*a1+t2*a2
    #     diameter_periodic = float(np.linalg.norm(Ch)/np.pi)
    #     diameter_periodic = "{:.2f}".format(diameter_periodic)
    #     length_periodic = np.linalg.norm(T) * repeat_units_periodic
    #     length_periodic = "{:.2f}".format(float(length_periodic))
    #     self.periodic.lineEdit_diameter_periodic.setText(str(diameter_periodic))
    #     self.periodic.lineEdit_length_periodic.setText(str(length_periodic))

    # def periodic_builder_callback(self):
    #     SM.number_of_walls = int(self.periodic.spinBox_num_walls_periodic.text())
    #     value_n_periodic = int(self.periodic.spinBox_chirality_N_periodic.text())
    #     value_m_periodic = int(self.periodic.spinBox_chirality_M_periodic.text())
    #     repeat_units_periodic = int(self.periodic.spinBox_repeat_units_periodic.text())
    #     periodic_type_1 = self.periodic.comboBox_type1_periodic.currentText()
    #     periodic_type_2 = self.periodic.comboBox_type2_periodic.currentText()
    #     i = 1
    #     universe = periodic_builder(value_n_periodic, value_m_periodic, repeat_units_periodic, length=None, a=SM.bond_length_periodic, species=(periodic_type_1, periodic_type_2), centered=True)
    #     for i in range(SM.number_of_walls):
    #         next_universe = periodic_builder(value_n_periodic + (6*i), value_m_periodic + (6*i), repeat_units_periodic, length=None, a=SM.bond_length_periodic, species=(periodic_type_1, periodic_type_2), centered=True)
    #         universe = MDAnalysis.Merge(universe.atoms, next_universe.atoms)
    #     file_name = 'fname.pdb'
    #     universe.atoms.write(file_name)
    #     self.win.process_load_file(fname=file_name)
