import numpy.testing as npt
import numpy as np
from fury import window, actor, molecular as mol
from furiousatoms.io import  load_files
from furiousatoms.viewer3d import Viewer3D
import os
from PySide2 import QtWidgets
###############################################################################
# Creating a `PeriodicTable()` object to obtain atomic numbers from names of
# elements
table = mol.PTable()
class VTK_Rep():
    def __init__(self, app_path=None, parent=None):
        self.viewer = Viewer3D()
        # self.viewer.scene.add(actor.axes())
        self.viewer.show()
        # sys.exit(app.exec_())
        # super(Ui_SWNT, self).__init__(parent)
        # self.SWNT = io.load_ui_widget("SWNT.ui")
        # self.v_layout = QtWidgets.QVBoxLayout()
        # self.v_layout.addWidget(self.SWNT)
        # self.setCentralWidget(self.SWNT)
        # self.setLayout(self.v_layout)
        # self.resize(247, 285)
        # self.scene = window.Scene()
        # self.showm = window.ShowManager(scene=self.scene, order_transparent=True)
        # self.init_settings()
        self.get_default_molecular_info()
    # def create_connections(self):
    #     self.radioButton_Ribbon.toggled.connect(self.get_default_molecular_info)
    ###############################################################################
    # Creating empty lists which will be filled with atomic information as we
    # parse the pdb file.
    NumberOfAtoms = 0
    points = []
    elements = []
    atom_names = []
    model = []
    sheets = []
    helix = []
    residue_seq = []
    chain = []
    is_hetatm = []
    current_model_number = 1

    def get_default_molecular_info(self):#, fine_name, all_info=False):
        all_info=True
        fine_name = 'C:/Users/nasim/Downloads/4ury'
        # downloadurl = "https://files.rcsb.org/download/"
        pdbfn = fine_name + ".pdb"
        flag = 0

        # if not os.path.isfile(pdbfn):
        #     flag = 1
        #     url = downloadurl + pdbfn
        #     outfnm = os.path.join(pdbfn)
            # try:
            #     urllib.request.urlretrieve(url, outfnm)
            # except Exception:
                # print("Error in downloading the file!")
        ###############################################################################
        # Creating empty lists which will be filled with atomic information as we
        # parse the pdb file.
        NumberOfAtoms = 0
        atom_coords = []
        atomic_numbers = []
        atom_types = []
        model = []
        sheets = []
        helix = []
        residue_seq = []
        chain = []
        is_hetatm = []
        current_model_number = 1
        ###############################################################################
        # Parsing the pdb file for information about coordinates and atoms

        pdbfile = open(pdbfn, 'r')
        pdb_lines = pdbfile.readlines()
        for line in pdb_lines:
            line = line.split()
            try:
                if line[0] == 'ATOM' or line[0] == 'HETATM':
                    if line[-1] != 'H':
                        coorX, coorY, coorZ = float(line[6]), float(line[7]), \
                                            float(line[8])
                        resi = line[5]
                        current_chain = ord(line[4])
                        atom_coords += [[coorX, coorY, coorZ]]
                        residue_seq += [resi]
                        chain += [current_chain]
                        atomic_numbers += [table.atomic_number(line[-1])]
                        atom_types += [line[2]]
                        model += [current_model_number]
                        NumberOfAtoms += 1
                        if(line[0] == 'HETATM'):
                            is_hetatm += [1]
                        else:
                            is_hetatm += [0]
                if line[0] == 'SHEET':
                    start_chain = ord(line[5])
                    start_resi = int(line[6])
                    end_chain = ord(line[8])
                    end_resi = int(line[9])
                    r = [start_chain, start_resi, end_chain, end_resi]
                    sheets += [r]
                if line[0] == 'HELIX':
                    start_chain = ord(line[4])
                    start_resi = int(line[5])
                    end_chain = ord(line[7])
                    end_resi = int(line[8])
                    r = [start_chain, start_resi, end_chain, end_resi]
                    helix += [r]
            except Exception:
                continue

        atom_coords = np.array(atom_coords)
        residue_seq = np.array(residue_seq, dtype=int)
        chain = np.array(chain)
        atomic_numbers = np.array(atomic_numbers)
        atom_types = np.array(atom_types)
        model = np.array(model)
        sheets = np.array(sheets)
        helix = np.array(helix)
        is_hetatm = np.array(is_hetatm)
        axes_actor = actor.axes()
        self.viewer.scene.add(axes_actor)

        molecule = mol.Molecule(atomic_numbers, atom_coords, atom_types, model,
                                residue_seq, chain, sheets, helix, is_hetatm)
        mol.compute_bonding(molecule)

        # stick representation

        self.viewer.scene.add(mol.stick(molecule, bond_thickness=0.2))

        # ribbon representation
        self.viewer.scene.add(mol.ribbon(molecule))

        # ball and stick representation
        # scene.add(mol.ball_stick(molecule, atom_scale_factor=0.3,
        #                          bond_thickness=0.2))

        # bounding box
        self.viewer.scene.add(mol.bounding_box(molecule, linewidth=0.4))




