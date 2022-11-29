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
# add_shader_callback(utah, shader_callback)

screen_x_dim = 800
screen_y_dim = 300
dims = (screen_x_dim, screen_y_dim)
table = mol.PTable()
# class VTK_Rep():
#     def __init__(self, app_path=None, parent=None):
#         self.get_default_molecular_info()

def make_aesthetic(self, molecule_rep):
    molecule_rep.GetProperty().SetAmbient(0.2)
    molecule_rep.GetProperty().SetDiffuse(1)
    molecule_rep.GetProperty().SetSpecular(1)
    molecule_rep.GetProperty().SetSpecularPower(100.0)

def get_default_molecular_info(self, pdbfn):#, fine_name, all_info=False):
    all_info = True
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
    # creating a ShowManager object
    showm = window.ShowManager(size=dims)#, title=pdb_code)
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
    atomic_numbers = np.array(atomic_numbers)
    atom_types = np.array(atom_types)
    # model = np.array(model)
    # sheets = np.array(sheets)
    # helix = np.array(helix)
    # is_hetatm = np.array(is_hetatm)


    if len(helix) > 0:
        residue_seq = np.array(residue_seq, dtype=int)
        chain = np.array(chain)
        model = np.array(model)
        sheets = np.array(sheets)
        helix = np.array(helix)
        is_hetatm = np.array(is_hetatm)
        molecule = mol.Molecule(atomic_numbers, atom_coords, atom_types, model,
                                residue_seq, chain, sheets, helix, is_hetatm)
        mol.compute_bonding(molecule)
        mol.get_all_bond_orders(molecule)
        ribbon = mol.ribbon(molecule)
        showm.scene.add(ribbon)
    else:
        molecule = mol.Molecule(atomic_numbers, atom_coords, atom_types)

    mol.compute_bonding(molecule)
    mol.get_all_bond_orders(molecule)
    # if molecule.total_num_bonds > 0:
    #     ball_stick_rep = mol.ball_stick(molecule, atom_scale_factor=0.3,
    #                                     bond_thickness=0.2)
    #     self.make_aesthetic(ball_stick_rep)
    #     showm.scene.add(ball_stick_rep)

    if molecule.total_num_bonds > 0:
        stick = mol.stick(molecule, bond_thickness=0.2)
        showm.scene.add(stick)


    # shader_to_actor(ball_stick_rep, "vertex", decl_code=vertex_shader_code_decl)
    # shader_to_actor(ball_stick_rep, "vertex", impl_code=vertex_shader_code_impl, block="light")

    # shader_to_actor(ball_stick_rep, "fragment", decl_code=fragment_shader_code_decl)

    # shader_to_actor(ball_stick_rep, "fragment", impl_code=fragment_shader_code_impl,
    #                 block="light", debug=False)

    # def shader_callback(_caller, _event, calldata=None):
    #     program = calldata
    #     global timer
    #     if program is not None:
    #         try:
    #             program.SetUniformf("time", timer)
    #         except ValueError:
    #             pass


    # bounding box
    b_box = mol.bounding_box(molecule, colors=(0, 0.8, 1), linewidth=0.1)


    showm.scene.add(b_box)
    axes_actor = actor.axes()
    showm.scene.add(axes_actor)


    # Delete the PDB file.
    flag = 0
    if flag:
        os.remove(pdbfn)
    showm.start()





    # # stick representation
    # ball_stick_rep = mol.stick(molecule, bond_thickness=0.2)
    # # self.make_aesthetic(ball_stick_rep)
    # self.viewer.scene.add(ball_stick_rep)

    # # ribbon representation
    # self.viewer.scene.add(mol.ribbon(molecule))

    # # ball and stick representation
    # # scene.add(mol.ball_stick(molecule, atom_scale_factor=0.3,
    # #                          bond_thickness=0.2))

    # # bounding box
    # self.viewer.scene.add(mol.bounding_box(molecule, linewidth=0.4))




