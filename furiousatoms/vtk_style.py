from fury import actor, utils, molecular as mol
from furiousatoms.molecular_info import get_default_molecular_info
from furiousatoms.structure import bbox



def box_vtk_style(SM, active_vtk_window):
    vtk_box_actor, _ = bbox(SM.box_lx, SM.box_ly, SM.box_lz, colors=(0, 0, 0), linewidth=2, fake_tube=True)
    utils.update_actor(vtk_box_actor)
    active_vtk_window.scene.add(vtk_box_actor)


def make_aesthetic(molecule_rep):
    molecule_rep.GetProperty().SetAmbient(0.2)
    molecule_rep.GetProperty().SetDiffuse(1)
    molecule_rep.GetProperty().SetSpecular(1)
    molecule_rep.GetProperty().SetSpecularPower(100.0)


def get_vtk_ribbon(self, SM, file_name, active_vtk_window):
    molecule, all_info = get_default_molecular_info(self, SM, file_name)
    if all_info is True:
        ribbon = mol.ribbon(molecule)
        box_vtk_style(SM, active_vtk_window)
        active_vtk_window.scene.add(ribbon)
        active_vtk_window.make_title()
        active_vtk_window.render()
        active_vtk_window.show()
    else:
        print('There is no ribbon structure')

def get_vtk_ball_stick(self, SM, file_name, active_vtk_window):
    molecule, _ = get_default_molecular_info(self, SM, file_name)
    if molecule.total_num_bonds > 0:
        ball_stick_rep = mol.ball_stick(molecule, atom_scale_factor=0.3,
                                bond_thickness=0.2)
        make_aesthetic(ball_stick_rep)
        box_vtk_style(SM, active_vtk_window)
        active_vtk_window.scene.add(ball_stick_rep)
        active_vtk_window.make_title()
        active_vtk_window.render()
        active_vtk_window.show()
    else:
        print('There is no bond in the structure')

def get_vtk_stick(self, SM, file_name, active_vtk_window):
    molecule, _ = get_default_molecular_info(self, SM, file_name)
    if molecule.total_num_bonds > 0:
        stick_rep = mol.stick(molecule)
        make_aesthetic(stick_rep)
        box_vtk_style(SM, active_vtk_window)
        active_vtk_window.scene.add(stick_rep)
        active_vtk_window.make_title()
        active_vtk_window.render()
        active_vtk_window.show()
    else:
        print('There is no bond in the structure')


def get_vtk_sphere(self, SM, file_name, active_vtk_window):
    molecule, _ = get_default_molecular_info(self, SM, file_name)
    if molecule.total_num_atoms > 0:
        sphere_cpk_rep = mol.sphere_cpk(molecule, colormode='discrete')
        make_aesthetic(sphere_cpk_rep)
        box_vtk_style(SM, active_vtk_window)
        active_vtk_window.scene.add(sphere_cpk_rep)
        active_vtk_window.render()
        active_vtk_window.show()
    else:
        print('There is no atom in the structure')


