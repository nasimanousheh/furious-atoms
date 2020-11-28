import numpy as np
from fury import utils


def delete_bonds(bond_actor, selected_bond, no_bonds):
    vertices_bond = utils.vertices_from_actor(bond_actor)
    no_vertices_all_bond = vertices_bond.shape[0]
    object_indices_bonds = np.where(selected_bond is True)[0]
    sec_bond = np.int(no_vertices_all_bond / no_bonds)
    color_add_bond = np.array([255, 0, 0, 0], dtype='uint8')
    vcolors_bond = utils.colors_from_actor(bond_actor, 'colors')
    for object_index_bond in object_indices_bonds:
        vcolors_bond[object_index_bond * sec_bond: object_index_bond * sec_bond + sec_bond] = color_add_bond
    utils.update_actor(bond_actor)
    bond_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()
    print('The bond is deleted')


def delete_particles(sphere_actor, no_atoms, selected):
    vertices = utils.vertices_from_actor(sphere_actor)
    no_vertices_all_sphere = vertices.shape[0]
    object_indices_spheres = np.where(selected is True)[0]
    sec = np.int(no_vertices_all_sphere / no_atoms)
    color_add = np.array([255, 0, 0, 0], dtype='uint8')
    vcolors = utils.colors_from_actor(sphere_actor, 'colors')
    for object_index in object_indices_spheres:
        vcolors[object_index * sec: object_index * sec + sec] = color_add
    utils.update_actor(sphere_actor)
    sphere_actor.GetMapper().GetInput().GetPointData().GetArray('colors').Modified()

    print('The particle is deleted')