import numpy as np


class SharedMemory(object):

    def __init__(self):

        self.cnt = 0
        self.no_vertices_per_particle = 0
        self.initial_vertices_particles = 0
        self.n_frames = 0
        self.pos = 0
        self.sphere_actor = None
        self.all_vertices_particles = 0
        self.box = 0
        self.MainWindow = 0
        self.no_atoms = 0
        self.vcolors_particle = 0
        self.colors_backup = 0
        self.selected_particle = 0
        self.sec_particle = 0
        self.no_vertices_all_particles = 0
        self.particle_color_add = 0
        self.colors = 0
        self.colors_unique_types = 0

        self.bonds = 0
        self.bond_actor = None
        self.no_bonds = 0
        self.selected_bond = 0
        self.colors_backup_bond = 0
        self.vcolors_bond = 0
        self.no_vertices_per_bond = 0
        self.initial_vertices_bonds = 0
        self.all_vertices_bonds = 0
        self.sec_bond = 0
        self.no_vertices_all_bonds = 0
