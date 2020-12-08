import numpy as np


class SharedMemory(object):

    def __init__(self):

        self.cnt = 0
        self.enable_timer = None
        self.load_file = None
        self.no_vertices_per_particle = 0
        self.initial_vertices_particles = np.array([0, 0, 0.])
        self.n_frames = 0
        self.pos = None
        self.sphere_actor = None
        self.all_vertices_particles = np.array([0, 0, 0.])
        self.box = None
        self.MainWindow = 0
        self.no_atoms = 0
        self.vcolors_particle = 0
        self.colors_backup_particles = 0
        self.selected_particle = 0
        self.sec_particle = 0
        self.no_vertices_all_particles = 0
        self.particle_color_add = 0
        self.colors = 0
        self.no_unique_types_particles = 0
        self.colors_unique_types = 0
        self.radii_spheres = 0
        self.radii_unique_types = 0
        self.colors_particles = 0
        self.set_value_radius = 0
        self.unique_types = None
        self.metallicCoefficient_particle = 0
        self.roughnessCoefficient_particle =0


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
        self.bond_pos = 0
        self.line_thickness = 0
        self.atom_type = None
        self.box_colors = 0
        self.box_actor = None
        self.vcolors_box = 0
        self.colors_backup_box = 0
        self.colors_line = 0
        self.vcolors_line_actor = 0
        self.colors_backup_line = 0
        self.line_color_add = 0
        self.line_actor = None
        self.file_directory = None
        self.extension = None
        self.file_name = None
        self.unique_types_bond = 0
        self.center_bond = None
        self.bond_type = None