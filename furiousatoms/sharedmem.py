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
        self.bond_center = 0



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
        self.bond_length_MWNT = 0

        self.play_factor = 0
        self.bonds_initial_structure = 0
        self.load_file_initial_structure = None

        #############File################
        self.universe = None
        self.H_termination_SWNT = None
        self.bond_length_SWNT = 0
        self.number_of_walls = 0
        ###########Electrolyte##########
        self.spacing_dia = 0
        self.water_diameter = 0
        self.spacing = 0
        self.box_lx = 0
        self.box_ly = 0
        self.box_lz = 0
        self.charge_density = 0.0
        self.all_information_element = 0
        self.valency_counterion = 0.0
        self.mass_counter = 0.0
        self.type_counter = None

        self.mass_cation_salt_1 = 0.0
        self.mass_anion_salt_1 = 0.0
        self.mass_cat_3 = 0.0
        self.mass_cat_4 = 0.0
        self.mass_cat_5 = 0.0
        self.mass_an_1 = 0.0
        self.mass_an_2 = 0.0
        self.mass_an_3 = 0.0
        self.mass_an_4 = 0.0
        self.mass_an_5 = 0.0
        self.charge_p_cation_salt_1 = 0.0
        self.charge_n_cation_salt_1 = 0.0
        self.charge_cat_3 = 0.0
        self.charge_cat_4 = 0.0
        self.charge_cat_5 = 0.0
        self.charge_an_1 = 0.0
        self.charge_an_2 = 0.0
        self.charge_an_3 = 0.0
        self.charge_an_4 = 0.0
        self.charge_an_5 = 0.0
        self.type_cat_1 =  None
        self.type_cat_2 =  None
        self.type_cat_3 =  None
        self.type_cat_4 =  None
        self.type_cat_5 =  None
        self.type_an_1 =  None
        self.type_an_2 =  None
        self.type_an_3 =  None
        self.type_an_4 =  None
        self.type_an_5 =  None
        self.type_cation_salt_1 = None
        self.type_anion_salt_1 = None
        self.total_ions_inside = 0
        self.total_ions_inside =0
        self.total_p_cation_salt_1 = 0
        self.total_n_cation_salt_1 = 0
        self.total_cation_salt_2 = 0
        self.total_anion_salt_2 = 0
        self.total_cation_salt_3 = 0
        self.total_anion_salt_3 = 0
        self.total_cation_salt_4 = 0
        self.total_anion_salt_4 = 0
        self.total_ions_inside = 0
        self.total_ions_concentration = 0
        self.wallR = None

        self.con_an_5 = 0

        self.total_surface_charge = 0
        self.counterions = 0
        self.info_element = None
        self.total_saltions_inside = 0
        self.charge_anion_salt_1 = 0
        self.charge_cation_salt_1 = 0
        self.total_anion_salt_1 = 0
        self.total_cation_salt_1 =0
        self.con_cation_salt_1 = 0
        self.con_cation_salt_2 = 0
        self.con_cation_salt_3 = 0
        self.con_cation_salt_4 = 0
        self.con_anion_salt_1 = 0
        self.con_anion_salt_2 = 0
        self.con_anion_salt_3 = 0
        self.con_anion_salt_4 = 0


