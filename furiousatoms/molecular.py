import numpy as np
from fury import utils, actor

from furiousatoms.structure import bbox
from fury import actor, molecular as mol

table = mol.PTable()

SPHERE_SIZE = 0.35
POS_DIM = 3
BOND_DIM = 2

class MolecularStructure:        
    def __init__(self, box_size: list, pos: np.array, bonds: np.array, atom_types: np.array):
        '''
        `box_size`: a 1x3 list with X, Y, and Z numerical dimensions. 
        `pos`: a Nx3 array with XYZ float positions of each atom.
        `bonds`: a Nx2 array where each row is a pair of integer atom IDs. 
            Each ID corresponds to a row in `pos`.
        `atom_types`: a Nx1 array of strings representing element symbols. 
            The length must match with `pos`.
        '''
        self.box_size = box_size
        self.pos = pos
        self.bonds = bonds
        self.atom_types = atom_types

    def create_empty():
        box_size = [0, 0, 0]
        pos = np.empty(shape=(0, 3))
        bonds = np.empty(shape=(0, 2), dtype=int)
        atom_types = np.empty(shape=(0), dtype=str)
        return MolecularStructure(box_size, pos, bonds, atom_types)

    def merge(s1, s2, offset_bonds=False):
        '''
        Merge a MolecularStructure with another, combining the positions, bonds, and atom_types.
        Only s1's box_size will be kept. Neither s1 nor s1 will be modified. 
        `offset_bonds`: If set to True, bond IDs will be corrected so that the two structures
            will not connect.
        '''
        pos_merged = np.zeros(shape=(len(s1.pos) + len(s2.pos), 3))
        i = 0
        for position in s1.pos:
            for j in range(3):
                pos_merged[i][j] = position[j]
            i += 1
        for position in s2.pos:
            for j in range(3):
                pos_merged[i][j] = position[j]
            i += 1

        offset = s1.pos.shape[0] if offset_bonds else 0
        bonds_merged = np.zeros(shape=(len(s1.bonds) + len(s2.bonds), 2), dtype=int)
        i = 0
        for bond in s1.bonds:
            for j in range(2):
                bonds_merged[i][j] = bond[j]
            i += 1
        for bond in s2.bonds:
            for j in range(2):
                bonds_merged[i][j] = bond[j] + offset
            i += 1

        types_merged = np.zeros(shape=(len(s1.atom_types) + len(s2.atom_types)), dtype=str)
        i = 0
        for typ in s1.atom_types:
            types_merged[i] = typ
            i += 1
        for typ in s2.atom_types:
            types_merged[i] = typ
            i += 1

        return MolecularStructure(s2.box_size, pos_merged, bonds_merged, types_merged)

    def center(self):
        for k in range(3):
            correction = (self.pos[:,k].max() + self.pos[:,k].min()) / 2
            self.pos[:, k] -= correction
        return self


class ViewerMemoryManager:

    def __init__(self, box, pos, bonds, atom_types):
        """
        """
        self.universe = self
        self.box_lx, self.box_ly, self.box_lz = box
        self.box_color = (0, 0, 0)
        self.bbox_actor, _ = bbox(self.box_lx, self.box_ly, self.box_lz,
                                  colors=self.box_color, linewidth=1, fake_tube=True)

        self.pos = pos
        self.no_atoms = pos.shape[0]
        self.bonds = bonds
        self.atom_types = atom_types
        self.no_bonds = bonds.shape[0]
        self.have_bonds = self.no_bonds > 0

        #Set colors
        self.colors = np.ones((self.no_atoms, 4))
        self.unique_types = np.unique(self.atom_types)
        self.colors_unique_types = np.random.rand(len(self.unique_types), 4)
        self.colors_unique_types[:, 3] = 1 #set opacity to 1 
        for i, typ in enumerate(self.atom_types):
            indexOfColor = np.where(self.unique_types == typ)
            self.colors[i] = self.colors_unique_types[indexOfColor]

        # set all radii to 1
        self.radii_spheres = np.ones((self.no_atoms)) * SPHERE_SIZE
        # but then switch to 0.2 for each type
        self.radii_unique_types = 0.4 + np.zeros(len(self.unique_types))
        self.selected_particle = np.zeros(self.no_atoms, dtype=np.bool)
        self.selected_bond = np.zeros(self.no_bonds, dtype=np.bool)

        self.sphere_actor = actor.sphere(centers=self.pos, colors=self.colors, radii=self.radii_spheres, phi=3, theta=6, use_primitive=False)
        self.bond_actor = self.generate_bond_actor() if self.have_bonds else None

        self.all_vertices_particles = utils.vertices_from_actor(self.sphere_actor)
        self.no_vertices_per_particle = len(self.all_vertices_particles) / self.no_atoms
        self.initial_vertices_particles = self.all_vertices_particles.copy() - \
            np.repeat(self.pos, self.no_vertices_per_particle, axis=0)
        vertices_particle = utils.vertices_from_actor(self.sphere_actor)
        self.no_vertices_all_particles = vertices_particle.shape[0]
        self.sec_particle = np.int(self.no_vertices_all_particles / self.no_atoms)
        self.vcolors_particle = utils.colors_from_actor(self.sphere_actor, 'colors')
        self.colors_backup_particles = self.vcolors_particle.copy()

        self.deleted_particles = np.zeros(self.no_atoms, dtype=np.bool)
        self.deleted_bonds = np.zeros(self.no_bonds, dtype=np.bool)

        # Animation Player
        self.cnt = 0
        # self.selected_value_radius = 0
        self.roughness = 0.01
        self.metallic = 0.5
        self.anisotropic = 0.01
        self.opacity = 1.0
        self.selected_value_radius = 0.01
        self.anisotropic_rot = 0.01
        self.anisotropic_X = 0.01
        self.anisotropic_Y = 0.01
        self.anisotropic_Z = 0.01
        self.coat_rough = 0.01
        self.coat_strength = 0.01
        self.pbr_params_atom = None
        self.pbr_params_bond = None
        self.particle_resolution = "Low"


    @property
    def n_frames(self):
        return 1

    @property
    def atom_type(self):
        return self.atom_types

    def get_bonds(self):
        if not self.have_bonds:
            return 0
        return self._bonds

    def generate_bond_actor(self):
        bonds_indices = self.bonds
        first_pos_bond = self.pos[(bonds_indices[:, 0])]
        second_pos_bond = self.pos[(bonds_indices[:, 1])]
        bonds = np.hstack((first_pos_bond, second_pos_bond))
        self._bonds = bonds.reshape(self.no_bonds, 2, 3)
        self._bonds_2 = np.zeros((self.no_bonds*2, 2, 3))
        self.bond_colors_2 = np.zeros(((self.no_bonds*2*2, 4)))
        # self.unique_types = np.unique(self.universe.atoms.types)
        self.unique_types = np.unique(self.atom_types)

        for i in range(self.no_bonds):
            p1, p2 = self._bonds[i]
            phalf = 0.5 * (p1 + p2)
            self._bonds_2[i * 2][0] = p1
            self._bonds_2[i * 2][1] = phalf
            self._bonds_2[i * 2 + 1][0] = phalf
            self._bonds_2[i * 2 + 1][1] = p2
            self.bond_colors_2[i * 4] = self.colors[bonds_indices[i][0]]
            self.bond_colors_2[i * 4 + 1] = self.colors[bonds_indices[i][0]]
            self.bond_colors_2[i * 4 + 2] = self.colors[bonds_indices[i][1]]
            self.bond_colors_2[i * 4 + 3] = self.colors[bonds_indices[i][1]]

        self.line_thickness = 0.2
        bond_actor = actor.streamtube(self._bonds_2, self.bond_colors_2, linewidth=self.line_thickness, tube_sides=4,
                                      lod=False, replace_strips=True)
        self.all_vertices_bonds = utils.vertices_from_actor(bond_actor)
        self.no_vertices_per_bond = len(self.all_vertices_bonds) / (2 * self.no_bonds)
        self.no_vertices_all_bonds = self.all_vertices_bonds.shape[0]
        self.sec_bond = np.int(self.no_vertices_all_bonds / (2 * self.no_bonds))
        self.vcolors_bond = utils.colors_from_actor(bond_actor, 'colors')
        self.colors_backup_bond = self.vcolors_bond.copy()
        return bond_actor

    def actors(self):
        l_actors = [self.sphere_actor, self.bbox_actor]
        if self.have_bonds:
            l_actors += [self.bond_actor, ]

        return l_actors
