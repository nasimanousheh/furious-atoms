import numpy as np
from fury import utils, actor, primitive

from furiousatoms.structure import bbox
import MDAnalysis


# TODO: I have to replace this class with the ViewerMemoryManager class
# that does not depend on universe (md analysis).
class UniverseManager:

    def __init__(self, universe, no_bonds=0):
        """
        """
        self.universe = universe
        self.bbox_actor, _ = bbox(self.box_lx, self.box_ly, self.box_lz,
                                  colors=(0, 0, 0), linewidth=1, fake_tube=True)
        # Anispmation Player
        self.cnt = 0
        if self.n_frames > 1:
            pos_R = self.universe.trajectory[self.cnt].positions.copy().astype('f8')
            self._pos = MDAnalysis.lib.distances.transform_RtoS(pos_R, self.box, backend='serial')
        else:
            self._pos = self.universe.trajectory[0].positions.copy().astype(float)

        try:
            self._bonds = self.universe.bonds.to_indices()
            no_bonds = len(self.universe.bonds)
            self.selected_bond = np.zeros(no_bonds, dtype=np.bool)
        except:
            pass
        self.have_bonds = no_bonds > 0
        self.bond_actor = self.generate_bond_actor() if self.have_bonds else None
        colors = np.ones((self.no_atoms, 4))
        self.unique_types = np.unique(self.universe.atoms.types)
        self.colors_unique_types = np.random.rand(len(self.unique_types), 4)
        self.colors_unique_types[:, 3] = 1
        for i, typ in enumerate(self.unique_types):
            colors[self.atom_type == typ] = self.colors_unique_types[i]

        # TODO: the two variables below take a fixed value. Maybe this should be provided by atom_type
        # self.radii_spheres = np.ones((self.no_atoms))
        self.radii_spheres = 0.4 * np.ones((self.no_atoms))
        # self.radii_unique_types = np.zeros(len(self.unique_types))
        self.radii_unique_types = 0.4 + np.zeros(len(self.unique_types))
        self.selected_particle = np.zeros(self.no_atoms, dtype=np.bool)

        # self.sphere_actor = actor.sphere(centers=self.pos, colors=colors,
        #                                  radii=self.radii_spheres, theta=32, phi=32)
        vertices, faces = primitive.prim_sphere(name='repulsion724', gen_faces=False)
        res = primitive.repeat_primitive(vertices, faces, centers=self.pos, colors=colors, scales=self.radii_spheres)#, dtype='uint8')
        big_verts, big_faces, big_colors, _ = res
        self.sphere_actor = utils.get_actor_from_primitive(big_verts, big_faces, big_colors)
        self.all_vertices_particles = utils.vertices_from_actor(self.sphere_actor)
        self.no_vertices_per_particle = len(self.all_vertices_particles) / self.no_atoms
        self.initial_vertices_particles = self.all_vertices_particles.copy() - np.repeat(self.pos, self.no_vertices_per_particle, axis=0)
        vertices_particle = utils.vertices_from_actor(self.sphere_actor)
        self.no_vertices_all_particles = vertices_particle.shape[0]
        self.sec_particle = np.int(self.no_vertices_all_particles / self.no_atoms)
        self.vcolors_particle = utils.colors_from_actor(self.sphere_actor, 'colors')
        self.colors_backup_particles = self.vcolors_particle.copy()
        # TODO: check if this is truly used
        self.set_value_radius = 0
        self.play_factor = 0
        self.enable_timer = True
        self.cnt = 0
        self.metallicCoefficient_particle = 0
        self.roughnessCoefficient_particle = 0
        self.degree_of_metallicity = 0

    @property
    def no_atoms(self):
        return len(self.universe.atoms)

    @property
    def n_frames(self):
        return self.universe.trajectory.n_frames

    @property
    def box(self):
        return self.universe.trajectory.ts.dimensions

    @property
    def box_lx(self):
        return self.universe.trajectory.ts.dimensions[0]

    @property
    def box_ly(self):
        return self.universe.trajectory.ts.dimensions[1]

    @property
    def box_lz(self):
        return self.universe.trajectory.ts.dimensions[2]

    @property
    def no_unique_types_particles(self):
        return len(np.unique(self.universe.atoms.types))

    @property
    def atom_type(self):
        return self.universe.atoms.types

    @property
    def pos(self):
        return self._pos

    @property
    def no_bonds(self):
        try:
            number_bonds = len(self.universe.bonds)
        except:
            number_bonds = 0
        return number_bonds

    def get_bonds(self):
        if not self.have_bonds:
            return 0
        return self._bonds

    def generate_bond_actor(self):
        bonds_indices = self.universe.bonds.to_indices()
        first_pos_bond = self.pos[(bonds_indices[:, 0])]
        second_pos_bond = self.pos[(bonds_indices[:, 1])]
        bonds = np.hstack((first_pos_bond, second_pos_bond))
        self._bonds = bonds.reshape((self.no_bonds), 2, 3)
        self.bond_colors = (0.8275, 0.8275, 0.8275, 1)
        self.line_thickness = 0.2
        bond_actor = actor.streamtube(self._bonds, self.bond_colors, linewidth=self.line_thickness)
        self.colors_backup_bond = utils.colors_from_actor(bond_actor, 'colors').copy()
        self.all_vertices_bonds = utils.vertices_from_actor(bond_actor)
        self.no_vertices_per_bond = len(self.all_vertices_bonds) / self.no_bonds
        self.no_vertices_all_bonds = self.all_vertices_bonds.shape[0]
        self.sec_bond = np.int(self.no_vertices_all_bonds / self.no_bonds)
        self.unique_types_bond = np.unique(self.universe.bonds.types)
        self.vcolors_bond = utils.colors_from_actor(bond_actor, 'colors')
        self.colors_backup_bond = self.vcolors_bond.copy()
        return bond_actor

    def actors(self):
        l_actors = [self.sphere_actor, self.bbox_actor]
        if self.have_bonds:
            l_actors += [self.bond_actor, ]

        return l_actors


class ViewerMemoryManager:

    def __init__(self, box, pos, bonds, atom_types):
        """
        """
        self.box_lx, self.box_ly, self.box_lz = box
        self.bbox_actor, _ = bbox(self.box_lx, self.box_ly, self.box_lz,
                                  colors=(0, 0, 0), linewidth=1, fake_tube=True)

        self.pos = pos
        self.no_atoms = pos.shape[0]
        self.bonds = bonds
        self.atom_types = atom_types
        self.no_bonds = bonds.shape[0]

        self.have_bonds = self.no_bonds > 0

        self.bond_actor = self.generate_bond_actor() if self.have_bonds else None

        colors = np.ones((self.no_atoms, 4))
        self.unique_types = np.unique(self.atom_types)
        self.colors_unique_types = np.random.rand(len(self.unique_types), 4)
        self.colors_unique_types[:, 3] = 1
        for i, typ in enumerate(self.unique_types):
            colors[self.atom_type == typ] = self.colors_unique_types[i]

        # set all radii to 1
        self.radii_spheres = np.ones((self.no_atoms))
        # but then switch to 0.2 for each type
        self.radii_unique_types = 0.4 + np.zeros(len(self.unique_types))
        self.selected_particle = np.zeros(self.no_atoms, dtype=np.bool)
        self.selected_bond = np.zeros(self.no_bonds, dtype=np.bool)

        # create the sphere particles actor using primitives
        vertices, faces = primitive.prim_sphere(name='repulsion724', gen_faces=False)
        res = primitive.repeat_primitive(vertices, faces, centers=self.pos, colors=colors, scales=self.radii_spheres)#, dtype='uint8')
        big_verts, big_faces, big_colors, _ = res
        self.sphere_actor = utils.get_actor_from_primitive(big_verts, big_faces, big_colors)


        self.all_vertices_particles = utils.vertices_from_actor(self.sphere_actor)
        self.no_vertices_per_particle = len(self.all_vertices_particles) / self.no_atoms
        self.initial_vertices_particles = self.all_vertices_particles.copy() - \
            np.repeat(self.pos, self.no_vertices_per_particle, axis=0)
        vertices_particle = utils.vertices_from_actor(self.sphere_actor)
        self.no_vertices_all_particles = vertices_particle.shape[0]
        self.sec_particle = np.int(self.no_vertices_all_particles / self.no_atoms)
        self.vcolors_particle = utils.colors_from_actor(self.sphere_actor, 'colors')
        self.colors_backup_particles = self.vcolors_particle.copy()
        # TODO: check if this is truly used
        self.set_value_radius = 0
        # Animation Player
        self.cnt = 0
        self.metallicCoefficient_particle = 0
        self.roughnessCoefficient_particle = 0


    @property
    def n_frames(self):
        return 1

    @property
    def atom_type(self):
        return ['C']

    def get_bonds(self):
        if not self.have_bonds:
            return 0
        return self._bonds

    def generate_bond_actor(self):
        # bonds_indices = self.universe.bonds.to_indices()
        bonds_indices = self.bonds
        first_pos_bond = self.pos[(bonds_indices[:, 0])]
        second_pos_bond = self.pos[(bonds_indices[:, 1])]
        bonds = np.hstack((first_pos_bond, second_pos_bond))
        self._bonds = bonds.reshape((self.no_bonds), 2, 3)
        self.bond_colors = (0.8275, 0.8275, 0.8275, 1)
        self.line_thickness = 0.2
        bond_actor = actor.streamtube(self._bonds, self.bond_colors, linewidth=self.line_thickness)
        self.colors_backup_bond = utils.colors_from_actor(bond_actor, 'colors').copy()
        self.all_vertices_bonds = utils.vertices_from_actor(bond_actor)
        self.no_vertices_per_bond = len(self.all_vertices_bonds) / self.no_bonds
        self.no_vertices_all_bonds = self.all_vertices_bonds.shape[0]
        self.sec_bond = np.int(self.no_vertices_all_bonds / self.no_bonds)
        self.unique_types_bond = np.unique(self.bond_types)
        return bond_actor

    def actors(self):
        l_actors = [self.sphere_actor, self.bbox_actor]
        if self.have_bonds:
            l_actors += [self.bond_actor, ]

        return l_actors
