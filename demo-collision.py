import sys
import math
import numpy as np
import time
import sys
import os
import shutil
import numpy as np
from fury import window, actor, ui
from vtk.util import numpy_support
import itertools
from pdb import set_trace

class Vector3D:
    def __init__(self, initial_x = 0.0, initial_y = 0.0, initial_z = 0.0):
        self.x = initial_x
        self.y = initial_y
        self.z = initial_z

    def magnitude(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def sqd_magnitude(self):
        return self.x**2 + self.y**2 + self.z**2

    # Operator overloading for adding two vecs
    def __add__(self, v):
        return Vector3D(self.x + v.x, self.y + v.y, self.z + v.z)

    def __mul__(self, scalar):
        return Vector3D(self.x*scalar, self.y*scalar, self.z*scalar)

    def __sub__(self, v):
        return Vector3D(self.x - v.x, self.y - v.y, self.z - v.z)

    def __eq__(self, other):
        if isinstance(other, Vector3D):
            return self.x == other.x and self.y == other.y and self.z == other.z

    #printing overloaded
    def __str__(self):
        return "x=" + str(self.x) + ", y=" + str(self.y) + ", z=" + str(self.z)

    def apply_periodic_boundries(self, lx, ly, lz):
        # periodic boundary
        if (self.x > lx/2.0):
            self.x = self.x - lx
        if (self.x < -lx/2.0):
            self.x = self.x + lx
        if (self.y > ly/2.0):
            self.y = self.y - ly
        if (self.y < -ly/2.0):
            self.y = self.y + ly
        if (self.z > lz/2.0):
            self.z = self.z - lz
        if (self.z < -lz/2.0):
            self.z = self.z + lz

"""# Define Particle class"""

import math

class Particle:

    def __init__(self, id_value=0, diameter = 0.0, initial_m = 0.0,  initial_position = Vector3D(0.0, 0.0, 0.0), initial_velocity = Vector3D(0.0, 0.0, 0.0), initial_lx = 0.0, initial_ly = 0.0, initial_lz = 0.0):
        self.id = id_value
        self.m = initial_m
        self.d = diameter
        self.position = initial_position
        self.velocity = initial_velocity
        self.lx = initial_lx
        self.ly = initial_ly
        self.lz = initial_lz
        self.pe = 0.0
        self.ke = 0.0
        self.force = Vector3D(0.0, 0.0, 0.0)

    def volume(self):
        self.volume = (4.0/3.0) * math.pi * ((self.d/2.0)**3)

    #position updated to a full time-step
    def update_position(self, dt):
        self.position = self.position + (self.velocity * dt)
        # periodic boundary
        self.position.apply_periodic_boundries(self.lx, self.ly, self.lz)

    #velocity updated to a half time-step
    def update_velocity(self, dt):
        self.velocity = self.velocity + (self.force * ( 0.5 * dt / self.m))


"""# Define Box class"""

class Box:

    def __init__(self, initial_position = Vector3D(0.0, 0.0, 0.0), m_lx = 10.0, m_ly = 10.0, m_lz = 10.0):
        self.position = initial_position #origin point of the simulation box
        self.lx = m_lx                        #length of the box in x direction
        self.ly = m_ly                        #length of the box in y direction
        self.lz = m_lz                        #length of the box in z direction


    def put_ljatoms_inside(self, number_ljatom=0, ljatom_diameter=0.0, ljatom_mass=0.0, density=0.0):
        self.ljatom_diameter=ljatom_diameter
       # total_molecule_water # = 500
       # total_positive_ion = total_negative_ion = number_ljatom
        total_salt_ion = number_ljatom

        total_positive_ion = total_negative_ion = total_salt_ion/2
        a = 1
        num_atoms_linear_in_lx = int(self.lx / a)
        num_atoms_linear_in_ly = int(self.ly / a)
        num_atoms_linear_in_lz = int(self.lz / a)
        self.ljatom = []
        for i in range(num_atoms_linear_in_lx):
            for j in range(num_atoms_linear_in_ly):
                for k in range(num_atoms_linear_in_lz):
                    if len(self.ljatom) < (total_salt_ion):
                        x = (-self.lx/2) + a + i*a
                        y = (-self.ly/2) + a + j*a
                        z = (-self.lz/2) + a + k*a
                        if (z > (self.lz/2.0 - a/2.0) or y > (self.ly/2.0 - a/2.0 ) or x > (self.lx/2.0 - a/2.0)):
                            continue
                        position = Vector3D(x,y,z)
                        velocity = Vector3D(0.0,0.0,0.0)
                        fresh_atom = Particle(len(self.ljatom)+1, self.ljatom_diameter, ljatom_mass, position, velocity, self.lx, self.ly, self.lz)
                        self.ljatom.append(fresh_atom)

        self.make_movie(num=0, file_name="final_ip.lammpstrj")

        return len(self.ljatom)


    def animation(self):
        lammps_file = open("many_particle_data/propagation.lammpstrj", "r")
        lines = lammps_file.readlines()
        no_atoms = int(lines[3])
        line_lx = list(map(float, lines[5].split('\n')[0].split('\t')))
        line_ly = list(map(float, lines[6].split('\n')[0].split('\t')))
        line_lz = list(map(float, lines[7].split('\n')[0].split('\t')))

        box_lx = (np.abs(line_lx[0])+np.abs(line_lx[1]))
        box_ly = (np.abs(line_ly[0])+np.abs(line_ly[1]))
        box_lz = (np.abs(line_lz[0])+np.abs(line_lz[1]))
        header = 9
        step = 0
        no_frames = 0

        # all_nframes = np.empty((0, 6), dtype=np.float)
        all_nframes = []

        while True:
            frame = lines[header + step: header + step + no_atoms]

            numerical_frame = np.zeros((len(frame), 5))
            for i in range(len(frame)):
                numerical_frame[i] = list(map(float, frame[i].split('\n')[0].split('\t')))
            all_nframes.append(numerical_frame)
            step += header + no_atoms
            no_frames += 1

            if step >= len(lines):
                break
        pos = all_nframes[0][:, 2:5]

        scene = window.Scene()
        box_centers = np.array([[0, 0, 0.]])
        box_directions = np.array([[0, 1, 0]])
        box_colors = np.array([[1, 0, 0, 0.2]])
        box_actor = actor.box(box_centers, box_directions, box_colors,
                            size=(box_lx, box_ly, box_lz),
                            heights=2, vertices=None, faces=None)
        box_actor.GetProperty().SetRepresentationToWireframe()
        box_actor.GetProperty().SetLineWidth(10)
        atom_types = all_nframes[0][:, 1]

        #print(atom_types)
        colors = np.ones((no_atoms, 3))
        colors[atom_types == 1] = np.array([1., 0, 0])
        colors[atom_types == -1] = np.array([0, 0, 1.])

        radii = 0.1 * np.ones(no_atoms)
        #radii[atom_types == 1] = 0.5
        #radii[atom_types == -1] = 0.66
        scene = window.Scene()

        sphere_actor = actor.sphere(centers=pos,
                                    colors=colors,
                                    radii=radii)
        #, theta=6, phi=6
        scene.add(sphere_actor)
        scene.add(box_actor)

        showm = window.ShowManager(scene,
                                size=(900, 768), reset_camera=False,
                                order_transparent=True)

        showm.initialize()

        tb = ui.TextBlock2D(bold=True)

        # use itertools to avoid global variables
        counter = itertools.count()


        def get_vertices(act):

            all_vertices = np.array(numpy_support.vtk_to_numpy(
                act.GetMapper().GetInput().GetPoints().GetData()))
            return all_vertices


        def set_vertices(act, num_arr):

            vtk_num_array = numpy_support.numpy_to_vtk(num_arr)
            act.GetMapper().GetInput().GetPoints().SetData(vtk_num_array)

        def modified(act):
            act.GetMapper().GetInput().GetPoints().GetData().Modified()
            act.GetMapper().GetInput().ComputeBounds()


        global all_vertices
        all_vertices = get_vertices(sphere_actor)
        initial_vertices = all_vertices.copy()
        no_vertices_per_sphere = len(all_vertices) / no_atoms

        # set_trace()
        pos = all_nframes[0][:, 2:5]


        def timer_callback(_obj, _event):
            cnt = next(counter)
            #cnt = next(counter)

            tb.message = "Let's count up to 100 and exit :" + str(cnt)

            pos = all_nframes[cnt][:, 2:5]
            all_vertices[:] = initial_vertices + \
                np.repeat(pos, no_vertices_per_sphere, axis=0)
            set_vertices(sphere_actor, all_vertices)
            modified(sphere_actor)

            showm.render()
            if cnt == no_frames - 1:
                showm.exit()


        scene.add(tb)

        # Run every 200 milliseconds
        showm.add_timer_callback(True, 500, timer_callback)

        showm.start()
    #make movie
    def make_movie(self, num=0, file_name = ''):
        total_positive_ion = 0.5 * (len(self.ljatom))
        if num==0:
            outdump = open("many_particle_data/"+file_name, "w")
        else:
            outdump = open("many_particle_data/"+file_name, "a")
        outdump.write("ITEM: TIMESTEP\n")
        outdump.write("{}\n".format(num - 1))
        outdump.write("ITEM: NUMBER OF ATOMS\n")
        outdump.write("{}\n".format(len(self.ljatom)))
        outdump.write("ITEM: BOX BOUNDS\n")
        outdump.write("{}\t{}\n".format(-0.5*self.lx, 0.5*self.lx))
        outdump.write("{}\t{}\n".format(-0.5*self.ly, 0.5*self.ly))
        outdump.write("{}\t{}\n".format(-0.5*self.lz, 0.5*self.lz))
        #outdump.write("ITEM: ATOMS index type x y z v\n")
        outdump.write("ITEM: ATOMS index type x y z \n")
        for i in range(len(self.ljatom)):
            if (i <= int(total_positive_ion)):
                outdump.write("{}\t{}\t{}\t{}\t{}\n".format(i+1, "1", self.ljatom[i].position.x, self.ljatom[i].position.y, self.ljatom[i].position.z))
            else:
                outdump.write("{}\t{}\t{}\t{}\t{}\n".format(i+1, "-1", self.ljatom[i].position.x, self.ljatom[i].position.y, self.ljatom[i].position.z))
        outdump.close()
       # print(frame[:])
    #Velocity init: Currenly zero
    def initialize_particle_velocities(self):
        print("Velocity init")

    #Force computation
    def update_forces(self, dcut=0.0):
        self.dcut = dcut
        elj = 1.0
        dcut2 = self.dcut*self.dcut
        d = 1 #recall that we are working in reduced units where the unit of length is the diameter of the particle
        d2 = d*d
        d6 = d2 * d2 * d2
        d12 = d6 * d6
        d_cut_off = dcut2 * d2

        for i in range(len(self.ljatom)):
            flj = Vector3D(0.0,0.0,0.0)
            for j in range(len(self.ljatom)):
                if (j == i):
                    continue
                r_vec = self.ljatom[i].position - self.ljatom[j].position
                #Add periodic boundries
                r_vec.apply_periodic_boundries(self.lx, self.ly, self.lz)
                r2 = r_vec.sqd_magnitude()
                #Cut off check
                if (r2 < d_cut_off):
                    r6 = r2 * r2 * r2
                    r12 = r6 * r6
                    flj += ( r_vec * ( 48 * elj * (  (d12 / r12)  - 0.5 *  (d6 / r6) ) * ( 1 / r2 ) ) )
            self.ljatom[i].force = flj

    #Potential Energy computation
    def update_p_energies(self, dcut=0.0):
        self.dcut = dcut
        elj = 1.0
        dcut2 = self.dcut  * self.dcut
        dcut6 = dcut2 * dcut2 * dcut2
        dcut12 = dcut6 * dcut6
        d = 1 #recall that we are working in reduced units where the unit of length is the diameter of the particle
        d2 = d * d
        d6 = d2 * d2 * d2
        d_cut_off = dcut2 * d2

        energy_shift = 4*elj*((1/dcut12) - (1/dcut6))
        total_lj_atom_atom = 0.0

        for i in range(len(self.ljatom)):
            uljpair = 0.0
            for j in range(len(self.ljatom)):
                if (j == i):
                    continue
                r_vec = self.ljatom[i].position - self.ljatom[j].position
                #Add periodic boundries
                r_vec.apply_periodic_boundries( self.lx, self.ly, self.lz)
                r2 = r_vec.sqd_magnitude()
                #Cut off check
                if (r2 < d_cut_off):
                    r6 = r2 * r2 * r2
                    uljpair +=  4 * elj * (d6 / r6) * ( ( d6 / r6 ) - 1 ) - energy_shift

            self.ljatom[i].pe = uljpair
            total_lj_atom_atom += uljpair
        total_lj_atom_atom = 0.5*total_lj_atom_atom
        return total_lj_atom_atom


    #update half time velocity
    def update_velocity(self, delta_t):
        for i in range(len(self.ljatom)):
            self.ljatom[i].update_velocity(delta_t)

    #update full time velocity
    def update_position(self, delta_t):
        for i in range(len(self.ljatom)):
            self.ljatom[i].update_position(delta_t)

"""# Other functions"""

#display progress
def progressBar(fraction_completed):
    val = int(fraction_completed * 100)
    lpad = int(fraction_completed * 60)
    sys.stdout.write("\r[%-60s] %d%%" % ('='*lpad, val))

"""# Velocity verlet Function"""

def velocity_verlet(ljatom_density=None, number_ljatom=None):
    bx=0.0                      #box edge lengths
    by=0.0
    bz=0.0
    dcut=2.5
    ljatom_diameter = 1.0         #in reduced units
    ljatom_mass = 1.0

    edge_length = math.pow(number_ljatom/ljatom_density,1.0/3.0)

    bx = edge_length
    by = edge_length
    bz = edge_length
    simulation_box = Box(Vector3D(0.0, 0.0, 0.0), bx, by, bz)
    ljatom_size = simulation_box.put_ljatoms_inside(number_ljatom, ljatom_diameter, ljatom_mass, ljatom_density)
    #Inital force computation
    simulation_box.update_forces(dcut)

    totaltime = 5 #20
    steps = 500# 20000
    delta_t = totaltime/steps
    movie_file_name = "propagation.lammpstrj"

    #Initial movie dump
    simulation_box.make_movie(num=0, file_name="propagation.lammpstrj")

    movie_frequency = 100
    average_pe = 0.0
    average_ke = 0.0
    data_collect_frequency = 1000
    samples = 0
    hit_eqm = 3000               #this is your choice of where you think the system hit equilibrium

    #Molecular Dynamics
    print("progress...")


    for num in range(1, steps+1):
        #velocity-Verlet
        #update velocity half timestep
        simulation_box.update_velocity(delta_t)
        #update position full timestep
        simulation_box.update_position(delta_t)
        #expensive step
        simulation_box.update_forces(dcut)
        totalpe = simulation_box.update_p_energies(dcut)
        #update other half velocity
        simulation_box.update_velocity(delta_t)
        #calling a movie function to get a movie of the simulation
       # if (num%movie_frequency == 0):
        simulation_box.make_movie(num=num, file_name=movie_file_name)
    simulation_box.animation()



"""# Run the code"""

# Run the program
def input_parameters(number_ljatom=None):

    ljatom_diameter=0.0
    ljatom_density=0.0
    ljatom_mass=0.0
    bx=0.0                      #box edge lengths
    by=0.0
    bz=0.0
    ljatom_diameter = 1.0         #in reduced units
    ljatom_mass = 1.0
    b = input("lx = ly ? ")
    edge_length_x = float(b)
    c = input("lz? ")
    edge_length_z = float(c)
    edge_length_y = edge_length_x
    bx = edge_length_x
    by = edge_length_y
    bz = edge_length_z
    simulation_box = Box(Vector3D(0.0, 0.0, 0.0), bx, by, bz)
    ljatom_size = simulation_box.put_ljatoms_inside(number_ljatom, ljatom_diameter, ljatom_mass, ljatom_density)

a = input("How many ions are available? ")
total_salt_ion = int(a)

input_parameters(number_ljatom=total_salt_ion)
velocity_verlet(ljatom_density=0.8442, number_ljatom=total_salt_ion)
