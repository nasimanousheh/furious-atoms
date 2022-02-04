# Standard package
from json import load
import os
import fnmatch

# Local package
from furiousatoms import io
from fury import disable_warnings

disable_warnings()

# 3rd Party package
import numpy as np
from fury import window, actor, utils, pick, ui, primitive, material
from PySide2 import QtCore
from PySide2 import QtGui
from PySide2.QtGui import QActionEvent, QIcon
from PySide2 import QtWidgets
from PySide2.QtCore import QTimer
import MDAnalysis


from furiousatoms.fullerenes_builder import load_CC1_file
from fury.utils import (get_actor_from_primitive, normals_from_actor,
                        tangents_to_actor, update_polydata_normals,
                        tangents_from_direction_of_anisotropy)

def save_PDB2LMP(SM, file_name):

    masses = MDAnalysis.topology.guessers.guess_masses(SM.atom_type)
    SM.universe.add_TopologyAttr('masses', masses)
    angles = MDAnalysis.topology.guessers.guess_angles(SM.universe.bonds)
    SM.universe.add_angles(angles)
    dihedrals = MDAnalysis.topology.guessers.guess_dihedrals(SM.universe.angles)
    SM.universe.add_dihedrals(dihedrals)
    improper_dihedrals = MDAnalysis.topology.guessers.guess_improper_dihedrals(SM.universe.angles)
    SM.universe.add_impropers(improper_dihedrals)
    mass_unique_types = np.unique(SM.universe.atoms.masses)
    bond_unique_types = np.unique(SM.universe.atoms.bonds.types())
    angle_unique_types = np.unique(SM.universe.atoms.angles.types())
    dihedral_unique_types = np.unique(SM.universe.atoms.dihedrals.types())
    improper_unique_types = np.unique(SM.universe.atoms.impropers.types())

    outdump = open(file_name, "w")
    outdump.write("LAMMPS data file\n\n")
    outdump.write("{}\t{} \n".format(SM.no_atoms, ' atoms'))
    outdump.write("{}\t{} \n".format(SM.no_bonds, ' bonds'))
    outdump.write("{}\t{} \n".format(len(angles), ' angles'))
    outdump.write("{}\t{} \n".format(len(dihedrals), ' dihedrals'))
    outdump.write("{}\t{} \n".format(len(improper_dihedrals), ' impropers'))
    outdump.write("{}\t{} \n".format(len(SM.unique_types), 'atom types'))
    if SM.no_bonds > 0 :
        outdump.write("{}\t{} \n".format(len(np.unique(SM.universe.bonds.types)), 'bond types'))
    if len(angles) > 0:
        outdump.write("{}\t{} \n".format(len(np.unique(SM.universe.angles.types)), 'angle types'))
    if len(dihedrals) > 0:
        outdump.write("{}\t{} \n\n".format(len(np.unique(SM.universe.dihedrals.types)), 'dihedral types'))
    if len(improper_dihedrals) > 0:
        outdump.write("{}\t{} \n\n".format(len(np.unique(SM.universe.impropers.types)), 'improper types'))
    outdump.write("{}\t{}\t{} \n".format(-0.5* SM.box_lx, 0.5* SM.box_lx, ' xlo xhi'))
    outdump.write("{}\t{}\t{} \n".format(-0.5* SM.box_ly, 0.5* SM.box_ly, ' ylo yhi'))
    outdump.write("{}\t{}\t{} \n\n".format((-0.5* SM.box_lz), (0.5* SM.box_lz), ' zlo zhi'))
    outdump.write("Masses\n\n")

    for i in range (len(SM.unique_types)):
        outdump.write("{}\t{} \n".format(SM.unique_types[i], mass_unique_types[i]))

    outdump.write("\n")
    outdump.write("Atoms          # full\n\n")
    num_molecules = 0
    for i in range(SM.no_atoms):
        outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 1), num_molecules, SM.universe.atoms.types[i], '0', SM.pos[i][0], SM.pos[i][1], SM.pos[i][2], '0   0   0 '))
    outdump.write("\n")
    if SM.no_bonds > 0:
        outdump.write("Bonds\n\n")
    #    num_molecules = 0
        for g in range(SM.no_bonds):
            outdump.write("{}\t{}\t{}\t{}\n".format((g + 1), '1', SM.universe.bonds.indices[g][0]+1, SM.universe.bonds.indices[g][1]+1))
        #    num_molecules = num_molecules + 1
        if len(angles) > 0:
            outdump.write("\n")
            outdump.write("Angles\n\n")
            num_molecules = 0
            for n in range(len(angles)):
                outdump.write("{}\t{}\t{}\t{}\t{}\n".format( n + 1, '1', SM.universe.angles.indices[n][0]+1, SM.universe.angles.indices[n][1]+1, SM.universe.angles.indices[n][2]+1))
            # num_molecules = num_molecules + 1

            if len(dihedrals) > 0:
                outdump.write("\n")
                outdump.write("Dihedrals\n\n")
                num_molecules = 0
                for s in range(len(dihedrals)):
                    outdump.write("{}\t{}\t{}\t{}\t{}\t{}\n".format( s + 1, '1', SM.universe.dihedrals.indices[s][0]+1, SM.universe.dihedrals.indices[s][1]+1, SM.universe.dihedrals.indices[s][2]+1,SM.universe.dihedrals.indices[s][3]+1))
                #    num_molecules = num_molecules + 1

                if len(improper_dihedrals) > 0:
                    outdump.write("\n")
                    outdump.write("Impropers\n\n")
                    num_molecules = 0
                    for t in range(len(improper_dihedrals)):
                        outdump.write("{}\t{}\t{}\t{}\t{}\t{}\n".format( t + 1, '1', SM.universe.impropers.indices[t][0], SM.universe.impropers.indices[t][1], SM.universe.impropers.indices[t][2],SM.universe.impropers.indices[t][3]))
                    #   num_molecules = num_molecules + 1