from fury import disable_warnings
disable_warnings()
import numpy as np
import MDAnalysis


def save_PDB2LMP(SM, file_name, universe):
    unique_types = np.unique(universe.atoms.types)
    no_bonds = len(universe.bonds)
    no_atoms = len(universe.atoms)

    masses = MDAnalysis.topology.guessers.guess_masses(universe.atoms.types)
    universe.add_TopologyAttr('masses', masses)
    angles = MDAnalysis.topology.guessers.guess_angles(universe.bonds)
    universe.add_angles(angles)
    dihedrals = MDAnalysis.topology.guessers.guess_dihedrals(universe.angles)
    universe.add_dihedrals(dihedrals)
    improper_dihedrals = MDAnalysis.topology.guessers.guess_improper_dihedrals(universe.angles)
    universe.add_impropers(improper_dihedrals)
    mass_unique_types = np.unique(universe.atoms.masses)
    bond_unique_types = np.unique(universe.atoms.bonds.types())
    angle_unique_types = np.unique(universe.atoms.angles.types())
    dihedral_unique_types = np.unique(universe.atoms.dihedrals.types())
    improper_unique_types = np.unique(universe.atoms.impropers.types())

    outdump = open(file_name, "w")
    outdump.write("LAMMPS data file\n\n")
    outdump.write("{}\t{} \n".format(no_atoms, ' atoms'))
    outdump.write("{}\t{} \n".format(no_bonds, ' bonds'))
    outdump.write("{}\t{} \n".format(len(angles), ' angles'))
    outdump.write("{}\t{} \n".format(len(dihedrals), ' dihedrals'))
    outdump.write("{}\t{} \n".format(len(improper_dihedrals), ' impropers'))
    outdump.write("{}\t{} \n".format(len(unique_types), 'atom types'))
    if no_bonds > 0:
        outdump.write("{}\t{} \n".format(len(np.unique(universe.bonds.types)), 'bond types'))
    if len(angles) > 0:
        outdump.write("{}\t{} \n".format(len(np.unique(universe.angles.types)), 'angle types'))
    if len(dihedrals) > 0:
        outdump.write("{}\t{} \n\n".format(len(np.unique(universe.dihedrals.types)), 'dihedral types'))
    if len(improper_dihedrals) > 0:
        outdump.write("{}\t{} \n\n".format(len(np.unique(universe.impropers.types)), 'improper types'))
    outdump.write("{}\t{}\t{} \n".format(-0.5* SM.box_lx, 0.5* SM.box_lx, ' xlo xhi'))
    outdump.write("{}\t{}\t{} \n".format(-0.5* SM.box_ly, 0.5* SM.box_ly, ' ylo yhi'))
    outdump.write("{}\t{}\t{} \n\n".format((-0.5* SM.box_lz), (0.5* SM.box_lz), ' zlo zhi'))
    # if len(unique_types) == len(mass_unique_types):
    outdump.write("Masses\n\n")
    for i in range(len(unique_types)):
        outdump.write("{}\t{} \n".format(unique_types[i], MDAnalysis.topology.guessers.guess_masses(unique_types[i])[0]))

    outdump.write("\n")
    outdump.write("Atoms          # full\n\n")
    num_molecules = 0
    for i in range(no_atoms):
        outdump.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format((i + 1), num_molecules, universe.atoms.types[i], '0', universe.atoms.positions[i][0], universe.atoms.positions[i][1], universe.atoms.positions[i][2], '0   0   0 '))
    outdump.write("\n")
    if no_bonds > 0:
        outdump.write("Bonds\n\n")
    #    num_molecules = 0
        for g in range(no_bonds):
            outdump.write("{}\t{}\t{}\t{}\n".format((g + 1), '1', universe.bonds.indices[g][0]+1, universe.bonds.indices[g][1]+1))
        #    num_molecules = num_molecules + 1
        if len(angles) > 0:
            outdump.write("\n")
            outdump.write("Angles\n\n")
            num_molecules = 0
            for n in range(len(angles)):
                outdump.write("{}\t{}\t{}\t{}\t{}\n".format( n + 1, '1', universe.angles.indices[n][0]+1, universe.angles.indices[n][1]+1, universe.angles.indices[n][2]+1))
            # num_molecules = num_molecules + 1

            if len(dihedrals) > 0:
                outdump.write("\n")
                outdump.write("Dihedrals\n\n")
                num_molecules = 0
                for s in range(len(dihedrals)):
                    outdump.write("{}\t{}\t{}\t{}\t{}\t{}\n".format( s + 1, '1', universe.dihedrals.indices[s][0]+1, universe.dihedrals.indices[s][1]+1, universe.dihedrals.indices[s][2]+1,universe.dihedrals.indices[s][3]+1))
                #    num_molecules = num_molecules + 1

                if len(improper_dihedrals) > 0:
                    outdump.write("\n")
                    outdump.write("Impropers\n\n")
                    num_molecules = 0
                    for t in range(len(improper_dihedrals)):
                        outdump.write("{}\t{}\t{}\t{}\t{}\t{}\n".format( t + 1, '1', universe.impropers.indices[t][0], universe.impropers.indices[t][1], universe.impropers.indices[t][2],universe.impropers.indices[t][3]))
                    #   num_molecules = num_molecules + 1