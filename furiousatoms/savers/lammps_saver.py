import numpy as np

HEADER = "LAMMPS data file. Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class LAMMPSSaver():
    def write_all_lines(self, fname, box_size: list, positions: np.array, bonds: np.array, atom_types: np.array):
        with open(fname, 'w') as fp:
            fp.write(HEADER)
            fp.write("%d atoms\n"%positions.shape[0]) #atom count
            fp.write("%d bonds\n"%bonds.shape[0]) #bond count

            #Unimplemented
            fp.write("0 angles\n")
            fp.write("0 dihedrals\n")
            fp.write("0 impropers\n")

            #Number of uniques
            fp.write("%d atom types\n"%(np.unique(atom_types).size))
            fp.write("1 bond types\n") #TODO Implement bond types

            #Unimplemented
            fp.write("0 angle types\n")
            fp.write("0 dihedral types\n")
            fp.write("0 improper types\n")
            fp.write("\n")

            #Box size: make lo and hi dimensions the same
            fp.write("%f %f xlo xhi\n"%(box_size[0], box_size[0]))
            fp.write("%f %f ylo yhi\n"%(box_size[1], box_size[1]))
            fp.write("%f %f zlo zhi\n"%(box_size[2], box_size[2]))
            fp.write("\n")

            #Masses
            for element in np.unique(atom_types):
                defaultMass = 0 #TODO write mass to LAMMPS file
                fp.write("%2s %f\n"%(element, defaultMass))
            fp.write("\n")

            #Atoms: always write in 'atomic' style
            #Syntax: `atom-ID atom-type x y z`
            fp.write("Atoms          #atomic\n")
            for i in range(len(positions)):
                x, y, z = positions[i]
                fp.write("%d %2s %8.6f %8.6f %8.6f\n"%(i, atom_types[i], x, y, z))
            fp.write("\n")

            #Bond syntax: `bondID bondType atom1 atom2`
            fp.write("Bonds\n")
            i = 1
            for bond in bonds:
                defaultBondType = 1 #TODO Implement bond types
                atom1, atom2 = bond + 1
                fp.write("%d %d %d %d\n"%(i, defaultBondType, atom1, atom2))
                i += 1
            fp.write("\n")

            #Omit angles and dihedrals sections.