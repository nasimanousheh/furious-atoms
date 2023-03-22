import numpy as np
from furiousatoms.molecular import MolecularStructure
from furiousatoms.savers.base_saver import BaseSaver
from furiousatoms.parsers.parser_util import float_or_zero
from furiousatoms.parsers.lammps_parser import LAMMPSParser #TODO remove

HEADER = "LAMMPS data file. Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

#TODO replace keep_line with simple bool flag

class LAMMPSSaver(BaseSaver):
    def keep_line(self):
        pass

    def write_all_atoms(self, new_fp, old_fp, structure):
        '''
        Always write in 'atomic' style.
        Syntax: `atom-ID atom-type x y z`
        '''
        i = 0
        while True:
            line = old_fp.readline()
            words = line.split()
            if len(words) < 5: #num. fields in 'atomic' style
                break
            atom_id = int(words[0]) - 1
            if not self.is_atom_deleted(atom_id):            
                x, y, z = structure.pos[i]
                typ = structure.atom_types[i]
                new_fp.write("%d %2s %8.6f %8.6f %8.6f\n"%(self.atom_serial, typ, x, y, z))
            i += 1

            self.atom_serial += 1
        new_fp.write("\n")


        # words = line.split()
        # if len(words < 2):
        #     return
        # atom_id = words[0]
        # if self.is_atom_deleted(atom_id):
        #     return
        
        # #Always write in 'atomic' style
        # #Syntax: `atom-ID atom-type x y z`
        # for i in range(len(structure.pos)):
        #     x, y, z = structure.pos[i]
        #     typ = structure.atom_types[i]
        #     fp.write("%d %2s %8.6f %8.6f %8.6f\n"%(self.atom_serial, typ, x, y, z))
        # new_fp.write("\n")

        # self.atom_serial += 1


    def write_all_bonds(self, new_fp, old_fp, structure):
        i = 0
        while True:
            line = old_fp.readline()
            words = line.split()
            if len(words) < 4 or i > len(structure.bonds):
                break
            
            bond = structure.bonds[i]
            if not self.is_bond_deleted(i):
                atom1, atom2 = bond + 1
                bond_type = words[1]
                new_fp.write("%d %s %d %d\n"%(i + 1, bond_type, atom1, atom2))
            i += 1
        new_fp.write("\n")
        

    def save_to_file(self, new_fname, old_fname, structure: MolecularStructure):
        '''
        Write all data to a file using the format specified in 
        https://www.smcm.iqfr.csic.es/docs/lammps/read_data.html

        Assuming the file referred to by old_fname was the source of the structure data,
        copy over any fields that were not read by our parser.
        '''
        self.saved_atom_ids = np.zeros(shape=(len(structure.pos)))
        self.atom_serial = 1
        saver_method = self.keep_line
        with open(new_fname, 'w') as new_fp:
            with open(old_fname, 'r') as old_fp:
                new_fp.write(old_fp.readline()) #file comment
                
                while True:
                    old_line = old_fp.readline()
                    if not old_line:
                        break
                    line = old_line.strip()
                    line_lower = line.lower()

                    #LAMMPS BODY ##############################
                    if line_lower.startswith("masses"): 
                        saver_method = self.keep_line
                    elif line_lower.startswith("atoms"): 
                        new_fp.write("Atoms          #atomic\n")
                        self.write_all_atoms(new_fp, old_fp, structure)
                        saver_method = None
                    elif line_lower.startswith("bonds"): 
                        new_fp.write("Bonds\n")
                        self.write_all_bonds(new_fp, old_fp, structure)
                        saver_method = None
                    #LAMMPS HEADER ##############################
                    elif line_lower.endswith("atoms"):
                        num_atoms = len(np.argwhere(self.deleted_particles == False))
                        new_fp.write("%d\tatoms\n"%(num_atoms))
                        saver_method = None
                    elif line_lower.endswith("bonds"):
                        num_bonds = len(np.argwhere(self.deleted_bonds == False))
                        new_fp.write("%d\tbonds\n"%(num_bonds))
                        saver_method = None
                    elif line_lower.endswith("atom types"):
                        num_types = len(np.unique(structure.atom_types))
                        new_fp.write("%d\tatom types\n"%(num_types))
                        saver_method = None
                    elif line_lower.endswith("hi") and "lo" in line_lower:
                        words = line.split()
                        if len(words) != 4:
                            new_fp.write(line) #handles invalid box size lines
                            continue
                        lo = float_or_zero(words[0])
                        if line_lower.endswith("xlo xhi"):
                            hi = structure.box_size[0]
                        elif line_lower.endswith("ylo yhi"):
                            hi = structure.box_size[1]
                        elif line_lower.endswith("zlo zhi"):
                            hi = structure.box_size[2]
                        else:
                            hi = 0.0
                        new_fp.write("%f\t%f\t%s %s\n"%(lo, hi, words[2], words[3]))
                        saver_method = None
                    else:
                        saver_method = self.keep_line

                    
                    if saver_method == self.keep_line:
                        new_fp.write(line)
                        if "\n" not in line:
                            new_fp.write("\n")


# class LAMMPSSaver():
#     def writeAllLines(self, fname, boxSize: list, positions: np.array, bonds: np.array, atomTypes: np.array):
#         with open(fname, 'w') as fp:
#             fp.write(HEADER)
#             fp.write("%d atoms\n"%positions.shape[0]) #atom count
#             fp.write("%d bonds\n"%bonds.shape[0]) #bond count

#             #Unimplemented
#             fp.write("0 angles\n")
#             fp.write("0 dihedrals\n")
#             fp.write("0 impropers\n")

#             #Number of uniques
#             fp.write("%d atom types\n"%(np.unique(atomTypes).size))
#             fp.write("1 bond types\n") #TODO Implement bond types

#             #Unimplemented
#             fp.write("0 angle types\n")
#             fp.write("0 dihedral types\n")
#             fp.write("0 improper types\n")
#             fp.write("\n")

#             #Box size: make lo and hi dimensions the same
#             fp.write("%f %f xlo xhi\n"%(boxSize[0], boxSize[0]))
#             fp.write("%f %f ylo yhi\n"%(boxSize[1], boxSize[1]))
#             fp.write("%f %f zlo zhi\n"%(boxSize[2], boxSize[2]))
#             fp.write("\n")

#             #Masses
#             for element in np.unique(atomTypes):
#                 defaultMass = 0 #TODO write mass to LAMMPS file
#                 fp.write("%2s %f\n"%(element, defaultMass))
#             fp.write("\n")

#             #Atoms: always write in 'atomic' style
#             #Syntax: `atom-ID atom-type x y z`
#             fp.write("Atoms          #atomic\n")
#             for i in range(len(positions)):
#                 x, y, z = positions[i]
#                 fp.write("%d %2s %8.6f %8.6f %8.6f\n"%(i, atomTypes[i], x, y, z))
#             fp.write("\n")

#             #Bond syntax: `bondID bondType atom1 atom2`
#             fp.write("Bonds\n")
#             i = 1
#             for bond in bonds:
#                 defaultBondType = 1 #TODO Implement bond types
#                 atom1, atom2 = bond + 1
#                 fp.write("%d %d %d %d\n"%(i, defaultBondType, atom1, atom2))
#                 i += 1
#             fp.write("\n")

#             #Omit angles and dihedrals sections.

            


#TODO remove
if __name__ == "__main__":
    INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\furiousatoms\\tests\\test_data\\CB_18\\CB_18.data"
    OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\saving_test.data"
    structure = LAMMPSParser().parse(INPUT_FNAME)
    
    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    # deleted_particles[5] = True
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = LAMMPSSaver(deleted_particles, deleted_bonds)
    saver.save_to_file(OUTPUT_FNAME, INPUT_FNAME, structure)

    with open(OUTPUT_FNAME, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            print(line, end='')