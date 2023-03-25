from furiousatoms.molecular import MolecularStructure
from furiousatoms.savers.base_saver import BaseSaver
from furiousatoms.parsers.gromacs_parser import GROMACSParser #TODO remove

DEFAULT_HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"
DEFAULT_RESIDUE_NUMBER = 2
DEFAULT_RESIDUE_NAME = "MOL"

class GROMACSSaver(BaseSaver):
    def save_to_file(self, new_fname, old_fname, structure: MolecularStructure):
        '''
        Write all data to a file using the format specified in 
        https://manual.gromacs.org/archive/5.0.3/online/gro.html 

        Assuming the file referred to by old_fname was the source of the structure data,
        copy over any fields that were not read by our parser.
        '''
        self.atom_serial = 1
        with open(new_fname, 'w') as fp:
            with open(old_fname, 'r') as old_fp:
                header = old_fp.readline()
                fp.write(header)
                old_fp.readline()
                fp.write("%d\n"%structure.pos.shape[0]) #atom count

                for i in range(len(structure.pos)):
                    old_line = old_fp.readline()
                    if self.is_atom_deleted(i):
                        continue
                    x, y, z = structure.pos[i]
                    if len(old_line) >= 10 and len(old_line.split()) > 3: #is it an atom line?
                        residue_number = old_line[:5]
                        residue_name = old_line[5:10]
                    else:
                        residue_number = DEFAULT_RESIDUE_NUMBER
                        residue_name = DEFAULT_RESIDUE_NAME
                    typ = structure.atom_types[i]
                    fp.write("%5s%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(
                        residue_number,
                        residue_name,
                        typ,
                        self.atom_serial,
                        x, y, z
                    ))
                    self.atom_serial += 1
                
                fp.write("%4.5f %4.5f %4.5f\n"%(structure.box_size[0], structure.box_size[1], structure.box_size[2]))

#TODO remove
if __name__ == "__main__":
    INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\furiousatoms\\tests\\test_data\\CB_18\\CB_18.gro"
    OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\saving_test.gro"
    structure = GROMACSParser().parse(INPUT_FNAME)

    import numpy as np
    structure.atom_types = np.hstack(( structure.atom_types, 'H'))
    structure.pos = np.vstack(( structure.pos, np.array([9.0,8.0,7.0]) ))
    structure.atom_types = np.hstack(( structure.atom_types, 'Ga'))
    structure.pos = np.vstack(( structure.pos, np.array([3.0,2.0,1.0]) ))

    from timeit import default_timer as timer

    start = timer()

    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    deleted_particles[5] = True
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = GROMACSSaver(deleted_particles, deleted_bonds)
    saver.save_to_file(OUTPUT_FNAME, INPUT_FNAME, structure)

    end = timer()
    print(end - start)

    with open(OUTPUT_FNAME, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            print(line, end='')
