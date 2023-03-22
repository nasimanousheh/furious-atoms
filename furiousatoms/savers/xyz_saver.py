from furiousatoms.molecular import MolecularStructure
from furiousatoms.savers.base_saver import BaseSaver
from furiousatoms.parsers.xyz_parser import XYZParser #TODO remove

HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class XYZSaver(BaseSaver):
    def save_to_file(self, new_fname, old_fname, structure: MolecularStructure):
        '''
        Write all data to a file using the format specified in 
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html

        Param `old_fname` is ignored.
        '''
        with open(new_fname, 'w') as fp:
            fp.write("%d\n"%structure.pos.shape[0])
            fp.write(HEADER)

            for i in range(len(structure.pos)):
                if not self.is_atom_deleted(i):
                    x, y, z = structure.pos[i]
                    fp.write("%2s %8.6f %8.6f %8.6f\n"%(structure.atom_types[i], x, y, z))


#TODO remove
if __name__ == "__main__":
    INPUT_FNAME = "C:\\Users\\Pete\\Desktop\\furious-atoms\\furiousatoms\\tests\\test_data\\CB_18\\CB_18.xyz"
    OUTPUT_FNAME = "C:\\Users\\Pete\\Desktop\\saving_test.xyz"
    structure = XYZParser().parse(INPUT_FNAME)

    import numpy as np
    structure.atom_types = np.hstack(( structure.atom_types, 'H'))
    structure.pos = np.vstack(( structure.pos, np.array([9.0,8.0,7.0]) ))
    structure.atom_types = np.hstack(( structure.atom_types, 'Ga'))
    structure.pos = np.vstack(( structure.pos, np.array([3.0,2.0,1.0]) ))

    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    deleted_particles[5] = True
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = XYZSaver(deleted_particles, deleted_bonds)
    saver.save_to_file(OUTPUT_FNAME, INPUT_FNAME, structure)

    with open(OUTPUT_FNAME, 'r') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            print(line, end='')