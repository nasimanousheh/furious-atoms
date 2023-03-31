from furiousatoms.molecular import MolecularStructure
from furiousatoms.savers.base_saver import BaseSaver
import numpy as np

DEFAULT_HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"
DEFAULT_RESIDUE_NUMBER = 2
DEFAULT_RESIDUE_NAME = "MOL"

ATOM_FORMAT = "%5s%-5s%5s%5d%8.3f%8.3f%8.3f\n"
BOX_SIZE_FORMAT = "%4.5f %4.5f %4.5f\n"

class GROMACSSaver(BaseSaver):
    def _save_to_file(self, new_fp, old_fp, structure: MolecularStructure):
        '''
        Write all data to a file using the format specified in 
        https://manual.gromacs.org/archive/5.0.3/online/gro.html 

        Assuming the file referred to by old_fname was the source of the structure data,
        copy over any fields that were not read by our parser.
        '''
        self.atom_serial = 1
        header = old_fp.readline()
        new_fp.write(header)

        old_fp.readline()
        num_atoms = len(np.argwhere(self.deleted_particles == False))
        new_fp.write("%d\n"%num_atoms)

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
            new_fp.write(ATOM_FORMAT%(
                residue_number,
                residue_name,
                structure.atom_types[i],
                self.atom_serial,
                x, y, z
            ))
            self.atom_serial += 1
        
        new_fp.write(BOX_SIZE_FORMAT%(structure.box_size[0], structure.box_size[1], structure.box_size[2]))


    def _save_to_file_use_defaults(self, new_fp, structure: MolecularStructure):
        self.atom_serial = 1
        new_fp.write(DEFAULT_HEADER)
        
        num_atoms = len(np.argwhere(self.deleted_particles == False))
        new_fp.write("%d\n"%num_atoms)

        for i in range(len(structure.pos)):
            if self.is_atom_deleted(i):
                continue
            x, y, z = structure.pos[i]
            new_fp.write("%5s%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(
                DEFAULT_RESIDUE_NUMBER,
                DEFAULT_RESIDUE_NAME,
                structure.atom_types[i],
                self.atom_serial,
                x, y, z
            ))
            self.atom_serial += 1
        
        new_fp.write(BOX_SIZE_FORMAT%(structure.box_size[0], structure.box_size[1], structure.box_size[2]))