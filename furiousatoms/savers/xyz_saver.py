from furiousatoms.molecular import MolecularStructure
from furiousatoms.savers.base_saver import BaseSaver

DEFAULT_HEADER = "Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"

class XYZSaver(BaseSaver):
    def _write_file_contents(self, new_fp, old_fp, structure: MolecularStructure):
        '''
        Write all data to a file using the format specified in 
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html

        Parameter `old_fp` is unused
        '''
        for i in range(len(structure.pos)):
            if not self.is_atom_deleted(i):
                x, y, z = structure.pos[i]
                new_fp.write("%2s %8.6f %8.6f %8.6f\n"%(structure.atom_types[i], x, y, z))


    def _save_to_file(self, new_fp, old_fp, structure: MolecularStructure):
        old_fp.readline()
        new_fp.write("%d\n"%structure.pos.shape[0])

        header = old_fp.readline()
        new_fp.write(header)
        
        self._write_file_contents(new_fp, old_fp, structure)


    def _save_to_file_use_defaults(self, new_fp, structure: MolecularStructure):
        '''
        Everything is the same as in _save_to_file(...), but the file's header
        is a default value.
        '''
        new_fp.write("%d\n"%structure.pos.shape[0])
        new_fp.write(DEFAULT_HEADER)
        self._write_file_contents(new_fp, None, structure)