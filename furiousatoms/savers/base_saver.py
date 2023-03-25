from abc import abstractmethod
from furiousatoms.molecular import MolecularStructure
import numpy as np 
import tempfile
import shutil
import os

class BaseSaver:
    '''
    An abstract class which other file format-specific savers can extend.
    Note: all savers must use `self.corrected_bonds` to account for bond/particle deletion
    rather than `structure.bonds`.
    '''
    def __init__(self, deleted_particles, deleted_bonds):
        '''
        Atoms or bonds listed as True in `deleted_particles` or `deleted_bonds` will be ignored
        during saving.
        '''
        self.atom_serial = 1 #the ID that an atom will be saved with. Must be incremented by subclass.
        self.deleted_particles = deleted_particles
        self.deleted_bonds = deleted_bonds

    def save_to_file(self, new_fname, old_fname, structure: MolecularStructure) -> None:
        self.corrected_bonds = self.get_corrected_bonds(structure)

        if not os.path.exists(old_fname):
            try:
                print("No file exists to infer data from during saving. Using defaults instead.")
                new_fp = open(new_fname, 'w')
                output = self._save_to_file_use_defaults(new_fp, structure)
            finally:
                new_fp.close()
            return output

        if new_fname == old_fname:
            temp_fname = None
            try:
                old_fp = open(old_fname, 'r')
                temp_fp = tempfile.NamedTemporaryFile(mode='w', delete=False)
                
                output = self._save_to_file(temp_fp, old_fp, structure)
                temp_fname = temp_fp.name
            finally:
                old_fp.close()
                temp_fp.close()
                if temp_fname and os.path.exists(temp_fname):
                    shutil.copy(temp_fname, new_fname)
                    os.remove(temp_fname)
        else:
            try:
                old_fp = open(old_fname, 'r')
                new_fp = open(new_fname, 'w')
                
                output = self._save_to_file(new_fp, old_fp, structure)
            finally:
                old_fp.close()
                new_fp.close()

        return output

    @abstractmethod
    def _save_to_file(self, new_fp, old_fp, structure: MolecularStructure) -> None:
        pass

    @abstractmethod
    def _save_to_file_use_defaults(self, old_fp, structure: MolecularStructure) -> None:
        pass

    def is_atom_deleted(self, atom_id):
        if atom_id >= len(self.deleted_particles):
            return False
        if atom_id < 0:
            return True
        return self.deleted_particles[atom_id]
        
    def is_bond_deleted(self, bond_id):
        if bond_id >= len(self.deleted_bonds):
            return False
        if bond_id < 0:
            return True
        return self.deleted_bonds[bond_id]
    
    def get_corrected_bonds(self, structure: MolecularStructure):
        deleted_atom_indices = np.argwhere(self.deleted_particles == True)
        new_bonds = np.empty(shape=(len(structure.bonds), 2), dtype=int)
        for i in range(len(structure.bonds)):
            for b in range(2):
                bondValue = structure.bonds[i][b]
                k = len(deleted_atom_indices) - 1
                while k > -1:
                    if bondValue >= deleted_atom_indices[k]:
                        bondValue -= (k + 1)
                        break #exit the loop over k
                    k -= 1
                new_bonds[i][b] = bondValue
        return new_bonds