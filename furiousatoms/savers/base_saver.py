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

    @abstractmethod
    def _save_to_file(self, new_fp, old_fp, structure: MolecularStructure) -> None:
        pass
   
    @abstractmethod
    def _save_to_file_use_defaults(self, old_fp, structure: MolecularStructure) -> None:
        pass

    def save_to_file(self, new_fpath, old_fpath, structure: MolecularStructure) -> None:
        '''
        Assuming the file referred to by old_fp was the source of the structure data,
        copy over any fields that were not read during parsing and save any changes to
        the new file. Additionally, exclude any deleted atoms or bonds.

        `new_fp` and `old_fp` are file pointers.
        '''
        self.corrected_bonds = self.get_corrected_bonds(structure)

        if old_fpath == None or not os.path.exists(old_fpath):
            print("No file exists to infer data from during saving. Using defaults instead.")
            try:
                new_fp = open(new_fpath, 'w')
                output = self._save_to_file_use_defaults(new_fp, structure)
            finally:
                new_fp.close()
            return output

        if new_fpath == old_fpath:
            temp_fpath = None
            try:
                old_fp = open(old_fpath, 'r')
                temp_fp = tempfile.NamedTemporaryFile(mode='w', delete=False)
                
                output = self._save_to_file(temp_fp, old_fp, structure)
                temp_fpath = temp_fp.name
            except:
                print("Something went wrong during saving. Trying again with defaults.")
                output = self.save_to_file_use_defaults(new_fpath, structure)
            finally:
                old_fp.close()
                temp_fp.close()
                if temp_fpath and os.path.exists(temp_fpath):
                    shutil.copy(temp_fpath, new_fpath)
                    os.remove(temp_fpath)
                    print("Saved successfully.")
        else:
            try:
                old_fp = open(old_fpath, 'r')
                new_fp = open(new_fpath, 'w')
                
                output = self._save_to_file(new_fp, old_fp, structure)
                print("Saved successfully.")
            except:
                print("Something went wrong during saving. Trying again with defaults.")
                output = self.save_to_file_use_defaults(new_fpath, structure)
            finally:
                old_fp.close()
                new_fp.close()

        return output
    
    def save_to_file_use_defaults(self, new_fpath, structure: MolecularStructure) -> None:
        try:
            new_fp = open(new_fpath, 'w')
            output = self._save_to_file_use_defaults(new_fp, structure)
            print("Saved with defaults successfully.")
        except:
            print("Failed to save with default settings.")
        finally:
            new_fp.close()

        return output

    def is_atom_deleted(self, atom_id):
        if atom_id >= len(self.deleted_particles):
            return False #probably a new atom
        if atom_id < 0:
            return True #don't save an invalid atom
        return self.deleted_particles[atom_id]
        
    def is_bond_deleted(self, bond_id):
        if bond_id >= len(self.deleted_bonds):
            return False #probably a new bond
        if bond_id < 0:
            return True #don't save an invalid bond
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