from abc import abstractmethod
from furiousatoms.molecular import MolecularStructure

class BaseSaver:
    '''
    An abstract class which other file format-specific savers can extend.
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
    def save_to_file(self, new_fname, old_fname, structure: MolecularStructure) -> None:
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