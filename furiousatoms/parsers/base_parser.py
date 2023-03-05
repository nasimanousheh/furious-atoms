from abc import abstractmethod
from typing import Tuple
import numpy as np
import warnings

from furiousatoms.molecular import MolecularStructure
from furiousatoms.parsers.parser_util import clean_bonds

'''
An abstract class which other file format-specific parsers can extend.
'''
class BaseParser:
    def __init__(self) -> None:
        self.box_size = [0, 0, 0]
        self.positions = []
        self.bonds = []
        self.atom_types = []
        self.errors = "" #error messages to help user debug their file
        self.lineId = 0 #index of line being processed
        self.parserMethod = None #optional; some parsers use this to switch between how they parse a line

    def parse(self, fname) -> MolecularStructure:
        '''
        Given a file name, extract as much information as possible then return
        a tuple of the following:
            - a 1x3 list indicating the box size
            - an array of 1x3 position arrays
            - an array of 1x2 bonds (each value is an atom index in the position array)
            - an string array of atom types as element symbols. Length should match position array.
        '''
        with open(fname, "r") as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                self.parseLine(line)
                self.lineId += 1

        # warnings.warn(self.errors) #TODO display popup

        for i in range(len(self.box_size)):
            self.box_size[i] = abs(self.box_size[i])
        if len(self.positions) == 0:
            warnings.warn("MolecularStructure is empty. No positions could be parsed.")
            empty = MolecularStructure.create_empty()
            empty.box_size = self.box_size
            return empty

        self.positions = np.array(self.positions)
        self.bonds = np.array(self.bonds)
        self.bonds = clean_bonds(self.bonds, len(self.positions))
        self.atom_types = np.array(self.atom_types)

        return MolecularStructure(self.box_size, self.positions, self.bonds, self.atom_types)
    
    @abstractmethod
    def parseLine(self, fname) -> Tuple[list, np.array, np.array, np.array]:
        pass