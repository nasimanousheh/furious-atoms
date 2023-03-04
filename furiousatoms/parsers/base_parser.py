from abc import abstractmethod
from typing import Tuple
import numpy as np

from furiousatoms.molecular import MolecularStructure

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

        # print(self.errors) #TODO display popup

        for i in range(len(self.box_size)):
            self.box_size[i] = abs(self.box_size[i])

        return MolecularStructure(self.box_size, np.array(self.positions), np.array(self.bonds), np.array(self.atom_types))
    
    @abstractmethod
    def parseLine(self, fname) -> Tuple[list, np.array, np.array, np.array]:
        pass