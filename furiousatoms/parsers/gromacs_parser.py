import numpy as np
from util import float_or_zero

class GROMACSParser:
    def __init__(self) -> None:
        self.box_size = [0, 0, 0]
        self.positions = []
        self.bonds = []
        self.atom_types = []
        self.errors = "" #error messages to help user debug their file
        self.lineId = 0 #index of line being processed

    def parse(self, fname):
        raise NotImplementedError()

#TODO remove
if __name__ == "__main__":
    parser = GROMACSParser()
    parser.parse("C:\\Users\\Pete\\Desktop\\furious-atoms-exercises\\one.xyz")