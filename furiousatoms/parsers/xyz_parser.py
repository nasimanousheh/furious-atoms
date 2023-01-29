import numpy as np
from util import float_or_zero

class XYZParser:
    def __init__(self) -> None:
        self.box_size = [0, 0, 0]
        self.positions = []
        self.bonds = []
        self.atom_types = []
        self.errors = "" #error messages to help user debug their file
        self.lineId = 0 #index of line being processed

    def parse(self, fname):
        def parseAtom(line):
            try:
                pos = np.zeros((3))
                pos = np.zeros((3))
                # print("|"+line[8:20]+"|")
                # print("|"+line[22:31]+"|")
                # print("|"+line[32:42]+"|")
                pos[0] = float_or_zero(line[8:20])
                pos[1] = float_or_zero(line[22:31])
                pos[2] = float_or_zero(line[32:42])
                self.positions.append(pos)
            except IndexError:
                self.errors += "Line #%d is too short to read an atom's position. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.lineId)
            except:
                self.errors += "Failure processing line #%d.\n"%(self.lineId)


        with open(fname, "r") as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                
                #Try to parse regardless of whether line is a comment or a position
                parseAtom(line)
                self.lineId += 1
        
        print(self.errors) #TODO display popup
    
        return self.box_size, np.array(self.positions), np.array(self.bonds), np.array(self.atom_types)

#TODO remove
if __name__ == "__main__":
    parser = XYZParser()
    parser.parse("C:\\Users\\Pete\\Desktop\\furious-atoms-exercises\\one.xyz")