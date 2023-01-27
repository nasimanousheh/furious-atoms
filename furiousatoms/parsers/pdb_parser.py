'''
"ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
            "{chainID:1s}{resSeq:4d}{iCode:1s}"
            "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
            "{tempFactor:6.2f}      {segID:<4s}{element:>2s}{charge:2s}\n"),
'''
import os
import numpy as np
import traceback

def float_or_zero(val):
    try:
        return float(val)
    except ValueError:
        return 0.0

class PDBParser:
    def parse(self, fname):
        positions = []
        bonds = []
        atom_types = []
        with open(fname, "r") as f:
            for line in f.read().split("\n"):
                if line.startswith("CONECT"): #Note the spelling. 
                    # return NotImplementedError()
                    pass
                elif line.startswith("ATOM"):
                    try:
                        pos = np.zeros((3))
                        # print(line[31:38])
                        # print(line[39:46])
                        # print(line[47:54])
                        pos[0] = float_or_zero(line[31:38])
                        pos[1] = float_or_zero(line[39:46])
                        pos[2] = float_or_zero(line[47:54])
                        positions.append(pos)
                    except:
                        print("Failure processing line `%s`"%line)
                        traceback.print_exc()
                elif line.startswith("LINK"):
                    # return NotImplementedError()
                    pass
                elif line.startswith("REMARK"):
                    pass
                else:
                    #TODO throw custom exception
                    print("Unable to process line `%s`"%line)
    
        return np.array(positions), np.array(bonds), np.array(atom_types)


if __name__ == "__main__":
    parser = PDBParser
    parser.parse(None, "C:\\Users\\Pete\\Desktop\\furious-atoms-exercises\\CB_18.pdb")