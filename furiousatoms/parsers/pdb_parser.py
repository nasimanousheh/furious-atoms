import numpy as np
from util import float_or_zero

#TODO maybe make error handlers less rigid?
#TODO atom types!

class PDBParser:
    def parse(self, fname):
        box_size = [0, 0, 0]
        positions = []
        bonds = []
        atom_types = []
        errors = "" #error messages to help user debug their file
        with open(fname, "r") as fp:
            lineId = 0
            while True:
                line = fp.readline()
                if not line:
                    break

                if line.startswith("CONECT"): #Note the spelling. 
                    for i, j in ((12,16), (17,21), (22,26)): #indices of atom IDs to connect to within the string
                        try:
                            otherId = int(line[i:j]) - 1 #ID of atom to connect to
                            bond = np.zeros((2), dtype='int')
                            bond[0] = abs(int(line[7:11]) - 1) #ID of first atom
                            bond[1] = abs(otherId)
                            if bond[0] != bond[1]:
                                bonds.append(bond)
                        except ValueError:
                            pass
                            # errors += "Refusing to connect a bond on line #%d.\n"%(lineId)
                        except IndexError:
                            if len(line) < 13:
                                errors += "Line #%d is too short. At least two atom IDs are needed for a bond.\n"%(lineId)
                elif line.startswith("ATOM"):
                    try:
                        pos = np.zeros((3))
                        pos[0] = float_or_zero(line[31:38])
                        pos[1] = float_or_zero(line[39:46])
                        pos[2] = float_or_zero(line[47:54])
                        positions.append(pos)
                    except:
                        errors += "Failure processing line #%d.\n"%(lineId)
                elif line.startswith("REMARK"):
                    pass
                elif line.startswith("CRYST1 "):
                    if len(line) >= 33:
                        box_size[0] = float_or_zero(line[7:15])
                        box_size[1] = float_or_zero(line[16:24])
                        box_size[2] = float_or_zero(line[25:33])
                    else:
                        errors += "Line #%d is too short. Three numbers are needed for the box size.\n"%(lineId)
                else:
                    errors += "Unable to process line #%d.\n"%(lineId)
                lineId += 1
        
        print(errors) #TODO display popup
    
        return box_size, np.array(positions), np.array(bonds), np.array(atom_types)