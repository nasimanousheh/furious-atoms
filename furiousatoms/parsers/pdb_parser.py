import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from furiousatoms.parsers.parser_util import float_or_zero, has_bond

class PDBParser(BaseParser):
    def parseLine(self, line):
        lineUpper = line.upper()
        if lineUpper.startswith("CONECT"): #Note the spelling. 
            for i, j in ((12,16), (17,21), (22,26)): #indices of atom IDs to connect to within the string
                try:
                    otherId = int(line[i:j]) - 1 #ID of atom to connect to
                    bond = np.zeros((2), dtype='int')
                    bond[0] = abs(int(line[7:11]) - 1) #ID of first atom
                    bond[1] = abs(otherId)
                    if bond[0] != bond[1] and not has_bond(self.bonds, bond):
                        self.bonds.append(bond)
                except ValueError:
                    self.errors += "Refusing to connect a bond on line #%d.\n"%(self.lineId)
                except IndexError:
                    if len(line) < 13:
                        self.errors += "Line #%d is too short. At least two atom IDs are needed for a bond.\n"%(self.lineId)
        elif lineUpper.startswith("ATOM") or lineUpper.startswith("HETATM"):
            try:
                self.atom_types.append(line[13:15].strip())

                pos = np.zeros((3))
                pos[0] = float_or_zero(line[31:38])
                pos[1] = float_or_zero(line[39:46])
                pos[2] = float_or_zero(line[47:54])
                self.positions.append(pos)
            except:
                self.errors += "Failure processing line #%d.\n"%(self.lineId)
        elif lineUpper.startswith("REMARK"):
            pass
        elif lineUpper.startswith("CRYST1 "):
            if len(line) >= 33:
                self.box_size[0] = float_or_zero(line[7:15])
                self.box_size[1] = float_or_zero(line[16:24])
                self.box_size[2] = float_or_zero(line[25:33])
            else:
                self.errors += "Line #%d is too short. Three numbers are needed for the box size.\n"%(self.lineId)
        else:
            self.errors += "Unable to process line #%d.\n"%(self.lineId)