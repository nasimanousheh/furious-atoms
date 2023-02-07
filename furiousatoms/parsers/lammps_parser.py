import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from furiousatoms.parsers.parser_util import float_or_zero

#TODO implement box size

class LAMMPSParser(BaseParser):
    def parseAtom(self, line):
        '''
        Atom lines come in many styles, but we assume the XYZ coordinates are the last 3 
        values as per the 'atomic' and 'full' styles. 
        '''
        try:
            pos = np.zeros((3))
            words = line.split()
            if len(words) >= 7: #'Full' style
                self.atom_types.append(words[2])
                pos[0] = float_or_zero(words[4])
                pos[1] = float_or_zero(words[5])
                pos[2] = float_or_zero(words[6])
            elif len(words) >= 5: #'Atomic' style
                self.atom_types.append(words[1])
                pos[0] = float_or_zero(words[2])
                pos[1] = float_or_zero(words[3])
                pos[2] = float_or_zero(words[4])
            else:
                raise IndexError("Line too short")
            self.positions.append(pos)
        except IndexError:
            self.errors += "Line #%d is too short to read an atom's position. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.lineId)
        except:
            self.errors += "Failure processing line #%d.\n"%(self.lineId)
    
    def parseBond(self, line):
        words = line.split()
        try:
            if len(words) >= 4:
                bond = np.zeros((2), dtype='int')
                #Record IDs of the two atoms to connect.
                bond[0] = int(words[2]) - 1
                bond[1] = int(words[3]) - 1
                if bond[0] < 0 or bond[1] < 0:
                    raise ValueError("Atom ID cannot be zero or negative.")
                self.bonds.append(bond)
            else: 
                raise IndexError("Line too short")
        except:
            if len(words) < 2:
                self.errors += "Line #%d is too short. At least two atom IDs are needed for a bond.\n"%(self.lineId)
            else:
                self.errors += "Could not read a bond from line #%d since it is not well-formatted. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.lineId)

    def parseLine(self, line):

        #Choose how to parse subsequent lines
        header = line.strip().lower()
        if len(header) == 0:
            pass
        elif header.startswith("masses"):
            self.parserMethod = None
        elif header.startswith("atoms"):
            self.parserMethod = self.parseAtom
        elif header.startswith("bonds"):
            self.parserMethod = self.parseBond
        elif len(header.split()) == 1:
            #For non-implemented headers like Angles and Dihedrals
            self.parserMethod = None
        
        #Perform the parse
        elif len(header.strip()) > 0:
            if self.parserMethod != None:
                self.parserMethod(line)
            else:
                self.errors += "Unable to process line #%d.\n"%(self.lineId)
        self.lineId += 1
        
#TODO remove
if __name__ == "__main__":
    parser = LAMMPSParser()
    parser.parse("C:\\Users\\Pete\\Desktop\\\Example_with_less_atoms\\graphdiyne_unitcell.data")