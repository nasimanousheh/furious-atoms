import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from furiousatoms.parsers.parser_util import float_or_zero, has_bond
from furiousatoms.element_lookup import lookup_element_by_mass


class LAMMPSParser(BaseParser):
    def __init__(self) -> None:
        super().__init__()
        self.elemInfoDict = {}
    
    def parseMass(self, line):
        words = line.split()
        key = words[0]
        mass = float(words[1])
        self.elemInfoDict[key] = lookup_element_by_mass(mass)
    
    def getAtomType(self, typeKey: int):
        #Assume element symbol has already been read in by parseMass(...).
        #If not, return typeKey.
        if typeKey in self.elemInfoDict:
            return self.elemInfoDict[typeKey]['symbol']
        return typeKey

    def parseAtom(self, line):
        '''
        Atom lines come in many styles, but we assume the XYZ coordinates are the last 3 
        values as per the 'atomic' and 'full' styles. 
        '''
        try:
            pos = np.zeros((3))
            typ = ''
            words = line.split()
            if len(words) >= 7: #'Full' style
                typ = self.getAtomType(words[2])
                pos[0] = float_or_zero(words[4])
                pos[1] = float_or_zero(words[5])
                pos[2] = float_or_zero(words[6])
            elif len(words) >= 5: #'Atomic' style
                typ = self.getAtomType(words[1])
                pos[0] = float_or_zero(words[2])
                pos[1] = float_or_zero(words[3])
                pos[2] = float_or_zero(words[4])
            else:
                raise IndexError("Line too short")
            self.atom_types.append(typ)
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
                if bond[0] != bond[1] and not has_bond(self.bonds, bond):
                    self.bonds.append(bond)
            else: 
                raise IndexError("Line too short")
        except:
            if len(words) < 2:
                self.errors += "Line #%d is too short. At least two atom IDs are needed for a bond.\n"%(self.lineId)
            else:
                self.errors += "Could not read a bond from line #%d since it is not well-formatted. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.lineId)

    def parseBoxSize(self, line):
        '''
        A box size dimension will look like this:
        ` -2.365000 7.095000  xlo xhi` 
        This function checks what label is at the end of the line and extracts the hi dimension.
        '''
        try:
            words = line.split()
            trimmed = line.strip()
            if len(words) >= 4: 
                #Always parse the "hi" dimension
                if trimmed.endswith("xlo xhi"):
                    self.box_size[0] = float_or_zero(words[1])
                elif trimmed.endswith("ylo yhi"):
                    self.box_size[1] = float_or_zero(words[1])
                elif trimmed.endswith("zlo zhi"):
                    self.box_size[2] = float_or_zero(words[1])
                else:
                    raise ValueError("Invalid format")
            else:
                self.errors += "Line #%d is too short to read an atom's position. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.lineId)
        except:
            self.errors += "Failure processing line #%d.\n"%(self.lineId)

    def skipLine(self, line):
        pass

    def parseLine(self, line):
        #Choose how to parse subsequent lines
        line = line.strip()
        header = line.lower()
        if len(header) == 0:
            pass
        elif header.startswith("masses"):
            self.parserMethod = self.parseMass
        elif header.startswith("atoms"):
            self.parserMethod = self.parseAtom
        elif header.startswith("bonds"):
            self.parserMethod = self.parseBond
        elif header.endswith("hi") and "lo" in header:
            #Special case: box size doesn't use a header
            self.parseBoxSize(line)
        elif "coeffs" in header:
            self.parserMethod = self.skipLine
        elif len(header.split()) == 1:
            #For non-implemented headers like Angles and Dihedrals
            self.parserMethod = self.skipLine
        
        #Perform the parse
        elif len(header.strip()) > 0:
            if self.parserMethod != None:
                self.parserMethod(line)
            else:
                self.errors += "Unable to process line #%d.\n"%(self.lineId)
        self.lineId += 1