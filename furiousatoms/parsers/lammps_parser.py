import numpy as np
from util import float_or_zero

class LAMMPSParser:
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
                #TODO remove
                # #Since some files have more spacing than necessary, search for 
                # #the fifth non-blank `word` (number) and record its index in `wordsRead`.
                # wordsRead = 0
                # words = line.split()
                # for word in words:
                #     if wordsRead == 4:
                #         break
                #     if len(word) > 0:
                #         wordsRead += 1
                # pos[0] = float_or_zero(words[wordsRead])
                # pos[1] = float_or_zero(words[wordsRead + 1])
                # pos[2] = float_or_zero(words[wordsRead + 2])

                words = line.split()
                if len(words) >= 7:
                    pos[0] = float_or_zero(words[4])
                    pos[1] = float_or_zero(words[5])
                    pos[2] = float_or_zero(words[6])
                    self.positions.append(pos)
                elif len(words) == 6: #If element is absent
                    pos[0] = float_or_zero(words[3])
                    pos[1] = float_or_zero(words[4])
                    pos[2] = float_or_zero(words[5])
                    self.positions.append(pos)
                elif len(words) == 5: 
                    pos[0] = float_or_zero(words[2])
                    pos[1] = float_or_zero(words[3])
                    pos[2] = float_or_zero(words[4])
                    self.positions.append(pos)
                else:
                    raise IndexError("Line too short")
            except IndexError:
                self.errors += "Line #%d is too short to read an atom's position. Please refer to a guide for LAMMPS format if you are unsure.\n"%(self.lineId)
            except:
                self.errors += "Failure processing line #%d.\n"%(self.lineId)
        
        def parseBond(line):
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


        with open(fname, "r") as fp:
            parserMethod = None
            while True:
                line = fp.readline()
                if not line:
                    break
                
                #Choose how to parse subsequent lines
                header = line.strip().lower()
                if len(header) == 0:
                    pass
                elif header.startswith("masses") or header.startswith("atoms"):
                    parserMethod = parseAtom
                elif header.startswith("bonds"):
                    parserMethod = parseBond
                elif header.startswith("angles"):
                    parserMethod = None #TODO: Implement Angles
                elif header.startswith("dihedrals"):
                    parserMethod = None #TODO: Implement Dihedrals
                
                #Perform the parse
                elif len(line.strip()) > 0 and parserMethod != None:
                    parserMethod(line)
                else:
                    self.errors += "Unable to process line #%d.\n"%(self.lineId)
                self.lineId += 1
        
        print(self.errors) #TODO display popup
    
        return self.box_size, np.array(self.positions), np.array(self.bonds), np.array(self.atom_types)