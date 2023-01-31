import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from util import float_or_zero

#TODO: Refactor to use iterator pattern for the words

class GROMACSParser(BaseParser):
    def __init__(self) -> None:
        super().__init__()
        self.numAtoms = 2 ** 31 #default is near max int

    def parseAtomCount(self, line):
            try:
                self.numAtoms = int(line)
            except:
                #Atom count is not truly necessary, but tell the user anyway.
                self.errors += "Line #2 must contain the number of atoms as per GROMACS format."

    def parseAtom(self, line):
        #Since GROMACS supports arbitrary decimal precision, 
        #search for non-blank words only.
        wordsRead = 0
        words = line.split()
        pos = np.zeros((3))
        for word in words:
            if len(word) > 0:
                wordsRead += 1
            
            #First (0th) word is residue number
            if wordsRead == 1: 
                pass
            #Residue name
            elif wordsRead == 2: 
                pass
            #Atom name
            elif wordsRead == 3: 
                pass
            #Word #4 is X, word #5 is Y, and word #6 is Z => store these as a position.
            elif wordsRead >= 4 and wordsRead <= 6: #Position
                k = wordsRead - 4
                pos[k] = float_or_zero(words[wordsRead - 1])
                if wordsRead == 6: #all 3 coords read
                    self.positions.append(pos)
            #Velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places
            elif wordsRead >= 7 and wordsRead <= 9: 
                pass
            else:
                continue

    def parseBoxSize(self, line):
        wordsRead = 0
        words = line.split()
        for word in words:
            if len(word) > 0:
                wordsRead += 1

            #First 3 words indicate XYZ box size
            if wordsRead >= 1 and wordsRead <= 3: 
                self.box_size[wordsRead - 1] = float_or_zero(word)
        if wordsRead < 3:
            self.errors += "At least 3 dimensions are needed for the box size on line #%d."%(self.lineId)

    def parseLine(self, line):
        #Choose how to parse the current line
        if self.lineId == 1:
            self.parserMethod = self.parseAtomCount
        elif self.lineId == self.numAtoms + 2: #"final" line
            self.parserMethod = self.parseBoxSize
        else:
            self.parserMethod = self.parseAtom
        
        #Perform the parse
        if len(line.strip()) > 0 and self.parserMethod != None:
            self.parserMethod(line)
        else:
            self.errors += "Unable to process line #%d.\n"%(self.lineId)
        self.lineId += 1

#TODO remove
if __name__ == "__main__":
    parser = GROMACSParser()
    parser.parse("C:\\Users\\Pete\\Desktop\\furious-atoms-exercises\\one.gro")