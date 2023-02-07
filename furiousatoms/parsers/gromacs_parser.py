import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from furiousatoms.parsers.parser_util import float_or_zero

class GROMACSParser(BaseParser):
    def __init__(self) -> None:
        super().__init__()
        self.numAtoms = 2 ** 31 #default is near max int
        self.nextWordIndex = 0

    def parseAtomCount(self, line):
            try:
                self.numAtoms = int(line)
            except:
                #Atom count is not truly necessary, but tell the user anyway.
                self.errors += "Line #2 must contain the number of atoms as per GROMACS format."

    def nextWord(self):
        for i in range(self.nextWordIndex, len(self.words)):
            if len(self.words[i]) > 0:
                self.nextWordIndex += 1
                return self.words[i]
        raise EOFError("No more words remain for the current line")

    def parseAtom(self, line):
        try:
            #Since GROMACS supports arbitrary decimal precision, 
            #search for non-blank words only.
            self.nextWordIndex = 0
            self.words = line.split()
            
            #First (0th) word is residue number
            self.nextWord() #skip the word
            #Residue name
            self.nextWord()
            #Atom name
            self.atom_types.append(self.nextWord())
            
            #Words #4-#6 form the XYZ position.
            pos = np.zeros((3))
            pos[0] = float_or_zero(self.nextWord())
            pos[1] = float_or_zero(self.nextWord())
            pos[2] = float_or_zero(self.nextWord())
            self.positions.append(pos)
            
            #Remaining words #7-9 form the velocity, but we skip these.
        except:
            self.errors += "Unable to parse atom on line #%d"%(self.lineId)

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