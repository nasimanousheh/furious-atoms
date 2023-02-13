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
            self.errors += "Line #2 must contain the number of atoms as per GROMACS format."

    def nextWord(self):
        for i in range(self.nextWordIndex, len(self.words)):
            if len(self.words[i]) > 0:
                self.nextWordIndex = i + 1
                return self.words[i]
        raise EOFError("No more words remain for the current line")

    def parseAtom(self, line):
        '''
        Extract the atom position and atom type of one line
        '''
        atomType = None
        pos = np.zeros((3))
        try:
            #Since GROMACS supports arbitrary decimal precision, 
            #search for non-blank words only.
            self.nextWordIndex = 0
            self.words = line.split()
            
            #First (0th) word is residue number
            self.nextWord() #skip the word
            #Atom name/abbreviation
            atomType = self.nextWord()
            #Atom ID (skipped)
            self.nextWord()
            
            #Words #4-#6 form the XYZ position.
            pos[0] = float_or_zero(self.nextWord())
            pos[1] = float_or_zero(self.nextWord())
            pos[2] = float_or_zero(self.nextWord())
            
            #Remaining words #7-9 form the velocity, but we skip these.
        except:
            self.errors += "Unable to parse atom on line #%d"%(self.lineId)
        else:
            self.atom_types.append(atomType)
            self.positions.append(pos)

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
        elif len(self.positions) == self.numAtoms: #"final" line: all atoms have been read
            self.parserMethod = self.parseBoxSize
        elif self.lineId != 0: #line 0 is header/comment; all other lines are atoms
            self.parserMethod = self.parseAtom
        
        #Perform the parse
        if len(line.strip()) > 0 and self.parserMethod != None:
            self.parserMethod(line)
        else:
            self.errors += "Unable to process line #%d.\n"%(self.lineId)