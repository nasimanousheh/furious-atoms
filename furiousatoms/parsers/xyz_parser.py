import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from util import float_or_zero

class XYZParser(BaseParser):
    def __init__(self) -> None:
        super().__init__()
        self.nextWordIndex = 0

    def nextWord(self):
        for i in range(self.nextWordIndex, len(self.words)):
            if len(self.words[i]) > 0:
                self.nextWordIndex += 1
                return self.words[i]
        raise EOFError("No more words remain for the current line")


    def parseLine(self, line):
        #Ignore 0th & 1st lines
        if len(line.strip()) > 0 and self.lineId >= 2:
            try:
                #Since XYZ supports arbitrary decimal precision, 
                #search for non-blank words only.
                self.nextWordIndex = 0
                self.words = line.split()
                
                #Element symbol or atomic number
                self.nextWord() #skip the word

                #Parse XYZ position
                pos = np.zeros((3))
                pos[0] = float_or_zero(self.nextWord())
                pos[1] = float_or_zero(self.nextWord())
                pos[2] = float_or_zero(self.nextWord())
                self.positions.append(pos)
            except:
                self.errors += "Unable to parse atom on line #%d"%(self.lineId)

        self.lineId += 1

#TODO remove
if __name__ == "__main__":
    parser = XYZParser()
    parser.parse("C:\\Users\\Pete\\Desktop\\\Example_with_less_atoms\\graphdiyne_unitcell.xyz")