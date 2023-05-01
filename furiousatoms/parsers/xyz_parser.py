import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from furiousatoms.parsers.parser_util import float_or_zero

class XYZParser(BaseParser):
    '''
    Based on https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html
    '''

    def __init__(self) -> None:
        super().__init__()
        self.next_word_index = 0

    def next_word(self):
        for i in range(self.next_word_index, len(self.words)):
            if len(self.words[i]) > 0:
                self.next_word_index = i + 1
                return self.words[i]
        raise EOFError("No more words remain for the current line")


    def parse_line(self, line):
        #Ignore 0th & 1st lines
        if len(line.strip()) > 0 and self.line_id >= 2:
            atomType = None
            pos = pos = np.zeros((3))
            try:
                #Since XYZ supports arbitrary decimal precision, 
                #search for non-blank words only.
                self.next_word_index = 0
                self.words = line.split()
                if len(self.words) < 2:
                    self.words = line.split(",")
                
                #Element symbol or atomic number
                atomType = self.next_word()

                #Parse XYZ position
                pos[0] = float_or_zero(self.next_word())
                pos[1] = float_or_zero(self.next_word())
                pos[2] = float_or_zero(self.next_word())
            except:
                self.errors += "Unable to parse atom on line #%d"%(self.line_id)
            else:
                self.atom_types.append(atomType)
                self.positions.append(pos)