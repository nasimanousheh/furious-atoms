import numpy as np
from furiousatoms.parsers.base_parser import BaseParser
from furiousatoms.parsers.parser_util import float_or_zero

class GROMACSParser(BaseParser):
    '''
    Based on https://manual.gromacs.org/archive/5.0.3/online/gro.html
    '''
    def __init__(self) -> None:
        super().__init__()
        self.num_atoms = 2 ** 31 #default is near max int
        self.next_word_index = 0

    def parse_atom_count(self, line):
        try:
            self.num_atoms = int(line)
        except:
            self.errors += "Line #2 must contain the number of atoms as per GROMACS format.\n"

    def next_word(self):
        for i in range(self.next_word_index, len(self.words)):
            if len(self.words[i]) > 0:
                self.next_word_index = i + 1
                return self.words[i]
        raise EOFError("No more words remain for the current line")

    def parse_atom(self, line):
        '''
        Extract the atom position and atom type of one line
        '''
        atomType = None
        pos = np.zeros((3))
        try:
            #Since GROMACS supports arbitrary decimal precision, 
            #search for non-blank words only.
            self.next_word_index = 0
            self.words = line.split()
            
            #First (0th) word is residue number
            self.next_word() #skip the word
            #Atom name/abbreviation
            atomType = self.next_word()
            #Atom ID (skipped)
            self.next_word()
            
            #Words #4-#6 form the XYZ position.
            pos[0] = float_or_zero(self.next_word())
            pos[1] = float_or_zero(self.next_word())
            pos[2] = float_or_zero(self.next_word())
            
            #Remaining words #7-9 form the velocity, but we skip these.
        except:
            self.errors += "Unable to parse atom on line #%d\n"%(self.line_id)
        else:
            self.atom_types.append(atomType)
            self.positions.append(pos)

    def parse_box_size(self, line):
        wordsRead = 0
        words = line.split()
        for word in words:
            if len(word) > 0:
                wordsRead += 1

            #First 3 words indicate XYZ box size
            if wordsRead >= 1 and wordsRead <= 3: 
                self.box_size[wordsRead - 1] = float_or_zero(word)
        if wordsRead < 3:
            self.errors += "At least 3 dimensions are needed for the box size on line #%d.\n"%(self.line_id)

    def parse_line(self, line):
        #Choose how to parse the current line
        if self.line_id == 1:
            self.parser_method = self.parse_atom_count
        elif len(self.positions) == self.num_atoms: #"final" line: all atoms have been read
            self.parser_method = self.parse_box_size
        elif self.line_id != 0: #line 0 is header/comment; all other lines are atoms
            self.parser_method = self.parse_atom
        
        #Perform the parse
        if len(line.strip()) > 0 and self.parser_method != None:
            self.parser_method(line)
        else:
            self.errors += "Unable to process line #%d.\n"%(self.line_id)