import unittest
import pytest
import os

from furiousatoms.parsers.gromacs_parser import GROMACSParser

@pytest.fixture
def newGROMACSParser():
    parser = GROMACSParser()
    return parser


def test_next_word(newGROMACSParser):
    parser = newGROMACSParser
    parser.words = "   9.00000   2.00050   9.00000".split(" ")
    assert parser.next_word() == "9.00000"
    assert parser.next_word() == "2.00050"
    assert parser.next_word() == "9.00000"
    with pytest.raises(EOFError) as err:
        parser.next_word()
        parser.next_word()

def test_parse_atom(newGROMACSParser):
    parser = newGROMACSParser
    LINE = "    1MOL      C    3   0.068  -0.000   0.062"

    parser.parse_atom(LINE)
    parser.parse_line(LINE)

    assert len(parser.atom_types) == 1
    assert parser.atom_types[0] == 'C'

    assert len(parser.positions) == 1
    assert parser.positions[0][0] == 0.068
    assert parser.positions[0][1] == 0.000
    assert parser.positions[0][2] == 0.062

def test_parse_line_parser_method(newGROMACSParser):
    parser = newGROMACSParser
    
    parser.line_id = 0
    assert parser.parser_method == None #header of file

    parser.line_id = 1
    LINE = "    4"
    parser.parse_line(LINE)
    assert parser.parser_method == parser.parse_atom_count

    parser.line_id = 2
    LINE = "    1MOL      C    4  -0.034   0.059   0.062"
    parser.parse_line(LINE)
    assert parser.parser_method == parser.parse_atom

    #Ensure the next line will be treated as the last one
    #i.e. if all 4 positions have been read in
    parser.line_id = 4 + 2
    LINE = "   9.00000   9.00000   9.00000"
    parser.positions = [[0,1,2], [0,1,2], [0,1,2], [0,1,2]]
    parser.parse_line(LINE)
    assert parser.parser_method == parser.parse_box_size

def test_parse_box_size(newGROMACSParser):
    parser = newGROMACSParser
    LINE = "   0.94600   0.81926   6.00000   0.00000   0.00000  -0.47300   0.00000   0.00000   0.00000"
    
    parser.parse_box_size(LINE)

    assert parser.box_size == [0.946, 0.81926, 6.0]
    assert len(parser.positions) == 0
    assert len(parser.bonds) == 0
    assert len(parser.atom_types) == 0

def test_parse(newGROMACSParser):
    parser = newGROMACSParser
    FPATH = os.path.join('furiousatoms','tests','test_data','CB_18','CB_18.gro')
    
    parser.parse(FPATH)

    assert parser.box_size == [9.0, 9.0, 9.0]

    assert len(parser.positions) == 18
    assert list(parser.positions[0]) == [-0.546, -0.000, 0.000]
    assert list(parser.positions[1]) == [-0.390, -0.000, 0.000]
    assert list(parser.positions[15]) == [0.312, 0.135, 0.000]
    assert list(parser.positions[17]) == [0.546, -0.000, 0.000]

    assert len(parser.bonds) == 0

    assert len(parser.atom_types) == 18
    assert parser.atom_types[0] == 'C'
    assert parser.atom_types[1] == 'B'
    assert parser.atom_types[6] == 'C'
    assert parser.atom_types[17] == 'B'


if __name__ == '__main__':
    unittest.main()