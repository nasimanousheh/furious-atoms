import unittest
import pytest
import os

from furiousatoms.parsers.xyz_parser import XYZParser

@pytest.fixture
def newXYZParser():
    parser = XYZParser()
    return parser


def test_next_word(newXYZParser):
    parser = newXYZParser
    parser.words = "This is a sentence.\n".split(" ")
    assert parser.next_word() == "This"
    assert parser.next_word() == "is"
    assert parser.next_word() == "a"
    assert parser.next_word() == "sentence.\n"
    with pytest.raises(EOFError) as err:
        parser.next_word()
        parser.next_word()

def test_parse_line(newXYZParser):
    parser = newXYZParser
    parser.line_id = 33
    LINE = "       C    -0.67848    0.00000   -0.61531"

    parser.parse_line(LINE)

    assert len(parser.atom_types) == 1
    assert parser.atom_types[0] == 'C'

    assert len(parser.positions) == 1
    assert parser.positions[0][0] == -0.67848
    assert parser.positions[0][1] ==  0.00000
    assert parser.positions[0][2] == -0.61531


def test_parse(newXYZParser):
    parser = newXYZParser
    FPATH = os.path.join('furiousatoms','tests','test_data','CB_18','CB_18.xyz')
    
    parser.parse(FPATH)

    assert parser.box_size == [0, 0, 0]

    assert len(parser.positions) == 18
    assert list(parser.positions[0]) == [-5.46, -0.0, 0.0]
    assert list(parser.positions[1]) == [-3.9, -0.0, 0.0]
    assert list(parser.positions[13]) == [0.78, 2.702, 0.0]
    assert list(parser.positions[17]) == [5.46, -0.0, 0.0]

    assert len(parser.bonds) == 0

    assert len(parser.atom_types) == 18
    assert parser.atom_types[0] == 'C'
    assert parser.atom_types[1] == 'B'
    assert parser.atom_types[13] == 'B'
    assert parser.atom_types[17] == 'B'


if __name__ == '__main__':
    unittest.main()