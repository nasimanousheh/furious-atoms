import unittest
import pytest
import os

from furiousatoms.parsers.xyz_parser import XYZParser

@pytest.fixture
def newXYZParser():
    parser = XYZParser()
    return parser


def test_nextWord(newXYZParser):
    parser = newXYZParser
    parser.words = "This is a sentence.\n".split(" ")
    assert parser.nextWord() == "This"
    assert parser.nextWord() == "is"
    assert parser.nextWord() == "a"
    assert parser.nextWord() == "sentence.\n"
    with pytest.raises(EOFError) as err:
        parser.nextWord()
        parser.nextWord()

def test_parseLine(newXYZParser):
    parser = newXYZParser
    parser.lineId = 33
    LINE = "       C    -0.67848    0.00000   -0.61531"

    parser.parseLine(LINE)

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