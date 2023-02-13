import unittest
import pytest

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
    FNAME = "./examples/WaterProtein/water_protein.xyz"
    
    parser.parse(FNAME)

    assert parser.box_size == [0, 0, 0]

    assert len(parser.positions) == 1800
    assert list(parser.positions[0]) == [17.047001, 14.099, 3.625]
    assert list(parser.positions[1]) == [16.966999, 12.784, 4.338]
    assert list(parser.positions[85]) == [5.929, 6.358, 5.055]
    assert list(parser.positions[1799]) == [12.951, 13.245, -2.112]

    assert len(parser.bonds) == 0

    assert len(parser.atom_types) == 1800
    assert parser.atom_types[0] == 'N'
    assert parser.atom_types[1] == 'C'
    assert parser.atom_types[42] == 'O'
    assert parser.atom_types[1799] == 'H'


if __name__ == '__main__':
    unittest.main()