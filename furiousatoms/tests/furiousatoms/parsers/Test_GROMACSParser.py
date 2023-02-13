import unittest
import pytest

from furiousatoms.parsers.gromacs_parser import GROMACSParser

@pytest.fixture
def newGROMACSParser():
    parser = GROMACSParser()
    return parser


def test_nextWord(newGROMACSParser):
    parser = newGROMACSParser
    parser.words = "   9.00000   2.00050   9.00000".split(" ")
    assert parser.nextWord() == "9.00000"
    assert parser.nextWord() == "2.00050"
    assert parser.nextWord() == "9.00000"
    with pytest.raises(EOFError) as err:
        parser.nextWord()
        parser.nextWord()

def test_parseAtom(newGROMACSParser):
    parser = newGROMACSParser
    LINE = "    1MOL      C    3   0.068  -0.000   0.062"

    parser.parseAtom(LINE)
    parser.parseLine(LINE)

    assert len(parser.atom_types) == 1
    assert parser.atom_types[0] == 'C'

    assert len(parser.positions) == 1
    assert parser.positions[0][0] == 0.068
    assert parser.positions[0][1] == 0.000
    assert parser.positions[0][2] == 0.062

def test_parseLine_parserMethod(newGROMACSParser):
    parser = newGROMACSParser
    
    parser.lineId = 0
    assert parser.parserMethod == None #header of file

    parser.lineId = 1
    LINE = "    4"
    parser.parseLine(LINE)
    assert parser.parserMethod == parser.parseAtomCount

    parser.lineId = 2
    LINE = "    1MOL      C    4  -0.034   0.059   0.062"
    parser.parseLine(LINE)
    assert parser.parserMethod == parser.parseAtom

    #Ensure the next line will be treated as the last one
    #i.e. if all 4 positions have been read in
    parser.lineId = 4 + 2
    LINE = "   9.00000   9.00000   9.00000"
    parser.positions = [[0,1,2], [0,1,2], [0,1,2], [0,1,2]]
    parser.parseLine(LINE)
    assert parser.parserMethod == parser.parseBoxSize

def test_parseBoxSize(newGROMACSParser):
    parser = newGROMACSParser
    LINE = "   0.94600   0.81926   6.00000   0.00000   0.00000  -0.47300   0.00000   0.00000   0.00000"
    
    parser.parseBoxSize(LINE)

    assert parser.box_size == [0.946, 0.81926, 6.0]
    assert len(parser.positions) == 0
    assert len(parser.bonds) == 0
    assert len(parser.atom_types) == 0

#TODO Test GROMACS parse after transform RtoS bug is fixed
# def test_parse(newGROMACSParser):
#     parser = newGROMACSParser
#     FNAME = "./examples/WaterProtein/water_protein.gro"
    
#     parser.parse(FNAME)

#     assert parser.box_size == [3.0, 3.0, 2.252]

#     assert len(parser.positions) == 1800
#     assert list(parser.positions[0]) == [17.047001, 14.099, 3.625]
#     assert list(parser.positions[1]) == [16.966999, 12.784, 4.338]
#     assert list(parser.positions[85]) == [5.929, 6.358, 5.055]
#     assert list(parser.positions[1799]) == [12.951, 13.245, -2.112]

#     assert len(parser.bonds) == 0

#     assert len(parser.atom_types) == 1800
#     assert parser.atom_types[0] == 'N'
#     assert parser.atom_types[1] == 'C'
#     assert parser.atom_types[42] == 'O'
#     assert parser.atom_types[1799] == 'H'


if __name__ == '__main__':
    unittest.main()