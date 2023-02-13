import unittest
import pytest

from furiousatoms.parsers.lammps_parser import LAMMPSParser

@pytest.fixture
def newLAMMPSParser():
    parser = LAMMPSParser()
    return parser


def test_parseAtom(newLAMMPSParser):
    parser = newLAMMPSParser
    LINE = "1 1 3 0.000000 17.047001 14.099000 3.625000 # N MOL"

    parser.parseAtom(LINE)

    assert len(parser.positions) == 1
    assert parser.positions[0][0] == 17.047001
    assert parser.positions[0][1] == 14.099000
    assert parser.positions[0][2] == 3.625000

    assert len(parser.atom_types) == 1
    assert parser.atom_types[0] == 'N' #FIXME

def test_parseBond(newLAMMPSParser):
    parser = newLAMMPSParser
    LINE = "791 1 880 882"

    parser.parseBond(LINE)

    assert len(parser.bonds) == 1
    assert parser.bonds[0][0] == 879
    assert parser.bonds[0][1] == 881

def test_parseBoxSize(newLAMMPSParser):
    parser = newLAMMPSParser
    LINE1 = " -8.651000 21.349000  xlo xhi"
    LINE2 = " -10.067500 19.932500  ylo yhi"
    LINE3 = " -5.881000 16.639000  zlo zhi"

    parser.parseBoxSize(LINE1)
    assert parser.box_size == [21.349, 0, 0]
    parser.parseBoxSize(LINE2)
    assert parser.box_size == [21.349, 19.9325, 0]
    parser.parseBoxSize(LINE3)
    assert parser.box_size == [21.349, 19.9325, 16.639]
    
    assert len(parser.positions) == 0
    assert len(parser.bonds) == 0
    assert len(parser.atom_types) == 0

def test_parseLine_parserMethod(newLAMMPSParser):
    parser = newLAMMPSParser
    
    parser.lineId = 0
    assert parser.parserMethod == None #header of file

    parser.lineId = 23
    LINE = " -10.067500 19.932500  ylo yhi"
    parser.parseLine(LINE)
    assert parser.parserMethod == None #Note: NOT parseBoxSize

    #TODO test parseMass
    # parser.lineId = 1
    # LINE = " Masses"
    # parser.parseLine(LINE)
    # assert parser.parserMethod == parser.parseMass

    parser.lineId = 42
    LINE = " Atoms # full"
    parser.parseLine(LINE)
    assert parser.parserMethod == parser.parseAtom

    parser.lineId = 86
    LINE = "Bonds"
    parser.parseLine(LINE)
    assert parser.parserMethod == parser.parseBond

def test_parse(newLAMMPSParser):
    parser = newLAMMPSParser
    FNAME = "./examples/WaterProtein/water_protein.data"
    
    parser.parse(FNAME)

    assert parser.box_size == [21.349, 19.9325, 16.639]

    assert len(parser.positions) == 1800
    assert list(parser.positions[0]) == [17.047001, 14.099, 3.625]
    assert list(parser.positions[1]) == [16.966999, 12.784, 4.338]
    assert list(parser.positions[1625]) == [9.849, 3.939, 0.99]
    assert list(parser.positions[1799]) == [12.951, 13.245, -2.112]

    assert len(parser.bonds) == 1403
    assert parser.bonds[0][0] == 0
    assert parser.bonds[0][1] == 1
    assert parser.bonds[1][0] == 1
    assert parser.bonds[1][1] == 2
    assert parser.bonds[783][0] == 870
    assert parser.bonds[783][1] == 872
    assert parser.bonds[1402][0] == 1797
    assert parser.bonds[1402][1] == 1799

    assert len(parser.atom_types) == 1800
    assert parser.atom_types[0] == 'N'
    assert parser.atom_types[1] == 'C'
    assert parser.atom_types[42] == 'O'
    assert parser.atom_types[1799] == 'H'


if __name__ == '__main__':
    unittest.main()
