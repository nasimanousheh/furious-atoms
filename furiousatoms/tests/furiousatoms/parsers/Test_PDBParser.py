import unittest
import pytest

from furiousatoms.parsers.pdb_parser import PDBParser

@pytest.fixture
def newPDBParser():
    parser = PDBParser()
    return parser


def test_parseLine_boxSize(newPDBParser):
    parser = newPDBParser
    LINE = "CRYST1   90.000   40.000   62.000  90.00  90.00  90.00 P 1           1"

    parser.parseLine(LINE)

    assert parser.box_size == [90.0, 40.0, 62.0]
    
    assert len(parser.positions) == 0
    assert len(parser.bonds) == 0
    assert len(parser.atom_types) == 0

def test_parseLine_positions(newPDBParser):
    parser = newPDBParser
    LINE = "ATOM      1  C   MOL X   1      -0.678   0.000  -0.615 "

    parser.parseLine(LINE)

    assert len(parser.positions) == 1
    assert parser.positions[0][0] == -0.678
    assert parser.positions[0][1] ==  0.000
    assert parser.positions[0][2] == -0.615

def test_parseLine_bonds(newPDBParser):
    parser = newPDBParser
    LINE = "CONECT    1    2    4"

    parser.parseLine(LINE)

    assert len(parser.bonds) == 2
    assert parser.bonds[0][0] == 0
    assert parser.bonds[0][1] == 1
    assert parser.bonds[1][0] == 0
    assert parser.bonds[1][1] == 3

def test_parseLine_atomTypes(newPDBParser):
    parser = newPDBParser
    LINE = "ATOM      3  C   MOL X   1       0.678  -0.000   0.615  1.00  0.00              "

    parser.parseLine(LINE)

    assert len(parser.atom_types) == 1
    assert parser.atom_types[0] == 'C'

def test_parse(newPDBParser):
    parser = newPDBParser
    FNAME = "./examples/WaterProtein/water_protein.pdb"
    
    parser.parse(FNAME)

    assert parser.box_size == [30.0, 30.0, 22.52]

    assert len(parser.positions) == 1800
    assert list(parser.positions[0]) == [17.047, 14.099, 3.625]
    assert list(parser.positions[1]) == [16.967, 12.784, 4.338]
    assert list(parser.positions[85]) == [5.929, 6.358, 5.055]
    assert list(parser.positions[1799]) == [12.951, 13.245, -2.112]

    assert len(parser.bonds) == 1970
    assert parser.bonds[0][0] == 19
    assert parser.bonds[0][1] == 281
    assert parser.bonds[96][0] == 394
    assert parser.bonds[96][1] == 393
    assert parser.bonds[1969][0] == 1799
    assert parser.bonds[1969][1] == 1797

    assert len(parser.atom_types) == 1800
    assert parser.atom_types[0] == 'N'
    assert parser.atom_types[1] == 'C'
    assert parser.atom_types[42] == 'O'
    assert parser.atom_types[1799] == 'H'


if __name__ == '__main__':
    unittest.main()