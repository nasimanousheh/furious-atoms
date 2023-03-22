import unittest
import pytest
import os

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
    FPATH = os.path.join('furiousatoms','tests','test_data','CB_18','CB_18.pdb')
    
    parser.parse(FPATH)

    assert parser.box_size == [90.0, 90.0, 90.0]

    assert len(parser.positions) == 18
    assert list(parser.positions[0]) == [-5.460, -0.000, 0.000]
    assert list(parser.positions[1]) == [-3.900, -0.000, 0.000]
    assert list(parser.positions[3]) == [-1.560, -1.351, 0.000]
    assert list(parser.positions[17]) == [5.460, -0.000, 0.000]

    assert len(parser.bonds) == 21
    assert parser.bonds[0][0] == 0
    assert parser.bonds[0][1] == 1
    assert parser.bonds[25][0] == 10
    assert parser.bonds[25][1] == 9
    assert parser.bonds[41][0] == 17
    assert parser.bonds[41][1] == 16

    assert len(parser.atom_types) == 18
    assert parser.atom_types[0] == 'C'
    assert parser.atom_types[1] == 'B'
    assert parser.atom_types[12] == 'C'
    assert parser.atom_types[17] == 'B'


if __name__ == '__main__':
    unittest.main()