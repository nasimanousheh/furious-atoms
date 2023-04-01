import unittest
import pytest
import os

from furiousatoms.parsers.lammps_parser import LAMMPSParser

@pytest.fixture
def newLAMMPSParser():
    parser = LAMMPSParser()
    return parser


def test_parse_atom(newLAMMPSParser):
    parser = newLAMMPSParser
    LINE = "1 1 3 0.000000 17.047001 14.099000 3.625000 # N MOL"

    parser.parse_atom(LINE)

    assert len(parser.positions) == 1
    assert parser.positions[0][0] == 17.047001
    assert parser.positions[0][1] == 14.099000
    assert parser.positions[0][2] == 3.625000

    assert len(parser.atom_types) == 1
    #fall back to number if element symbol is unknown
    assert parser.atom_types[0] == '3' 

def test_parse_bond(newLAMMPSParser):
    parser = newLAMMPSParser
    LINE = "791 1 880 882"

    parser.parse_bond(LINE)

    assert len(parser.bonds) == 1
    assert parser.bonds[0][0] == 879
    assert parser.bonds[0][1] == 881

def test_parse_box_size(newLAMMPSParser):
    parser = newLAMMPSParser
    LINE1 = " -8.651000 21.349000  xlo xhi"
    LINE2 = " -10.067500 19.932500  ylo yhi"
    LINE3 = " -5.881000 16.639000  zlo zhi"

    parser.parse_box_size(LINE1)
    assert parser.box_size == [21.349 + 8.651000, 0, 0]
    parser.parse_box_size(LINE2)
    assert parser.box_size == [21.349 + 8.651000, 19.932500 + 10.067500, 0]
    parser.parse_box_size(LINE3)
    assert parser.box_size == [21.349 + 8.651000, 19.932500 + 10.067500, 16.639000 + 5.881000]
    
    assert len(parser.positions) == 0
    assert len(parser.bonds) == 0
    assert len(parser.atom_types) == 0

def test_parse_line_parser_method(newLAMMPSParser):
    parser = newLAMMPSParser
    
    parser.line_id = 0
    assert parser.parser_method == None #header of file

    parser.line_id = 23
    LINE = " -10.067500 19.932500  ylo yhi"
    parser.parse_line(LINE)
    assert parser.parser_method == None #Note: NOT parse_box_size

    parser.line_id = 77
    LINE = " Masses"
    parser.parse_line(LINE)
    assert parser.parser_method == parser.parse_mass

    parser.line_id = 42
    LINE = " Atoms # full"
    parser.parse_line(LINE)
    assert parser.parser_method == parser.parse_atom

    parser.line_id = 86
    LINE = "Bonds"
    parser.parse_line(LINE)
    assert parser.parser_method == parser.parse_bond

def test_parse(newLAMMPSParser):
    parser = newLAMMPSParser
    FPATH = os.path.join('furiousatoms','tests','test_data','CB_18','CB_18.data')
    
    parser.parse(FPATH)

    assert parser.box_size == [21.349, 19.9325, 16.639]

    assert len(parser.positions) == 18
    assert list(parser.positions[0]) == [-5.460000038146973, -2.6490953430879927e-08, 0.0]
    assert list(parser.positions[1]) == [-3.8999998569488525, -2.6490953430879927e-08, 0.0]
    assert list(parser.positions[14]) == [1.5599998235702515, 1.3509995937347412, 0.0]
    assert list(parser.positions[17]) == [5.460000038146973, -2.6490953430879927e-08, 0.0]

    assert len(parser.bonds) == 21
    assert parser.bonds[0][0] == 0
    assert parser.bonds[0][1] == 1
    assert parser.bonds[1][0] == 1
    assert parser.bonds[1][1] == 2
    assert parser.bonds[9][0] == 7
    assert parser.bonds[9][1] == 8
    assert parser.bonds[20][0] == 16
    assert parser.bonds[20][1] == 17

    assert len(parser.atom_types) == 18
    assert parser.atom_types[0] == 'C'
    assert parser.atom_types[1] == 'B'
    assert parser.atom_types[14] == 'C'
    assert parser.atom_types[17] == 'B'


if __name__ == '__main__':
    unittest.main()
