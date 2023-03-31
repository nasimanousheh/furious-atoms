import os
import unittest
import pytest

from furiousatoms.savers.gromacs_saver import GROMACSSaver
from furiousatoms.parsers.gromacs_parser import GROMACSParser
from furiousatoms.tests.conftest import IN_FILE_PATH, NON_EXISTENT_FILE_PATH, SAVE_FILE_PATH
from saver_fixtures import *

IN_FILE_PATH = "%s.gro" % (IN_FILE_PATH)
SAVE_FILE_PATH = "%s.gro" % (SAVE_FILE_PATH)

@pytest.fixture
def newGROStructure():
    structure = GROMACSParser().parse(IN_FILE_PATH)
    assert len(structure.pos) == 18

    return structure

@pytest.fixture
def newGROSaver(newGROStructure):
    structure = newGROStructure
    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = GROMACSSaver(deleted_particles, deleted_bonds)
    
    yield saver
    #Teardown
    if os.path.exists(SAVE_FILE_PATH):
        os.remove(SAVE_FILE_PATH)



def test_save_to_file(newGROSaver, newGROStructure):
    saver = newGROSaver
    structure = newGROStructure

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        assert fp.readline() == "GROMACS file intended for test cases\n"
        assert fp.readline() == "18\n"
        assert fp.readline() == "    1MOL      C    1  -0.546  -0.000   0.000\n"
        assert fp.readline() == "    1MOL      B    2  -0.390  -0.000   0.000\n"
        for i in range(14): #Skip atoms lines 3-16
            fp.readline()
        assert fp.readline() == "    1MOL      C   17   0.390  -0.000   0.000\n"
        assert fp.readline() == "    1MOL      B   18   0.546  -0.000   0.000\n"
        assert fp.readline().strip() == "9.00000 9.00000 9.00000"
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()

    
def test_save_to_file_atom_deletion(newGROSaver, newGROStructure):
    saver = newGROSaver
    structure = newGROStructure
    saver.deleted_particles[0] = True
    saver.deleted_particles[1] = True
    saver.deleted_particles[3] = True

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        assert fp.readline() == "GROMACS file intended for test cases\n"
        assert fp.readline() == "15\n"
        assert fp.readline() == "    1MOL      C    1  -0.312  -0.135   0.000\n"
        assert fp.readline() == "    1MOL      C    2  -0.078  -0.270   0.000\n"
        assert fp.readline() == "    1MOL      B    3   0.078  -0.270   0.000\n"
        for i in range(10):
            fp.readline()
        assert fp.readline() == "    1MOL      C   14   0.390  -0.000   0.000\n"
        assert fp.readline() == "    1MOL      B   15   0.546  -0.000   0.000\n"
        assert fp.readline().strip() == "9.00000 9.00000 9.00000"
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()


if __name__ == '__main__':
    unittest.main()