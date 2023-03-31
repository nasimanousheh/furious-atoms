import os
import unittest
import pytest

from furiousatoms.savers.xyz_saver import XYZSaver
from furiousatoms.parsers.xyz_parser import XYZParser
from furiousatoms.tests.conftest import IN_FILE_PATH, NON_EXISTENT_FILE_PATH, SAVE_FILE_PATH
from saver_fixtures import *

IN_FILE_PATH = "%s.xyz" % (IN_FILE_PATH)
SAVE_FILE_PATH = "%s.xyz" % (SAVE_FILE_PATH)

@pytest.fixture
def newXYZStructure():
    structure = XYZParser().parse(IN_FILE_PATH)
    assert len(structure.pos) == 18

    return structure

@pytest.fixture
def newXYZSaver(newXYZStructure):
    structure = newXYZStructure
    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = XYZSaver(deleted_particles, deleted_bonds)
    
    yield saver
    #Teardown
    if os.path.exists(SAVE_FILE_PATH):
        os.remove(SAVE_FILE_PATH)



def test_save_to_file(newXYZSaver, newXYZStructure):
    saver = newXYZSaver
    structure = newXYZStructure

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        assert fp.readline() == "18\n"
        assert fp.readline() == "XYZ file intended for test cases\n"
        assert fp.readline() == " C -5.460000 -0.000000 0.000000\n"
        assert fp.readline() == " B -3.900000 -0.000000 0.000000\n"
        assert fp.readline() == " C -3.120000 -1.351000 0.000000\n"
        for i in range(13): #Skip atoms
            fp.readline()
        assert fp.readline() == " C 3.900000 -0.000000 0.000000\n"
        assert fp.readline() == " B 5.460000 -0.000000 0.000000\n"
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()


def test_save_to_file_atom_deletion(newXYZSaver, newXYZStructure):
    saver = newXYZSaver
    structure = newXYZStructure
    saver.deleted_particles[0] = True
    saver.deleted_particles[3] = True
    saver.deleted_particles[17] = True

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        assert fp.readline() == "18\n"
        assert fp.readline() == "XYZ file intended for test cases\n"
        assert fp.readline() == " B -3.900000 -0.000000 0.000000\n"
        assert fp.readline() == " C -3.120000 -1.351000 0.000000\n"
        assert fp.readline() == " C -0.780000 -2.702000 0.000000\n"
        assert fp.readline() == " B 0.780000 -2.702000 0.000000\n"
        for i in range(10): #Skip atoms
            fp.readline()
        assert fp.readline() == " C 3.900000 -0.000000 0.000000\n"
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()


def test_save_to_file_use_defaults(newXYZSaver, newXYZStructure):
    saver = newXYZSaver
    structure = newXYZStructure

    saver.save_to_file(SAVE_FILE_PATH, NON_EXISTENT_FILE_PATH, structure)
    from furiousatoms.savers.xyz_saver import DEFAULT_HEADER

    with open(SAVE_FILE_PATH, "r") as fp:
        assert fp.readline() == "18\n"
        assert fp.readline().strip() == DEFAULT_HEADER.strip() #line under test
        for i in range(17): #Skip atoms
            fp.readline()
        assert fp.readline() == " B 5.460000 -0.000000 0.000000\n" #last atom
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()


if __name__ == '__main__':
    unittest.main()