import os
import unittest
import pytest

from furiousatoms.savers.pdb_saver import PDBSaver
from furiousatoms.parsers.pdb_parser import PDBParser
from furiousatoms.tests.conftest import IN_FILE_PATH, NON_EXISTENT_FILE_PATH, SAVE_FILE_PATH
from saver_fixtures import *

IN_FILE_PATH = "%s.pdb" % (IN_FILE_PATH)
SAVE_FILE_PATH = "%s.pdb" % (SAVE_FILE_PATH)

@pytest.fixture
def newPDBStructure():
    structure = PDBParser().parse(IN_FILE_PATH)
    assert len(structure.pos) == 18
    assert len(structure.bonds) == 21

    return structure

@pytest.fixture
def newPDBSaver(newPDBStructure):
    structure = newPDBStructure
    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = PDBSaver(deleted_particles, deleted_bonds)
    
    yield saver
    #Teardown
    if os.path.exists(SAVE_FILE_PATH):
        os.remove(SAVE_FILE_PATH)



def test_save_to_file(newPDBSaver, newPDBStructure):
    saver = newPDBSaver
    structure = newPDBStructure
    structure.box_size = [45.7, 60.6, 70.5]

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        assert fp.readline() == "REMARK Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom \n"
        assert fp.readline() == "CRYST1   45.700   60.600   70.500  90.00  90.00  90.00 P 1           1\n"
        assert fp.readline() == "ATOM      1  C   MOL X   1      -5.460  -0.000   0.000  1.00  0.00           C  \n"
        assert fp.readline() == "ATOM      2  B   MOL X   1      -3.900  -0.000   0.000  1.00  0.00           B  \n"
        for i in range(14): #Skip atoms lines 3-16
            fp.readline()
        assert fp.readline() == "ATOM     17  C   MOL X   1       3.900  -0.000   0.000  1.00  0.00           C  \n"
        assert fp.readline() == "ATOM     18  B   MOL X   1       5.460  -0.000   0.000  1.00  0.00           B  \n"
        assert fp.readline() == "CONECT    1    2\n"
        assert fp.readline() == "CONECT    2    1    3    7\n"
        for i in range(14): #Skip some bonds
            fp.readline()
        assert fp.readline() == "CONECT   17   12   16   18\n"
        assert fp.readline() == "CONECT   18   17\n"
        assert fp.readline().strip() == "END"
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()

    
def test_save_to_file_atom_deletion(newPDBSaver, newPDBStructure):
    saver = newPDBSaver
    structure = newPDBStructure
    saver.deleted_particles[5] = True

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        for i in range(5):
            fp.readline()
        assert fp.readline() == "ATOM      4  B   MOL X   1      -1.560  -1.351   0.000  1.00  0.00           B  \n"
        assert fp.readline() == "ATOM      5  C   MOL X   1      -0.780  -2.702   0.000  1.00  0.00           C  \n"
        assert fp.readline() == "ATOM      6  C   MOL X   1      -3.120   1.351   0.000  1.00  0.00           C  \n"
        for i in range(10):
            fp.readline()
        assert fp.readline() == "ATOM     17  B   MOL X   1       5.460  -0.000   0.000  1.00  0.00           B  \n"
        assert fp.readline() == "CONECT    1    2\n"
        assert fp.readline() == "CONECT    2    1    3    6\n"
        assert fp.readline() == "CONECT    3    2    4\n"
        for i in range(12):
            fp.readline()
        assert fp.readline() == "CONECT   16   11   15   17\n"
        assert fp.readline() == "CONECT   17   16\n"
        assert fp.readline().strip() == "END"


def test_save_to_file_bond_deletion(newPDBSaver, newPDBStructure):
    saver = newPDBSaver
    structure = newPDBStructure
    saver.deleted_bonds[1] = True

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        for i in range(19):
            fp.readline()
        assert fp.readline() == "ATOM     18  B   MOL X   1       5.460  -0.000   0.000  1.00  0.00           B  \n"
        assert fp.readline() == "CONECT    1    2\n"
        assert fp.readline() == "CONECT    2    1    7\n" #bond connecting 2-3 was deleted
        assert fp.readline() == "CONECT    3    4\n"
        for i in range(13):
            fp.readline()
        assert fp.readline() == "CONECT   17   12   16   18\n"
        assert fp.readline() == "CONECT   18   17\n"
        assert fp.readline().strip() == "END"


def test_save_to_file_use_defaults(newPDBSaver, newPDBStructure):
    saver = newPDBSaver
    structure = newPDBStructure

    saver.save_to_file(SAVE_FILE_PATH, NON_EXISTENT_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        assert fp.readline() == "REMARK Created by FURIOUS ATOMS. By default the box dimensions are 90x90x90 cubic angstrom\n"
        assert fp.readline() == "CRYST1   90.000   90.000   90.000   0.00   0.00   0.00 P             1\n"
        assert fp.readline() == "ATOM      1  C    MOLX   11     -5.460  -0.000   0.000  1.00  0.00           C  \n"
        assert fp.readline() == "ATOM      2  B    MOLX   11     -3.900  -0.000   0.000  1.00  0.00           B  \n"
        for i in range(14): #Skip atoms lines 3-16
            fp.readline()
        assert fp.readline() == "ATOM     17  C    MOLX   11      3.900  -0.000   0.000  1.00  0.00           C  \n"
        assert fp.readline() == "ATOM     18  B    MOLX   11      5.460  -0.000   0.000  1.00  0.00           B  \n"
        assert fp.readline() == "CONECT    1    2\n"
        assert fp.readline() == "CONECT    2    1    3    7\n"
        for i in range(14): #Skip some bonds
            fp.readline()
        assert fp.readline() == "CONECT   17   12   16   18\n"
        assert fp.readline() == "CONECT   18   17\n"
        assert fp.readline().strip() == "END"
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()


if __name__ == '__main__':
    unittest.main()