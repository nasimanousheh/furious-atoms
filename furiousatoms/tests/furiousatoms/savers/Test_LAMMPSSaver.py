import os
import unittest
import pytest

from furiousatoms.savers.lammps_saver import LAMMPSSaver
from furiousatoms.parsers.lammps_parser import LAMMPSParser
from furiousatoms.tests.conftest import IN_FILE_PATH, NON_EXISTENT_FILE_PATH, SAVE_FILE_PATH
from saver_fixtures import *

IN_FILE_PATH = "%s.data" % (IN_FILE_PATH)
SAVE_FILE_PATH = "%s.data" % (SAVE_FILE_PATH)

@pytest.fixture
def newLAMMPSStructure():
    structure = LAMMPSParser().parse(IN_FILE_PATH)
    assert len(structure.pos) == 18
    assert len(structure.bonds) == 21

    return structure

@pytest.fixture
def newLAMMPSSaver(newLAMMPSStructure):
    structure = newLAMMPSStructure
    deleted_particles = np.zeros(len(structure.atom_types), dtype=bool)
    deleted_bonds = np.zeros(len(structure.bonds), dtype=bool)
    saver = LAMMPSSaver(deleted_particles, deleted_bonds)
    
    yield saver
    #Teardown
    if os.path.exists(SAVE_FILE_PATH):
        os.remove(SAVE_FILE_PATH)



def test_save_to_file(newLAMMPSSaver, newLAMMPSStructure):
    saver = newLAMMPSSaver
    structure = newLAMMPSStructure

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        #HEADER SECTION ####################################
        fp.readline() #header
        fp.readline() #blank
        assert fp.readline() == "18	atoms\n"
        assert fp.readline() == "21	bonds\n"
        assert fp.readline() == "32	 angles\n"
        assert fp.readline() == "46	 dihedrals\n"
        assert fp.readline() == "24	 impropers\n"
        assert fp.readline() == "2	atom types\n"
        assert fp.readline() == "1	bond types\n"
        assert fp.readline() == "1	angle types\n"
        assert fp.readline() == "1	dihedral types\n"
        fp.readline() #blank
        assert fp.readline() == "1	improper types\n"
        fp.readline() #blank
        assert fp.readline() == "-0.000000	21.349000	xlo xhi\n"
        assert fp.readline() == "-0.000000	19.932500	ylo yhi\n"
        assert fp.readline() == "-0.000000	16.639000	zlo zhi\n"
        fp.readline() #blank
        #BODY SECTION ####################################
        assert fp.readline() == "Masses\n"
        fp.readline() #blank
        assert fp.readline() == "B	10.811\n"
        assert fp.readline() == "C	12.011\n"
        fp.readline() #blank
        assert fp.readline().startswith("Atoms")
        fp.readline() #blank
        assert fp.readline() == "1  C -5.460000 -0.000000 0.000000\n"
        assert fp.readline() == "2  B -3.900000 -0.000000 0.000000\n"
        for i in range(14):
            fp.readline() #skip some atoms
        assert fp.readline() == "17  C 3.900000 -0.000000 0.000000\n"
        assert fp.readline() == "18  B 5.460000 -0.000000 0.000000\n"
        fp.readline() #blank
        assert fp.readline() == "Bonds\n"
        fp.readline() #blank
        assert fp.readline() == "1 1 1 2\n"
        assert fp.readline() == "2 1 2 3\n"
        assert fp.readline() == "3 1 2 7\n"
        for i in range(16):
            fp.readline() #skip some bonds
        assert fp.readline() == "20 1 16 17\n"
        assert fp.readline() == "21 1 17 18\n"
        fp.readline() #blank
        #UNUSED (as of 3/25/2023) DATA SECTION ####################################
        assert fp.readline() == "Angles\n"
        fp.readline() #blank
        assert fp.readline() == "1	1	1	2	3\n"
        for i in range(32):
            fp.readline()
        assert fp.readline() == "Dihedrals\n"
        fp.readline() #blank
        assert fp.readline() == "1	1	1	2	3	4\n"
        for i in range(46):
            fp.readline()
        assert fp.readline() == "Impropers\n"
        fp.readline() #blank
        assert fp.readline() == "1	1	0	2	6	1\n"
        for i in range(25):
            fp.readline()

        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()


def test_save_to_file_atom_deletion(newLAMMPSSaver, newLAMMPSStructure):
    saver = newLAMMPSSaver
    structure = newLAMMPSStructure
    saver.deleted_particles[2] = True

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        #HEADER SECTION ####################################
        fp.readline() #header
        fp.readline() #blank
        assert fp.readline() == "17	atoms\n"
        for _ in range(22):
            fp.readline()
        assert fp.readline() == "1  C -5.460000 -0.000000 0.000000\n"
        assert fp.readline() == "2  B -3.900000 -0.000000 0.000000\n"
        #atom deleted here
        assert fp.readline() == "3  B -1.560000 -1.351000 0.000000\n"
        assert fp.readline() == "4  C -0.780000 -2.701999 0.000000\n"
        for _ in range(11):
            fp.readline()
        assert fp.readline() == "16  C 3.900000 -0.000000 0.000000\n"
        assert fp.readline() == "17  B 5.460000 -0.000000 0.000000\n"
        assert len(fp.readline().strip()) == 0
        assert fp.readline() == "Bonds\n"
        fp.readline() #blank
        assert fp.readline() == "1 1 1 2\n"
        assert fp.readline() == "2 1 2 2\n"
        assert fp.readline() == "3 1 2 6\n"
        for i in range(16):
            fp.readline() #skip some bonds
        assert fp.readline() == "20 1 15 16\n"
        assert fp.readline() == "21 1 16 17\n"
        fp.readline() #blank


def test_save_to_file_bond_deletion(newLAMMPSSaver, newLAMMPSStructure):
    saver = newLAMMPSSaver
    structure = newLAMMPSStructure
    saver.deleted_bonds[8] = True #expect bond from IDs 6-to-7 to be deleted (or 7-to-8 in the file)

    saver.save_to_file(SAVE_FILE_PATH, IN_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        #HEADER SECTION ####################################
        fp.readline() #header
        fp.readline() #blank
        assert fp.readline() == "18	atoms\n"
        assert fp.readline() == "20	bonds\n"
        for _ in range(21):
            fp.readline()
        assert fp.readline() == "1  C -5.460000 -0.000000 0.000000\n"
        assert fp.readline() == "2  B -3.900000 -0.000000 0.000000\n"
        for _ in range(4):
            fp.readline()
        assert fp.readline() == "7  C -3.120000 1.351000 0.000000\n"
        assert fp.readline() == "8  B -1.560000 1.351000 0.000000\n"
        assert fp.readline() == "9  C -0.780000 -0.000000 0.000000\n"
        for _ in range(8):
            fp.readline()
        assert fp.readline() == "18  B 5.460000 -0.000000 0.000000\n"
        fp.readline() #blank
        assert fp.readline() == "Bonds\n"
        fp.readline() #blank
        assert fp.readline() == "1 1 1 2\n"
        assert fp.readline() == "2 1 2 3\n"
        assert fp.readline() == "3 1 2 7\n"
        for i in range(3):
            fp.readline() #skip some bonds
        assert fp.readline() == "7 1 5 6\n"
        assert fp.readline() == "8 1 6 11\n"
        #bond is deleted here
        assert fp.readline() == "9 1 8 9\n"
        assert fp.readline() == "10 1 8 13\n"
        for i in range(8):
            fp.readline() #skip some bonds
        assert fp.readline() == "19 1 16 17\n"
        assert fp.readline() == "20 1 17 18\n"
        assert len(fp.readline().strip()) == 0


def test_save_to_file_use_defaults(newLAMMPSSaver, newLAMMPSStructure):
    saver = newLAMMPSSaver
    structure = newLAMMPSStructure

    saver.save_to_file(SAVE_FILE_PATH, NON_EXISTENT_FILE_PATH, structure)

    with open(SAVE_FILE_PATH, "r") as fp:
        #HEADER SECTION ####################################
        fp.readline() #header
        fp.readline() #blank
        assert fp.readline() == "18	atoms\n"
        assert fp.readline() == "21	bonds\n"
        assert fp.readline() == "2	atom types\n"
        assert fp.readline() == "1	bond types\n"
        fp.readline() #blank
        assert fp.readline() == "21.349000	21.349000	xlo xhi\n"
        assert fp.readline() == "19.932500	19.932500	ylo yhi\n"
        assert fp.readline() == "16.639000	16.639000	zlo zhi\n"
        fp.readline() #blank
        #BODY SECTION ####################################
        assert fp.readline() == "Masses\n"
        fp.readline() #blank
        assert fp.readline() == "B 10.811000\n"
        assert fp.readline() == "C 12.010700\n"
        fp.readline() #blank
        assert fp.readline().startswith("Atoms")
        fp.readline() #blank
        assert fp.readline() == "1  C -5.460000 -0.000000 0.000000\n"
        assert fp.readline() == "2  B -3.900000 -0.000000 0.000000\n"
        for i in range(14):
            fp.readline() #skip some atoms
        assert fp.readline() == "17  C 3.900000 -0.000000 0.000000\n"
        assert fp.readline() == "18  B 5.460000 -0.000000 0.000000\n"
        fp.readline() #blank
        assert fp.readline() == "Bonds\n"
        fp.readline() #blank
        assert fp.readline() == "1 1 1 2\n"
        assert fp.readline() == "2 1 2 3\n"
        assert fp.readline() == "3 1 2 7\n"
        for i in range(16):
            fp.readline() #skip some bonds
        assert fp.readline() == "20 1 16 17\n"
        assert fp.readline() == "21 1 17 18\n"
        #Ensure this the end of file
        with pytest.raises(Exception):
            assert fp.readline()
            assert fp.readline()


if __name__ == '__main__':
    unittest.main()