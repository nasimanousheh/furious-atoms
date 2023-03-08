import os
import unittest
import pytest

from furiousatoms.savers.xyz_saver import XYZSaver, HEADER
from furiousatoms.tests.conftest import OUT_FILE_NAME
from saver_fixtures import *

OUT_FILE_PATH = "%s.xyz"%OUT_FILE_NAME


@pytest.fixture
def newXYZSaver():
    saver = XYZSaver()
    yield saver
    #Teardown
    if os.path.exists(OUT_FILE_PATH):
        os.remove(OUT_FILE_PATH)


def test_write_all_lines(newXYZSaver, box_size, positions, bonds, atom_types):
    saver = newXYZSaver
    
    saver.write_all_lines(OUT_FILE_PATH, box_size, positions, bonds, atom_types)

    with open(OUT_FILE_PATH, "r") as fp:
        assert fp.readline() == "4\n"
        assert fp.readline() == HEADER
        assert fp.readline() == " N 17.047001 14.099000 3.625000\n"
        assert fp.readline() == " C 16.966999 12.784000 4.338000\n"
        assert fp.readline() == " O 5.929000 6.358000 5.055000\n"
        assert fp.readline() == " H 12.951000 13.245000 -2.112000\n"
        assert fp.readline() == ''


if __name__ == '__main__':
    unittest.main()