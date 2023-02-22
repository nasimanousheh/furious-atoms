import pytest
import numpy as np

@pytest.fixture
def boxSize():
    return [0, 0, 0]

@pytest.fixture
def positions():
    positions = np.zeros(shape=(4, 3))
    positions[0] = [17.047001, 14.099, 3.625]
    positions[1] = [16.966999, 12.784, 4.338]
    positions[2] = [5.929, 6.358, 5.055]
    positions[3] = [12.951, 13.245, -2.112]
    return positions

@pytest.fixture
def bonds():
    return np.zeros(shape=(0, 0))

@pytest.fixture
def atomTypes():
    atom_types = np.zeros(shape=(4), dtype=str)
    atom_types[0] = 'N'
    atom_types[1] = 'C'
    atom_types[2] = 'O'
    atom_types[3] = 'H'
    return atom_types