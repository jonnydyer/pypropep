import pytest
import os

@pytest.fixture
def pypropep():
    import pypropep
    pypropep.init()
    return pypropep

def test_propellant_str(pypropep):
    assert 'Propellant' in pypropep.PROPELLANTS[pypropep.PROPELLANTS.keys()[0]].__str__()

def test_atoms_of(pypropep):
    p = pypropep.PROPELLANTS['METHANE']
    assert p.atoms_of('C') == 1
    assert p.atoms_of('H') == 4
    assert p.atoms_of('N') == 0
    assert p.atoms_of('X') == 0

