import pytest
import os

@pytest.fixture
def pypropep():
    import pypropep
    pypropep.init()
    return pypropep


def test_propellants_dict(pypropep):
    assert len(pypropep.PROPELLANTS) == 1030


def test_thermo_dict(pypropep):
    assert len(pypropep.SPECIES) == 1921


def test_prop_file_override(pypropep):
    prop_file = os.path.dirname(pypropep.__file__) + '/data/propellant.dat'
    pypropep.init(propellant_file=prop_file)


def test_thermo_file_override(pypropep):
    therm_file = os.path.dirname(pypropep.__file__) + '/data/thermo.dat'
    pypropep.init(thermo_file=therm_file)


def test_find_propellant(pypropep):
    assert len(pypropep.find_propellant('oxygen')) > 1
    assert len(pypropep.find_propellant('OXYGEN')) > 1
