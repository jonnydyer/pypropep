import pytest
import os

@pytest.fixture
def pypropep():
    import pypropep
    pypropep.init()
    return pypropep
    
def test_propellant_str(pypropep):
    assert 'Propellant' in pypropep.PROPELLANTS[pypropep.PROPELLANTS.keys()[0]].__str__()
