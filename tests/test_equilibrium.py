import pytest


@pytest.fixture
def pypropep():
    import pypropep
    pypropep.init()
    return pypropep


def test_equillibrium(pypropep):
    e = pypropep.Equillibrium()
    return e.__str__()


def test_add_propellant(pypropep):
    e = pypropep.Equillibrium()
    o2 = pypropep.PROPELLANTS['OXYGEN (GAS)']
    ch4 = pypropep.PROPELLANTS['METHANE']
    e.add_propellants([(o2, 1.), (ch4, 1.)])
    print e
