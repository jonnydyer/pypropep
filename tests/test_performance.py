import pytest


@pytest.fixture
def pypropep():
    import pypropep
    pypropep.init()
    return pypropep


def test_simp_equilibrium(pypropep):
    p = pypropep.RocketPerformance()
    return p.__str__()
