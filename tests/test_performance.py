import pytest


@pytest.fixture
def pypropep():
    import pypropep
    pypropep.init()
    return pypropep


def test_simp_performance(pypropep):
    p = pypropep.RocketPerformance()
    return p.__str__()

def test_shifting_performance(pypropep):
    p = pypropep.ShiftingPerformance()
    lh2 = pypropep.PROPELLANTS['HYDROGEN (CRYOGENIC)']
    lox = pypropep.PROPELLANTS['OXYGEN (LIQUID)']
    # The below is from Sutton pg. 181
    OF = 5.551
    N_lh2 = 1.0 / lh2.mw
    N_lox = OF / lox.mw
    p.add_propellants([(lh2, N_lh2), (lox, N_lox)])
    p.set_state(P=53.317*0.986923, Ae_At=25.)
    assert p.performance.Isp == pytest.approx(4124, 1e-2)
    assert p.performance.Ivac == pytest.approx(4348, 1e-2)
    assert p.performance.cstar == pytest.approx(2332.1, 1e-3)
    assert p.properties[0].T == pytest.approx(3389, 1e-2)
    assert p.properties[2].T == pytest.approx(1468, 1e-2)


def test_frozen_performance(pypropep):
    p = pypropep.FrozenPerformance()
    lh2 = pypropep.PROPELLANTS['HYDROGEN (CRYOGENIC)']
    lox = pypropep.PROPELLANTS['OXYGEN (LIQUID)']
    # The below is from Sutton pg. 181
    OF = 5.551
    N_lh2 = 1.0 / lh2.mw
    N_lox = OF / lox.mw
    p.add_propellants([(lh2, N_lh2), (lox, N_lox)])
    p.set_state(P=53.317*0.986923, Ae_At=25.)
    #assert p.performance.cstar == pytest.approx(2332.1, 1e-2)
    assert p.properties[0].T == pytest.approx(3389, 1e-2)
