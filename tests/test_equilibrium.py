import pytest


@pytest.fixture
def pypropep():
    import pypropep
    pypropep.init()
    return pypropep


def test_simp_equilibrium(pypropep):
    e = pypropep.Equilibrium()
    return e.__str__()


def test_add_propellant(pypropep):
    e = pypropep.Equilibrium()
    o2 = pypropep.PROPELLANTS['OXYGEN (GAS)']
    ch4 = pypropep.PROPELLANTS['METHANE']
    e.add_propellants([(o2, 1.), (ch4, 1.)])
    print e

def test_equil_modes(pypropep):
    e = pypropep.Equilibrium()
    with pytest.raises(ValueError):
        e.set_state(P=1., type='TP')

    with pytest.raises(ValueError):
        e.set_state(P=1., T=300., type='HP')

    with pytest.raises(ValueError):
        e.set_state(P=1., T=300., type='SP')

def test_simp_equil(pypropep):
    e = pypropep.Equilibrium()
    o2 = pypropep.PROPELLANTS['OXYGEN (GAS)']
    ch4 = pypropep.PROPELLANTS['METHANE']
    e.add_propellants([(o2, 1.), (ch4, 1.)])
    e.set_state(P=1.0, T=3000., type='TP')
    assert e.equilibrated is True
    assert e.properties_computed is True

def test_TP_composition(pypropep):
    e = pypropep.Equilibrium()
    n2 = pypropep.PROPELLANTS['NITROGEN (GASEOUS)']
    e.add_propellant(n2, 1.0)
    with pytest.raises(RuntimeError):
        print e._compute_product_composition()

    # Compostion should be None prior to equilibration
    assert e.composition is None

    # equilibrate
    e.set_state(P=1.0, T=273., type='TP')
    assert e.equilibrated is True
    assert e.properties_computed is True
    assert 'N2' in e.composition
    for k,v in e.composition.items():
        if k == 'N2':
            assert v == pytest.approx(1.0, 1e-6)
        else:
            assert v == pytest.approx(0.0, 1e-6)

def test_HP_equil(pypropep):
    e = pypropep.Equilibrium()
    o2 = pypropep.PROPELLANTS['OXYGEN (GAS)']
    ch4 = pypropep.PROPELLANTS['METHANE']
    e.add_propellants([(o2, 1.), (ch4, 1.)])
    e.set_state(P=1.0, type='HP')
    assert e.equilibrated is True
    assert e.properties_computed is True
    assert e.properties.T > 300.

# def test_SP_equil(pypropep):
#     e = pypropep.Equilibrium()
#     o2 = pypropep.PROPELLANTS['OXYGEN (GAS)']
#     ch4 = pypropep.PROPELLANTS['METHANE']
#     e.add_propellants([(o2, 1.), (ch4, 1.)])
#     e.set_state(P=1.0, type='HP')
#     T = e.properties.T
#     e.set_state(P=1, type='SP')
#     assert e.equilibrated is True
#     assert e.properties_computed is True
#     assert e.properties.T < 0.5*T
