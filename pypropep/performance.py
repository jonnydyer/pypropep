import re
from cpropep._cpropep import ffi, lib
from pypropep.equilibrium import Equilibrium

__all__ = ['RocketPerformance']

class RocketPerformance(object):
    '''
    A generic container class for cpropep case's.
    The GenericCase is equivalent to the "TP" option in cpropep where
    temperature and pressure are specified for the equillibrium calculation.
    '''
    def __init__(self, T=300., P=1.):
        super(RocketPerformance, self).__init__()
        self._equil_structs = ffi.new("equilibrium_t[3]")
        self._equil_objs = list()
        for i in xrange(3):
            e = ffi.addressof(self._equil_structs[i])
            self._equil_objs.append(Equilibrium(e))

    @property
    def equilibrated(self):
        equilibrated = False
        for i in xrange(3):
            equilibrated = equilibrated and self._equil_objs[i].equilibrated
        return equilibrated

    @property
    def properties_computed(self):
        properties_computed = False
        for i in xrange(3):
            properties_computed = properties_computed and \
                self._equil_objs[i].properties_computed
        return properties_computed

    @property
    def performance_computed(self):
        computed = False
        for i in xrange(3):
            computed = computed and self._equil_objs[i].performance_computed
        return computed

    @property
    def properties(self):
        return [p.properties for p in self._equil_structs]


    def add_propellant(self, propellant, mol):
        self._equil_objs[0].add_propellant(propellant, mol)

    def add_propellants(self, propellant_list):
        self._equil_objs[0].add_propellants(propellant_list)

    def __str__(self):
        s = "Status:\n"
        s += "\tEquillibrium Computed: {}\n".format(str(self.equilibrated))
        s += "\tProperties Computed: {}\n".format(str(self.properties_computed))
        s += "\tPerformance Computed: {}\n".format(str(self.performance_computed))
        s += "Composition:\n"
        for i in xrange(self._equil_structs[0].propellant.ncomp):
            ind = self._equil_structs[0].propellant.molecule[i]
            name = ffi.string(lib.propellant_list[ind].name)
            s += "\t{} - {:.3f} mol\n".format(name,
                                              self._equil_structs[0].propellant.coef[i])

        for i,c in enumerate([
                            '======= Chamber =======',
                            '======= Throat =======',
                            '======= Exit =======']):
            s += "{0}: \n".format(c)
            s += "\t" + re.sub(r"(\n)", r"\1\t", self._equil_objs[i].state_str)

        return s

class FrozenPerformance(RocketPerformance):
    def __init__(self, *args):
        super(FrozenPerformance, self).__init__(*args)

    def set_state(P, Pe=None, Ae_At=None):
        if (Pe is not None) and (Ae_At is not None):
            raise RuntimeError("Only one of Pe or At_Ae may be set at a time")

        if (Pe is None) and (At_Ae is None):
            raise RuntimeError("At least one of Pe or Ae_at
        elif Pe is not None:
            pass
        elif At_Ae is not None:
            pass
        else:
            raise RuntimeError("Shouldn't have gotten here")

