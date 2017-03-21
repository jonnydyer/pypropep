from cpropep._cpropep import ffi, lib
from equilibrum import Equilibrium

__all__ = ['GenericCase']

class GenericCase(object):
    '''
    A generic container class for cpropep case's.
    The GenericCase is equivalent to the "TP" option in cpropep where
    temperature and pressure are specified for the equillibrium calculation.
    '''
    def __init__(self, T=300., P=1.):
        super(GenericCase, self).__init__()
        self._equils = ffi.new("equilibrium_t[3]")
        self._equil_objs = list()
        for i in xrange(3):
            self._equil_objs.append(
                Equilibrium(self._equils[i])
            )

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

    def __str__(self):
        s = "Status:\n"
        s += "\tEquillibrium Computed: {}\n".format(str(self.equilibrated))
        s += "\tProperties Computed: {}\n".format(str(self.properties_computed))
        s += "\tPerformance Computed: {}\n".format(str(self.performance_computed))
        s += "Composition:\n"
        for i in xrange(self._equil.propellant.ncomp):
            ind = self._equil.propellant.molecule[i]
            name = ffi.string(lib.propellant_list[ind].name)
            s += "\t{} - {:.3f} mol\n".format(name,
                                              self._equil.propellant.coef[i])

        s += "State:\n"
        s += "\tPressure: {:.3f} atm \n".format(self.properties.P)
        s += "\tTemperature: {:.1f} K \n".format(self.properties.T)
        s += "\tEnthalpy: {:.3f} kJ/kg \n".format(self.properties.H)
        s += "\tInt. Energy: {:.3f} kJ/kg \n".format(self.properties.U)
        s += "\tGibbs Free Energy: {:.3f} kJ/kg \n".format(self.properties.G)
        s += "\tEntropy: {:.3f} kJ/kg-K \n".format(self.properties.S)
        s += "\tMolar Mass: {:.3f} g/mol \n".format(self.properties.M)
        s += "\tdV_P: {:.3f}\n".format(self.properties.dV_P)
        s += "\tdV_T: {:.3f}\n".format(self.properties.dV_T)
        s += "\tCp: {:.3f} kJ/kg-K\n".format(self.properties.Cp)
        s += "\tCv: {:.3f} kJ/kg-K\n".format(self.properties.Cv)
        s += "\tgamma: {:.3f}\n".format(self.properties.Isex)
        s += "\tSound Speed: {:.1f} m/s\n".format(self.properties.Vson)
