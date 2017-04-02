import re
from .cpropep._cpropep import ffi, lib
from pypropep.equilibrium import Equilibrium
from pypropep.error import RET_ERRORS

__all__ = ['RocketPerformance', 'FrozenPerformance', 'ShiftingPerformance']

Ge = 9.80665

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
        for i in range(3):
            e = ffi.addressof(self._equil_structs[i])
            self._equil_objs.append(Equilibrium(e))

    @property
    def equilibrated(self):
        equilibrated = False
        for i in range(3):
            equilibrated = equilibrated and self._equil_objs[i].equilibrated
        return equilibrated

    @property
    def properties_computed(self):
        properties_computed = False
        for i in range(3):
            properties_computed = properties_computed and \
                self._equil_objs[i].properties_computed
        return properties_computed

    @property
    def performance_computed(self):
        computed = False
        for i in range(3):
            computed = computed and self._equil_objs[i].performance_computed
        return computed

    @property
    def properties(self):
        return [p.properties for p in self._equil_structs]

    @property
    def composition(self):
        return {
            "chamber" : self._equil_objs[0].composition_sorted,
            "throat" : self._equil_objs[1].composition_sorted,
            "exit" : self._equil_objs[2].composition_sorted
        }

    @property
    def composition_condensed(self):
        return {
            "chamber" : self._equil_objs[0].composition_condensed,
            "throat" : self._equil_objs[1].composition_condensed,
            "exit" : self._equil_objs[2].composition_condensed
        }

    @property
    def performance(self):
        return self._equil_structs[2].performance


    def add_propellant(self, propellant, mol):
        '''
        Addes propellant to the pre-equilibrium mixture.
        propellant must be a Propellant instance and mol is the number of
        mols.
        '''
        self._equil_objs[0].add_propellant(propellant, mol)

    def add_propellants(self, propellant_list):
        '''
        Addes a list of propellants to the pre-equilibrium mixture.
        propellant_list is a list of (propellant, mol) tuples as are passed
        in to add_propellant().  Example:
            add_propellants([
                            (pypropep.PROPELLANTS['METHANE'], 1.0),
                            (pypropep.PROPELLANTS['OXYGEN (GAS)'], 1.0)])
        '''
        self._equil_objs[0].add_propellants(propellant_list)

    def add_propellants_by_mass(self, propellant_list):
        '''
        Addes a list of propellants to the pre-equilibrium mixture by mass.
        This is equivalent to add_propellants except that quantity is specified
        by mass not mols.
        propellant_list is a list of (propellant, mol) tuples as are passed
        in to add_propellant().  Example:
            add_propellants([
                            (pypropep.PROPELLANTS['METHANE'], 1.0),
                            (pypropep.PROPELLANTS['OXYGEN (GAS)'], 1.0)])
        '''
        propellant_list_by_mol = []
        for p, w in propellant_list:
            N = w / p.mw
            propellant_list_by_mol.append((p, N))
        self.add_propellants(propellant_list_by_mol)

    def set_state(self):
        for e in self._equil_objs:
            e._compute_product_composition()

    def __str__(self):
        s = "Status:\n"
        s += "\tEquillibrium Computed: {}\n".format(str(self.equilibrated))
        s += "\tProperties Computed: {}\n".format(str(self.properties_computed))
        s += "\tPerformance Computed: {}\n".format(str(self.performance_computed))
        s += "Composition:\n"
        for i in range(self._equil_structs[0].propellant.ncomp):
            ind = self._equil_structs[0].propellant.molecule[i]
            name = ffi.string(lib.propellant_list[ind].name).decode('utf-8')
            s += "\t{} - {:.3f} mol\n".format(name,
                                              self._equil_structs[0].propellant.coef[i])

        for i,c in enumerate([
                            '======= Chamber =======',
                            '======= Throat =======',
                            '======= Exit =======']):
            s += "{0}: \n".format(c)
            s += "\t" + re.sub(r"(\n)", r"\1\t", self._equil_objs[i].state_str)
            s += "Ae/At: {:.5f}\n".format(self._equil_structs[i].performance.ae_at)
            s += "\tA/dotm: {:.5f} m/s/atm\n".format(self._equil_structs[i].performance.a_dotm)
            s += "\tC*: {:.5f} m/s\n".format(self._equil_structs[i].performance.cstar)
            s += "\tCf: {:.5f}\n".format(self._equil_structs[i].performance.cf)
            s += "\tIvac (m/s): {:.5f}\n".format(self._equil_structs[i].performance.Ivac)
            s += "\tIsp (m/s): {:.5f}\n".format(self._equil_structs[i].performance.Isp)
            s += "\tIsp/g (s): {:.5f}\n\n".format(self._equil_structs[i].performance.Isp/Ge)

        return s

class FrozenPerformance(RocketPerformance):
    def __init__(self, *args):
        super(FrozenPerformance, self).__init__(*args)

    def set_state(self, P, Pe=None, Ae_At=None):
        if (Pe is not None) and (Ae_At is not None):
            raise RuntimeError("Only one of Pe or At_Ae may be set at a time")

        self._equil_structs[0].properties.P = P

        if (Pe is None) and (Ae_At is None):
            raise RuntimeError("At least one of Pe or Ae_At must be specified")
        elif Pe is not None:
            err = lib.frozen_performance(ffi.addressof(self._equil_structs[0]),
                                         lib.PRESSURE, Pe)
        elif Ae_At is not None:
            err = lib.frozen_performance(ffi.addressof(self._equil_structs[0]),
                                         lib.SUPERSONIC_AREA_RATIO, Ae_At)
        else:
            raise RuntimeError("Shouldn't have gotten here")

        if err < 0:
            raise RuntimeError("Frozen performance failed with {}".format(
                RET_ERRORS[err]))

        super(FrozenPerformance, self).set_state()


class ShiftingPerformance(RocketPerformance):
    def __init__(self, *args):
        super(ShiftingPerformance, self).__init__(*args)

    def set_state(self, P, Pe=None, Ae_At=None):
        if (Pe is not None) and (Ae_At is not None):
            raise RuntimeError("Only one of Pe or At_Ae may be set at a time")

        self._equil_structs[0].properties.P = P

        if (Pe is None) and (Ae_At is None):
            raise RuntimeError("At least one of Pe or Ae_At must be specified")
        elif Pe is not None:
            err = lib.shifting_performance(ffi.addressof(self._equil_structs[0]),
                                         lib.PRESSURE, Pe)
        elif Ae_At is not None:
            err = lib.shifting_performance(ffi.addressof(self._equil_structs[0]),
                                         lib.SUPERSONIC_AREA_RATIO, Ae_At)
        else:
            raise RuntimeError("Shouldn't have gotten here")

        if err < 0:
            raise RuntimeError("Frozen performance failed with {}".format(
                RET_ERRORS[err]))

        super(ShiftingPerformance, self).set_state()
