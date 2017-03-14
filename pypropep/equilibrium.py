import numpy as np
import operator
from cpropep._cpropep import ffi, lib
from error import RET_ERRORS

__all__ = ['Equilibrium']

class Equilibrium(object):
    def __init__(self):
        super(Equilibrium, self).__init__()
        self._equil = ffi.new("equilibrium_t *")
        lib.initialize_equilibrium(self._equil)
        self.reset()

    def reset(self):
        lib.reset_equilibrium(self._equil)
        self._composition = dict()
        self.propellants = []
        self.equilibrated = False

    def __del__(self):
        del self._equil

    def add_propellant(self, propellant, mol):
        try:
            lib.add_in_propellant(self._equil, propellant['id'], mol)
            self.propellants.append(propellant)
        except:
            # TODO: Write tests for this case!
            TypeError("Problem adding propellant.  \
                Did you pass in a Propellant object?")

    def add_propellants(self, propellant_list):
        '''
        Propellant list must be a list of tuples with:
        [(Propellant_1, mol_1), ..., (Propellant_n, mol_n)]
        '''
        for p, m in propellant_list:
            self.add_propellant(p, m)

    @property
    def properties_computed(self):
        return bool(self._equil.properties_ok)

    @property
    def performance_computed(self):
        return bool(self._equil.performance_ok)

    @property
    def properties(self):
        return self._equil.properties

    @property
    def composition(self):
        if self.equilibrated is False:
            return None
        return self._composition

    @property
    def composotion_sorted(self):
        if self.equilibrated is False:
            return None
        return sorted(self._composition.items(), key=operator.itemgetter(1),
                      reverse=True)


    def _compute_product_composition(self):
        if self.equilibrated is False:
            raise RuntimeError("Can't compute product composition until \
                               equilibrum is computed")
        mol_g = self._equil.itn.n
        for i in xrange(self._equil.product.n[lib.CONDENSED]):
            mol_g += self._equil.product.coef[lib.CONDENSED][i]

        for i in xrange(self._equil.product.n[lib.GAS]):
            ind = self._equil.product.species[lib.GAS][i]
            name = ffi.string(lib.thermo_list[ind].name)
            self._composition[name] = self._equil.product.coef[lib.GAS][i] / \
                mol_g

    def set_state(self, P, T=None, type='HP'):
        '''
        Set state for equillibrium calculation.  Note that type
        must be in ('TP', 'SP', 'HP').  If 'TP' temperature must
        be specified; otherwise it must not be.
        '''
        # TODO: Write tests for these cases!
        self._equil.product.n[lib.CONDENSED] = 0
        if type == 'TP':
            if T is None:
                raise ValueError("set_state: Temperature must be specified \
                    for mode 'TP'")
            self.properties.T = T
            self.properties.P = P
            eq_type = lib.TP

        elif type == 'HP' or type == 'SP':
            if T is not None:
                raise ValueError("set_state: Temperature must not be specified \
                    for mode 'HP' or 'SP'")
            self.properties.P = P
            if type == 'HP':
                eq_type = lib.HP
            else:
                eq_type = lib.SP

        else:
            raise ValueError("set_state: type must be one of ('TP', 'SP', 'HP')!")

        err = lib.equilibrium(self._equil, eq_type)
        if err != 0:
            raise RuntimeError("Equilibrium failed with error: '{}'".format(
                                RET_ERRORS[err]))

        self.equilibrated = True
        self._compute_product_composition()

    def __str__(self):
        s = "Status:\n"
        s += "\tEquillibrium Computed: {}\n".format(str(self.equilibrated))
        s += "\tProperties Computed: {}\n".format(str(self.properties_computed))
        s += "\tPerformance Computed: {}\n".format(str(self.performance_computed))
        s += "Composition:\n"
        for i in xrange(self._equil.propellant.ncomp):
            ind = self._equil.propellant.molecule[i]
            name = ffi.string(lib.propellant_list[ind].name)
            s += "\t{} - {:.3f} mol\n".format(name, self._equil.propellant.coef[i])

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

        return s

    def __repr__(self):
        return self.__str__()

