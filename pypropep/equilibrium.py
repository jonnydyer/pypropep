import re
import numpy as np
import operator
from .cpropep._cpropep import ffi, lib
from .error import RET_ERRORS

__all__ = ['Equilibrium']

class Equilibrium(object):
    def __init__(self, equilibrium_t_ptr=None):
        super(Equilibrium, self).__init__()
        if equilibrium_t_ptr is not None:
            self._equil = equilibrium_t_ptr
        else:
            self._equil = ffi.new("equilibrium_t *")
        lib.initialize_equilibrium(self._equil)
        self.reset()

    def reset(self):
        lib.reset_equilibrium(self._equil)
        self._composition = dict()
        self._composition_condensed = dict()
        self.propellants = []

    def __del__(self):
        del self._equil

    def add_propellant(self, propellant, mol):
        '''
        Addes propellant to the pre-equilibrium mixture.
        propellant must be a Propellant instance and mol is the number of
        mols.
        '''
        try:
            lib.add_in_propellant(self._equil, propellant['id'], mol)
            self.propellants.append(propellant)
        except:
            # TODO: Write tests for this case!
            TypeError("Problem adding propellant.  \
                Did you pass in a Propellant object?")

    def add_propellants(self, propellant_list):
        '''
        Addes a list of propellants to the pre-equilibrium mixture.
        propellant_list is a list of (propellant, mol) tuples as are passed
        in to add_propellant().  Example:
            add_propellants([
                            (pypropep.PROPELLANTS['METHANE'], 1.0),
                            (pypropep.PROPELLANTS['OXYGEN (GAS)'], 1.0)])
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
    def equilibrated(self):
        return bool(self._equil.equilibrium_ok)

    @property
    def properties(self):
        return self._equil.properties

    @property
    def composition(self):
        if self.equilibrated is False:
            return None
        return self._composition

    @property
    def composition_sorted(self):
        if self.equilibrated is False:
            return None
        return sorted(list(self._composition.items()), key=operator.itemgetter(1),
                      reverse=True)

    @property
    def composition_condensed(self):
        if self.equilibrated is False:
            return None
        return self._composition_condensed


    def _compute_product_composition(self):
        if self.equilibrated is False:
            raise RuntimeError("Can't compute product composition until \
                               equilibrum is computed")
        mol_g = self._equil.itn.n
        for i in range(self._equil.product.n[lib.CONDENSED]):
            mol_g += self._equil.product.coef[lib.CONDENSED][i]

        for i in range(self._equil.product.n[lib.GAS]):
            ind = self._equil.product.species[lib.GAS][i]
            name = ffi.string(lib.thermo_list[ind].name).decode('utf-8')
            self._composition[name] = \
                self._equil.product.coef[lib.GAS][i] / mol_g

        for i in range(self._equil.product.n[lib.CONDENSED]):
            ind = self._equil.product.species[lib.CONDENSED][i]
            name = ffi.string(lib.thermo_list[ind].name).decode('utf-8')
            self._composition_condensed[name] = \
                self._equil.product.coef[lib.CONDENSED][i] / mol_g

    def set_state(self, P, T=None, type='HP'):
        '''
        Set state for equillibrium calculation.  Note that type
        must be in ('TP', 'SP', 'HP').  If 'TP' temperature must
        be specified; otherwise it must not be.
        '''
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

        self._compute_product_composition()

    @property
    def state_str(self):
        s = "Pressure: {:.3f} atm \n".format(self.properties.P)
        s += "Temperature: {:.1f} K \n".format(self.properties.T)
        s += "Enthalpy: {:.3f} kJ/kg \n".format(self.properties.H)
        s += "Int. Energy: {:.3f} kJ/kg \n".format(self.properties.U)
        s += "Gibbs Free Energy: {:.3f} kJ/kg \n".format(self.properties.G)
        s += "Entropy: {:.3f} kJ/kg-K \n".format(self.properties.S)
        s += "Molar Mass: {:.3f} g/mol \n".format(self.properties.M)
        s += "dV_P: {:.3f}\n".format(self.properties.dV_P)
        s += "dV_T: {:.3f}\n".format(self.properties.dV_T)
        s += "Cp: {:.3f} kJ/kg-K\n".format(self.properties.Cp)
        s += "Cv: {:.3f} kJ/kg-K\n".format(self.properties.Cv)
        s += "gamma: {:.3f}\n".format(self.properties.Isex)
        s += "Sound Speed: {:.1f} m/s\n\n".format(self.properties.Vson)
        return s

    def __str__(self):
        s = "Status:\n"
        s += "\tEquillibrium Computed: {}\n".format(str(self.equilibrated))
        s += "\tProperties Computed: {}\n".format(str(self.properties_computed))
        s += "\tPerformance Computed: {}\n".format(str(self.performance_computed))
        s += "Composition:\n"
        for i in range(self._equil.propellant.ncomp):
            ind = self._equil.propellant.molecule[i]
            name = ffi.string(lib.propellant_list[ind].name).decode('utf-8')
            s += "\t{} - {:.3f} mol\n".format(name,
                                              self._equil.propellant.coef[i])
        s += "State:\n"
        s += "\t" + re.sub(r"(\n)", r"\1\t", self.state_str)
        return s

    def __repr__(self):
        return self.__str__()
