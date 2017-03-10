from cpropep._cpropep import ffi, lib

__all__ = ['Equillibrium']

class Equillibrium(object):
    def __init__(self):
        super(Equillibrium,self).__init__()
        self._equil = ffi.new("equilibrium_t *")
        lib.initialize_equilibrium(self._equil)
        self.reset()
    
    def reset(self):
        lib.reset_equilibrium(self._equil)
        self.propellants = []
        
    def __del__(self):
        del self._equil
    
    def add_propellant(self, propellant, mol):
        try:
            lib.add_in_propellant(self._equil, propellant['id'], mol)
            propellants.append(propellant)
        except:
            # TODO: Write tests for this case!
            TypeError("Problem adding propellant.  \
                Did you pass in a Propellant object?")
    
    def add_propellants(self, propellant_list):
        '''
        Propellant list must be a list of tuples with:
        [(Propellant_1, mol_1), ..., (Propellant_n, mol_n)]
        '''
        for p, m in propellants:
            self.add_propellant(p, m)
    
    def set_state(self, P, T=None, type='HP'):
        '''
        Set state for equillibrium calculation.  Note that type
        must be in ('TP', 'SP', 'HP').  If 'TP' temperature must
        be specified; otherwise it must not be.
        '''
        # TODO: Write tests for these cases!
        if type == 'TP':
            if T == None:
                raise ValueError("set_state: Temperature must be specified \
                    for mode 'TP'")
                    
        elif type == 'HP':
            if T != None:
                raise ValueError("set_state: Temperature must not be specified \
                    for mode 'HP'")
                    
        elif type == 'SP':
            if T != None:
                raise ValueError("set_state: Temperature must not be specified \
                    for mode 'SP'")
                    
        else:
            raise ValueError("set_state: type must be one of ('TP', 'SP', 'HP')!")
    
    def __str__(self):
        s = "cpropep composition:\n"
        for i in xrange(self._equil.propellant.ncomp):
            ind = self._equil.propellant.molecule[i]
            name = ffi.string(lib.propellant_list[ind].name)
            s += "\t%s - %0.3f mol\n" % (name, self._equil.propellant.coeff[i])
        s += "python composition list: %s" % (str(self.propellants))
    
    def __repr__(self):
        return self.__str__()
    
    