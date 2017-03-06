from cpropep._cpropep import ffi, lib

class Case(object):
    def __init__(self):
        super(Case, self).__init__()
        
        self._equil = ffi.new("equilibrium_t *")
        