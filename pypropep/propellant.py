from cpropep._cpropep import ffi, lib

__all__ = ['Propellant']

class Propellant(dict):
    '''
    Thin shell around attrdict for propellant type.
    Purpose of a new class is so we can do instance check
    elsewhere.
    '''
    def __init__(self, *args, **kwargs):
        super(Propellant, self).__init__(*args, **kwargs)
    
    def __str__(self):
        return "Propellant: %s [%d]" % (self['name'], self['id'])
    
    def __repr__(self):
        return super(Propellant,self).__repr__()