import os
from attrdict import AttrDict
from cpropep._cpropep import ffi, lib

FILE_PATH = os.path.abspath(__file__)
THERMO_FILE = os.path.dirname(FILE_PATH) + '/data/thermo.dat'
PROPELLANT_FILE = os.path.dirname(FILE_PATH) + '/data/propellant.dat'

def __convert_struct_field( s, fields ):
    for field,fieldtype in fields:
        if fieldtype.type.kind == 'primitive':
            yield (field,getattr( s, field ))
        else:
            yield (field, convert_to_python( getattr( s, field ) ))

def convert_to_python(s):
    '''
    Given a cdata struct, returns a dict with all of the parsable
    members.  Borrowed with mods from :
    http://stackoverflow.com/questions/20444546/python-cffi-convert-structure-to-dictionary
    '''
    if isinstance(s, ffi.CData) == False:
        return s
    
    type=ffi.typeof(s)
    if type.kind == 'struct':
        return dict(__convert_struct_field( s, type.fields ) )
    elif type.kind == 'array':
        if type.item.kind == 'primitive':
            if type.item.cname == 'char':
                return ffi.string(s)
            else:
                return [ s[i] for i in range(type.length) ]
        else:
            return [ convert_to_python(s[i]) for i in range(type.length) ]
    elif type.kind == 'primitive':
        return int(s)
        
class Propellant(AttrDict):
    '''
    Thin shell around attrdict for propellant type.
    Purpose of a new class is so we can do instance check
    elsewhere.
    '''
    def __init__(self, *args, **kwargs):
        super(Propellant, self).__init__(*args, **kwargs)
    
    def __str__(self):
        print super(Propellant, self)
        lib.print_propellant_info(self.id)

class Equillibrium(object):
    def __init__(self):
        super(Equillibrium,self).__init__()
        self._equil = ffi.new("equilibrium_t *")
        lib.initialize_equillibrium(self._equil)
        self.reset()
    
    def reset(self):
        lib.reset_equillibrium(self._equil)
        self.propellants = []
        
    def __del__(self):
        del self._equil
    
    def add_propellant(self, propellant, mol):
        try:
            lib.add_in_propellant(self._equil, propellant['id'], mol)
        except:
            TypeError("Problem adding propellant.  Did you pass in a Propellant object?")
    
    def add_propellants(self, propellant_list):
        '''
        Propellant list must be a list of tuples with:
        [(Propellant_1, mol_1), ..., (Propellant_n, mol_n)]
        '''
        for p, m in propellants:
            self.add_propellant(p, m)

class GenericCase(object):
    '''
    A generic container class for cpropep case's.
    The GenericCase is equivalent to the "TP" option in cpropep where 
    temperature and pressure are specified for the equillibrium calculation.
    '''
    def __init__(self, T=300., P=1.):
        super(Case, self).__init__()
        self.equil = Equillibrium()


r = lib.load_thermo(THERMO_FILE)
if r > 0:
    print 'Loaded %d thermo species' % (r)
else:
    print 'Failed to load thermo file %s' % (THERMO_FILE)

r = lib.load_propellant(PROPELLANT_FILE)

if r > 0:
    print 'Loaded %d propellants' % (r)
else:
    print 'Failed to load propellant file %s' % (PROPELLANT_FILE)

species = dict()
propellants = dict()

# Build species dict
for i in xrange(lib.num_thermo):
    s = lib.thermo_list[i]
    name = ffi.string(s.name)
    species[name] = AttrDict(convert_to_python(s))
    species[name]['id'] = i
    
# Build propellants dict
for i in xrange(lib.num_propellant):
    p = lib.propellant_list[i]
    name = ffi.string(p.name)
    propellants[name] = Propellant(convert_to_python(p))
    propellants[name]['id'] = i
