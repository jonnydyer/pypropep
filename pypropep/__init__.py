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
    members.  Motivated by :
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

species = AttrDict()
propellants = AttrDict()

# Build species dict
for i in xrange(lib.num_thermo):
    s = lib.thermo_list[i]
    species[ffi.string(s.name)] = AttrDict(convert_to_python(s))
    
# Build propellants dict
for i in xrange(lib.num_propellant):
    p = lib.propellant_list[i]
    propellants[ffi.string(p.name)] = AttrDict(convert_to_python(p))
