import os
from attrdict import AttrDict
from .cpropep._cpropep import ffi, lib

from pypropep.propellant import Propellant
from pypropep.equilibrium import Equilibrium
from pypropep.performance import RocketPerformance, FrozenPerformance, \
                                 ShiftingPerformance

__all__ = ['Propellant', 'Equilibrium', 'RocketPerformance',
           'FrozenPerformance', 'ShiftingPerformance', 'init']

FILE_PATH = os.path.abspath(__file__)
THERMO_FILE = os.path.dirname(FILE_PATH) + '/data/thermo.dat'
PROPELLANT_FILE = os.path.dirname(FILE_PATH) + '/data/propellant.dat'

def __convert_struct_field(s, fields):
    for field, fieldtype in fields:
        if fieldtype.type.kind == 'primitive':
            yield (field, getattr(s, field))
        else:
            yield (field, convert_to_python(getattr(s, field)))


def convert_to_python(s):
    '''
    Given a cdata struct, returns a dict with all of the parsable
    members.  Borrowed with mods from :
    http://stackoverflow.com/questions/20444546/python-cffi-convert-structure-to-dictionary
    '''
    if isinstance(s, ffi.CData) == False:
        return s

    type = ffi.typeof(s)
    if type.kind == 'struct':
        return dict(__convert_struct_field(s, type.fields))
    elif type.kind == 'array':
        if type.item.kind == 'primitive':
            if type.item.cname == 'char':
                return ffi.string(s).decode('utf-8')
            else:
                return [s[i] for i in range(type.length)]
        else:
            return [convert_to_python(s[i]) for i in range(type.length)]
    elif type.kind == 'primitive':
        return int(s)


def init(thermo_file=None, propellant_file=None):
    global THERMO_FILE, PROPELLANT_FILE, SPECIES, PROPELLANTS
    if thermo_file is not None:
        THERMO_FILE = thermo_file
    if propellant_file is not None:
        PROPELLANT_FILE = propellant_file

    r = lib.load_thermo(THERMO_FILE.encode('utf-8'))
    if r > 0:
        print("Loaded {} thermo species".format(r))
    else:
        print("Failed to load thermo file {}".format(THERMO_FILE))

    r = lib.load_propellant(PROPELLANT_FILE.encode('utf-8'))

    if r > 0:
        print("Loaded {} propellants".format(r))
    else:
        print("Failed to load propellant file {}".format(PROPELLANT_FILE))

    SPECIES = dict()
    PROPELLANTS = dict()

    # Build species dict
    for i in range(lib.num_thermo):
        s = lib.thermo_list[i]
        l = len(SPECIES)
        name = ffi.string(s.name).decode('utf-8')
        while name in SPECIES:
            name += "'"
        SPECIES[name] = AttrDict(convert_to_python(s))
        SPECIES[name]['id'] = i
        if len(SPECIES) <= l:
            raise RuntimeWarning("Species {}, {}:{} dropped".format(i, name,
                                                            str(SPECIES[name])))

    # Build propellants dict
    for i in range(lib.num_propellant):
        p = lib.propellant_list[i]
        name = ffi.string(p.name).decode('utf-8')
        l = len(PROPELLANTS)
        while name in PROPELLANTS:
            name += "'"
        PROPELLANTS[name] = Propellant(convert_to_python(p))
        PROPELLANTS[name]['id'] = i
        if len(PROPELLANTS) <= l:
            raise RuntimeWarning("Propellant {}, {}:{} dropped".format(i,
                                                name, str(PROPELLANTS[name])))

def find_propellant(substr):
    return [value for key, value in list(PROPELLANTS.items())
                        if substr.lower() in key.lower()]
