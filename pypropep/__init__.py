import os
from _cpropep import ffi, lib

FILE_PATH = os.path.abspath(__file__)
THERMO_FILE = os.path.dirname(FILE_PATH) + '/data/thermo.dat'
PROPELLANT_FILE = os.path.dirname(FILE_PATH) + '/data/propellant.dat'

r = lib.load_thermo(THERMO_FILE)
if r > 0:
    print 'Loaded %d thermo species' % (r)
else:
    print 'Failed to load thermo file %s' % (THERMO_FILE)

#r = lib.load_propellant(PROPELLANT_FILE)

if r > 0:
    print 'Loaded %d propellants' % (r)
else:
    print 'Failed to load propellant file %s' % (PROPELLANT_FILE)
