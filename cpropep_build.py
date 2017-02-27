# file "cpropep_build.py"

from cffi import FFI
from glob import glob

ffibuilder = FFI()

cpropep_libs = ['libnum', 'libthermo', 'libcpropep', 'libcompat']
inc_dir = [('pypropep/cpropep/' + d + '/include/') for d in cpropep_libs]

src_files = []
for l in cpropep_libs:
    src_files += glob('pypropep/cpropep/' + l + '/src/*.c')

ffibuilder.set_source("pypropep.cpropep._cpropep",
    r"""
    // blah
    """,
    sources=src_files,
    include_dirs=inc_dir)

ffibuilder.cdef("""
// FROM libthermo/load.h
int load_thermo(char *filename);
int load_propellant(char *filename);

//FROM libthermo/thermo.h
    """)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
