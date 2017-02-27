'''
Python interface to cpropep
'''
from distutils.core import setup
import cpropep_build

setup(name='pypropep',
      version='0.1',
      description='Python wrapper for cpropep rocket performance tool',
      license='GPLv3',
      author='Jonny Dyer',
      author_email='jonny.dyer@gmail.com',
      packages=['pypropep'],
      package_data={"pypropep" : ["data/*.dat"]},
      ext_modules=[cpropep_build.ffibuilder.distutils_extension()],
      )
