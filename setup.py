'''
Python interface to cpropep
'''
from setuptools import setup

setup(
    name="pypropep",
    version='0.1.1',
    description='Python wrapper for cpropep rocket performance tool',
    license='GPLv3',
    author='Jonny Dyer',
    author_email='jonny.dyer@gmail.com',
    url='https://github.com/jonnydyer/pypropep',
    packages=['pypropep', 'pypropep.cpropep'],
    package_data={"pypropep" : ["data/*.dat"]},
    setup_requires=["cffi>=1.0.0"],
    cffi_modules=["cpropep_build.py:ffibuilder"],
    install_requires=["cffi>=1.0.0", "attrdict"],
    py_modules=['cpropep_build'],
)

# from distutils.core import setup
# import cpropep_build
#
# setup(name='pypropep',
#       version='0.1',
#       description='Python wrapper for cpropep rocket performance tool',
#       license='GPLv3',
#       author='Jonny Dyer',
#       author_email='jonny.dyer@gmail.com',
#       packages=['pypropep', 'pypropep.cpropep'],
#       package_data={"pypropep" : ["data/*.dat"]},
#       py_modules=['cpropep_build'],
#       ext_modules=[cpropep_build.ffibuilder.distutils_extension()]
#       )

