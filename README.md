# pypropep
[![Build Status](https://travis-ci.org/jonnydyer/pypropep.svg?branch=master)](https://travis-ci.org/jonnydyer/pypropep)
[![Coverage Status](https://coveralls.io/repos/github/jonnydyer/pypropep/badge.svg?branch=master)](https://coveralls.io/github/jonnydyer/pypropep?branch=master)

Python interfaces to [cpropep](https://sourceforge.net/projects/rocketworkbench/?source=navbar), a project started ~15 years ago by Antoine Lefebvre and Ray Calkins to implement the [classic Gordon and McBride](https://www.grc.nasa.gov/WWW/CEAWeb/RP-1311.pdf) chemical equillibrium code in C.  

The RocketWorkBench (and hence cpropep) project hasn't seen any activity in close to 15 years and yet I still found myself prefering cpropep over [CEA](https://www.grc.nasa.gov/WWW/CEAWeb/) for a number of reasons including accessibility, it's C-heritage (rather than Fortran) and the extensive propellant library distributed with all flavors of PROPEP.  Calling it manually from the command line or writing and compiling new C executables for every analysis task is cumbersome to say the least.  The goal of this project is to bring cpropep into the 21st century with a clean, useful Python interface.

## Version 0.1
Today, this module is at v0.1.  It is only tested and guaranteed to work on Python 2.7 (something I hope to fix for the next release) and requires the users machine to have a compiler available for installation (also hope to fix for 0.2).  Otherwise, it is fairly functional.

# Installation

Currently the two methods for installing pypropep are pip and from source using setuptools.  Pip is recommended.

## Pip
Not much to it -

    pip install pypropep

## From source

    git clone https://github.com/jonnydyer/pypropep.git
    cd pypropep
    python setup.py install

## Conda

    Coming soon...
    
# Usage

## Basic Usage
Here is a brief example of how to use pypropep::

    >>> import pypropep as ppp
    >>> ppp.init()
    Loaded 1921 thermo species
    Loaded 1030 propellants
    
    >>> o2 = ppp.PROPELLANTS['OXYGEN (GAS)']
    >>> ch4 = ppp.PROPELLANTS['METHANE']
    >>> sp = ppp.ShiftingPerformance()
    
    >>> OF = 2.8
    >>> sp.add_propellants_by_mass([(ch4, 1.0), (o2, OF)])
    >>> sp.set_state(P=50., Pe=1.)     # Pressure in atm
    
    >>> print sp.performance.cstar      # in m/s
    1892.82959658
    
    >>> print sp.performance.cf
    1.57123484882
    
    >>> print sp.performance.Isp/9.8      # in seconds
    303.477533166
    
    >>> print sp.performance.cstar * sp.performance.cf / 9.8     # in seconds
    303.477533166

## iPython examples
More detailed examples demonstrating the utility of the library are given in the form of two Jupyter notebooks (kindly rendered here by Git!)

- [Basic Usage and Background](ipython_doc/BasicUsage.ipynb)
- [Rocket Performance Examples](ipython_doc/BasicRocketPerformance.ipynb)

# Roadmap

## v0.2
There are several things I'd like to add to this module for the v0.2 release including:

- Set up [Python Wheels](http://pythonwheels.com/) distribution so that installation doesn't require local compiling
- Set up an [Anaconda](https://www.continuum.io/anaconda-overview) distribution for the module
- Support Python versions other than 2.7 in both test and deployment

## v0.3
- Add Finite Area Contraction-ratio (FAC) support to cpropep library and python interface
