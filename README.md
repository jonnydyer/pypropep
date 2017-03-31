# pypropep
[![Build Status](https://travis-ci.org/jonnydyer/pypropep.svg?branch=master)](https://travis-ci.org/jonnydyer/pypropep)
[![Coverage Status](https://coveralls.io/repos/github/jonnydyer/pypropep/badge.svg?branch=master)](https://coveralls.io/github/jonnydyer/pypropep?branch=master)

Python interfaces to [cpropep](https://sourceforge.net/projects/rocketworkbench/?source=navbar), a project started ~15 years ago by Antoine Lefebvre and Ray Calkins to implement the [classic Gordon and McBride](https://www.grc.nasa.gov/WWW/CEAWeb/RP-1311.pdf) chemical equillibrium code in C.  

The RocketWorkBench (and hence cpropep) project hasn't seen any activity in close to 15 years and yet I still found myself prefering cpropep over [CEA](https://www.grc.nasa.gov/WWW/CEAWeb/) for a number of reasons.  Calling it manually from the command line or writing new C scripts for every analysis task is cumbersome to say the least.  The goal of this project is to bring cpropep into the 21st century with a clean, useful Python interface.

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
    >>> sp.set_state(P=50., Pe=1.)
    >>> print sp.performance.cstar
    1892.82959658
    >>> print sp.performance.cf
    1.57123484882
    >>> print sp.performance.Isp/9.8
    303.477533166
    >>> print sp.performance.cstar * sp.performance.cf / 9.8
    303.477533166

## iPython examples
The best intro to this library is an [iPython Notebook included in the repo](ipython_doc/BasicUsage.ipynb).  Git will kindly render it for you if you don't want to run it yourself.
