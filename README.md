# pypropep
[![Build Status](https://travis-ci.org/jonnydyer/pypropep.svg?branch=master)](https://travis-ci.org/jonnydyer/pypropep)
Python interfaces to [cpropep](https://sourceforge.net/projects/rocketworkbench/?source=navbar), a project started ~15 years ago by Antoine Lefebvre and Ray Calkins to implement the [classic Gordon and McBride](https://www.grc.nasa.gov/WWW/CEAWeb/RP-1311.pdf) chemical equillibrium code in C.  

The RocketWorkBench (and hence cpropep) project hasn't seen any activity in close to 15 years and yet I still found myself prefering cpropep over [CEA](https://www.grc.nasa.gov/WWW/CEAWeb/) for a number of reasons.  Calling it manually from the command line or writing new C scripts for every analysis task is cumbersome to say the least.  The goal of this project is to bring cpropep into the 21st century with a clean, useful Python interface.

##Basic usage

