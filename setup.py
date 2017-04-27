#!/usr/bin/env python

from distutils.core import setup

setup(name='Pizza.py',
      version='1.0',
      description='Pizza.py is a loosely integrated collection of tools written in Python, many of which provide pre- '
                  'and post-processing capability for the LAMMPS molecular dynamics package. There are tools to create '
                  'input files, convert between file formats, process log and dump files, create plots, and visualize '
                  'and animate simulation snapshots.',
      author='Steve Plimpton',
      author_email='sjplimp@sandia.gov',
      url='http://pizza.sandia.gov/',
      packages=['pizza'],
     )