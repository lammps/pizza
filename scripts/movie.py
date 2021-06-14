#!/usr/bin/python

# Script:  movie.py
# Purpose: create images from LAMMPS dump snapshots
# Syntax:  movie.py raster/svg theta phi dump.1 dump.2 ...
#          raster/svg = style of image to create
#	   theta/phi = vertical (z) and azimuthal angle to view from
#          files = one or more dump files
# Example: movie.py svg 60 130 dump.*
# Author:  Steve Plimpton (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from dump import dump
from raster import raster
from svg import svg
if not globals().has_key("argv"): argv = sys.argv

# main script

if len(argv) < 5:
  raise StandardError, "Syntax: movie.py raster/svg theta phi dump.1 ..."

style = argv[1]
theta = float(argv[2])
phi = float(argv[3])
files = ' '.join(argv[4:])

d = dump(files)
exec "viz = %s(d)" % style
viz.rotate(theta,phi)
viz.all()
