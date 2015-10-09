#!/usr/bin/python

# Script:  dview.py
# Purpose: launch vcr tool on LAMMPS dump files
# Syntax:  dview.py --n 512 dump.1 dump.2 ...
#          files = one or more dump files
# Example: dview.py dump.*
# Author:  Steve Plimpton (Sandia)

# main script

if len(argv) < 2:
  raise StandardError, "Syntax: dview.py --n 512 dump.1 ..."

if argv[1] == "--n":
  n = int(argv[2])
  files = ' '.join(argv[3:])
else:
  n = 512
  files = ' '.join(argv[1:])

d = dump(files)
g = gl(d,n)
v = vcr(g)
