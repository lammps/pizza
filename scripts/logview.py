#!/usr/bin/python

# Script:  logview.py
# Purpose: plots of LAMMPS log-file thermodynamic data
# Syntax:  logview.py gnu/matlab files ...
#          gnu/matlab = style of plots to create
#          files = one or more log files
# Example: logview.py gnu log.*
# Author:  Steve Plimpton (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from log import log
from plotview import plotview
from gnu import gnu
from matlab import matlab
if not globals().has_key("argv"): argv = sys.argv

# main script

if len(argv) < 3:
  raise StandardError, "Syntax: logview.py gnu/matlab files  ..."

style = argv[1]
files = ' '.join(argv[2:])

lg = log(files)
exec "plot = %s()" % style
p = plotview(lg,plot)

