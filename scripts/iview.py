#!/usr/bin/python

# Script:  iview.py
# Purpose: launch animate tool on series of image files
# Syntax:  iview.py files ...
#          files = one or more image files
# Example: iview.py image*png
# Author:  Steve Plimpton (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from animate import animate
if not globals().has_key("argv"): argv = sys.argv

# main script
# this could be done with one-line alias in shell start-up file

if len(argv) < 2: raise StandardError, "Syntax: iview.py files  ..."
a = animate(' '.join(argv[1:]))
