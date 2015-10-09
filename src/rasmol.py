# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# rasmol tool

oneline = "3d visualization via RasMol program"

docstr = """
r = rasmol(p)           create RasMol wrapper for pdb object p

r.file = "image"        file prefix for created images (def = "image")

r.show(N)               show snapshot at timestep N with default script
r.show(N,"my.rasmol")   use file as RasMol script

r.all()                 make images of all selected snapshots with def script
r.all("my.rasmol")      use file as RasMol script
  
r.run(N)                run RasMol interactivly on snapshot N
r.run(N,"new.rasmol")                 adjust via mouse or RasMol commands
r.run(N,"new.rasmol","old.rasmol")    type quit to save RasMol script file

  if 2 args, 2nd arg is new script file, else save to "tmp.rasmol"
  if 3 args, 3rd arg is initial script file, else use default script
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list
#   allow zoom, view, movie, pan settings like other viz tools
#     so can specify rock & make fly-by movies of an image

# Variables

# Imports and external programs

import sys, os, commands, re, types
from pdbfile import pdbfile

try: from DEFAULTS import PIZZA_RASMOL
except: PIZZA_RASMOL = "rasmol"
try: from DEFAULTS import PIZZA_DISPLAY
except: PIZZA_DISPLAY = "display"

# Class definition

class rasmol:

  # --------------------------------------------------------------------

  def __init__(self,pdb):
    self.pdb = pdb
    self.file = "image"
    
  # --------------------------------------------------------------------

  def start(self):
    self.RASMOL = os.popen(PIZZA_RASMOL,'w')

  # --------------------------------------------------------------------

  def enter(self):
    while 1:
      command = raw_input("rasmol> ")
      if command == "quit" or command == "exit": return
      self.__call__(command)

  # --------------------------------------------------------------------

  def __call__(self,command):
    self.RASMOL.write(command + '\n')
    self.RASMOL.flush()

  # --------------------------------------------------------------------

  def stop(self):
    self.__call__("quit")
    del self.RASMOL

  # --------------------------------------------------------------------

  def show(self,*list):

    # create tmp.pdb with atom data
    
    n = list[0]
    self.pdb.single(n,"tmp.pdb")

    # if RasMol input script specified, read it
    # replace load pdb "file" with load pdb "%s"
    # if no RasMol input script specified, use rasmol_template

    if len(list) == 2:
      rasmol_text = open(list[1],"r").read()
      rasmol_text = re.sub('load pdb ".*"','load pdb "%s"',rasmol_text)
    else:
      rasmol_text = rasmol_template

    # write rasmol_text to tmp.rasmol, substituting tmp.pdb for filename
    
    f = open("tmp.rasmol","w")
    text = rasmol_text % "tmp.pdb"
    print >>f,text
    f.close()

    # run RasMol to create image in tmp.gif
    
    self.start()
    self.__call__("source tmp.rasmol")
    self.__call__("write tmp.gif")
    self.stop()

    # display the image
    
    cmd = "%s tmp.gif" % (PIZZA_DISPLAY)
    commands.getoutput(cmd)

  # --------------------------------------------------------------------

  def all(self,*list):

    # if RasMol script specified, read it
    # replace load pdb "file" with load pdb "%s"
    # if no RasMol script specified, just use rasmol_template

    if len(list) == 1:
      rasmol_text = open(list[0],"r").read()
      rasmol_text = re.sub('load pdb ".*"','load pdb "%s"',rasmol_text)
    else:
      rasmol_text = rasmol_template
      
    # iterate over all timesteps
    # write snapshot to tmpN.pdb
    # write RasMol input script to tmpN.rasmol

    n = flag = 0
    while 1:
      which,time,flag = self.pdb.iterator(flag)

      if flag == -1: break

      if n < 10:
        ncount = "000" + str(n)
      elif n < 100:
        ncount = "00" + str(n)
      elif n < 1000:
        ncount = "0" + str(n)
      else:
        ncount = str(n)

      file_pdb = "tmp%s.pdb" % ncount
      self.pdb.single(time,file_pdb)

      text = rasmol_text % file_pdb
      file_rasmol = "tmp%s.rasmol" % ncount
      f = open(file_rasmol,"w")
      print >>f,text
      f.close()

      print time,
      sys.stdout.flush()
      n += 1

    # run RasMol on each pair of RasMol scripts and PDB files

    self.start()

    loop = n
    for n in xrange(loop):
      if n < 10:
        ncount = "000" + str(n)
      elif n < 100:
        ncount = "00" + str(n)
      elif n < 1000:
        ncount = "0" + str(n)
      else:
        ncount = str(n)

      source_str = "source tmp%s.rasmol" % ncount
      self.__call__(source_str)
      write_str = "write %s%s.gif" % (self.file,ncount)
      self.__call__(write_str)

    self.stop()

    # clean up

    commands.getoutput("rm tmp*.pdb")
    commands.getoutput("rm tmp*.rasmol")
    
  # --------------------------------------------------------------------

  def run(self,*list):

    # create tmp.pdb with atom data

    n = list[0]
    self.pdb.single(n,"tmp.pdb")

    # if RasMol script specified, read it
    # replace load pdb "file" with load pdb "%s"
    # if no RasMol script specified, just use rasmol_template

    if len(list) == 3:
      rasmol_text = open(list[2],"r").read()
      rasmol_text = re.sub('load pdb ".*"','load pdb "%s"',rasmol_text)
    else:
      rasmol_text = rasmol_template

    # write rasmol_text to tmp.rasmol
    
    f = open("tmp.rasmol","w")
    text = rasmol_template % "tmp.pdb"
    print >>f,text
    f.close()

    # run RasMol to create image in tmp.gif
    
    self.start()
    self.__call__("source tmp.rasmol")
    self.enter()
    
    if len(list) > 1: newfile = list[1]
    else: newfile = "tmp.rasmol"
    self.__call__("write script %s" % newfile)
    self.stop()

# --------------------------------------------------------------------
# generic Rasmol script with spacefill option
# PDB filename must be filled in for load command

rasmol_template = """
# Creator: RasMol Version 2.7.1

zap
load pdb "%s"
background [0,0,0]
set ambient 40
set specular off

reset
slab off

set axes off
set boundingbox off
set unitcell off
set bondmode and
dots off

# Avoid Colour Problems!
select all
colour bonds none
colour backbone none
colour hbonds none
colour ssbonds none
colour ribbons none
colour white

# Atoms
colour atoms cpk
spacefill on
set shadow off

# Bonds
wireframe off

# Ribbons
ribbons off

# Backbone
backbone off

# Labels
labels off

# Monitors
monitors off

ssbonds off
hbonds off
"""
