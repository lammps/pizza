# Pizza.py toolkit, https://pizza.sandia.gov/
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# tdump tool

oneline = "Read dump files with triangle info"

docstr = """
t = tdump("dump.one")             read in one or more dump files
t = tdump("dump.1 dump.2.gz")	  can be gzipped
t = tdump("dump.*")		  wildcard expands to multiple files
t = tdump("dump.*",0)		  two args = store filenames, but don't read

  incomplete and duplicate snapshots are deleted
  no column name assignment is performed

time = t.next()             	  read next snapshot from dump files

  used with 2-argument constructor to allow reading snapshots one-at-a-time
  snapshot will be skipped only if another snapshot has same time stamp
  return time stamp of snapshot read
  return -1 if no snapshots left or last snapshot is incomplete
  no column name assignment is performed

t.map(1,"id",3,"x")               assign names to atom columns (1-N)

  must assign id,type,corner1x,corner1y,corner1z,corner2x,corner2y,corner2z,corner3x,corner3y,corner3z

time,box,atoms,bonds,tris,lines = t.viz(index)   return list of viz objects

  viz() returns line info for specified timestep index
    can also call as viz(time,1) and will find index of preceding snapshot
    time = timestep value
    box = \[xlo,ylo,zlo,xhi,yhi,zhi\]
    atoms = NULL
    bonds = NULL
    tris = id,type,x1,y1,z1,x2,y2,z2,x3,y3,z3 for each tri as 2d array
      id,type are from associated atom
    lines = NULL

t.owrap(...)		          wrap tris to same image as their atoms

  owrap() is called by dump tool's owrap()
  useful for wrapping all molecule's atoms/tris the same so it is contiguous
"""

# History
#   4/11, Steve Plimpton (SNL): original version

# Variables
#   flist = list of dump file names
#   increment = 1 if reading snapshots one-at-a-time
#   nextfile = which file to read from via next()
#   eof = ptr into current file for where to read via next()
#   nsnaps = # of snapshots
#   snaps = list of snapshots
#   names = dictionary of column names:
#     key = "id", value = column # (0 to M-1)
#   Snap = one snapshot
#     time = time stamp
#     natoms = # of atoms
#     xlo,xhi,ylo,yhi,zlo,zhi = box bounds (float)
#     atoms[i][j] = 2d array of floats, i = 0 to natoms-1, j = 0 to ncols-1

# Imports and external programs

import sys, commands, re, glob, types
from math import sqrt
from os import popen

try:
    import numpy as np
    oldnumeric = False
except:
    import Numeric as np
    oldnumeric = True

try: from DEFAULTS import PIZZA_GUNZIP
except: PIZZA_GUNZIP = "gunzip"

# Class definition

class tdump:

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.snaps = []
    self.nsnaps = 0
    self.names = {}

    # flist = list of all dump file names

    words = list[0].split()
    self.flist = []
    for word in words: self.flist += glob.glob(word)
    if len(self.flist) == 0 and len(list) == 1:
      raise StandardError,"no ldump file specified"
    
    if len(list) == 1:
      self.increment = 0
      self.read_all()
    else:
      self.increment = 1
      self.nextfile = 0
      self.eof = 0

  # --------------------------------------------------------------------

  def read_all(self):

    # read all snapshots from each file
    # test for gzipped files

    for file in self.flist:
      if file[-3:] == ".gz":
        f = popen("%s -c %s" % (PIZZA_GUNZIP,file),'r')
      else: f = open(file)

      snap = self.read_snapshot(f)
      while snap:
        self.snaps.append(snap)
        print snap.time,
        sys.stdout.flush()
        snap = self.read_snapshot(f)

      f.close()
    print

    # sort entries by timestep, cull duplicates

    self.snaps.sort(self.compare_time)
    self.cull()
    self.nsnaps = len(self.snaps)
    print "read %d snapshots" % self.nsnaps

  # --------------------------------------------------------------------
  # read next snapshot from list of files

  def next(self):

    if not self.increment: raise StandardError,"cannot read incrementally"

    # read next snapshot in current file using eof as pointer
    # if fail, try next file
    # if new snapshot time stamp already exists, read next snapshot

    while 1:
      f = open(self.flist[self.nextfile],'rb')
      f.seek(self.eof)
      snap = self.read_snapshot(f)
      if not snap:
        self.nextfile += 1
	if self.nextfile == len(self.flist): return -1
        f.close()
	self.eof = 0
	continue
      self.eof = f.tell()
      f.close()
      try:
        self.findtime(snap.time)
	continue
      except: break

    self.snaps.append(snap)
    snap = self.snaps[self.nsnaps]
    self.nsnaps += 1

    return snap.time

  # --------------------------------------------------------------------
  # read a single snapshot from file f
  # return snapshot or 0 if failed
  
  def read_snapshot(self,f):
    try:
      snap = Snap()
      item = f.readline()
      snap.time = int(f.readline().split()[0])    # just grab 1st field
      item = f.readline()
      snap.natoms = int(f.readline())

      item = f.readline()
      words = f.readline().split()
      snap.xlo,snap.xhi = float(words[0]),float(words[1])
      words = f.readline().split()
      snap.ylo,snap.yhi = float(words[0]),float(words[1])
      words = f.readline().split()
      snap.zlo,snap.zhi = float(words[0]),float(words[1])

      item = f.readline()

      if snap.natoms:
        words = f.readline().split()
        ncol = len(words)
        for i in xrange(1,snap.natoms):
          words += f.readline().split()
        floats = map(float,words)
        if oldnumeric: atoms = np.zeros((snap.natoms,ncol),np.Float)
        else: atoms = np.zeros((snap.natoms,ncol),np.float)
        start = 0
        stop = ncol
        for i in xrange(snap.natoms):
          atoms[i] = floats[start:stop]
          start = stop
          stop += ncol
      else: atoms = None
      snap.atoms = atoms
      return snap
    except:
      return 0

  # --------------------------------------------------------------------
  # map atom column names
  
  def map(self,*pairs):
    if len(pairs) % 2 != 0:
      raise StandardError, "tdump map() requires pairs of mappings"
    for i in range(0,len(pairs),2):
      j = i + 1
      self.names[pairs[j]] = pairs[i]-1

  # --------------------------------------------------------------------
  # return vector of snapshot time stamps

  def time(self):
    vec = self.nsnaps * [0]
    i = 0
    for snap in self.snaps:
      vec[i] = snap.time
      i += 1
    return vec

  # --------------------------------------------------------------------
  # sort snapshots on time stamp

  def compare_time(self,a,b):
    if a.time < b.time:
      return -1
    elif a.time > b.time:
      return 1
    else:
      return 0

  # --------------------------------------------------------------------
    
  def findtime(self,n):
    for i in xrange(self.nsnaps):
      if self.snaps[i].time == n: return i
    raise StandardError, "no step %d exists" % n

  # --------------------------------------------------------------------
  # delete successive snapshots with duplicate time stamp

  def cull(self):
    i = 1
    while i < len(self.snaps):
      if self.snaps[i].time == self.snaps[i-1].time:
        del self.snaps[i]
      else:
        i += 1
  
  # --------------------------------------------------------------------
  # return list of lines to viz for snapshot isnap
  # if called with flag, then index is timestep, so convert to snapshot index

  def viz(self,index,flag=0):
    if not flag: isnap = index
    else:
      times = self.time()
      n = len(times)
      i = 0
      while i < n:
        if times[i] > index: break
        i += 1
      isnap = i - 1
    snap = self.snaps[isnap]

    time = snap.time
    box = [snap.xlo,snap.ylo,snap.zlo,snap.xhi,snap.yhi,snap.zhi]
    id = self.names["id"]
    type = self.names["type"]
    corner1x = self.names["corner1x"]
    corner1y = self.names["corner1y"]
    corner1z = self.names["corner1z"]
    corner2x = self.names["corner2x"]
    corner2y = self.names["corner2y"]
    corner2z = self.names["corner2z"]
    corner3x = self.names["corner3x"]
    corner3y = self.names["corner3y"]
    corner3z = self.names["corner3z"]

    # create trid list from id,type,corner1x,...corner3z
    # don't add tri if all 4 values are 0 since not a line
    
    tris = []
    for i in xrange(snap.natoms):
      atom = snap.atoms[i]
      c1 = [atom[corner1x],atom[corner1y],atom[corner1z]]
      c2 = [atom[corner2x],atom[corner2y],atom[corner2z]]
      c3 = [atom[corner3x],atom[corner3y],atom[corner3z]]
      n = normal(c1,c2,c3)
      if c1[0] == 0.0 and c1[1] == 0.0 and c1[2] == 0.0 and \
              c2[0] == 0.0 and c2[1] == 0.0 and c2[2] == 0.0: continue
      tris.append([atom[id],atom[type]] + c1 + c2 + c3 + n)

    return time,box,None,None,tris,None

  # --------------------------------------------------------------------
  # wrap tri corner points associated with atoms thru periodic boundaries
  # invoked by dump() when it does an owrap() on its atoms

  def owrap(self,time,xprd,yprd,zprd,idsdump,atomsdump,iother,ix,iy,iz):
    id = self.names["id"]
    corner1x = self.names["corner1x"]
    corner1y = self.names["corner1y"]
    corner1z = self.names["corner1z"]
    corner2x = self.names["corner2x"]
    corner2y = self.names["corner2y"]
    corner2z = self.names["corner2z"]
    corner3x = self.names["corner3x"]
    corner3y = self.names["corner3y"]
    corner3z = self.names["corner3z"]

    isnap = self.findtime(time)
    snap = self.snaps[isnap]
    atoms = snap.atoms

    # idump = index of my line I in dump's atoms
    # jdump = atom J in dump's atoms that atom I was owrapped on
    # delx,dely = offset applied to atom I and thus to line I
    
    for i in xrange(snap.natoms):
      tag = atoms[i][id]
      idump = idsdump[tag]
      jdump = idsdump[atomsdump[idump][iother]]
      delx = (atomsdump[idump][ix]-atomsdump[jdump][ix])*xprd
      dely = (atomsdump[idump][iy]-atomsdump[jdump][iy])*yprd
      delz = (atomsdump[idump][iz]-atomsdump[jdump][iz])*zprd
      atoms[i][corner1x] += delx
      atoms[i][corner1y] += dely
      atoms[i][corner1z] += delz
      atoms[i][corner2x] += delx
      atoms[i][corner2y] += dely
      atoms[i][corner2z] += delz
      atoms[i][corner3x] += delx
      atoms[i][corner3y] += dely
      atoms[i][corner3z] += delz

# --------------------------------------------------------------------
# one snapshot

class Snap:
  pass

# --------------------------------------------------------------------
# compute normal for a triangle with 3 vertices

def normal(x,y,z):
  v1 = 3*[0]
  v1[0] = y[0] - x[0]
  v1[1] = y[1] - x[1]
  v1[2] = y[2] - x[2]
  
  v2 = 3*[0]
  v2[0] = z[0] - y[0]
  v2[1] = z[1] - y[1]
  v2[2] = z[2] - y[2]
      
  n = 3*[0]
  n[0] = v1[1]*v2[2] - v1[2]*v2[1]
  n[1] = v1[2]*v2[0] - v1[0]*v2[2]
  n[2] = v1[0]*v2[1] - v1[1]*v2[0]
  
  length = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
  n[0] /= length
  n[1] /= length
  n[2] /= length
  
  return n
