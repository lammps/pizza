# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# bdump tool

oneline = "Read dump files with bond info"

docstr = """
b = bdump("dump.one")             read in one or more dump files
b = bdump("dump.1 dump.2.gz")	  can be gzipped
b = bdump("dump.*")		  wildcard expands to multiple files
b = bdump("dump.*",0)		  two args = store filenames, but don't read

  incomplete and duplicate snapshots are deleted
  no column name assignment is performed

time = b.next()             	  read next snapshot from dump files

  used with 2-argument constructor to allow reading snapshots one-at-a-time
  snapshot will be skipped only if another snapshot has same time stamp
  return time stamp of snapshot read
  return -1 if no snapshots left or last snapshot is incomplete
  no column name assignment is performed

b.map(1,"id",3,"x")               assign names to atom columns (1-N)

  must assign id,type,atom1,atom2

time,box,atoms,bonds,tris,lines = b.viz(index)   return list of viz objects

  viz() returns line info for specified timestep index
    can also call as viz(time,1) and will find index of preceding snapshot
    time = timestep value
    box = NULL
    atoms = NULL
    bonds = id,type,atom1,atom2 for each line as 2d array
    tris = NULL
    lines = NULL
"""

# History
#   11/10, Steve Plimpton (SNL): original version

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
#     atoms[i][j] = 2d array of floats, i = 0 to natoms-1, j = 0 to ncols-1

# Imports and external programs

import sys, commands, re, glob, types
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

class bdump:

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
      raise StandardError,"no bdump file specified"
    
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
  # assign column names if not already done and file is self-describing
  # convert xs,xu to x
  
  def read_snapshot(self,f):
    try:
      snap = Snap()
      item = f.readline()
      snap.time = int(f.readline().split()[0])    # just grab 1st field
      item = f.readline()
      snap.natoms = int(f.readline())
      item = f.readline()

      f.readline()    # read past BOX BOUNDS
      f.readline()
      f.readline()
      f.readline()

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
      raise StandardError, "bdump map() requires pairs of mappings"
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
  # delete successive snapshots with duplicate time stamp

  def cull(self):
    i = 1
    while i < len(self.snaps):
      if self.snaps[i].time == self.snaps[i-1].time:
        del self.snaps[i]
      else:
        i += 1
  
  # --------------------------------------------------------------------
  # return list of bonds to viz for snapshot isnap
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
    id = self.names["id"]
    type = self.names["type"]
    atom1 = self.names["atom1"]
    atom2 = self.names["atom2"]

    # create line list from id,type,atom1,atom2
    # abs() of type since could be negative
    
    bonds = []
    for i in xrange(snap.natoms):
      atom = snap.atoms[i]
      bonds.append([int(atom[id]),abs(int(atom[type])),
                    int(atom[atom1]),int(atom[atom2])])

    return time,None,None,bonds,None,None

# --------------------------------------------------------------------
# one snapshot

class Snap:
  pass
