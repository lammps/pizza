# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# mdump tool

oneline = "Read, write, manipulate mesh dump files"

docstr = """
m = mdump("mesh.one")             read in one or more mesh dump files
m = mdump("mesh.1 mesh.2.gz")	  can be gzipped
m = mdump("mesh.*")		  wildcard expands to multiple files
m = mdump("mesh.*",0)		  two args = store filenames, but don't read

  incomplete and duplicate snapshots are deleted

time = m.next()             	  read next snapshot from dump files

  used with 2-argument constructor to allow reading snapshots one-at-a-time
  snapshot will be skipped only if another snapshot has same time stamp
  return time stamp of snapshot read
  return -1 if no snapshots left or last snapshot is incomplete
  no column name assignment or unscaling is performed

m.map(2,"temperature")            assign names to element value columns (1-N)

m.tselect.all()			  select all timesteps
m.tselect.one(N)		  select only timestep N
m.tselect.none()		  deselect all timesteps
m.tselect.skip(M)		  select every Mth step
m.tselect.test("$t >= 100 and $t < 10000")      select matching timesteps
m.delete()	      	      	  delete non-selected timesteps

  selecting a timestep also selects all elements in the timestep
  skip() and test() only select from currently selected timesteps
  test() uses a Python Boolean expression with $t for timestep value
    Python comparison syntax: == != < > <= >= and or

m.eselect.all()	      	                      select all elems in all steps
m.eselect.all(N)      	                      select all elems in one step
m.eselect.test("$id > 100 and $type == 2")    select match elems in all steps
m.eselect.test("$id > 100 and $type == 2",N)  select matching elems in one step

  all() with no args selects atoms from currently selected timesteps
  test() with one arg selects atoms from currently selected timesteps
  test() sub-selects from currently selected elements
  test() uses a Python Boolean expression with $ for atom attributes
    Python comparison syntax: == != < > <= >= and or
    $name must end with a space

t = m.time()  	     	       	   return vector of selected timestep values
fx,fy,... = m.vecs(1000,"fx","fy",...)  return vector(s) for timestep N

  vecs() returns vectors with one value for each selected elem in the timestep

index,time,flag = m.iterator(0/1)          loop over mesh dump snapshots
time,box,atoms,bonds,tris,lines = m.viz(index)  return list of viz objects
nodes,elements,nvalues,evalues = m.mviz(index)  return list of mesh viz objects
m.etype = "color"                          set column returned as "type" by viz

  iterator() loops over selected timesteps
  iterator() called with arg = 0 first time, with arg = 1 on subsequent calls
    index = index within dump object (0 to # of snapshots)
    time = timestep value
    flag = -1 when iteration is done, 1 otherwise
  viz() returns info for selected elements for specified timestep index
    can also call as viz(time,1) and will find index of preceding snapshot
    time = timestep value
    box = \[xlo,ylo,zlo,xhi,yhi,zhi\]
    atoms = NULL
    bonds = NULL
    tris = id,type,x1,y1,z1,x2,y2,z2,x3,y3,z3,nx,ny,nz for each tri as 2d array
      each element is decomposed into tris
    lines = NULL
  mviz() returns info for all elements for specified timestep index
    can also call as mviz(time,1) and will find index of preceding snapshot
    time = timestep value
    box = \[xlo,ylo,zlo,xhi,yhi,zhi\]
    nodes = list of nodes = id,type,x,y,z
    elements = list of elements = id,type,node1,node2,...
    nvalues = list of node values = id,type,value1,value2,...
    evalues = list of element values = id,type,value1,value2,...
  etype is column name viz() will return as element type (def = "" = elem type)
"""

# History
#   11/06, Steve Plimpton (SNL): original version
#   12/09, David Hart (SNL): allow use of NumPy or Numeric

# Variables
#   flist = list of dump file names
#   increment = 1 if reading snapshots one-at-a-time
#   nextfile = which file to read from via next()
#   eof = ptr into current file for where to read via next()
#   nsnaps = # of snapshots
#   nselect = # of selected snapshots
#   snaps = list of snapshots
#   names = dictionary of column names:
#     key = "id", value = column # (0 to M-1)
#   tselect = class for time selection
#   eselect = class for element selection
#   etype = name of vector used as element type by viz extract (def = "")
#   Snap = one snapshot
#     time = time stamp
#     tselect = 0/1 if this snapshot selected
#     nflag = 0/1 if this snapshot has nodal coords (may be copies)
#     eflag = 0/n if this snapshot has elements (may be copies)
#             1 = tri, 2 = tet, 3 = square, 4 = cube
#     nvalueflag = 0/1 if this snapshot has nodal values
#     evalueflag = 0/1 if this snapshot has element values
#     nselect = # of selected elements in this snapshot
#     eselect[i] = 0/1 for each element
#     xlo,xhi,ylo,yhi,zlo,zhi = box bounds (float)
#     nnodes = # of nodes
#     nodes[i][j] = 2d array of floats, i = 0 to Nnod-1, j = 0 to Ncol
#     nelements = # of elements
#     elements[i][j] = 2d array of floats, i = 0 to Nel-1, j = 0 to Ncol
#     nnvalues = # of node values
#     nvalues[i][j] = 2d array of floats, i = 0 to Nnod-1, j = 0 to Ncol
#     nevalues = # of node values
#     evalues[i][j] = 2d array of floats, i = 0 to Nel-1, j = 0 to Ncol

# Imports and external programs

import sys, commands, re, glob, types
from os import popen
from math import *             # any function could be used by set()

try:
    import numpy as np
    oldnumeric = False
except:
    import Numeric as np
    oldnumeric = True

try: from DEFAULTS import PIZZA_GUNZIP
except: PIZZA_GUNZIP = "gunzip"

# Class definition

class mdump:

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.snaps = []
    self.nsnaps = self.nselect = 0
    self.names = {}
    self.tselect = tselect(self)
    self.eselect = eselect(self)
    self.etype = ""

    # flist = list of all dump file names

    words = list[0].split()
    self.flist = []
    for word in words: self.flist += glob.glob(word)
    if len(self.flist) == 0 and len(list) == 1:
      raise StandardError,"no dump file specified"
    
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

    # sort entries by timestep, cull and combine duplicates

    self.snaps.sort(self.compare_time)
    self.cull()

    # sort all node, element, nvalue, evalue arrays by ID

    for snap in self.snaps:
      if snap.nflag:
        array = snap.nodes
        ids = array[:,0]
        ordering = np.argsort(ids)
        for i in xrange(len(array[0])):
          array[:,i] = np.take(array[:,i],ordering)
      if snap.eflag:
        array = snap.elements
        ids = array[:,0]
        ordering = np.argsort(ids)
        for i in xrange(len(array[0])):
          array[:,i] = np.take(array[:,i],ordering)
      if snap.nvalueflag:
        array = snap.nvalues
        ids = array[:,0]
        ordering = np.argsort(ids)
        for i in xrange(len(array[0])):
          array[:,i] = np.take(array[:,i],ordering)
      if snap.evalueflag:
        array = snap.evalues
        ids = array[:,0]
        ordering = np.argsort(ids)
        for i in xrange(len(array[0])):
          array[:,i] = np.take(array[:,i],ordering)

    # reference definitions of nodes and elements in previous timesteps

    self.reference()
    
    self.nsnaps = len(self.snaps)
    print "read %d snapshots" % self.nsnaps

    # select all timesteps and elements

    self.tselect.all()

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

    # select the new snapshot with all its elements

    self.snaps.append(snap)
    snap = self.snaps[self.nsnaps]
    snap.tselect = 1
    snap.nselect = snap.nelements
    if snap.eflag:
      for i in xrange(snap.nelements): snap.eselect[i] = 1
    self.nsnaps += 1
    self.nselect += 1

    return snap.time

  # --------------------------------------------------------------------
  # read a single snapshot from file f
  # return snapshot or 0 if failed

  def read_snapshot(self,f):
    try:
      snap = Snap()
      item = f.readline()
      snap.time = int(f.readline())
      snap.nflag = snap.eflag = snap.nvalueflag = snap.evalueflag = 0
      str = f.readline()
      if "NUMBER OF NODES" in str: snap.nflag = 1
      elif "NUMBER OF TRIANGLES" in str: snap.eflag = 1
      elif "NUMBER OF TETS" in str: snap.eflag = 2
      elif "NUMBER OF SQUARES" in str: snap.eflag = 3
      elif "NUMBER OF CUBES" in str: snap.eflag = 4
      elif "NUMBER OF NODE VALUES" in str: snap.nvalueflag = 1
      elif "NUMBER OF ELEMENT VALUES" in str: snap.evalueflag = 1
      else: raise StandardError,"unrecognized snapshot in dump file"
      n = int(f.readline())

      if snap.eflag: snap.eselect = np.zeros(n)

      if snap.nflag:
        item = f.readline()
        words = f.readline().split()
        snap.xlo,snap.xhi = float(words[0]),float(words[1])
        words = f.readline().split()
        snap.ylo,snap.yhi = float(words[0]),float(words[1])
        words = f.readline().split()
        snap.zlo,snap.zhi = float(words[0]),float(words[1])
        
      item = f.readline()
      if n:
        words = f.readline().split()
        ncol = len(words)
        for i in xrange(1,n):
          words += f.readline().split()
        floats = map(float,words)
        if oldnumeric: values = np.zeros((n,ncol),np.Float)
        else: values = np.zeros((n,ncol),np.float)
        start = 0
        stop = ncol
        for i in xrange(n):
          values[i] = floats[start:stop]
          start = stop
          stop += ncol
      else: values = None

      if snap.nflag:
        snap.nodes = values; snap.nnodes = n
      elif snap.eflag:
        snap.elements = values; snap.nelements = n
      elif snap.nvalueflag:
        snap.nvalues = values; snap.nnvalues = n
      elif snap.evalueflag:
        snap.evalues = values; snap.nevalues = n
      return snap
    except:
      return 0

  # --------------------------------------------------------------------
  # map atom column names
  
  def map(self,*pairs):
    if len(pairs) % 2 != 0:
      raise StandardError, "mdump map() requires pairs of mappings"
    for i in range(0,len(pairs),2):
      j = i + 1
      self.names[pairs[j]] = pairs[i]-1

  # delete unselected snapshots

  # --------------------------------------------------------------------

  def delete(self):
    ndel = i = 0
    while i < self.nsnaps:
      if not self.snaps[i].tselect:
        del self.snaps[i]
        self.nsnaps -= 1
        ndel += 1
      else: i += 1
    print "%d snapshots deleted" % ndel
    print "%d snapshots remaining" % self.nsnaps

  # --------------------------------------------------------------------
  # sort elements by ID in all selected timesteps or one timestep

  def sort(self,*list):
    if len(list) == 0:
      print "Sorting selected snapshots ..."
      id = self.names["id"]
      for snap in self.snaps:
        if snap.tselect: self.sort_one(snap,id)
    else:
      i = self.findtime(list[0])
      id = self.names["id"]
      self.sort_one(self.snaps[i],id)


  # --------------------------------------------------------------------
  # sort a single snapshot by ID column

  def sort_one(self,snap,id):
    atoms = snap.atoms
    ids = atoms[:,id]
    ordering = np.argsort(ids)
    for i in xrange(len(atoms[0])):
      atoms[:,i] = np.take(atoms[:,i],ordering)

  # --------------------------------------------------------------------
  # return vector of selected snapshot time stamps

  def time(self):
    vec = self.nselect * [0]
    i = 0
    for snap in self.snaps:
      if not snap.tselect: continue
      vec[i] = snap.time
      i += 1
    return vec

  # --------------------------------------------------------------------
  # extract vector(s) of values for selected elements at chosen timestep

  def vecs(self,n,*list):
    snap = self.snaps[self.findtime(n)]
    if not snap.evalues:
      raise StandardError, "snapshot has no element values"
      
    if len(list) == 0:
      raise StandardError, "no columns specified"
    columns = []
    values = []
    for name in list:
      columns.append(self.names[name])
      values.append(snap.nselect * [0])
    ncol = len(columns)

    m = 0
    for i in xrange(len(snap.evalues)):
      if not snap.eselect[i]: continue
      for j in xrange(ncol):
        values[j][m] = snap.evalues[i][columns[j]]
      m += 1

    if len(list) == 1: return values[0]
    else: return values

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
  # if have same timestamp, combine them if internal flags are different
  
  def cull(self):
    i = 1
    while i < len(self.snaps):
      if self.snaps[i].time == self.snaps[i-1].time:
        if self.snaps[i].nflag:
          if not self.snaps[i-1].nflag:
            self.snaps[i-1].nflag = 1
            self.snaps[i-1].nnodes = self.snaps[i].nnodes
            self.snaps[i-1].nodes = self.snaps[i].nodes
            self.snaps[i-1].xlo = self.snaps[i].xlo
            self.snaps[i-1].xhi = self.snaps[i].xhi
            self.snaps[i-1].ylo = self.snaps[i].ylo
            self.snaps[i-1].yhi = self.snaps[i].yhi
            self.snaps[i-1].zlo = self.snaps[i].zlo
            self.snaps[i-1].zhi = self.snaps[i].zhi
        elif self.snaps[i].eflag:
          if not self.snaps[i-1].eflag:
            self.snaps[i-1].eflag = self.snaps[i].eflag
            self.snaps[i-1].nelements = self.snaps[i].nelements
            self.snaps[i-1].elements = self.snaps[i].elements
            self.snaps[i-1].eselect = self.snaps[i].eselect
        elif self.snaps[i].nvalueflag:
          if not self.snaps[i-1].nvalueflag:
            self.snaps[i-1].nvalueflag = 1
            self.snaps[i-1].nnvalues = self.snaps[i].nnvalues
            self.snaps[i-1].nvalues = self.snaps[i].nvalues
        elif self.snaps[i].evalueflag:
          if not self.snaps[i-1].evalueflag:
            self.snaps[i-1].evalueflag = 1
            self.snaps[i-1].nevalues = self.snaps[i].nevalues
            self.snaps[i-1].evalues = self.snaps[i].evalues
        del self.snaps[i]
      else:
        i += 1

  # --------------------------------------------------------------------
  # insure every snapshot has node and element connectivity info
  # if not, point it at most recent shapshot that does
  
  def reference(self):
    for i in xrange(len(self.snaps)):
      if not self.snaps[i].nflag:
        for j in xrange(i,-1,-1):
          if self.snaps[j].nflag:
            self.snaps[i].nflag = self.snaps[j].nflag
            self.snaps[i].nnodes = self.snaps[j].nnodes
            self.snaps[i].nodes = self.snaps[j].nodes
            self.snaps[i].xlo = self.snaps[j].xlo
            self.snaps[i].xhi = self.snaps[j].xhi
            self.snaps[i].ylo = self.snaps[j].ylo
            self.snaps[i].yhi = self.snaps[j].yhi
            self.snaps[i].zlo = self.snaps[j].zlo
            self.snaps[i].zhi = self.snaps[j].zhi
            break
      if not self.snaps[i].nflag:
        raise StandardError,"no nodal coords found in previous snapshots"
      if not self.snaps[i].eflag:
        for j in xrange(i,-1,-1):
          if self.snaps[j].eflag:
            self.snaps[i].eflag = self.snaps[j].eflag
            self.snaps[i].nelements = self.snaps[j].nelements
            self.snaps[i].elements = self.snaps[j].elements
            self.snaps[i].eselect = self.snaps[j].eselect
            break
      if not self.snaps[i].eflag:
        raise StandardError,"no elem connections found in previous snapshots"

  # --------------------------------------------------------------------
  # iterate over selected snapshots

  def iterator(self,flag):
    start = 0
    if flag: start = self.iterate + 1
    for i in xrange(start,self.nsnaps):
      if self.snaps[i].tselect:
        self.iterate = i
        return i,self.snaps[i].time,1
    return 0,0,-1
  
  # --------------------------------------------------------------------
  # return list of triangles to viz for snapshot isnap
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
    if self.etype == "": type = -1
    else: type = self.names[self.etype]

    atoms = []
    bonds = []

    # create triangle list from all elements
    # for type, either use element type (-1) or user-defined column in evalues
    
    tris = []
    nodes = snap.nodes
    for i in xrange(snap.nelements):
      if not snap.eselect[i]: continue
      element = snap.elements[i]
      if snap.evalueflag: evalue = snap.evalues[i]
      else: evalues = []

      # single tri, normal = up
      
      if snap.eflag == 1:
        v1 = nodes[int(element[2])-1][2:5].tolist()
        v2 = nodes[int(element[3])-1][2:5].tolist()
        v3 = nodes[int(element[4])-1][2:5].tolist()
        list = v1 + v2 + v3
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)

      # single tet, convert to 4 tris, normals = out
      
      elif snap.eflag == 2:
        v1 = nodes[int(element[2])-1][2:5].tolist()
        v2 = nodes[int(element[3])-1][2:5].tolist()
        v3 = nodes[int(element[4])-1][2:5].tolist()
        v4 = nodes[int(element[5])-1][2:5].tolist()
        list = v1 + v2 + v4
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v2 + v3 + v4
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v4 + v3
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v3 + v2
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)

      # single square, convert to 2 tris, normals = up
      
      elif snap.eflag == 3:
        v1 = nodes[int(element[2])-1][2:5].tolist()
        v2 = nodes[int(element[3])-1][2:5].tolist()
        v3 = nodes[int(element[4])-1][2:5].tolist()
        v4 = nodes[int(element[5])-1][2:5].tolist()
        list = v1 + v2 + v3
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v3 + v4
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)

      # single cube, convert to 12 tris, normals = out
        
      elif snap.eflag == 4:
        v1 = nodes[int(element[2])-1][2:5].tolist()
        v2 = nodes[int(element[3])-1][2:5].tolist()
        v3 = nodes[int(element[4])-1][2:5].tolist()
        v4 = nodes[int(element[5])-1][2:5].tolist()
        v5 = nodes[int(element[6])-1][2:5].tolist()
        v6 = nodes[int(element[7])-1][2:5].tolist()
        v7 = nodes[int(element[8])-1][2:5].tolist()
        v8 = nodes[int(element[9])-1][2:5].tolist()
        list = v1 + v3 + v2   # lower z face
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v4 + v3
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v5 + v6 + v7   # upper z face
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v5 + v7 + v8
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v2 + v6  # lower y face
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v6 + v5
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v4 + v7 + v3  # upper y face
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v4 + v8 + v7
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v8 + v4  # lower x face
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v1 + v5 + v8
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v2 + v3 + v7  # upper x face
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
        list = v2 + v7 + v6
        n = normal(list[0:3],list[3:6],list[6:9])
        if type == -1: tris.append([element[0],element[1]] + list + n)
        else: tris.append([element[0],evalue[type]] + list + n)
      
    lines = []

    return time,box,atoms,bonds,tris,lines

  # --------------------------------------------------------------------
  # return lists of node/element info for snapshot isnap
  # if called with flag, then index is timestep, so convert to snapshot index
  
  def mviz(self,index,flag=0):
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
    nvalues = []
    if snap.nvalueflag: nvalues = snap.nvalues
    evalues = []
    if snap.nvalueflag: evalues = snap.evalues
    
    return time,box,snap.nodes,snap.elements,nvalues,evalues

  # --------------------------------------------------------------------

  def findtime(self,n):
    for i in xrange(self.nsnaps):
      if self.snaps[i].time == n: return i
    raise StandardError, "no step %d exists" % n

  # --------------------------------------------------------------------
  # return maximum box size across all selected snapshots

  def maxbox(self):
    xlo = ylo = zlo = None
    xhi = yhi = zhi = None
    for snap in self.snaps:
      if not snap.tselect: continue
      if xlo == None or snap.xlo < xlo: xlo = snap.xlo
      if xhi == None or snap.xhi > xhi: xhi = snap.xhi
      if ylo == None or snap.ylo < ylo: ylo = snap.ylo
      if yhi == None or snap.yhi > yhi: yhi = snap.yhi
      if zlo == None or snap.zlo < zlo: zlo = snap.zlo
      if zhi == None or snap.zhi > zhi: zhi = snap.zhi
    return [xlo,ylo,zlo,xhi,yhi,zhi]
    
  # --------------------------------------------------------------------

  def compare_atom(self,a,b):
    if a[0] < b[0]:
      return -1
    elif a[0] > b[0]:
      return 1
    else:
      return 0  

# --------------------------------------------------------------------
# one snapshot

class Snap:
  pass

# --------------------------------------------------------------------
# time selection class

class tselect:

  def __init__(self,data):
    self.data = data
    
  # --------------------------------------------------------------------

  def all(self):
    data = self.data
    for snap in data.snaps:
      snap.tselect = 1
    data.nselect = len(data.snaps)
    data.eselect.all()
    print "%d snapshots selected out of %d" % (data.nselect,data.nsnaps)

  # --------------------------------------------------------------------

  def one(self,n):
    data = self.data
    for snap in data.snaps:
      snap.tselect = 0
    i = data.findtime(n)
    data.snaps[i].tselect = 1
    data.nselect = 1
    data.eselect.all()
    print "%d snapshots selected out of %d" % (data.nselect,data.nsnaps)

  # --------------------------------------------------------------------

  def none(self):
    data = self.data
    for snap in data.snaps:
      snap.tselect = 0
    data.nselect = 0
    print "%d snapshots selected out of %d" % (data.nselect,data.nsnaps)

  # --------------------------------------------------------------------

  def skip(self,n):
    data = self.data
    count = n-1
    for snap in data.snaps:
      if not snap.tselect: continue
      count += 1
      if count == n:
        count = 0
        continue
      snap.tselect = 0
      data.nselect -= 1
    data.eselect.all()
    print "%d snapshots selected out of %d" % (data.nselect,data.nsnaps)
  
  # --------------------------------------------------------------------

  def test(self,teststr):
    data = self.data
    snaps = data.snaps
    cmd = "flag = " + teststr.replace("$t","snaps[i].time")
    ccmd = compile(cmd,'','single')
    for i in xrange(data.nsnaps):
      if not snaps[i].tselect: continue
      exec ccmd
      if not flag:
        snaps[i].tselect = 0
        data.nselect -= 1
    data.eselect.all()
    print "%d snapshots selected out of %d" % (data.nselect,data.nsnaps)

# --------------------------------------------------------------------
# element selection class

class eselect:

  def __init__(self,data):
    self.data = data

  # --------------------------------------------------------------------

  def all(self,*args):
    data = self.data
    if len(args) == 0:                           # all selected timesteps
      for snap in data.snaps:
        if not snap.tselect: continue
        for i in xrange(snap.nelements): snap.eselect[i] = 1
        snap.nselect = snap.nelements
    else:                                        # one timestep
      n = data.findtime(args[0])
      snap = data.snaps[n]
      for i in xrange(snap.nelements): snap.eselect[i] = 1
      snap.nselect = snap.nelements

  # --------------------------------------------------------------------

  def test(self,teststr,*args):
    data = self.data

    # replace all $var with snap.atoms references and compile test string
    
    pattern = "\$\w*"
    list = re.findall(pattern,teststr)
    for item in list:
      name = item[1:]
      column = data.names[name]
      insert = "snap.atoms[i][%d]" % column
      teststr = teststr.replace(item,insert)
    cmd = "flag = " + teststr
    ccmd = compile(cmd,'','single')

    if len(args) == 0:                           # all selected timesteps
      for snap in data.snaps:
        if not snap.tselect: continue
        for i in xrange(snap.nelements):
          if not snap.eselect[i]: continue
          exec ccmd
          if not flag:
            snap.eselect[i] = 0
            snap.nselect -= 1
      for i in xrange(data.nsnaps):
        if data.snaps[i].tselect:
          print "%d atoms of %d selected in first step %d" % \
                (data.snaps[i].nselect,data.snaps[i].nelements,data.snaps[i].time)
          break
      for i in xrange(data.nsnaps-1,-1,-1):
        if data.snaps[i].tselect:
          print "%d atoms of %d selected in last step %d" % \
                (data.snaps[i].nselect,data.snaps[i].nelements,data.snaps[i].time)
          break

    else:                                        # one timestep
      n = data.findtime(args[0])
      snap = data.snaps[n]
      for i in xrange(snap.nelements):
        if not snap.eselect[i]: continue
        exec ccmd
        if not flag:
          snap.eselect[i] = 0
          snap.nselect -= 1

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
