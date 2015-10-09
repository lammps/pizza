# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# histo tool

oneline = "Particle density histogram from a dump"

docstr = """
h = histo(d)                        d = dump/cdump object

x,y = h.compute('x',N,lo,hi)        compute histogram in dim with N bins

  lo/hi are optional, if not used histo will be over entire box
"""

# History
#   12/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   data = dump object

# Imports and external programs

# Class definition

class histo:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data
    
  # --------------------------------------------------------------------

  def compute(self,dim,nbins,lo=None,hi=None):
    if dim == 'x': idim = 2
    elif dim == 'y': idim = 3
    elif dim == 'z': idim = 4
    else: raise StandardError,"illegal dim value"

    y = nbins*[0]
    
    count = 0
    n = flag = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)

      if not lo:
        if dim == 'x':
          lo = box[0]
          hi = box[3]
        elif dim == 'y':
          lo = box[1]
          hi = box[4]
        elif dim == 'z':
          lo = box[2]
          hi = box[5]

      delta = (hi-lo) / nbins;
      invdelta = 1.0/delta
      
      for atom in atoms:
        coord = atom[idim]
        ibin = int((coord-lo) * invdelta)
        if ibin < 0 or ibin >= nbins: continue
        y[ibin] += 1
        count += 1
        
      n += 1

    x = nbins*[0]
    for i in xrange(nbins): x[i] = (i+0.5)*delta
      
    print "histogram snapshots = ",n
    print "histogram counts (per snap) = %d (%g)" % (count,float(count)/n)
    print "histogram bounds = ",lo,hi
    return x,y
