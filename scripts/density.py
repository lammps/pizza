#!/usr/bin/python 

# Script:  density.py
# Purpose: binned atom density by atom type
# Syntax:  density.py x/y/z nbin outfile files ...
#          x/y/z = get density distribution along this axis
#          nbin = # of bins in desired direction
#          outfile = file to write flux stats to
#          files = series of dump files
# Example: density.py z 100 dens.out dump.*
# Author:  Paul Crozier (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from dump import dump
if not globals().has_key("argv"): argv = sys.argv

# main script

if len(argv) < 5:
  raise StandardError, "Syntax: density.py x/y/z nbin outfile files ..."

direction = argv[1]
nbins = int(argv[2])
outfile = argv[3]
files = ' '.join(argv[4:])

# read snapshots one-at-a-time

d = dump(files,0)
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

first = 1
nsnaps = 0
while 1:
  time = d.next()
  if time == -1: break

  if first:
    tmp,ntypes = d.minmax("type")
    ntypes = int(ntypes)
    bin = nbins * [0]
    for i in xrange(nbins): bin[i] = ntypes * [0]
    first = 0
    
  box = (d.snaps[-1].xlo,d.snaps[-1].ylo,d.snaps[-1].zlo,
         d.snaps[-1].xhi,d.snaps[-1].yhi,d.snaps[-1].zhi)
  vol = (box[3] - box[0]) * (box[4] - box[1]) * (box[5] - box[2]) 

  if direction == "x": type,x = d.vecs(time,"type","x")
  elif direction == "y": type,x = d.vecs(time,"type","y")
  elif direction == "z": type,x = d.vecs(time,"type","z")
  
  type = map(int,type)
  natoms = len(type)
  for i in xrange(natoms): type[i] -= 1
  
  for i in xrange(natoms):
    ibin = int(nbins*x[i] + 0.5)
    if (ibin < 0): ibin += nbins
    if (ibin > nbins-1): ibin -= nbins
    bin[ibin][type[i]] += nbins/vol
  nsnaps += 1
  print time,
  
print
print "Printing ", direction, "-directional density distribution in mol/L to",outfile
conversion = 1660.53873              # convert from atoms/Angs^3 to mol/L
    
fp = open(outfile,"w")
for i in xrange(nbins):
  print >>fp, float(i)/float(nbins), 
  for j in xrange(ntypes):
    print >>fp, conversion*bin[i][j]/nsnaps,
  print >>fp 
fp.close()
