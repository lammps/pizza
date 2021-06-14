#!/usr/bin/python

# Script:  cluster.py
# Purpose: histogram of cluster size of type2 atoms near type1
# Syntax:  cluster.py type1 type2 cutoff nbin dump.1 dump.2 ...
#          type1,type2 = atom types (1-N)
#          cutoff = only consider atom pairs within distance cutoff
#          nbin = # of bins for histogram
#          files = series of dump files
# Example: cluster.py 1 5 2.0 10 dump.*
# Author:  Steve Plimpton (Sandia)

# for all snapshots, for each type1 atom, count # of type2 atoms within cutoff
# will be very slow (N^2) if are too many type1 atoms

# enable script to run from Python directly w/out Pizza.py

import sys
from dump import dump
from gnu import gnu
if not globals().has_key("argv"): argv = sys.argv

# main script

# function to compute distance sq between 2 atoms with PBC

def distance(atom1,atom2,box):
  x1 = atom1[2]
  y1 = atom1[3]
  z1 = atom1[4]
  x2 = atom2[2]
  y2 = atom2[3]
  z2 = atom2[4]
  
  delx = x2 - x1
  dely = y2 - y1
  delz = z2 - z1
  
  xprd = box[3] - box[0]
  yprd = box[4] - box[1]
  zprd = box[5] - box[2]

  if abs(delx) > 0.5*xprd:
    if delx < 0.0:
      delx += xprd
    else:
      delx -= xprd
  if abs(dely) > 0.5*yprd:
    if dely < 0.0:
      dely += yprd
    else:
      dely -= yprd
  if abs(delz) > 0.5*zprd:
    if delz < 0.0:
      delz += zprd
    else:
      delz -= zprd
      
  distsq = delx*delx + dely*dely + delz*delz
  return distsq

# main script

if len(argv) < 6:
  raise StandardError,"cluster.py type1 type2 cutoff nbin dump.1 dump.2 ..."

type1 = int(argv[1])
type2 = int(argv[2])
cutoff = float(argv[3])
nbin = int(argv[4])
files = ' '.join(argv[5:])

# read in dump file(s) and select only type1 or type2 atoms

d = dump(files)
d.aselect.test("$type == %d or $type == %d" % (type1,type2))

# loop over snapshots
# viz returns list of atoms in one snapshot

cutsq = cutoff*cutoff
cluster = nbin*[0]

print "Clustering ..."

flag = 0
while 1:
  which,time,flag = d.iterator(flag)
  if flag == -1: break
  time,box,atoms,bonds,tris = d.viz(which)
  print time,
  sys.stdout.flush()

  # loop over all type1 atoms
  
  n = len(atoms)
  for i in xrange(n):
    itype = atoms[i][1]
    if itype != type1: continue
    ncount = 0

    # loop over all type2 atoms
    # increment cluster count if distance is within cutoff
    
    for j in xrange(n):
      jtype = atoms[j][1]
      if jtype != type2 or i == j: continue 
      distsq = distance(atoms[i],atoms[j],box)
      if distsq < cutsq: ncount += 1

    # increment histogram count
    
    if ncount >= nbin: cluster[nbin-1] += 1
    else: cluster[ncount] += 1

print
print "Cluster size and count:"
for i in range(nbin): print i,cluster[i]

# comment out if don't want plot

#g = gnu()
#g.plot(range(nbin),cluster)
