#!/usr/bin/python 

# Script:  angle_distribute.py 
# Purpose: binned angle distributions by angle type
# Syntax:  angle_distribute.py datafile nbin theta_min theta_max outfile files ...
#          datafile = lammps data file
#          nbin = # of bins per angle type
#          theta_min = min expected angle
#          theta_max = max expected angle length
#          outfile = file to write stats to
#          files = series of dump files
# Example: angle_distribute.py pore.data 1000 110. 120. angles.out pore.dump.1
# Author:  Paul Crozier (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from dump import dump
from math import sqrt,acos,atan
if not globals().has_key("argv"): argv = sys.argv

# main script

if len(argv) < 7:
  raise StandardError, \
  "Syntax: angle_distribute.py datafile nbin theta_min theta_max outfile files ..."

dt = data(argv[1])	
nbins = int(argv[2])
theta_min = float(argv[3])
theta_max = float(argv[4])
outfile = argv[5]
files = ' '.join(argv[6:])

# get the angles from the data file

angle = dt.get("Angles")
nangles = len(angle)
atype = nangles * [0]
iatom = nangles * [0]
jatom = nangles * [0]
katom = nangles * [0]
for i in xrange(nangles):
  atype[i] = int(angle[i][1] - 1)
  iatom[i] = int(angle[i][2] - 1)
  jatom[i] = int(angle[i][3] - 1)
  katom[i] = int(angle[i][4] - 1)

ntypes = 0
for i in xrange(nangles): ntypes = max(angle[i][1],ntypes)
ntypes = int(ntypes)
ncount = ntypes * [0]
bin = nbins * [0]
for i in xrange(nbins): 
  bin[i] = ntypes * [0] 

# read snapshots one-at-a-time

d = dump(files,0)
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

PI = 4.0*atan(1.0)

while 1:
  time = d.next()
  if time == -1: break
   
  box = (d.snaps[-1].xlo,d.snaps[-1].ylo,d.snaps[-1].zlo,
         d.snaps[-1].xhi,d.snaps[-1].yhi,d.snaps[-1].zhi)
         
  xprd = box[3] - box[0]
  yprd = box[4] - box[1] 
  zprd = box[5] - box[2]
  
  d.unscale() 
  d.sort()
  x,y,z = d.vecs(time,"x","y","z")
   
  for i in xrange(nangles):
  
    delx1 = x[iatom[i]] - x[jatom[i]]
    dely1 = y[iatom[i]] - y[jatom[i]]
    delz1 = z[iatom[i]] - z[jatom[i]]
        
    if abs(delx1) > 0.5*xprd:
      if delx1 < 0.0:
        delx1 += xprd
      else:
        delx1 -= xprd
    if abs(dely1) > 0.5*yprd:
      if dely1 < 0.0:
        dely1 += yprd
      else:
        dely1 -= yprd
    if abs(delz1) > 0.5*zprd:
      if delz1 < 0.0:
        delz1 += zprd
      else:
        delz1 -= zprd
        
    r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1)
    
    delx2 = x[katom[i]] - x[jatom[i]]
    dely2 = y[katom[i]] - y[jatom[i]]
    delz2 = z[katom[i]] - z[jatom[i]]
        
    if abs(delx2) > 0.5*xprd:
      if delx2 < 0.0:
        delx2 += xprd
      else:
        delx2 -= xprd
    if abs(dely2) > 0.5*yprd:
      if dely2 < 0.0:
        dely2 += yprd
      else:
        dely2 -= yprd
    if abs(delz2) > 0.5*zprd:
      if delz2 < 0.0:
        delz2 += zprd
      else:
        delz2 -= zprd
        
    r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2)
    
    c = delx1*delx2 + dely1*dely2 + delz1*delz2
    c /= r1*r2
        
    if (c > 1.0): c = 1.0
    if (c < -1.0): c = -1.0
        
    theta = 180.0*acos(c)/PI
    
    ibin = int(nbins*(theta - theta_min)/(theta_max - theta_min) + 0.5)
    if ((ibin >= 0) and (ibin <= nbins-1)): 
      bin[ibin][atype[i]] += nbins
      ncount[atype[i]] += 1
    else:
      print "Warning: angle outside specified range"
      print "angle type:", atype[i]+1
      print "angle number:", i
  print time,    
      
print
print "Printing normalized angle distributions to",outfile
    
fp = open(outfile,"w")
theta_range = theta_max - theta_min
for i in xrange(nbins):
  print >>fp, theta_min + theta_range*float(i)/float(nbins), 
  for j in xrange(ntypes):
    if (ncount[j] > 0):
      print >>fp, float(bin[i][j])/float(ncount[j])/theta_range,
    else:
      print >>fp, 0.0,    
  print >>fp 
fp.close()
