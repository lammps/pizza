#!/usr/bin/python

# Script:  distance.py 
# Purpose: check if any atom pairs are closer than specified distance
# Syntax:  distance.py maxcut dump.file1 dump.file2 ...
#          maxcut = flag atoms which are less than this distance apart
# Example: distance.py 0.95 dump.file1
# Author:  Paul Crozier (Sandia)

# print out 2 atoms less than maxcut apart (with PBC)

from math import sqrt

if len(argv) < 3:
  raise StandardError,"distance.py maxcut dump.file1 dump.file2 ..."
  
maxcut = float(argv[1])
maxcut_sq = maxcut*maxcut
  
files = ' '.join(argv[2:])		        # dump files
d = dump(files,0)
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

while 1:
  time = d.next()
  if time < 0: break
  d.unscale(time)

  box = (d.snaps[-1].xlo,d.snaps[-1].ylo,d.snaps[-1].zlo,
         d.snaps[-1].xhi,d.snaps[-1].yhi,d.snaps[-1].zhi)
  d.aselect.all(time)
  id,type,x,y,z = d.vecs(time,"id","type","x","y","z")
  n = len(x)

  xprd = box[3] - box[0]
  yprd = box[4] - box[1]
  zprd = box[5] - box[2]
  
  for i in xrange(n):
    for j in xrange(i+1,n):
  
      delx = x[j] - x[i]
      if abs(delx) > 0.5*xprd:
        if delx < 0.0:
          delx += xprd
        else:
          delx -= xprd
      if (delx*delx < maxcut_sq):
          
        dely = y[j] - y[i]
        if abs(dely) > 0.5*yprd:
          if dely < 0.0:
            dely += yprd
          else:
            dely -= yprd
        if ((dely*dely + delx*delx) < maxcut_sq):
            
          delz = z[j] - z[i]
          if abs(delz) > 0.5*zprd:
            if delz < 0.0:
              delz += zprd
            else:
              delz -= zprd
              
          rsq = delx*delx + dely*dely + delz*delz
          
          if rsq < maxcut_sq: 
            print "time = %d, id[i] = %d, id[j] = %d," \
               " type[i] = %d, type[j] = %d, distance = %g" % \
              (time, id[i], id[j], type[i], type[j], sqrt(rsq))

  d.tselect.none()
  d.tselect.one(time)
  print "timestep = ", time
  d.delete()
