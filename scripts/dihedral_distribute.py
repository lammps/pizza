#!/usr/bin/python 

# Script:  dihedral_distribute.py 
# Purpose: binned signed dihedral distributions by dihedral type
# Syntax:  dihedral_distribute.py datafile nbin theta_min theta_max outfile files ...
#          datafile = lammps data file
#          nbin = # of bins per dihedral type
#          theta_min = min expected angle
#          theta_max = max expected angle length
#          outfile = file to write stats to
#          files = series of dump files
# Example: dihedral_distribute.py pore.data 1000 110. 120. angles.out pore.dump.1
# Author:  Evangelos Voyiatzis (TU Darmstadt)

# enable script to run from Python directly w/out Pizza.py

import sys
import math 
from dump import dump

if not globals().has_key("argv"): argv = sys.argv

def PBC(distance,boxlength): # function to compute the minimum image 

  if abs(distance) > 0.5*xprd:
    if distance < 0.0:
      distance += boxlength
    else:
      distance -= boxlength

  return distance

def InnerProduct(Vector1,Vector2): # function to compute the inner product of two vectors

  return (Vector1[0]*Vector2[0] + Vector1[1]*Vector2[1] + Vector1[2]*Vector2[2])

def CrossProduct(Vector1,Vector2): # function to compute the Cross product of two vectors

  Vector3 = [0] *3
  Vector3[0] = Vector1[1]*Vector2[2] - Vector1[2]*Vector2[1]
  Vector3[1] = Vector1[2]*Vector2[0] - Vector1[0]*Vector2[2]
  Vector3[2] = Vector1[0]*Vector2[1] - Vector1[1]*Vector2[0]

  return Vector3

# main script

if len(argv) < 8:
  raise StandardError, \
  "Syntax: dihedral_distribute.py dihedrals/impropers datafile nbin theta_min theta_max outfile files ..."

dt = data(argv[2])	
nbins = int(argv[3])
theta_min = float(argv[4])
theta_max = float(argv[5])
outfile = argv[6]
files = ' '.join(argv[7:])

# get the angles from the data file
if == "dihedrals":
 angle = dt.get("Dihedrals")
else if == "impropers":
 angle = dt.get("Impropers") 
else: 
 raise StandardError, "The second keyword is neither 'dihedrals' nor 'impropers' "
 sys.exit()

nangles = len(angle)
atype = nangles * [0]
iatom = nangles * [0]
jatom = nangles * [0]
katom = nangles * [0]
latom = nangles * [0]
for i in xrange(nangles):
  atype[i] = int(angle[i][1] - 1)
  iatom[i] = int(angle[i][2] - 1)
  jatom[i] = int(angle[i][3] - 1)
  katom[i] = int(angle[i][4] - 1)
  latom[i] = int(angle[i][5] - 1)

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

  VectorA = [0] * 3
  VectorB = [0] * 3
  VectorC = [0] * 3
  
  for i in xrange(nangles):
    # first vector 
    VectorA[0] = PBC(x[jatom[i]] - x[iatom[i]],xprd)    
    VectorA[1] = PBC(y[jatom[i]] - y[iatom[i]],yprd)
    VectorA[1] = PBC(z[jatom[i]] - z[iatom[i]],zprd)

    # second vector          
    VectorB[0] = PBC(x[katom[i]] - x[jatom[i]],xprd)
    VectorB[1] = PBC(y[katom[i]] - y[jatom[i]],yprd)
    VectorB[2] = PBC(z[katom[i]] - z[jatom[i]],zprd)
    normVectorB = math.sqrt(InnerProduct(VectorB,VectorB))

    # third vector
    VectorC[0] = PBC(x[latom[i]] - x[katom[i]],xprd)
    VectorC[1] = PBC(y[latom[i]] - y[katom[i]],yprd)
    VectorC[2] = PBC(z[latom[i]] - z[katom[i]],zprd)
               
    # compute the signed dihedral angle
    HelpVector = CrossProduct(CrossProduct(VectorA,VectorB),CrossProduct(VectorB,VectorC))
    Argument1 = InnerProduct(HelpVector,VectorB/normVectorB)
    Argument2 = InnerProduct(CrossProduct(VectorA,VectorB),CrossProduct(VectorB,VectorC))
    theta = 180.0*math.atan2(Argument1,Argument2)/math.pi
           
    ibin = int(nbins*(theta - theta_min)/(theta_max - theta_min) + 0.5)
    if ((ibin >= 0) and (ibin <= nbins-1)): 
      bin[ibin][atype[i]] += nbins
      ncount[atype[i]] += 1
    else:
      print "Warning: dihedral outside specified range"
      print "dihedral type:", atype[i]+1
      print "dihedral number:", i
  print time,    
      
print
print "Printing normalized dihedral distributions to",outfile
    
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
