#!/usr/bin/python 

# Script:  bond_distribute.py 
# Purpose: binned bond length distributions by bond type
# Syntax:  bond_distribute.py datafile nbin rmin rmax outfile files ...
#          datafile = lammps data file
#          nbin = # of bins per bond type
#          rmin = min expected bond length
#          rmax = max expected bond length
#          outfile = file to write stats to
#          files = series of dump files
# Example: bond_distribute.py pore.data 1000 1.5 1.85 bonds.out pore.dump.1
# Author:  Paul Crozier (Sandia)

# enable script to run from Python directly w/out Pizza.py

from argparse import ArgumentParser
import numpy
from dump import dump
from data import data
from math import sqrt

# main script
def compute_bond_distribution():
 parser = ArgumentParser(description='A python script to bond distribution')

 parser.add_argument("input_data_file", help="The name of the lammps input data file")
 parser.add_argument("nbins", type=int, help="# of bins per bond type")
 parser.add_argument("rmin", type=float, help="min expected bond length")
 parser.add_argument("rmax", type=float, help="max expected bond length")
 parser.add_argument("output_file", help="The name of the file to write stats")
 parser.add_argument("dump_files", nargs='+', help="series of dump files")

 args = parser.parse_args()

 dt = data(args.input_data_file)	
 nbins = args.nbins
 rmin = args.rmin
 rmax = args.rmax
 files = ' '.join(args.dump_files[:])

 # get the bonds from the data file

 bond = dt.get("Bonds")
 nbonds = len(bond)
 btype = numpy.asarray(bond[:][1], dtype=int) - 1
 iatom = numpy.asarray(bond[:][2], dtype=int) - 1
 jatom = numpy.asarray(bond[:][3], dtype=int) - 1 

 ntypes = 0
 for i in xrange(nbonds): ntypes = max(bond[i][1],ntypes)
 ntypes = int(ntypes)
 ncount = numpy.zeros(ntypes, dtype=int)
 bin = nbins * [0]
 for i in xrange(nbins): 
   bin[i] = ntypes * [0] 

 # read snapshots one-at-a-time
 d = dump(files,0)
 d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

 while True:
   time = d.next()
   if time == -1:
     break
   
   xprd = d.snaps[-1].xhi - d.snaps[-1].xlo
   yprd = d.snaps[-1].yhi - d.snaps[-1].ylo
   zprd = d.snaps[-1].zhi - d.snaps[-1].zlo
  
   d.unscale() 
   d.sort()
   x,y,z = d.vecs(time,"x","y","z")
  
   x = numpy.asarray(x)
   y = numpy.asarray(y)
   z = numpy.asarray(z)
  
   delx = x[jatom[:]] - x[iatom[:]]
   dely = y[jatom[:]] - y[iatom[:]]
   delz = z[jatom[:]] - z[iatom[:]]
  
   delx -= xprd*numpy.rint((x[jatom[:]] - x[iatom[:]])/xprd)
   dely -= yprd*numpy.rint((y[jatom[:]] - y[iatom[:]])/yprd)
   delz -= zprd*numpy.rint((z[jatom[:]] - z[iatom[:]])/zprd)
  
   r = numpy.sqrt(delx*delx + dely*dely + delz*delz)

   for i in xrange(nbonds):   
    
     ibin = int(nbins*(r[i] - rmin)/(rmax - rmin) + 0.5)
     if ((ibin >= 0) and (ibin <= nbins-1)): 
       bin[ibin][btype[i]] += nbins
       ncount[btype[i]] += 1
     else:
       print "Warning: bond distance outside specified range"
       print "Bond type:", btype[i]+1
       print "Bond number:", i
   print time,    
      
 print
 print("Printing bond distance normalized distribution to {}".format(args.output_file))
    
 with open(args.output_file,"w") as fp:
  rrange = rmax - rmin
  for i in xrange(nbins):
    print >>fp, rmin + rrange*float(i)/float(nbins),
    for j in xrange(ntypes):
      if (ncount[j] > 0):
        print >>fp, float(bin[i][j])/float(ncount[j])/rrange,
      else:
        print >>fp, 0.0,
    print >>fp

if __name__ == "__main__":
  compute_bond_distribution()
