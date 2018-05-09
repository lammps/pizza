#!/usr/bin/python 
"""
Script:  bond_distribute.py 
Purpose: binned bond length distributions by bond type
Syntax:  bond_distribute.py datafile nbin rmin rmax outfile files ...
          datafile = lammps data file
          nbin = # of bins per bond type
          rmin = min expected bond length
          rmax = max expected bond length
          outfile = file to write stats to
          files = series of dump files
Example: bond_distribute.py pore.data 1000 1.5 1.85 bonds.out pore.dump.1
Author:  Paul Crozier (Sandia)

enable script to run from Python directly w/out Pizza.py
"""

from argparse import ArgumentParser
import numpy
from dump import dump
from data import data

# main script
def compute_bond_distribution():
 """The function doing the actual calculation of the bond distribution"""
 
 parser = ArgumentParser(description='A python script to compute bond distribution using pizza.py')

 parser.add_argument("input_data_file", help="The name of the lammps input data file")
 parser.add_argument("nbins", type=int, help="# of bins per bond type")
 parser.add_argument("rmin", type=float, help="min expected bond length")
 parser.add_argument("rmax", type=float, help="max expected bond length")
 parser.add_argument("output_file", help="The name of the file to write stats")
 parser.add_argument("dump_files", nargs='+', help="series of dump files")

 args = parser.parse_args()

 dt = data(args.input_data_file)	
 files = ' '.join(args.dump_files[:])

 # get the bonds from the data file
 bond = dt.get("Bonds")
 nbonds = len(bond)
 btype = numpy.asarray(bond[:][1], dtype=int) - 1
 iatom = numpy.asarray(bond[:][2], dtype=int) - 1
 jatom = numpy.asarray(bond[:][3], dtype=int) - 1 

 ntypes = int(numpy.max(btype)) + 1
 bin = numpy.zeros((args.nbins, ntypes), dtype=int)
 
 ncount = numpy.zeros(ntypes, dtype=int)
 for itype in range(0, ntypes):
  ncount[itype] = numpy.sum(btype==itype)
 
 # read snapshots one-at-a-time
 d = dump(files,0)
 d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

 nconfs = 0 
 while True:
   time = d.next()
   if time == -1:
     break
     
   nconfs += 1
   
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

   for itype in range(0, ntypes):
    xx, hist_edge = numpy.histogram(r[btype == itype],bins=args.nbins, range=(args.rmin, args.rmax))
    bin[:][itype] += xx
    
   print("{} ".format(time))    
      
 print("Printing bond distance normalized distribution to {}".format(args.output_file))
    
 with open(args.output_file,"w") as output_file:
  rrange = args.rmax - args.rmin
  for i in range(0, args.nbins):
    output_file.write("{} ".format(args.rmin + rrange*float(i)/float(args.nbins)))
    for j in range(0, ntypes):
      if (ncount[j] > 0):
        output_file.write("{} ".format(float(bin[i][j])/float(ncount[j]*nconfs)/rrange))
      else:
        output_file.write("{} ".format(0.0))
    output_file.write("\n")

if __name__ == "__main__":
  compute_bond_distribution()
