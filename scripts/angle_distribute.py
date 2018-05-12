#!/usr/bin/python 
"""
Script:  angle_distribute.py 
Purpose: binned angle distributions by angle type
Syntax:  angle_distribute.py datafile nbin theta_min theta_max outfile files ...
         datafile = lammps data file
         nbin = # of bins per angle type
         theta_min = min expected angle
         theta_max = max expected angle length
         outfile = file to write stats to
         files = series of dump files
Example: angle_distribute.py pore.data 1000 110. 120. angles.out pore.dump.1
Author:  Paul Crozier (Sandia)

enable script to run from Python directly w/out Pizza.py
"""

from argparse import ArgumentParser
import numpy
from dump import dump
from data import data

# main script
def compute_angle_distribution():
 """The function doing the actual calculation of the bond distribution"""
 
 parser = ArgumentParser(description='A python script to compute bond distribution using pizza.py')

 parser.add_argument("input_data_file", help="The name of the lammps input data file")
 parser.add_argument("nbins", type=int, help="# of bins per bond type")
 parser.add_argument("theta_min", type=float, help="min expected angle")
 parser.add_argument("theta_max", type=float, help="max expected angle")
 parser.add_argument("output_file", help="The name of the file to write stats")
 parser.add_argument("dump_files", nargs='+', help="series of dump files")

 args = parser.parse_args() 
  
 dt = data(args.input_data_file)	
 files = ' '.join(args.dump_files[:])

 # get the angles from the data file
 angle = dt.get("Angles")
 nangles = len(angle)
 atype = numpy.asarray(angle[:][1], dtype=int) - 1
 iatom = numpy.asarray(angle[:][2], dtype=int) - 1
 jatom = numpy.asarray(angle[:][3], dtype=int) - 1 
 katom = numpy.asarray(angle[:][4], dtype=int) - 1
 
 ntypes = int(numpy.max(btype)) + 1
 bin = numpy.zeros((args.nbins, ntypes), dtype=int)
 
 ncount = numpy.zeros(ntypes, dtype=int)
 for itype in range(0, ntypes):
  ncount[itype] = numpy.sum(atype==itype)

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
    
  delx = x[iatom[:]] - x[jatom[:]]
  dely = y[iatom[:]] - y[jatom[:]]
  delz = z[iatom[:]] - z[jatom[:]]
  
  delx -= xprd*numpy.rint((x[iatom[:]] - x[jatom[:]])/xprd)
  dely -= yprd*numpy.rint((y[iatom[:]] - y[jatom[:]])/yprd)
  delz -= zprd*numpy.rint((z[iatom[:]] - z[jatom[:]])/zprd)
        
  r1 = numpy.sqrt(delx*delx + dely*dely + delz*delz)
    
  delx2 = x[katom[:]] - x[jatom[:]]
  dely2 = y[katom[:]] - y[jatom[:]]
  delz2 = z[katom[:]] - z[jatom[:]]

  delx2 -= xprd*numpy.rint((x[katom[:]] - x[jatom[:]])/xprd)
  dely2 -= yprd*numpy.rint((y[katom[:]] - y[jatom[:]])/yprd)
  delz2 -= zprd*numpy.rint((z[katom[:]] - z[jatom[:]])/zprd)
        
  r2 = numpy.sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2)
    
  c = delx1*delx2 + dely1*dely2 + delz1*delz2
  c /= r1*r2
        
  if (c > 1.0): c = 1.0
  if (c < -1.0): c = -1.0
        
  theta = 180.0*numpy.arccos(c)/numpy.pi

  for itype in range(0, ntypes):
    xx, hist_edge = numpy.histogram(theta[atype == itype],bins=args.nbins, range=(args.theta_min, args.theta_max))
    if numpy.sum(xx) != ncount[itype]:
      print("Warning: angle distance outside specified range ")
      print("Angle type: {}".format(itype+1))
    bin[:][itype] += xx
  
  print("{} ".format(time))    
      
 print
 print("Printing normalized angle distributions to {}".format(args.output_file))
    
 with open(args.output_file,"w") as output_file:
  theta_range = args.theta_max - args.theta_min
  for i in range(0, args.nbins):
    output_file.write("{} ".format(theta_min + theta_range*float(i)/float(nbins)))
    for j in range(0, ntypes):
      if (ncount[j] > 0):
        output_file.write("{} ".format(float(bin[i][j])/float(ncount[j]*nconfs)/theta_range))
      else:
        output_file.write("{} ".format(0.0))    
    output_file.write("\n") 

if __name__ == "__main__":
  compute_angle_distribution()
