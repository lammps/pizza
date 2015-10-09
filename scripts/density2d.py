#!/usr/bin/python 

# Script:    density2d.py
# Purpose:   binned atom density by atom type
# Syntax:    density.py x/y/z nbin vmin vmax outfile files ...
#            x/y/z = get density distribution along this (vertical) axis,
#                    bining along the other two (planar) axes
#            nbin = # of bins along each of planar axes
#            vmin, vmax = min and max along the vertical axis,
#                    use 0 for both for the complete box
#            outfile = file to write flux stats to
#            files = series of dump files
# Example: density2d.py z 100 dens.out dump.*
# Modified from density.py:  David Hart (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
import numpy as np
from dump import dump
if not globals().has_key("argv"): argv = sys.argv

# main script

if len(argv) < 5:
    raise StandardError, "Syntax: density.py x/y/z nbin vmin vmax outfile files ..."

direction = argv[1]
nbins = int(argv[2])
if nbins < 1:
    nbins = 100
try:
    zmin = float(argv[3])
except ValueError:
    zmin = -np.inf
try:
    zmax = float(argv[4])
except ValueError:
    zmax = np.inf
if zmax == zmin:
    zmin = -np.inf
    zmax = np.inf
outfile = argv[5]
files = ' '.join(argv[6:])

# read snapshots one-at-a-time

d = dump(files,0)
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")
bidirect = 'x/y'
first = 1
nsnaps = 0
dx = 0
dy = 0
x0 = 0
y0 = 0
vol = 0

while 1:
    time = d.next()
    if time == -1: break

    if first:
        tmp,ntypes = d.minmax("type")
        ntypes = int(ntypes)
        bin = np.zeros(shape=(nbins,nbins,ntypes))
        first = 0
        
    box = (d.snaps[-1].xlo,d.snaps[-1].ylo,d.snaps[-1].zlo,
                 d.snaps[-1].xhi,d.snaps[-1].yhi,d.snaps[-1].zhi)
    vol = (box[3] - box[0]) * (box[4] - box[1]) * (box[5] - box[2]) 

    if direction == "z": 
        type,x,y,z = d.vecs(time,"type","x","y","z")
        bidirect = 'x/y'
        dx = box[3] - box[0]
        dy = box[4] - box[1]
        dz = box[5] - box[2]
        x0 = box[0] + float(dx)/float(nbins)/2.0
        y0 = box[1] + float(dy)/float(nbins)/2.0
        zmax = min(zmax,box[5])
        zmin = max(zmin,box[2])
    elif direction == "y":
        type,x,y,z = d.vecs(time,"type","x","z","y")
        bidirect = 'x/z'
        dx = box[3] - box[0]
        dy = box[5] - box[2]
        dz = box[4] - box[1]
        x0 = box[0] + float(dx)/float(nbins)/2.0
        y0 = box[2] + float(dy)/float(nbins)/2.0
        zmax = min(zmax,box[4])
        zmin = max(zmin,box[1])
    elif direction == "x": 
        type,x,y,z = d.vecs(time,"type","y","z","x")
        bidirect = 'y/z'
        dx = box[4] - box[1]
        dy = box[5] - box[2]
        dz = box[3] - box[0]
        x0 = box[1] + float(dx)/float(nbins)/2.0
        y0 = box[2] + float(dy)/float(nbins)/2.0
        zmax = min(zmax,box[3])
        zmin = max(zmin,box[0])
    vol = dx * dy * float(zmax - zmin)

    type = map(int,type)
    natoms = len(type)
    for i in xrange(natoms): type[i] -= 1
    
    for i in xrange(natoms):
        ibin = int(nbins*x[i])
        jbin = int(nbins*y[i])
        zloc = float(z[i])*float(dz)
        if zloc < zmin or zloc > zmax: continue
        if (ibin < 0): ibin = 0
        if (ibin > nbins-1): ibin = nbins - 1
        if (jbin < 0): jbin = 0
        if (jbin > nbins-1): jbin = nbins - 1
        bin[jbin][ibin][type[i]] += nbins*nbins/vol
    nsnaps += 1
    print time,

print 
print "Printing %s-mapped density distribution for %s-slice [%.2f,%.2f] in mol/L to %s" %(bidirect, direction, zmin, zmax, outfile)
conversion = 1660.53873  # convert from atoms/Angs^3 to mol/L
        
fp = open(outfile,"w")
# '''Uncomment for column headers. Commented for consistency with density.py'''
# print >>fp, " %8s  %8s " %('ra', 'rb'),
# for k in xrange(ntypes):
#     print >>fp, " %8s " %(k),
# print >>fp
for i in xrange(nbins):
    for j in xrange(nbins):
        print >>fp, " %8.3f  %8.3f " %(float(i)/float(nbins)*float(dx)+float(x0), float(j)/float(nbins)*float(dy)+float(y0)),
        for k in xrange(ntypes):
            print >>fp, " %8.3f " % (conversion*bin[j][i][k]/nsnaps),
        print >>fp 
fp.close()
