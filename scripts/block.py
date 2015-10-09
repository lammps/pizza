#!/usr/bin/python

# Script:  block.py
# Purpose: thermodynamic block averages from LAMMPS log files
# Syntax:  block.py nblocks nskip log.1 log.2 ...
#          nblocks = number of blocks for block averages 
#          nskip = skip first nskip samples in the log file(s)
#          files = series of log files (LAMMPS thermo one style only)
# Example: block.py 10 0 log.*
# Author:  Paul Crozier (Sandia)

# enable script to run from Python directly w/out Pizza.py

import sys
from log import log
if not globals().has_key("argv"): argv = sys.argv

# main script

if len(argv) < 4:
  raise StandardError, "Syntax: block.py nblocks nskip log.1 log.2 ..."

nblocks = int(argv[1])
nskip = int(argv[2])
files = ' '.join(argv[3:])

l = log(files)

# assumes LAMMPS thermo style one

e_bond,e_pair,e_total,pressure,step,temperature = \
  l.get("E_bond","E_pair", "E_total", "Pressure", "Step", "Temperature")
if "Volume" in l.names:
  vflag = 1
  volume = l.get("Volume")
else: vflag = 0

print "Computing %g block averages" % nblocks
print "Skipping first %g samples" % nskip

n = len(step)
n_per_block = (n-nskip)/nblocks
k = nskip
temperature_ave = []
e_bond_ave = []
e_pair_ave = [] 
e_total_ave = []
pressure_ave = []
volume_ave = []

print
print " Block    Samples Temperature      E_bond      E_pair",  \
  "    E_total    Pressure      Volume"

for i in xrange(nblocks):
  temperature_sum = 0
  e_bond_sum = 0
  e_pair_sum = 0
  e_total_sum = 0
  pressure_sum = 0
  volume_sum = 0
  samples = 0
  for j in xrange(n_per_block): 
    temperature_sum += temperature[k]
    e_bond_sum += e_bond[k]
    e_pair_sum += e_pair[k]
    e_total_sum += e_total[k]
    pressure_sum += pressure[k]
    if vflag: volume_sum += volume[k]
    samples += 1
    k += 1
  temperature_ave.append(temperature_sum/samples)
  e_bond_ave.append(e_bond_sum/samples)
  e_pair_ave.append(e_pair_sum/samples)
  e_total_ave.append(e_total_sum/samples)
  pressure_ave.append(pressure_sum/samples)
  volume_ave.append(volume_sum/samples)
  print " %5i %10i %11.2f %11.2f %11.2f %11.2f %11.2f %11.2f" %  \
    (i+1, samples, temperature_ave[i], e_bond_ave[i], e_pair_ave[i], \
    e_total_ave[i], pressure_ave[i], volume_ave[i])

grand_ave_temperature = 0.0
grand_ave_e_bond = 0.0
grand_ave_e_pair = 0.0
grand_ave_e_total = 0.0
grand_ave_pressure = 0.0
grand_ave_volume = 0.0 
stdev_temperature = 0.0
stdev_e_bond = 0.0
stdev_e_pair = 0.0
stdev_e_total = 0.0
stdev_pressure = 0.0
stdev_volume = 0.0 
for i in xrange(nblocks):
  grand_ave_temperature += temperature_ave[i]
  grand_ave_e_bond += e_bond_ave[i]
  grand_ave_e_pair += e_pair_ave[i]
  grand_ave_e_total += e_total_ave[i]
  grand_ave_pressure += pressure_ave[i]
  grand_ave_volume += volume_ave[i] 
grand_ave_temperature /= nblocks
grand_ave_e_bond /= nblocks
grand_ave_e_pair /= nblocks
grand_ave_e_total /= nblocks
grand_ave_pressure /= nblocks
grand_ave_volume /= nblocks 
for i in xrange(nblocks):
  stdev_temperature += (temperature_ave[i] - grand_ave_temperature)*  \
                       (temperature_ave[i] - grand_ave_temperature)
  stdev_e_bond += (e_bond_ave[i] - grand_ave_e_bond)*  \
                  (e_bond_ave[i] - grand_ave_e_bond)
  stdev_e_pair += (e_pair_ave[i] - grand_ave_e_pair)*  \
                  (e_pair_ave[i] - grand_ave_e_pair) 
  stdev_e_total += (e_total_ave[i] - grand_ave_e_total)*  \
                   (e_total_ave[i] - grand_ave_e_total)
  stdev_pressure += (pressure_ave[i] - grand_ave_pressure)*  \
                    (pressure_ave[i] - grand_ave_pressure)
  stdev_volume += (volume_ave[i] - grand_ave_volume)*  \
                  (volume_ave[i] - grand_ave_volume)
from math import *
stdev_temperature = sqrt(stdev_temperature/(nblocks-1))
stdev_e_bond = sqrt(stdev_e_bond/(nblocks-1))
stdev_e_pair = sqrt(stdev_e_pair/(nblocks-1))
stdev_e_total = sqrt(stdev_e_total/(nblocks-1))
stdev_pressure = sqrt(stdev_pressure/(nblocks-1))
stdev_volume = sqrt(stdev_volume/(nblocks-1)) 

print " ====================================================",  \
  "==================================="
  
print " Ave.             %11.2f %11.2f %11.2f %11.2f %11.2f %11.2f" % \
    (grand_ave_temperature, grand_ave_e_bond, grand_ave_e_pair, \
    grand_ave_e_total, grand_ave_pressure, grand_ave_volume)
    
print " Stdev            %11.2f %11.2f %11.2f %11.2f %11.2f %11.2f" % \
    (stdev_temperature, stdev_e_bond, stdev_e_pair, \
    stdev_e_total, stdev_pressure, stdev_volume)
