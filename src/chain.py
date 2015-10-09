# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# chain tool

oneline = "Create bead-spring chains for LAMMPS input"

docstr = """
c = chain(N,rho)            setup box with N monomers at reduced density rho
c = chain(N,rho,1,1,2)	    x,y,z = aspect ratio of box (def = 1,1,1)

c.seed = 48379              set random # seed (def = 12345)
c.mtype = 2    		    set type of monomers (def = 1)
c.btype = 1           	    set type of bonds (def = 1)
c.blen = 0.97               set length of bonds (def = 0.97)
c.dmin = 1.02               set min dist from i-1 to i+1 site (def = 1.02)

c.id = "chain"              set molecule ID to chain # (default)
c.id = "end1"               set molecule ID to count from one end of chain
c.id = "end2"               set molecule ID to count from either end of chain

c.build(100,10)		    create 100 chains, each of length 10

  can be invoked multiple times interleaved with different settings
  must fill box with total of N monomers
  
c.write("data.file")        write out all built chains to LAMMPS data file
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   n = number of monomers
#   rhostar = reduced density
#   seed = 12345
#   mtype = type of monomers
#   btype = type of bonds
#   blen = length of bonds
#   dmin = minimum distance from i-1 to i+1
#   id = "chain","end1",or "end2"
#   atoms = list of atoms
#   bonds = list of bonds
#   xprd,yprd,zprd = x,y,z box size
#   xlo,ylo,zlo = -xyz prd / 2
#   xhi,yhi,zhi = x,y,zprd /2

# Imports and external programs

import math
from data import data

# Class definition

class chain:
  
  # --------------------------------------------------------------------

  def __init__(self,n,rhostar,*list):
    self.n = n
    self.rhostar = rhostar
    xaspect = yaspect = zaspect = 1.0
    if len(list):
      xaspect = list[0]
      yaspect = list[1]
      zaspect = list[2]
    self.seed = 12345
    self.mtype = 1
    self.btype = 1
    self.blen = 0.97
    self.dmin = 1.02
    self.id = "chain"
    self.atoms = []
    self.bonds = []

    volume = n/rhostar
    prd = pow(volume/xaspect/yaspect/zaspect,1.0/3.0)
    self.xprd = xaspect * prd
    self.xlo = -self.xprd/2.0
    self.xhi = self.xprd/2.0
    self.yprd = yaspect * prd
    self.ylo = -self.yprd/2.0
    self.yhi = self.yprd/2.0
    self.zprd = zaspect * prd
    self.zlo = -self.zprd/2.0
    self.zhi = self.zprd/2.0

    print "Simulation box: %g by %g by %g" % (self.xprd,self.yprd,self.zprd)

  # --------------------------------------------------------------------

  def build(self,n,nper):
    for ichain in xrange(n):
      atoms = []
      bonds = []
      id_atom_prev = id_mol_prev = id_bond_prev = 0
      if len(self.atoms):
        id_atom_prev = self.atoms[-1][0]
        id_mol_prev = self.atoms[-1][1]
      if len(self.bonds):
        id_bond_prev = self.bonds[-1][0]

      for imonomer in xrange(nper):
        if imonomer == 0:
          x = self.xlo + self.random()*self.xprd
          y = self.ylo + self.random()*self.yprd
          z = self.zlo + self.random()*self.zprd
	  ix = iy = iz = 0
        else:
          restriction = True
          while restriction:
            rsq = 2.0
            while rsq > 1.0:
              dx = 2.0*self.random() - 1.0
              dy = 2.0*self.random() - 1.0
              dz = 2.0*self.random() - 1.0
              rsq = dx*dx + dy*dy + dz*dz
            r = math.sqrt(rsq)
            dx,dy,dz = dx/r,dy/r,dz/r
            x = atoms[-1][3] + dx*self.blen
            y = atoms[-1][4] + dy*self.blen
            z = atoms[-1][5] + dz*self.blen
            restriction = False
            if imonomer >= 2:
              dx = x - atoms[-2][3]
              dy = y - atoms[-2][4]
              dz = z - atoms[-2][5]
              if math.sqrt(dx*dx + dy*dy + dz*dz) <= self.dmin:
                restriction = True
              
        x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
        idatom = id_atom_prev + imonomer + 1
        if self.id == "chain":
          idmol = id_mol_prev + 1
        elif self.id == "end1":
          idmol = imonomer + 1
        elif self.id == "end2":
          idmol = imonomer + 1
          if idmol > nper/2:
            idmol = nper - imonomer
        else:
          raise StandardError,"chain ID is not a valid value"
	        
        atoms.append([idatom,idmol,self.mtype,x,y,z,ix,iy,iz])
        if imonomer:
	  bondid = id_bond_prev + imonomer
          bonds.append([bondid,self.btype,idatom-1,idatom])
        
      self.atoms += atoms
      self.bonds += bonds

  # --------------------------------------------------------------------

  def write(self,file):
    if len(self.atoms) != self.n:
      raise StandardError,"%d monomers instead of requested %d" % \
                           (len(self.atoms),self.n)

    list = [atom[2] for atom in self.atoms]
    atypes = max(list)

    btypes = 0
    if len(self.bonds):
      list = [bond[1] for bond in self.bonds]
      btypes = max(list)

    # create the data file

    d = data()
    d.title = "LAMMPS FENE chain data file"
    d.headers["atoms"] = len(self.atoms)
    d.headers["bonds"] = len(self.bonds)
    d.headers["atom types"] = atypes
    d.headers["bond types"] = btypes
    d.headers["xlo xhi"] = (self.xlo,self.xhi)
    d.headers["ylo yhi"] = (self.ylo,self.yhi)
    d.headers["zlo zhi"] = (self.zlo,self.zhi)

    lines = []
    for i in range(atypes): lines.append("%d 1.0\n" % (i+1))
    d.sections["Masses"] = lines
    
    lines = []
    for atom in self.atoms:
      line = "%d %d %d %g %g %g %d %d %d\n" % \
             (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
              atom[6], atom[7], atom[8])
      lines.append(line)
    d.sections["Atoms"] = lines
    
    lines = []
    for bond in self.bonds:
      line = "%d %d %d %d\n" % (bond[0], bond[1], bond[2], bond[3])
      lines.append(line)
    d.sections["Bonds"] = lines

    d.write(file)

  # --------------------------------------------------------------------

  def pbc(self,x,y,z,ix,iy,iz):
    if x < self.xlo:
      x += self.xprd
      ix -= 1
    elif x >= self.xhi:
      x -= self.xprd
      ix += 1
    if y < self.ylo:
      y += self.yprd
      iy -= 1
    elif y >= self.yhi:
      y -= self.yprd
      iy += 1
    if z < self.zlo:
      z += self.zprd
      iz -= 1
    elif z >= self.zhi:
      z -= self.zprd
      iz += 1
    return x,y,z,ix,iy,iz

  # --------------------------------------------------------------------

  def random(self):
    k = self.seed/IQ
    self.seed = IA*(self.seed-k*IQ) - IR*k
    if self.seed < 0:
      self.seed += IM
    return AM*self.seed

# --------------------------------------------------------------------
# random # generator constants

IM = 2147483647
AM = 1.0/IM
IA = 16807
IQ = 127773
IR = 2836
