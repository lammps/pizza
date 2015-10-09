# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# patch tool

oneline = "Create patchy or rigid particles for LAMMPS input"

docstr = """
p = patch(vfrac)           setup box with a specified volume fraction
p = patch(vfrac,1,1,2)     x,y,z = aspect ratio of box (def = 1,1,1)
    
p.seed = 48379		   set random # seed (def = 12345)
p.randomized = 0	   1 = choose next mol randomly (def), 0 = as generated
p.dim = 2		   set dimension of created box (def = 3)
p.blen = 0.97              set length of tether bonds (def = 0.97)
p.dmin = 1.02              set min r from i-1 to i+1 tether site (def = 1.02)
p.lattice = [Nx,Ny,Nz]     generate Nx by Ny by Nz lattice of particles
p.displace = [Dx,Dy,Dz]    displace particles randomly by +/- Dx,Dy,Dz
p.style = "sphere"         atom-style of data file, molecular or sphere
p.extra = "Molecules"      add extra Molecules section to data file
p.extratype = 1            add extra atom types when write data file

  randomized means choose molecules in random order when creating output
  if lattice is set, Nx*Ny*Nz must equal N for build (Nz = 1 for 2d)
  lattice = [0,0,0] = generate N particles randomly = default
  displace = [0,0,0] = default
  displacement applied when writing molecule to data file
  style = molecular by default
  style is auto-set to line,tri,box by corresponding keywords
  extratype = 0 by default

p.build(100,"hex2",1,2,3)  create 100 "hex2" particles with params 1,2,3
  
  can be invoked multiple times
  keywords:
    c60hex2: diam,1,2,3 = C-60 with 2 hex patches and ctr part, types 1,2,3
    hex2: diam,1,2 = one large particle with 2 7-mer hex patches, types 1,2
    hex4: diam,1,2 = one large particle with 4 7-mer hex patches, types 1,2
    ring: diam,N,1,2 = one large part with equatorial ring of N, types 1,2
    ball: diam,m1,m2,1,2,3 = large ball with m12-len tethers, types 1,2,3
    tri5: 1,2 = 3-layer 5-size hollow tri, types 1,2
    rod: N,m1,m2,1,2,3 = N-length rod with m12-len tethers, types 1,2,3
    tri: N,m1,m2,m3,1,2,3,4 = N-size tri with m123-len tethers, types 1-4
    trid2d: N,r,1 = 3d equilateral tri, N beads r apart, type 1, no bonds
    hex: m1,m2,m3,m4,m5,m6,1,2,3,4,5,6,7 = 7-atom hex with m-len tethers, t 1-7
    dimer: r,1 = two particles r apart, type 1, no bond
    star2d: N,r,1 = 2d star of length N (odd), beads r apart, type 1, no bonds
    box2d: N,M,r,1 = 2d NxM hollow box, beads r apart, type 1, no bonds
    pgon2d: Nlo,Nhi,m = 2d hollow polygons with random N beads from Nlo to Nhi
    sphere3d: Nlo,Nhi,m = 3d hollow spheres with random N beads/cube-edge 
                          from Nlo to Nhi
    tritet: A,m = 4-tri tet with edge length A, tri type m
    tribox: Alo,Ahi,Blo,Bhi,Clo,Chi,m = 12-tri box with side lengths A,B,C & m
    linebox: Alo,Ahi,Blo,Bhi,m = 4-line 2d rectangle with random side lengths
                                 from Alo to Ahi and Blo to Bhi, line type m
                                 built of line particles
    linetri: Alo,Ahi,Blo,Bhi,m = 3-line 2d triangle with random base
                                 from Alo to Ahi and height Blo to Bhi, type m
                                 built of triangle particles
    bodypgon: Nlo,Nhi,m = 2d polygons with random N particles from Nlo to Nhi
                                 built of body particles
    
p.write("data.patch")      write out system to LAMMPS data file
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   style = atom-style of output data file (e.g. molecular)
#   vfrac = desired volume fraction
#   x,y,z = aspect ratio of box (def = 1,1,1)
#   seed = random seed
#   molecules = list of atoms, grouped by molecule

# Imports and external programs

from math import pi,sqrt,cos,sin
from data import data

# Class definition

class patch:
  
  # --------------------------------------------------------------------

  def __init__(self,vfrac,*list):
    self.vfrac = vfrac
    self.xaspect = self.yaspect = self.zaspect = 1.0
    if len(list):
      self.xaspect = list[0]
      self.yaspect = list[1]
      self.zaspect = list[2]
    self.seed = 12345
    self.randomized = 1
    self.dim = 3
    self.molecules = []
    self.volume = 0
    self.blen = 0.97
    self.dmin = 1.02
    self.lattice = [0,0,0]
    self.displace = [0.0,0.0,0.0]
    self.style = "molecular"
    self.extra = ""
    self.extratype = 0
    
  # --------------------------------------------------------------------
  # call style method with extra args
  # adds to volume and atom list
  # reset self.style for lines and triangles

  def build(self,n,style,*types):
    cmd = "atoms,bonds,tris,segments,bodies,volume = self.%s(*types)" % style
    for i in xrange(n):
      exec cmd
      self.molecules.append([atoms,bonds,tris,segments,bodies])
      self.volume += volume

  # --------------------------------------------------------------------
  # create the atom coords in a scaled-size box of correct dimension
  # write them to LAMMPS data file

  def write(self,file):
    if self.dim == 3: self.write3d(file)
    else: self.write2d(file)

  # --------------------------------------------------------------------
  # write a 3d simulation to data file

  def write3d(self,file):
    volume = self.volume/self.vfrac
    prd = pow(volume/self.xaspect/self.yaspect/self.zaspect,1.0/3.0)
    self.xprd = self.xaspect * prd
    self.xlo = -self.xprd/2.0
    self.xhi = self.xprd/2.0
    self.yprd = self.yaspect * prd
    self.ylo = -self.yprd/2.0
    self.yhi = self.yprd/2.0
    self.zprd = self.zaspect * prd
    self.zlo = -self.zprd/2.0
    self.zhi = self.zprd/2.0

    if self.lattice[0] or self.lattice[1] or self.lattice[2]:
      latflag = 1
      if self.lattice[0]*self.lattice[1]*self.lattice[2] != \
             len(self.molecules):
        raise StandardError,"lattice inconsistent with # of molecules"
    else: latflag = 0

    idatom = idbond = idtri = idmol = 0
    atoms = []
    bonds = []
    tris = []
    mols = []
    xp = 3*[0]
    yp = 3*[0]
    zp = 3*[0]
    maxtypes = 0

    while self.molecules:
      if self.randomized: i = int(self.random()*len(self.molecules))
      else: i = 0
      molecule = self.molecules.pop(i)

      idmol += 1
      triples = []

      # xp[3],yp[3],zp[3] = randomly oriented, normalized basis vectors
      # xp is in random direction
      # yp is random dir crossed into xp
      # zp is xp crossed into yp

      xp[0] = self.random() - 0.5
      xp[1] = self.random() - 0.5
      xp[2] = self.random() - 0.5
      r = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2])
      xp[0],xp[1],xp[2] = xp[0]/r,xp[1]/r,xp[2]/r

      r0 = self.random() - 0.5
      r1 = self.random() - 0.5
      r2 = self.random() - 0.5
      yp[0] = r1*xp[2] - r2*xp[1]
      yp[1] = r2*xp[0] - r0*xp[2]
      yp[2] = r0*xp[1] - r1*xp[0]
      r = sqrt(yp[0]*yp[0] + yp[1]*yp[1] + yp[2]*yp[2])
      yp[0],yp[1],yp[2] = yp[0]/r,yp[1]/r,yp[2]/r

      zp[0] = xp[1]*yp[2] - xp[2]*yp[1]
      zp[1] = xp[2]*yp[0] - xp[0]*yp[2]
      zp[2] = xp[0]*yp[1] - xp[1]*yp[0]
      r = sqrt(zp[0]*zp[0] + zp[1]*zp[1] + zp[2]*zp[2])
      zp[0],zp[1],zp[2] = zp[0]/r,zp[1]/r,zp[2]/r

      #xp[0] = 1; xp[1] = 0; xp[2] = 0
      #yp[0] = 0; yp[1] = 1; yp[2] = 0
      #zp[0] = 0; zp[1] = 0; zp[2] = 1
      
      # random origin or lattice site for new particle

      if latflag == 0:
        xorig = self.xlo + self.random()*self.xprd
        yorig = self.ylo + self.random()*self.yprd
        zorig = self.zlo + self.random()*self.zprd
      else:
        ix = (idmol-1) % self.lattice[0]
        iy = (idmol-1)/self.lattice[0] % self.lattice[1]
        iz = (idmol-1) / (self.lattice[0]*self.lattice[1])
        xorig = self.xlo + ix*self.xprd/self.lattice[0]
        yorig = self.ylo + iy*self.yprd/self.lattice[1]
        zorig = self.zlo + iz*self.zprd/self.lattice[2]

      #xorig = 0; yorig = 0; zorig = 0
        
      # unpack bonds in molecule before atoms so idatom = all previous atoms
      
      for bond in molecule[1]:
        idbond += 1
        bonds.append([idbond,bond[0],bond[1]+idatom+1,bond[2]+idatom+1])

      # unpack triples in molecule as displacements from associated atom
        
      for triple in molecule[2]: triples.append(triple)

      # unpack atoms in molecule
      # xnew,ynew,xnew = coeffs in new rotated basis vectors
      # x,y,z = coeffs in original xyz axes
      # apply random displace to final x,y,z for molecular style
      # format data file for moleular or tri atom style

      if self.style == "molecular":
        for atom in molecule[0]:
          idatom += 1
          xnew = atom[1]
          ynew = atom[2]
          znew = atom[3]
          x = xorig + xnew*xp[0] + ynew*yp[0] + znew*zp[0]
          y = yorig + xnew*xp[1] + ynew*yp[1] + znew*zp[1]
          z = zorig + xnew*xp[2] + ynew*yp[2] + znew*zp[2]
          x += (self.random()-0.5)*2*self.displace[0]
          y += (self.random()-0.5)*2*self.displace[1]
          z += (self.random()-0.5)*2*self.displace[2]
          ix = iy = iz = 0
          x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
          atoms.append([idatom,idmol,atom[0],x,y,z,ix,iy,iz])
          mols.append((idatom,idmol))
          maxtypes = max(maxtypes,atom[0])
          
      elif self.style == "sphere":
        for atom in molecule[0]:
          idatom += 1
          xnew = atom[1]
          ynew = atom[2]
          znew = atom[3]
          x = xorig + xnew*xp[0] + ynew*yp[0] + znew*zp[0]
          y = yorig + xnew*xp[1] + ynew*yp[1] + znew*zp[1]
          z = zorig + xnew*xp[2] + ynew*yp[2] + znew*zp[2]
          x += (self.random()-0.5)*2*self.displace[0]
          y += (self.random()-0.5)*2*self.displace[1]
          z += (self.random()-0.5)*2*self.displace[2]
          ix = iy = iz = 0
          x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
          atoms.append([idatom,atom[0],1.0,1.0,x,y,z,ix,iy,iz])
          mols.append((idatom,idmol))
          maxtypes = max(maxtypes,atom[0])
          
      elif self.style == "tri":
        for i,atom in enumerate(molecule[0]):
          idatom += 1
          xnew = atom[1]
          ynew = atom[2]
          znew = atom[3]
          x = xorig + xnew*xp[0] + ynew*yp[0] + znew*zp[0]
          y = yorig + xnew*xp[1] + ynew*yp[1] + znew*zp[1]
          z = zorig + xnew*xp[2] + ynew*yp[2] + znew*zp[2]
          mass = 1.0
          ix = iy = iz = 0
          x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
          if not triples: triflag = 0
          else: triflag = 1
          atoms.append([idatom,idmol,atom[0],triflag,mass,x,y,z,ix,iy,iz])
          mols.append((idatom,idmol))
          maxtypes = max(maxtypes,atom[0])
          
          if triflag:
            triple = triples[i]
            xtri = triple[0]
            ytri = triple[1]
            ztri = triple[2]
            triple[0] = xorig + xtri*xp[0] + ytri*yp[0] + ztri*zp[0]
            triple[1] = yorig + xtri*xp[1] + ytri*yp[1] + ztri*zp[1]
            triple[2] = zorig + xtri*xp[2] + ytri*yp[2] + ztri*zp[2]
            triple[0],triple[1],triple[2] = \
                self.pbc_near(triple[0],triple[1],triple[2],x,y,z)
            xtri = triple[3]
            ytri = triple[4]
            ztri = triple[5]
            triple[3] = xorig + xtri*xp[0] + ytri*yp[0] + ztri*zp[0]
            triple[4] = yorig + xtri*xp[1] + ytri*yp[1] + ztri*zp[1]
            triple[5] = zorig + xtri*xp[2] + ytri*yp[2] + ztri*zp[2]
            triple[3],triple[4],triple[5] = \
                self.pbc_near(triple[3],triple[4],triple[5],x,y,z)
            xtri = triple[6]
            ytri = triple[7]
            ztri = triple[8]
            triple[6] = xorig + xtri*xp[0] + ytri*yp[0] + ztri*zp[0]
            triple[7] = yorig + xtri*xp[1] + ytri*yp[1] + ztri*zp[1]
            triple[8] = zorig + xtri*xp[2] + ytri*yp[2] + ztri*zp[2]
            triple[6],triple[7],triple[8] = \
                self.pbc_near(triple[6],triple[7],triple[8],x,y,z)
            tris.append([idatom] + triple)

    # create the data file

    d = data()
    d.title = "LAMMPS data file for Nanoparticles"
    d.headers["atoms"] = len(atoms)
    d.headers["atom types"] = maxtypes + self.extratype
    if bonds:
      d.headers["bonds"] = len(bonds)
      d.headers["bond types"] = 1
    if tris: d.headers["triangles"] = len(tris)

    if tris: print "TRIS",len(tris),d.headers["triangles"]
    d.headers["xlo xhi"] = (self.xlo,self.xhi)
    d.headers["ylo yhi"] = (self.ylo,self.yhi)
    d.headers["zlo zhi"] = (self.zlo,self.zhi)

    # atoms section of data file

    records = []
    if self.style == "molecular":
      for atom in atoms:
        one = "%d %d %d %g %g %g %d %d %d\n" % \
            (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
             atom[6], atom[7], atom[8])
        records.append(one)
    elif self.style == "sphere":
      for atom in atoms:
        one = "%d %d %g %g %g %g %g %d %d %d\n" % \
            (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
             atom[6], atom[7], atom[8], atom[9])
        records.append(one)
    elif self.style == "tri":
      for atom in atoms:
        one = "%d %d %d %d %g %g %g %g %d %d %d\n" % \
            (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
             atom[6], atom[7], atom[8], atom[9], atom[10])
        records.append(one)
    d.sections["Atoms"] = records

    # bonds section of data file

    if bonds:
      lines = []
      for bond in bonds:
        line = "%d %d %d %d\n" % (bond[0], bond[1], bond[2], bond[3])
        lines.append(line)
      d.sections["Bonds"] = lines

     # triangles section of data file

    if tris:
      records = []
      for tri in tris:
        one = "%d %g %g %g %g %g %g %g %g %g\n" % \
            (tri[0],tri[1],tri[2],tri[3],tri[4],tri[5],
             tri[6],tri[7],tri[8],tri[9])
        records.append(one)
      d.sections["Triangles"] = records

    # extra Molecules section if requested
      
    if self.extra == "Molecules":
      records = []
      for mol in mols:
        one = "%d %d\n" % (mol[0],mol[1])
        records.append(one)
      d.sections["Molecules"] = records

    d.write(file)

  # --------------------------------------------------------------------
  # write a 2d simulation to data file

  def write2d(self,file):
    volume = self.volume/self.vfrac
    prd = pow(volume/self.xaspect/self.yaspect,1.0/2.0)
    self.xprd = self.xaspect * prd
    self.xlo = -self.xprd/2.0
    self.xhi = self.xprd/2.0
    self.yprd = self.yaspect * prd
    self.ylo = -self.yprd/2.0
    self.yhi = self.yprd/2.0
    self.zprd = 1.0
    self.zlo = -0.5
    self.zhi = 0.5

    if self.lattice[0] or self.lattice[1]:
      latflag = 1
      if self.lattice[0]*self.lattice[1] != len(self.molecules):
        raise StandardError,"lattice inconsistent with # of molecules"
    else: latflag = 0
    
    idatom = idbond = idmol = 0
    atoms = []
    bonds = []
    lines = []
    bodies = []
    mols = []
    xp = 3*[0]
    yp = 3*[0]
    maxtypes = 0

    while self.molecules:
      if self.randomized: i = int(self.random()*len(self.molecules))
      else: i = 0
      molecule = self.molecules.pop(i)
        
      idmol += 1
      segments = []
      subs = []
      
      # xp[2],yp[2] = randomly oriented, normalized basis vectors
      # xp is in random direction
      # yp is (0,0,1) crossed into xp

      xp[0] = self.random() - 0.5
      xp[1] = self.random() - 0.5
      r = sqrt(xp[0]*xp[0] + xp[1]*xp[1])
      xp[0],xp[1] = xp[0]/r,xp[1]/r

      yp[0] = -xp[1]
      yp[1] = xp[0]
      r = sqrt(yp[0]*yp[0] + yp[1]*yp[1])
      yp[0],yp[1] = yp[0]/r,yp[1]/r

      # random origin or lattice site for new particle

      if latflag == 0:
        xorig = self.xlo + self.random()*self.xprd
        yorig = self.ylo + self.random()*self.yprd
        zorig = 0.0
      else:
        ix = (idmol-1) % self.lattice[0]
        iy = (idmol-1) / self.lattice[0]
        xorig = self.xlo + ix*self.xprd/self.lattice[0]
        yorig = self.ylo + iy*self.yprd/self.lattice[1]
        zorig = 0.0
        
      # unpack bonds in molecule before atoms so idatom = all previous atoms
      # segments = molecule[3] field = displacement from associated atom
      # subs = molecule[4] field = sub-particle displacements from body atom
      
      for bond in molecule[1]:
        idbond += 1
        bonds.append([idbond,bond[0],bond[1]+idatom+1,bond[2]+idatom+1])
      segments = molecule[3]
      subs = molecule[4]

      # unpack atoms in molecule
      # xnew,ynew,xnew = coeffs in new rotated basis vectors
      # x,y,z = coeffs in original xyz axes
      # apply random displace to final x,y for molecular style
      # format data file for moleular or line atom style

      if self.style == "molecular":
        for i,atom in enumerate(molecule[0]):
          idatom += 1
          xnew = atom[1]
          ynew = atom[2]
          x = xorig + xnew*xp[0] + ynew*yp[0]
          y = yorig + xnew*xp[1] + ynew*yp[1]
          z = atom[3]
          x += (self.random()-0.5)*2*self.displace[0]
          y += (self.random()-0.5)*2*self.displace[1]
          ix = iy = iz = 0
          x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
          atoms.append([idatom,idmol,atom[0],x,y,z,ix,iy,iz])
          mols.append((idatom,idmol))
          maxtypes = max(maxtypes,atom[0])

      elif self.style == "sphere":
        for i,atom in enumerate(molecule[0]):
          idatom += 1
          xnew = atom[1]
          ynew = atom[2]
          x = xorig + xnew*xp[0] + ynew*yp[0]
          y = yorig + xnew*xp[1] + ynew*yp[1]
          z = atom[3]
          x += (self.random()-0.5)*2*self.displace[0]
          y += (self.random()-0.5)*2*self.displace[1]
          ix = iy = iz = 0
          x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
          atoms.append([idatom,atom[0],1.0,1.0,x,y,z,ix,iy,iz])
          mols.append((idatom,idmol))
          maxtypes = max(maxtypes,atom[0])

      elif self.style == "line":
        for i,atom in enumerate(molecule[0]):
          idatom += 1
          xnew = atom[1]
          ynew = atom[2]
          x = xorig + xnew*xp[0] + ynew*yp[0]
          y = yorig + xnew*xp[1] + ynew*yp[1]
          z = atom[3]
          mass = 1.0
          ix = iy = iz = 0
          x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
          if not segments: lineflag = 0
          else: lineflag = 1
          atoms.append([idatom,idmol,atom[0],lineflag,mass,x,y,z,ix,iy,iz])
          mols.append((idatom,idmol))
          maxtypes = max(maxtypes,atom[0])

          if lineflag:
            segment = segments[i]
            xseg = xnew + segment[0]
            yseg = ynew + segment[1]
            segment[0] = xorig + xseg*xp[0] + yseg*yp[0]
            segment[1] = yorig + xseg*xp[1] + yseg*yp[1]
            segment[0],segment[1],tmp = \
                self.pbc_near(segment[0],segment[1],0,x,y,z)
            xseg = xnew + segment[2]
            yseg = ynew + segment[3]
            segment[2] = xorig + xseg*xp[0] + yseg*yp[0]
            segment[3] = yorig + xseg*xp[1] + yseg*yp[1]
            segment[2],segment[3],tmp = \
                self.pbc_near(segment[2],segment[3],0,x,y,z)
            lines.append([idatom] + segment)
            
      elif self.style == "body":
        atom = molecule[0][0]
        idatom += 1
        xnew = atom[1]
        ynew = atom[2]
        x = xorig + xnew*xp[0] + ynew*yp[0]
        y = yorig + xnew*xp[1] + ynew*yp[1]
        z = atom[3]
        mass = atom[4]
        ix = iy = iz = 0
        x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
        if not subs: bodyflag = 0
        else: bodyflag = 1
        atoms.append([idatom,atom[0],bodyflag,mass,x,y,z,ix,iy,iz])
        mols.append((idatom,idmol))
        maxtypes = max(maxtypes,atom[0])
        
        if bodyflag:
          ivalues = [mass]
          for one in subs:
            x = one[0]
            y = one[1]
            one[0] = x*xp[0] + y*yp[0]
            one[1] = x*xp[1] + y*yp[1]
          inertia = [0,0,0,0,0,0]
          for one in subs:
            inertia[0] += one[1]*one[1] + one[2]*one[2]
            inertia[1] += one[0]*one[0] + one[2]*one[2]
            inertia[2] += one[0]*one[0] + one[1]*one[1]
            inertia[3] -= one[0]*one[1]
            inertia[4] -= one[0]*one[2]
            inertia[5] -= one[1]*one[2]
          dvalues = inertia
          for one in subs:
            dvalues += one
          bodies.append([idatom,ivalues,dvalues])
            
    # create the data file

    d = data()
    d.title = "LAMMPS data file for Nanoparticles"
    d.headers["atoms"] = len(atoms)
    d.headers["atom types"] = maxtypes + self.extratype
    if bonds:
      d.headers["bonds"] = len(bonds)
      d.headers["bond types"] = 1
    if lines: d.headers["lines"] = len(lines)
    if bodies: d.headers["bodies"] = len(bodies)
    d.headers["xlo xhi"] = (self.xlo,self.xhi)
    d.headers["ylo yhi"] = (self.ylo,self.yhi)
    d.headers["zlo zhi"] = (self.zlo,self.zhi)

    # atoms section of data file
    
    records = []
    if self.style == "molecular":
      for atom in atoms:
        one = "%d %d %d %g %g %g %d %d %d\n" % \
            (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
             atom[6], atom[7], atom[8])
        records.append(one)
    elif self.style == "sphere":
      for atom in atoms:
        one = "%d %d %g %g %g %g %g %d %d %d\n" % \
            (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
             atom[6], atom[7], atom[8], atom[9])
        records.append(one)
    elif self.style == "line":
      for atom in atoms:
        one = "%d %d %d %d %g %g %g %g %d %d %d\n" % \
            (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
             atom[6], atom[7], atom[8], atom[9], atom[10])
        records.append(one)
    elif self.style == "body":
      for atom in atoms:
        one = "%d %d %d %d %g %g %g %g %d %d\n" % \
            (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
             atom[6], atom[7], atom[8], atom[9])
        records.append(one)
    d.sections["Atoms"] = records

    # bonds section of data file

    if bonds:
      records = []
      for bond in bonds:
        one = "%d %d %d %d\n" % (bond[0], bond[1], bond[2], bond[3])
        records.append(one)
      d.sections["Bonds"] = records

    # lines section of data file

    if lines:
      records = []
      for line in lines:
        one = "%d %g %g %g %g\n" % (line[0],line[1],line[2],line[3],line[4])
        records.append(one)
      d.sections["Lines"] = records

    # bodies section of data file
    # print ivalues and dvalues in sets of 10 per line

    if bodies:
      records = []
      for body in bodies:
        ivalues = body[1]
        dvalues = body[2]
        one = "%d %d %d\n" % (body[0],len(ivalues),len(dvalues))
        records.append(one)
        i = 0
        while i < len(ivalues):
          list = ivalues[i:min(i+10,len(ivalues)+1)]
          one = ""
          for value in list: one += "%g " % int(value)
          one += "\n"
          records.append(one)
          i += len(list)
        i = 0
        while i < len(dvalues):
          list = dvalues[i:min(i+10,len(dvalues)+1)]
          one = ""
          for value in list: one += "%g " % float(value)
          one += "\n"
          records.append(one)
          i += len(list)
      d.sections["Bodies"] = records

    # extra Molecules section if requested

    if self.extra == "Molecules":
      records = []
      for mol in mols:
        one = "%d %d\n" % (mol[0],mol[1])
        records.append(one)
      d.sections["Molecules"] = records
      
    d.write(file)

  # --------------------------------------------------------------------
  # adjust x,y,z to be inside periodic box
    
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
  # adjust xnew,ynew,znew to be near x,y,z in periodic sense
    
  def pbc_near(self,xnew,ynew,znew,x,y,z):
    if x-xnew > 0.5*self.xprd: xnew += self.xprd
    elif xnew-x > 0.5*self.xprd: xnew -= self.xprd
    if y-ynew > 0.5*self.yprd: ynew += self.yprd
    elif ynew-y > 0.5*self.yprd: ynew -= self.yprd
    if z-znew > 0.5*self.zprd: znew += self.zprd
    elif znew-z > 0.5*self.zprd: znew -= self.zprd
    return xnew,ynew,znew

  # --------------------------------------------------------------------
  # params = diam,type1,type2,type3
  # type1 = type of non-patch atoms, type2 = type of patch atoms
  # type3 = type of center-of-sphere atom
  
  def c60hex2(self,*params):
    template = BUCKY_60
    diam = params[0]
    patches = [[params[2],[1,2,3,4,5,6,34,35,46,47,48,60]],
               [params[3],[61]]]
    atoms = make_sphere(template,diam,params[1],patches)
    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],[],[],[],volume
  
  # --------------------------------------------------------------------
  # params = diam,type1,type2
  # type1 = type of large center atom, type2 = type of hex patch atoms
  
  def hex2(self,*params):
    diam = params[0]
    type1 = params[1]
    type2 = params[2]
    
    atoms = []
    atoms.append([type1,0.0,0.0,0.0])
    
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,-sqrt(3.0)/2))

    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],[],[],[],volume
  
  # --------------------------------------------------------------------
  # params = diam,type1,type2
  # type1 = type of large center atom, type2 = type of hex patch atoms
  
  def hex4(self,*params):
    diam = params[0]
    type1 = params[1]
    type2 = params[2]
    
    atoms = []
    atoms.append([type1,0.0,0.0,0.0])
    
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,0.0,0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,1.0,0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,-1.0,0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5,0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5,0.5*diam,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,0.5*diam,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,0.0,-0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,1.0,-0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,-1.0,-0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5,-0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,-0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5,-0.5*diam,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,-0.5*diam,-sqrt(3.0)/2))

    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = diam,nring,type1,type2
  # nring = # of particles in ring
  # type1 = type of large center atom, type2 = type of ring atoms
  
  def ring(self,*params):
    diam = params[0]
    nring = params[1]
    type1 = params[2]
    type2 = params[3]
    
    atoms = []
    atoms.append([type1,0.0,0.0,0.0])

    for i in range(nring):
      atoms.append([type2,0.5*diam*cos(i * 2*pi/nring),
                    0.5*diam*sin(i * 2*pi/nring),0.0])

    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = diam,m1,m2,ntype,m1type,m2type
  # ntype = type of big central ball
  # m12,m12type = length of tethers on each side of ball (m12 = 0 = no tether)
  # set three types of bonds:
  #   1 = big to small, 2 = small to small, 3 = across two tethers
  
  def ball(self,*params):
    diam = params[0]
    m1 = params[1]
    m2 = params[2]
    ntype = params[3]
    m1type = params[4]
    m2type = params[5]

    atoms = []
    atoms.append([ntype,0.0,0.0,0.0])

    if m1:
      atoms.append([m1type,0.5*diam+0.5,0.0,0.0])
      atoms += tether(m1-1,m1type,self.blen,self.dmin,
                      [ntype,0.5*diam-0.5,0.0,0.0],atoms[1],self.random)
    if m2:
      atoms.append([m2type,-0.5*diam-0.5,0.0,0.0])
      atoms += tether(m2-1,m2type,self.blen,self.dmin,
                      [ntype,-0.5*diam+0.5,0.0,0.0],atoms[1+m1],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,0,1])
      else: bonds.append([2,1+i-1,1+i])
    for i in range(m2):
      if i == 0: bonds.append([1,0,1+m1])
      else: bonds.append([2,1+m1+i-1,1+m1+i])
    if m1 and m2: bonds.append([3,1,m1+1])
      
    volume = 4.0/3.0 * pi * diam*diam*diam/8 + (m1+m2)*pi/6.0
    return atoms,bonds,[],[],[],volume

  # --------------------------------------------------------------------
  # params = type1,type2
  # type1 = type of 1/3 layers, type2 = type of middle layer
  
  def tri5(self,*params):
    template = TRI5_HOLLOW
    nlayer = 3
    patches = [[params[1],[13,14,15,16,17,18,19,20,21,22,23,24]]]

    atoms = []
    n = 0
    for iz in range(nlayer):
      for atom in template:
        n += 1
        type = params[0]
        for pair in patches:
          patch = pair[1]
          if n in patch: type = pair[0]
        x,y,z = atom[0],atom[1],atom[2] + iz
        atoms.append([type,x,y,z])

    volume = nlayer * 0.5 * 5.0 * 5.0*sqrt(3)/2
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = n,m1,m2,ntype,m1type,m2type
  # n,ntype = length/type of rod
  # m12,m12type = length of tethers on each end (m12 = 0 = no tether)

  def rod(self,*params):
    n = params[0]
    m1 = params[1]
    m2 = params[2]
    ntype = params[3]
    m1type = params[4]
    m2type = params[5]
    
    atoms = []
    for i in range(n):
      x,y,z = i*self.blen,0.0,0.0
      atoms.append([ntype,x,y,z])
      
    if m1: atoms += tether(m1,m1type,self.blen,self.dmin,
                           atoms[n-2],atoms[n-1],self.random)
    if m2: atoms += tether(m2,m2type,self.blen,self.dmin,
                           atoms[1],atoms[0],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,n-1,n])
      else: bonds.append([1,n+i-1,n+i])
    for i in range(m2):
      if i == 0: bonds.append([1,0,n+m1])
      else: bonds.append([1,n+m1+i-1,n+m1+i])
  
    volume = (n+m1+m2) * pi / 6.0
    return atoms,bonds,[],[],[],volume

  # --------------------------------------------------------------------
  # params = nsize,m1,m2,m3,ntype,m1type,m2type,m3type
  # nsize,ntype = size,type of triangle center
  # m123,m123type = length of tethers on each corner (m123 = 0 = no tether)
  
  def tri(self,*params):
    nsize = params[0]
    m1 = params[1]
    m2 = params[2]
    m3 = params[3]
    ntype = params[4]
    m1type = params[5]
    m2type = params[6]
    m3type = params[7]
  
    atoms = []
    for i in range(nsize):
      n = nsize - i
      for j in range(n):
        x,y,z = 0.5*i*self.blen + j*self.blen, i*self.blen*sqrt(3.0)/2.0, 0.0
        atoms.append([ntype,x,y,z])

    n = len(atoms)
    if m1: atoms += tether(m1,m1type,self.blen,self.dmin,
                           atoms[1],atoms[0],self.random)
    if m2: atoms += tether(m2,m2type,self.blen,self.dmin,
                           atoms[nsize-2],atoms[nsize-1],self.random)
    if m3: atoms += tether(m3,m3type,self.blen,self.dmin,
                           atoms[n-2],atoms[n-1],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,0,n])
      else: bonds.append([1,n+i-1,n+i])
    for i in range(m2):
      if i == 0: bonds.append([1,nsize-1,n+m1])
      else: bonds.append([1,n+m1+i-1,n+m1+i])
    for i in range(m3):
      if i == 0: bonds.append([1,n-1,n+m1+m2])
      else: bonds.append([1,n+m1+m2+i-1,n+m1+m2+i])

    volume = (nsize*(nsize+1)/2 + m1+m2+m3) * pi / 6.0
    return atoms,bonds,[],[],[],volume

  # --------------------------------------------------------------------
  # params = N,sep,type
  # N = length of each side of equilateral tri
  # sep = separation distance of consecutive particles
  # type = type of all particles

  def tri2d(self,*params):
    n = params[0]
    sep = params[1]
    type = params[2]

    length = (n-1) * sep
    
    atoms = []
    for i in range(n):
      x,y,z = i*sep,0.0,0.0
      atoms.append([type,x,y,z])
    for i in range(1,n):
      x,y,z = i*sep*cos(pi/3.0),i*sep*sin(pi/3.0),0.0
      atoms.append([type,x,y,z])
    for i in range(1,n-1):
      x,y,z = length-i*sep*cos(pi/3.0),i*sep*sin(pi/3.0),0.0
      atoms.append([type,x,y,z])

    volume = (3*n-2) * pi / 6.0     # need to account for overlap and 2d/3d
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = m1,m2,m3,m4,m5,m6,ntype,m1type,m2type,m3type,m4type,m5type,m6type
  # ntype = type of hex center
  # m123456,m123456type = length of tethers on each corner (m = 0 = no tether)
  
  def hex(self,*params):
    m1 = params[0]
    m2 = params[1]
    m3 = params[2]
    m4 = params[3]
    m5 = params[4]
    m6 = params[5]
    ntype = params[6]
    m1type = params[7]
    m2type = params[8]
    m3type = params[9]
    m4type = params[10]
    m5type = params[11]
    m6type = params[12]
  
    atoms = []
    atoms.append([ntype,0.0,0.0,0.0])
    atoms.append([ntype,self.blen,0.0,0.0])
    atoms.append([ntype,-self.blen,0.0,0.0])
    atoms.append([ntype,self.blen/2.0,self.blen*sqrt(3.0)/2.0,0.0])
    atoms.append([ntype,-self.blen/2.0,self.blen*sqrt(3.0)/2.0,0.0])
    atoms.append([ntype,self.blen/2.0,-self.blen*sqrt(3.0)/2.0,0.0])
    atoms.append([ntype,-self.blen/2.0,-self.blen*sqrt(3.0)/2.0,0.0])
    
    n = len(atoms)
    if m1: atoms += tether(m1,m1type,self.blen,self.dmin,
                           atoms[0],atoms[1],self.random)
    if m2: atoms += tether(m2,m2type,self.blen,self.dmin,
                           atoms[0],atoms[2],self.random)
    if m3: atoms += tether(m3,m3type,self.blen,self.dmin,
                           atoms[0],atoms[3],self.random)
    if m4: atoms += tether(m4,m4type,self.blen,self.dmin,
                           atoms[0],atoms[4],self.random)
    if m5: atoms += tether(m5,m5type,self.blen,self.dmin,
                           atoms[0],atoms[5],self.random)
    if m6: atoms += tether(m6,m6type,self.blen,self.dmin,
                           atoms[0],atoms[6],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,1,n])
      else: bonds.append([1,n+i-1,n+i])
    for i in range(m2):
      if i == 0: bonds.append([1,2,n+m1])
      else: bonds.append([1,n+m1+i-1,n+m1+i])
    for i in range(m3):
      if i == 0: bonds.append([1,3,n+m1+m2])
      else: bonds.append([1,n+m1+m2+i-1,n+m1+m2+i])
    for i in range(m4):
      if i == 0: bonds.append([1,4,n+m1+m2+m3])
      else: bonds.append([1,n+m1+m2+m3+i-1,n+m1+m2+m3+i])
    for i in range(m5):
      if i == 0: bonds.append([1,5,n+m1+m2+m3+m4])
      else: bonds.append([1,n+m1+m2+m3+m4+i-1,n+m1+m2+m3+m4+i])
    for i in range(m6):
      if i == 0: bonds.append([1,6,n+m1+m2+m3+m4+m5])
      else: bonds.append([1,n+m1+m2+m3+m4+m5+i-1,n+m1+m2+m3+m4+m5+i])

    volume = (7 + m1+m2+m3+m4+m5+m6) * pi / 6.0
    return atoms,bonds,[],[],[],volume

  # --------------------------------------------------------------------
  # params = sep,type
  # sep = separation distance of two particles
  # type = type of 2 particles

  def dimer(self,*params):
    sep = params[0]
    type = params[1]
    
    atoms = []
    x,y,z = 0.0,0.0,0.0
    atoms.append([type,x,y,z])
    x,y,z = sep,0.0,0.0
    atoms.append([type,x,y,z])

    volume = 2 * pi / 6.0     # need to account for overlap and 2d/3d
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = N,sep,type
  # N = lengths of each arm of star (must be odd)
  # sep = separation distance of consecutive particles
  # type = type of all particles

  def star2d(self,*params):
    n = params[0]
    sep = params[1]
    type = params[2]
    if n % 2 == 0:
      raise StandardError, "N in patch::star2d is not odd"
    middle = n/2
    
    atoms = []
    x,y,z = 0.0,0.0,0.0
    atoms.append([type,x,y,z])
    for i in range(n):
      i -= middle
      if i == 0: continue
      x,y,z = i*sep,0.0,0.0
      atoms.append([type,x,y,z])
    for i in range(n):
      i -= middle
      if i == 0: continue
      x,y,z = 0.0,i*sep,0.0
      atoms.append([type,x,y,z])

    volume = (2*n-1) * pi / 6.0     # need to account for overlap and 2d/3d
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = N,M,sep,type
  # N,M = lengths of each side of box
  # sep = separation distance of consecutive particles
  # type = type of all particles

  def box2d(self,*params):
    n = params[0]
    m = params[1]
    sep = params[2]
    type = params[3]

    height = (m-1) * sep
    width = (n-1) * sep
    
    atoms = []
    for i in range(n):
      x,y,z = i*sep,0.0,0.0
      atoms.append([type,x,y,z])
    for i in range(n):
      x,y,z = i*sep,height,0.0
      atoms.append([type,x,y,z])
    for i in range(1,m-1,1):
      x,y,z = 0.0,i*sep,0.0
      atoms.append([type,x,y,z])
    for i in range(1,m-1,1):
      x,y,z = width,i*sep,0.0
      atoms.append([type,x,y,z])

    volume = (2*(n+m-2)) * pi / 6.0     # need to account for overlap and 2d/3d
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = Nlo,Nhi,type
  # polygon with random N particles from Nlo to Nhi
  # type = type of all particles

  def pgon2d(self,*params):
    nlo = int(params[0])
    nhi = int(params[1])
    type = params[2]

    n = int(nlo + self.random()*(nhi+1-nlo))

    # N sub-particles around perimeter of polygon in xy plane
    # rad set so that successive sub-particles are distance 1.0 apart
    # adjust sub-particle coords so COM is at (0,0,0)

    list = []
    rad = sqrt(0.5/(1-cos(2.0*pi/n)))
    for i in xrange(n):
      theta = 2.0*pi * i/n
      list.append([rad*cos(theta),rad*sin(theta),0])

    com = [0,0]
    for one in list:
      com[0] += one[0]
      com[1] += one[1]
    com[0] /= n
    com[1] /= n
    for one in list:
      one[0] -= com[0]
      one[1] -= com[1]
    
    atoms = []
    for one in list:
      atoms.append([type,one[0],one[1],0.0])

    volume = n*pi / 6.0
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = Nlo,Nhi,type
  # sphere made from cube faces projected to sphere
  # random N particles on each cube edge from Nlo to Nhi
  # type = type of all particles

  def sphere3d(self,*params):
    nlo = int(params[0])
    nhi = int(params[1])
    type = params[2]

    n = int(nlo + self.random()*(nhi+1-nlo))

    # N sub-particles around edges of cube faces
    # rad set so that cube corners are on surface of sphere

    list = []
    rad = sqrt(3.0)/2.0 * (n-1)

    atoms = []
    for i in range(n):
      for j in range(n):
        c = [-0.5 + 1.0*i/n, -0.5 + 1.0*j/n, -0.5]
        scale = rad/sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
        atoms.append([type,scale*c[0],scale*c[1],scale*c[2]])
        atoms.append([type,scale*c[0],scale*c[1],-scale*c[2]])
    for i in range(n):
      for k in range(1,n-1):
        c = [-0.5 + 1.0*i/n, -0.5, -0.5 + 1.0*k/n]
        scale = rad/sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
        atoms.append([type,scale*c[0],scale*c[1],scale*c[2]])
        atoms.append([type,scale*c[0],-scale*c[1],scale*c[2]])
    for j in range(1,n-1):
      for k in range(1,n-1):
        c = [-0.5, -0.5 + 1.0*i/n, -0.5 + 1.0*k/n]
        scale = rad/sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
        atoms.append([type,scale*c[0],scale*c[1],scale*c[2]])
        atoms.append([type,-scale*c[0],scale*c[1],scale*c[2]])

    ntotal = 2 * (n*n + n*(n-2) + (n-2)*(n-2))
    volume = ntotal*pi / 6.0
    return atoms,[],[],[],[],volume

  # --------------------------------------------------------------------
  # params = a,type
  # a = edge length of tet
  # type = type of each vertex in tet
  
  def tritet(self,*params):
    self.style = "tri"
    a = params[0]
    type = params[1]
    cubeside = a/sqrt(2)
    v = 0.5*cubeside

    atoms = []
    tris = []

    # 4 corner pts of tet

    c1 = [+v/2,+v/2,+v/2]
    c2 = [-v/2,-v/2,+v/2]
    c3 = [-v/2,+v/2,-v/2]
    c4 = [+v/2,-v/2,-v/2]

    # 4 triangles on faces of tet, outward normals via rh rule

    tris.append(c1+c3+c2)
    tris.append(c1+c2+c4)
    tris.append(c1+c3+c4)
    tris.append(c2+c3+c4)

    # 4 atoms at centroids of 4 triangles

    atoms.append([type,(c1[0]+c3[0]+c2[0])/3.0,(c1[1]+c3[1]+c2[1])/3.0,
                  (c1[2]+c3[2]+c2[2])/3.0])
    atoms.append([type,(c1[0]+c2[0]+c4[0])/3.0,(c1[1]+c2[1]+c4[1])/3.0,
                  (c1[2]+c2[2]+c4[2])/3.0])
    atoms.append([type,(c1[0]+c3[0]+c4[0])/3.0,(c1[1]+c3[1]+c4[1])/3.0,
                  (c1[2]+c3[2]+c4[2])/3.0])
    atoms.append([type,(c2[0]+c3[0]+c4[0])/3.0,(c2[1]+c3[1]+c4[1])/3.0,
                  (c2[2]+c3[2]+c4[2])/3.0])

    volume = sqrt(2)/12 * a*a*a
    return atoms,[],tris,[],[],volume

  # --------------------------------------------------------------------
  # params = Alo,Ahi,Blo,Bhi,Clo,Chi,type
  # Alo to Ahi = bounds of edge length in x of box
  # Blo to Bhi = bounds of edge length in y of box
  # Clo to Chi = bounds of edge length in y of box
  # type = type of each vertex in rectangular box
  
  def tribox(self,*params):
    self.style = "tri"
    alo = params[0]
    ahi = params[1]
    blo = params[2]
    bhi = params[3]
    clo = params[4]
    chi = params[5]
    type = params[6]
    
    a = alo + self.random()*(ahi-alo)
    b = blo + self.random()*(bhi-blo)
    c = clo + self.random()*(chi-clo)

    atoms = []
    tris = []

    # 8 corner pts of rectangular box

    c1 = [-a/2,-b/2,-c/2]
    c2 = [+a/2,-b/2,-c/2]
    c3 = [-a/2,+b/2,-c/2]
    c4 = [+a/2,+b/2,-c/2]
    c5 = [-a/2,-b/2,+c/2]
    c6 = [+a/2,-b/2,+c/2]
    c7 = [-a/2,+b/2,+c/2]
    c8 = [+a/2,+b/2,+c/2]

    # 12 triangles on faces of cube, outward normals via rh rule

    tris.append(c1+c7+c3)
    tris.append(c1+c5+c7)
    tris.append(c2+c4+c8)
    tris.append(c2+c8+c6)
    tris.append(c1+c2+c6)
    tris.append(c1+c6+c5)
    tris.append(c3+c8+c4)
    tris.append(c3+c7+c8)
    tris.append(c1+c4+c2)
    tris.append(c1+c3+c4)
    tris.append(c5+c6+c8)
    tris.append(c5+c8+c7)

    # 12 atoms at centroids of 12 triangles

    atoms.append([type,(c1[0]+c7[0]+c3[0])/3.0,(c1[1]+c7[1]+c3[1])/3.0,
                  (c1[2]+c7[2]+c3[2])/3.0])
    atoms.append([type,(c1[0]+c5[0]+c7[0])/3.0,(c1[1]+c5[1]+c7[1])/3.0,
                  (c1[2]+c5[2]+c7[2])/3.0])
    atoms.append([type,(c2[0]+c4[0]+c8[0])/3.0,(c2[1]+c4[1]+c8[1])/3.0,
                  (c2[2]+c4[2]+c8[2])/3.0])
    atoms.append([type,(c2[0]+c8[0]+c6[0])/3.0,(c2[1]+c8[1]+c6[1])/3.0,
                  (c2[2]+c8[2]+c6[2])/3.0])
    atoms.append([type,(c1[0]+c2[0]+c6[0])/3.0,(c1[1]+c2[1]+c6[1])/3.0,
                  (c1[2]+c2[2]+c6[2])/3.0])
    atoms.append([type,(c1[0]+c6[0]+c5[0])/3.0,(c1[1]+c6[1]+c5[1])/3.0,
                  (c1[2]+c6[2]+c5[2])/3.0])
    atoms.append([type,(c3[0]+c8[0]+c4[0])/3.0,(c3[1]+c8[1]+c4[1])/3.0,
                  (c3[2]+c8[2]+c4[2])/3.0])
    atoms.append([type,(c3[0]+c7[0]+c8[0])/3.0,(c3[1]+c7[1]+c8[1])/3.0,
                  (c3[2]+c7[2]+c8[2])/3.0])
    atoms.append([type,(c1[0]+c4[0]+c2[0])/3.0,(c1[1]+c4[1]+c2[1])/3.0,
                  (c1[2]+c4[2]+c2[2])/3.0])
    atoms.append([type,(c1[0]+c3[0]+c4[0])/3.0,(c1[1]+c3[1]+c4[1])/3.0,
                  (c1[2]+c3[2]+c4[2])/3.0])
    atoms.append([type,(c5[0]+c6[0]+c8[0])/3.0,(c5[1]+c6[1]+c8[1])/3.0,
                  (c5[2]+c6[2]+c8[2])/3.0])
    atoms.append([type,(c5[0]+c8[0]+c7[0])/3.0,(c5[1]+c8[1]+c7[1])/3.0,
                  (c5[2]+c8[2]+c7[2])/3.0])

    volume = a*b*c
    return atoms,[],tris,[],[],volume

  # --------------------------------------------------------------------
  # params = Alo,Ahi,Blo,Bhi,type
  # Alo to Ahi = bounds of edge length in x of rectangle
  # Blo to Bhi = bounds of edge length in y of rectangle
  # type = type of each line segment in rectangle
  
  def linebox(self,*params):
    self.style = "line"
    alo = float(params[0])
    ahi = float(params[1])
    blo = float(params[2])
    bhi = float(params[3])
    type = params[4]

    a = alo + self.random()*(ahi-alo)
    b = blo + self.random()*(bhi-blo)
    
    atoms = []
    segments = []

    # 4 atoms at center points of line segments in 2d

    atoms.append([type,0,-b/2,0])
    atoms.append([type,a/2,0,0])
    atoms.append([type,0,b/2,0])
    atoms.append([type,-a/2,0,0])

    # 4 line segments as displacements from atom

    segments.append([-a/2,0,a/2,0])
    segments.append([0,-b/2,0,b/2])
    segments.append([a/2,0,-a/2,0])
    segments.append([0,b/2,0,-b/2])

    volume = a*b
    return atoms,[],[],segments,[],volume

  # --------------------------------------------------------------------
  # params = Alo,Ahi,Blo,Bhi,type
  # Alo to Ahi = bounds of base length in x of triangle
  # Blo to Bhi = bounds of heigth in y of isosceles triangle
  # type = type of each line segment in triangle
  
  def linetri(self,*params):
    self.style = "line"
    alo = float(params[0])
    ahi = float(params[1])
    blo = float(params[2])
    bhi = float(params[3])
    type = params[4]

    base = alo + self.random()*(ahi-alo)
    ht = blo + self.random()*(bhi-blo)
    
    atoms = []
    segments = []

    # 3 atoms at center points of line segments in 2d

    atoms.append([type,0,-ht/2,0])
    atoms.append([type,base/4,0,0])
    atoms.append([type,-base/4,0,0])

    # 3 line segments as displacements from atom

    segments.append([-base/2,0,base/2,0])
    segments.append([base/4,-ht/2,-base/4,ht/2])
    segments.append([base/4,ht/2,-base/4,-ht/2])

    volume = 0.5*base*ht
    return atoms,[],[],segments,[],volume

  # --------------------------------------------------------------------
  # params = Nlo,Nhi,type
  # polygon with random N sub-particles from Nlo to Nhi
  # type = type of each entire body particle
  
  def bodypgon(self,*params):
    self.style = "body"
    nlo = int(params[0])
    nhi = int(params[1])
    type = params[2]

    n = int(nlo + self.random()*(nhi+1-nlo))

    atoms = []
    bodies = []

    # 1 body atom at COM = (0,0,0)
    # append sub-particle count as mass

    atoms.append([type,0,0,0,n])

    # N sub-particles around perimeter of polygon in xy plane
    # rad set so that successive sub-particles are distance 1.0 apart
    # adjust sub-particle coords so COM is at (0,0,0)

    rad = sqrt(0.5/(1-cos(2.0*pi/n)))
    for i in xrange(n):
      theta = 2.0*pi * i/n
      bodies.append([rad*cos(theta),rad*sin(theta),0])

    com = [0,0]
    for body in bodies:
      com[0] += body[0]
      com[1] += body[1]
    com[0] /= n
    com[1] /= n
    for body in bodies:
      body[0] -= com[0]
      body[1] -= com[1]

    volume = n
    return atoms,[],[],[],bodies,volume

  # --------------------------------------------------------------------

  def random(self):
    k = self.seed/IQ
    self.seed = IA*(self.seed-k*IQ) - IR*k
    if self.seed < 0:
      self.seed += IM
    return AM*self.seed

# --------------------------------------------------------------------
# random number generator class

IM = 2147483647
AM = 1.0/IM
IA = 16807
IQ = 127773
IR = 2836

# --------------------------------------------------------------------
# push atom onto sphere surface of diam and return [type,x,y,z]

def atom_on_sphere(diam,type,x,y,z):
  if x == 0.0 and y == 0.0 and z == 0.0: scale = 0.0
  else: scale = 0.5*diam / sqrt(x*x + y*y + z*z)
  return [type,scale*x,scale*y,scale*z]

# --------------------------------------------------------------------
# make a sphere from a template with some atoms in patches w/ different types

def make_sphere(template,diam,type0,patches):
  atoms = []
  n = 0
  for atom in template:
    n += 1
    type = type0
    for pair in patches:
      patch = pair[1]
      if n in patch:
        type = pair[0]
    atoms.append(atom_on_sphere(diam,type,atom[0],atom[1],atom[2]))
  return atoms

# --------------------------------------------------------------------
# build a tether of length M of mtype, connected to atom1 with prev atom0
# blen = bond length between successive tether atoms
# dmin = length restriction between atoms i-1,i+1

def tether(m,mtype,blen,dmin,atom0,atom1,random):
  atoms = [atom0,atom1]
  for i in range(m):
    imonomer = len(atoms)
    restriction = True
    while restriction:
      rsq = 2.0
      while rsq > 1.0:
        dx = 2.0*random() - 1.0
        dy = 2.0*random() - 1.0
        dz = 2.0*random() - 1.0
        rsq = dx*dx + dy*dy + dz*dz
      r = sqrt(rsq)
      dx,dy,dz = dx/r,dy/r,dz/r
      x = atoms[imonomer-1][1] + dx*blen
      y = atoms[imonomer-1][2] + dy*blen
      z = atoms[imonomer-1][3] + dz*blen
      dx = x - atoms[imonomer-2][1]
      dy = y - atoms[imonomer-2][2]
      dz = z - atoms[imonomer-2][3]
      restriction = False
      if sqrt(dx*dx + dy*dy + dz*dz) <= dmin: restriction = True
      if not restriction: atoms.append([mtype,x,y,z])
  return atoms[2:]

# --------------------------------------------------------------------
# templates

TRI5_HOLLOW = ((0,0,0),(1,0,0),(2,0,0),(3,0,0),(4,0,0),
               (0.5,sqrt(3)/2,0),(3.5,sqrt(3)/2,0),
               (1.0,2*sqrt(3)/2,0),(3.0,2*sqrt(3)/2,0),
               (1.5,3*sqrt(3)/2,0),(2.5,3*sqrt(3)/2,0),
               (2.0,4*sqrt(3)/2,0))
               
SIMPLE_7 = ((0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))

# C60 with added center point at end

BUCKY_60 = ((0.014416695,    1.232517,      -3.260650),
            (1.236417,       1.719517,      -2.768650),
            (2.292417,       0.8265171,     -2.492650),
            (2.122416,      -0.5504832,     -2.709650),
            (0.8954167,     -1.040483,      -3.204650),
            (-0.1565833,     -0.1504831,     -3.479650),
            (-1.018583,       1.984517,      -2.678650),
            (-0.4355836,      2.936517,      -1.826650),
            (0.9584165,      2.773517,      -1.881650),
            (1.734416,       2.937517,      -0.7156501),
            (3.065417,       0.9905167,     -1.331650),
            (2.790417,      -1.238483,      -1.682650),
            (2.234416,      -2.418483,      -1.146650),
            (1.012417,      -2.905483,      -1.638650),
            (0.3414168,     -2.215483,      -2.669650),
            (-1.052583,      -2.051483,      -2.614650),
            (-1.360583,      -0.7754831,     -3.114650),
            (-2.397583,      -0.020483017,   -2.529650),
            (-2.227583,       1.356517,      -2.312650),
            (-1.059583,       3.265517,      -0.6046500),
            (-2.262583,       2.640517,      -0.2406499),
            (-2.848583,       1.684517,      -1.095650),
            (-3.402583,       0.5095167,     -0.5616500),
            (-3.124583,      -0.5444832,     -1.447650),
            (-2.815583,      -1.825483,      -0.9456501),
            (-1.781583,      -2.577483,      -1.527650),
            (-1.113584,      -3.265483,      -0.5006499),
            (0.2864165,     -3.429483,      -0.5566499),
            (3.373417,      -0.2854834,     -0.8306501),
            (1.018416,      -1.984483,       2.678350),
            (2.227417,      -1.356483,       2.312350),
            (2.398417,       0.020516872,    2.530350),
            (1.360417,       0.7755170,      3.114350),
            (0.1564164,      0.1505170,      3.479350),
            (-1.236583,      -1.719483,       2.768350),
            (-0.9585834,     -2.773483,       1.882350),
            (0.4354167,     -2.936483,       1.826350),
            (1.059417,      -3.265483,       0.6053500),
            (2.263417,      -2.640483,       0.2403500),
            (2.848417,      -1.684483,       1.096350),
            (3.124417,       0.5445170,      1.447350),
            (2.815417,       1.825517,       0.9453500),
            (1.781417,       2.577517,       1.527350),
            (1.052417,       2.051517,       2.614350),
            (-0.3415833,      2.215517,       2.669350),
            (-0.8955836,      1.040517,       3.204350),
            (-2.122583,       0.5505171,      2.710350),
            (-2.292583,      -0.8264832,      2.492350),
            (-1.734583,      -2.937483,       0.7163500),
            (-2.786583,      -2.048483,       0.4413500),
            (-3.065583,      -0.9904833,      1.331350),
            (-3.373584,       0.2855167,      0.8313500),
            (-2.790584,       1.238517,       1.683350),
            (-2.233583,       2.418517,       1.146350),
            (-1.011583,       2.905517,       1.638350),
            (-0.2855835,      3.429517,       0.5563500),
            (1.113417,       3.265517,       0.5003500),
            (3.402417,      -0.5094833,      0.5613500),
            (2.786417,       2.047517,      -0.4416499),
            (-0.014583588,   -1.232483,       3.261350),
            (0.0,            0.0,             0.0))

# --------------------------------------------------------------------
# C80

BUCKY_80 = ((-1.243762,     -0.7016125,       3.734700),
            (-1.248762,      0.6973875,       3.749700),
            (-0.054762483,   1.429388,       3.749700),
            (1.189238,      0.7893875,       3.735700),
            (1.194237,     -0.6086125,       3.783700),
            (0.00023752451,  -1.340613,       3.782700),
            (-2.204762,       1.428387,       3.034700),
            (-3.191762,      0.7883875,       2.276700),
            (-3.223763,     -0.6096125,       2.311700),
            (-2.267762,      -1.340613,       3.027700),
            (0.2282375,      -2.540613,       3.098700),
            (-0.7787625,      -3.147613,       2.340700),
            (-2.039762,      -2.541613,       2.343700),
            (2.160238,      -1.356613,       3.098700),
            (1.563237,      -2.550612,       2.675700),
            (-0.2727625,       2.612388,       3.034700),
            (-1.601763,       2.612388,       2.591700),
            (-3.586762,      -1.357612,       1.184700),
            (-2.854763,      -2.551612,       1.204700),
            (2.223238,       1.412387,       3.028700),
            (2.006238,       2.595387,       2.312700),
            (0.7452375,       3.202387,       2.276700),
            (0.3742375,       3.844388,       1.090700),
            (-0.9547625,       3.844388,      0.6487002),
            (-1.962763,       3.201387,       1.374700),
            (-2.992763,       2.595387,      0.6477003),
            (-3.595763,       1.411387,       1.090700),
            (-3.958763,      0.6633875,     -0.036299706),
            (-3.930763,     -0.7356125,     -0.020299911),
            (-3.568762,      -1.385612,      -1.205300),
            (-2.837763,      -2.579613,      -1.186300),
            (-2.439763,      -3.168612,      0.019700050),
            (-1.205762,      -3.828612,      0.035700321),
            (-0.3907624,      -3.817612,       1.174700),
            (0.9442375,      -3.827612,      0.7516999),
            (1.941237,      -3.168612,       1.478700),
            (3.190238,      0.6643875,       2.343700),
            (3.158237,     -0.7346125,       2.340700),
            (3.579237,      -1.384613,       1.175700),
            (2.982238,      -2.578612,      0.7527003),
            (1.244238,      0.7023875,      -3.735300),
            (1.249238,     -0.6976125,      -3.750300),
            (0.055237532,  -1.428612,      -3.750300),
            (-1.189762,     -0.7896125,      -3.735300),
            (-1.193763,      0.6093875,      -3.783300),
            (0.00023752451,   1.340387,      -3.783300),
            (2.205238,      -1.428612,      -3.034300),
            (3.192237,     -0.7886125,      -2.276300),
            (3.224237,      0.6093875,      -2.312300),
            (2.268238,       1.341388,      -3.027300),
            (-0.2277625,       2.541388,      -3.098300),
            (0.7792375,       3.147388,      -2.340300),
            (2.040237,       2.541388,      -2.343300),
            (-2.159763,       1.356387,      -3.098300),
            (-1.562762,       2.550387,      -2.675300),
            (0.2732375,      -2.612613,      -3.034300),
            (1.601238,      -2.612613,      -2.592300),
            (3.586237,       1.357388,      -1.185300),
            (2.854238,       2.551387,      -1.204300),
            (-2.223763,      -1.411613,      -3.028300),
            (-2.005763,      -2.595613,      -2.312300),
            (-0.7457625,      -3.201612,      -2.277300),
            (-0.3747625,      -3.844613,      -1.090300),
            (0.9542375,      -3.844613,     -0.6482999),
            (1.962238,      -3.201612,      -1.375300),
            (2.992238,      -2.595613,     -0.6482999),
            (3.596237,      -1.411613,      -1.090300),
            (3.958237,     -0.6636125,      0.036700249),
            (3.931237,      0.7353874,      0.019700050),
            (3.569237,       1.385387,       1.205700),
            (2.837237,       2.579388,       1.185700),
            (2.439238,       3.168387,     -0.019299984),
            (1.205238,       3.828387,     -0.036299706),
            (0.3902375,       3.818388,      -1.175300),
            (-0.9447625,       3.828387,     -0.7522998),
            (-1.941762,       3.168387,      -1.478300),
            (-3.189763,     -0.6646125,      -2.344300),
            (-3.157763,      0.7343875,      -2.340300),
            (-3.579762,       1.384387,      -1.175300),
            (-2.982763,       2.578387,     -0.7522998))
