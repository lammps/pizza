# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# how to specify and print 2d particles

# cdata tool

oneline = "Read, create, manipulate ChemCell data files"

docstr = """
c = cdata()			   create a datafile object
c = cdata("mem.surf")              read in one or more ChemCell data files
c = cdata("mem.part.gz mem.surf")  can be gzipped
c = cdata("mem.*")		   wildcard expands to multiple files
c.read("mem.surf")		   read in one or more data files

  read() has same argument options as constructor
  files contain the following kinds of entries, each of which becomes an object
    particles, triangles, region, facets
    particles is a list of particles -> becomes a group
    triangles is 3 lists of vertices, triangles, connections -> becomes a surf
    region is a ChemCell command defining a region -> becomes a region
    facets is a CUBIT format of vertices and triangles -> becomes a surf
  each object is assigned an ID = name in file
  ID can be any number or string, must be unique

c.box(ID,xlo,ylo,zlo,xhi,yhi,zhi)  create a box region
c.sphere(ID,x,y,z,r)		   create a sphere region
c.shell(ID,x,y,z,r,rinner)	   create a shell region
c.cyl(ID,'x',c1,c2,r,lo,hi)	   create a axis-aligned cylinder region
c.cap(ID,'x',c1,c2,r,lo,hi)	   create a axis-aligned capped-cylinder region
c.q(ID,q1,q2,...)                  set region triangulation quality factors

  box() can create an axis-aligned plane, line, or point if lo=hi
  cyl() can create an axis-aligned circle if lo=hi
  for cyl() and cap(): 'x' c1,c2 = y,z; 'y' c1,c2 = x,z; 'z' c,c2 = x,y
  q's are size factors for region triangulation
    for box, q1,q2,q3 = # of divisions per xyz of box
    for sphere or shell, q1 = # of divisions per face edge of embedded cube
    for cyl or cap, q1 = # of divisions per face edge of end cap, must be even
                    q2 = # of divisions along length of cylinder

c.line(ID,x1,y1,z1,x2,y2,z2)       create a line object with one line
c.lbox(ID,xlo,ylo,zlo,xhi,yhi,zhi) create a line object with 12 box lines

c.surf(ID,id-region)               create a triangulated surf from a region
c.surftri(ID,id-surf,t1,t2,...)    create a tri surf from list of id-surf tris
c.surfselect(ID,id-surf,test)      create a tri surf from test on id-surf tris
c.bins(ID,nx,ny)                   set binning parameters for a surf

  triangulation of a shell is just done for the outer sphere
  for surftri(), one or more tri indices (1-N) must be listed
  for surfselect(), test is string like "$x < 2.0 and $y > 0.0"
  bins are used when particles are created inside/outside a surf
  
c.part(ID,n,id_in)  	           create N particles inside object id_in
c.part(ID,n,id_in,id_out)	   particles are also outside object id_out
c.part2d(ID,n,id_on)               create 2d particles on object id_on
c.partarray(ID,nx,nz,nz,x,y,z,dx,dy,dz)   create 3d grid of particles
c.partring(ID,n,x,y,z,r,'x')              create ring of particles
c.partsurf(ID,id_on)               change surf of existing 2d particle group
c.seed(43284)			   set random # seed (def = 12345)

  generate particle positions randomly (unless otherwise noted)
  for part(), id_in and id_out must be IDs of a surf, region, or union object
    inside a union object means inside any of the lower-level objects
    outside a union object means outside all of the lower-level objects
  for part2d(), id_on must be ID of a surf, region, or union object
  for part2d(), particles will be written as 2d assigned to surf id_on
  for partring(), ring axis is in 'x','y', or 'z' direction
  partsurf() changes surf id_on for an existing 2d particle group

x,n = c.random(ID)                 pick a random pt on surf of object ID
c.project(ID,ID2,dx,dy,dz,eps,fg)  project particles in ID to surf of obj ID2

  random() returns pt = [x,y,z] and normal vec n [nx,ny,nz]
  for random(), ID can be surf or region obj
  project() remaps particle coords in group ID
    moves each particle along dir until they are within eps of surface
    if no fg arg, dir = (dx,dy,dz)
    if fg arg, dir = line from particle coord to (dx,dy,dz)
    ID2 can be surf or region obj
    particles are converted to 2d assigned to surf ID2

c.center(ID,x,y,z)                 set center point of object
c.trans(ID,dx,dy,dz)   	 	   translate an object
c.rotate(ID,'x',1,1,0,'z',-1,1,0)  rotate an object
c.scale(ID,sx,sy,sz)		   scale an object

  objects must be surface or particle group, regions cannot be changed
  for center(), default is middle of bounding box (set when obj is created)
  for rotate(), set any 2 axes, must be orthogonal, 3rd is inferred
    object is rotated so that it's current xyz axes point along new ones
  rotation and scaling occur relative to center point

c.union(ID,id1,id2,...)		   create a new union object from id1,id2,etc
c.join(ID,id1,id2,...)             create a new object by joining id1,id2,etc
c.delete(id1,id2,...)              delete one or more objects
c.rename(ID,IDnew)                 rename an object
c.copy(ID,IDnew) 	           create a new object as copy of old object

  for union, all lower-level objects must be of surface, region, or union style
  for join, all joined objects must be of same style: group, surf, line
    new object is the same style

c.select(id1,id2,...)              select one or more objects
c.select()                         select all objects
c.unselect(id1,id2,...)            unselect one or more objects
c.unselect()                       unselect all objects

  selection applies to write() and viz()
  
c.write("file")			   write all selected objs to ChemCell file
c.write("file",id1,id2,...)	   write only listed & selected objects to file
c.append("file")		   append all selected objs to ChemCell file
c.append("file",id1,id2,...)	   append only listed & selected objects

  union objects are skipped, not written to file
  
index,time,flag = c.iterator(0/1)          loop over single snapshot
time,box,atoms,bonds,tris,lines = c.viz(index)   return list of viz objects

  iterator() and viz() are compatible with equivalent dump calls
  iterator() called with arg = 0 first time, with arg = 1 on subsequent calls
    index = timestep index within dump object (only 0 for data file)
    time = timestep value (only 0 for data file)
    flag = -1 when iteration is done, 1 otherwise
  viz() returns info for selected objs for specified timestep index (must be 0)
    time = 0
    box = [xlo,ylo,zlo,xhi,yhi,zhi]
    atoms = id,type,x,y,z for each atom as 2d array
      NULL if atoms do not exist
    bonds = NULL
    tris = id,type,x1,y1,z1,x2,y2,z2,x3,y3,z3,nx,ny,nz for each tri as 2d array
      regions are triangulated according to q() settings by viz()
      NULL if surfaces do not exist
    lines = id,type,x1,y1,z1,x2,y2,z2 for each line as 2d array
      NULL if lines do not exist
    types are assigned to each object of same style in ascending order
"""

# History
#   11/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   nselect = 1 = # of snapshots
#   ids = dictionary of IDs that points to object index
#   objs = list of objects, style = REGION, SURFACE, GROUP, UNION

# Imports and external programs

import sys, glob
from os import popen
from math import sqrt,pi,cos,sin,fabs
from copy import deepcopy

try: from DEFAULTS import PIZZA_GUNZIP
except: PIZZA_GUNZIP = "gunzip"

# Class definition

class cdata:

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.nselect = 1
    self.ids = {}
    self.objs = []
    self.random = Random(12345)

    if len(list): self.read(*list)

  # --------------------------------------------------------------------

  def read(self,*list):

    # flist = list of all data file names

    words = list[0].split()
    flist = []
    for word in words: flist += glob.glob(word)
    if len(flist) == 0 and len(list) == 1:
      raise StandardError,"no data file specified"

    for file in flist:

      # test for gzipped file

      if file[-3:] == ".gz":
        f = popen("%s -c %s" % (PIZZA_GUNZIP,file),'r')
      else: f = open(file)

      # read all entries in file
      
      while 1:
        line = f.readline()
        if not line: break
        line = line.strip()
        if not line: continue
        elif line.find("triangles") == 0: flag = "triangles"
        elif line.find("particles") == 0: flag = "particles"
        elif line.find("facets") == 0: flag = "facets"
        elif line.find("region") == 0: flag = "region"
        else:
          print "unknown line:",line
          raise StandardError, "unrecognized ChemCell data file"

        # create a surface object from set of triangles or facets
        
        if flag == "triangles" or flag == "facets":
          tmp,id,nvert,ntri = line.split()
          nvert = int(nvert)
          ntri = int(ntri)

          if self.ids.has_key(id):
            raise StandardError,"ID %s is already in use" % id

          f.readline()
          vertices = []
          for i in xrange(nvert):
            list = f.readline().split()
            vertices.append([float(value) for value in list[1:]])
          f.readline()
          triangles = []
          for i in xrange(ntri):
            list = f.readline().split()
            triangles.append([int(value) for value in list[1:]])

          if flag == "triangles":
            f.readline()
            connections = []
            for i in xrange(ntri):
              list = f.readline().split()
              connections.append([int(value) for value in list[1:]])
          else:
            connections = connect(nvert,ntri,triangles)
            
          obj = Surface()
          obj.select = 1
          self.ids[id] = len(self.objs)
          self.objs.append(obj)
          obj.id = id
          obj.style = SURFACE
          obj.nvert = nvert
          obj.ntri = ntri
          obj.vertices = vertices
          obj.triangles = triangles
          obj.connections = connections
          obj.center()
	  
          print id,
          sys.stdout.flush()

        # create a group object from list of particles

        if flag == "particles":
          words = line.split()
          id = words[1]
          npart = int(words[2])

          if self.ids.has_key(id):
            raise StandardError,"ID %s is already in use" % id

          f.readline()
          xyz = []
          for i in xrange(npart):
            list = f.readline().split()
            xyz.append([float(value) for value in list[1:]])

          obj = Group()
          obj.select = 1
          self.ids[id] = len(self.objs)
          self.objs.append(obj)
          obj.id = id
          obj.style = GROUP
          if len(words) == 3: obj.on_id = ""
          else: obj.on_id = words[3]
          obj.npart = npart
          obj.xyz = xyz
          obj.center()
            
          print id,
          sys.stdout.flush()

        # create a region object from ChemCell region command

        if flag == "region":
          words = line.split()
          id = words[1]
          style = words[2]
          args = words[3:]
          
          if style == "box":
            obj = Box(*args)
            obj.substyle = BOX
          elif style == "sphere":
            obj = Sphere(*args)
            obj.substyle = SPHERE
          obj.select = 1
          self.ids[id] = len(self.objs)
          self.objs.append(obj)
          obj.id = id
          obj.style = REGION
            
          print id,
          sys.stdout.flush()

      f.close()
    print

  # --------------------------------------------------------------------
  # create box region

  def box(self,id,*args):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Box(*args)
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = REGION
    obj.substyle = BOX

  # --------------------------------------------------------------------
  # create sphere region

  def sphere(self,id,*args):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Sphere(*args)
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = REGION
    obj.substyle = SPHERE

  # --------------------------------------------------------------------
  # create shell region

  def shell(self,id,*args):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Shell(*args)
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = REGION
    obj.substyle = SHELL

  # --------------------------------------------------------------------
  # create cylinder region

  def cyl(self,id,*args):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Cylinder(*args)
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = REGION
    obj.substyle = CYLINDER

  # --------------------------------------------------------------------
  # create a capped-cylinder region

  def cap(self,id,*args):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Capped(*args)
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = REGION
    obj.substyle = CAPPED

  # --------------------------------------------------------------------
  # set quality factors for a region's triangulation routine

  def q(self,id,*args):
    obj = self.objs[self.ids[id]]
    if obj.style != REGION:
      raise StandardError,"Can only use q() on a region object"
    n = 1
    for arg in args:
      cmd = "obj.q%d = arg" % n
      exec cmd
      n += 1
      
  # --------------------------------------------------------------------
  # create a line object with single line

  def line(self,id,*args):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Line()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = LINE
    obj.nline = 0
    obj.pairs = []
    
    obj.addline(args)

  # --------------------------------------------------------------------
  # create a line object with 12 lines for a box

  def lbox(self,id,*args):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Line()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = LINE
    obj.nline = 0
    obj.pairs = []
    
    xlo,ylo,zlo,xhi,yhi,zhi = args
    obj.addline([xlo,ylo,zlo,xhi,ylo,zlo])
    obj.addline([xlo,yhi,zlo,xhi,yhi,zlo])
    obj.addline([xlo,yhi,zhi,xhi,yhi,zhi])
    obj.addline([xlo,ylo,zhi,xhi,ylo,zhi])
    obj.addline([xlo,ylo,zlo,xlo,yhi,zlo])
    obj.addline([xhi,ylo,zlo,xhi,yhi,zlo])
    obj.addline([xhi,ylo,zhi,xhi,yhi,zhi])
    obj.addline([xlo,ylo,zhi,xlo,yhi,zhi])
    obj.addline([xlo,ylo,zlo,xlo,ylo,zhi])
    obj.addline([xhi,ylo,zlo,xhi,ylo,zhi])
    obj.addline([xhi,yhi,zlo,xhi,yhi,zhi])
    obj.addline([xlo,yhi,zlo,xlo,yhi,zhi])
    
  # --------------------------------------------------------------------
  # create a triangulated surface object from a region object

  def surf(self,id,id_region):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id

    region = self.objs[self.ids[id_region]]
    region.triangulate()
    
    obj = Surface()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = SURFACE
    obj.nvert = region.nvert
    obj.ntri = region.ntri
    obj.vertices = deepcopy(region.vertices)
    obj.triangles = deepcopy(region.triangles)
    obj.connections = deepcopy(region.connections)
    obj.center()

  # --------------------------------------------------------------------
  # create a triangulated surface object
  # from list of tri indices (1-N) in id_surf obj

  def surftri(self,id,id_surf,*list):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id

    o = self.objs[self.ids[id_surf]]
    
    obj = Surface()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = SURFACE
    obj.nvert = 0
    obj.ntri = 0
    obj.vertices = []
    obj.triangles = []

    # subtract 1 from tri and vert to convert to C indexing from (1-N)
    
    for i in list:
      v1 = o.triangles[i-1][0]
      v2 = o.triangles[i-1][1]
      v3 = o.triangles[i-1][2]
      obj.vertices.append(o.vertices[v1-1][:])
      obj.vertices.append(o.vertices[v2-1][:])
      obj.vertices.append(o.vertices[v3-1][:])
      obj.triangles.append([obj.nvert+1,obj.nvert+2,obj.nvert+3])
      obj.nvert += 3
      obj.ntri += 1

    # make any connections in new set of triangles
    
    obj.connections = connect(obj.nvert,obj.ntri,obj.triangles)
    obj.center()

  # --------------------------------------------------------------------
  # create a triangulated surface object
  # from subset of tris in id_surf obj that meet test string

  def surfselect(self,id,id_surf,teststr):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id

    o = self.objs[self.ids[id_surf]]
    
    obj = Surface()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = SURFACE
    obj.nvert = 0
    obj.ntri = 0
    obj.vertices = []
    obj.triangles = []

    # replace $var with o.vertices reference and compile test string

    cmd1 = teststr.replace("$x","o.vertices[v1][0]")
    cmd1 = cmd1.replace("$y","o.vertices[v1][1]")
    cmd1 = "flag1 = " + cmd1.replace("$z","o.vertices[v1][2]")
    ccmd1 = compile(cmd1,'','single')

    cmd2 = teststr.replace("$x","o.vertices[v2][0]")
    cmd2 = cmd2.replace("$y","o.vertices[v2][1]")
    cmd2 = "flag2 = " + cmd2.replace("$z","o.vertices[v2][2]")
    ccmd2 = compile(cmd2,'','single')

    cmd3 = teststr.replace("$x","o.vertices[v3][0]")
    cmd3 = cmd3.replace("$y","o.vertices[v3][1]")
    cmd3 = "flag3 = " + cmd3.replace("$z","o.vertices[v3][2]")
    ccmd3 = compile(cmd3,'','single')

    # loop over triangles in id_surf
    # 3 vertices must satisfy all 3 tests for tri's inclusion in new surf obj
    
    for tri in o.triangles:
      v1 = tri[0] - 1
      v2 = tri[1] - 1
      v3 = tri[2] - 1
      exec ccmd1
      exec ccmd2
      exec ccmd3
      if flag1 and flag2 and flag3:
        obj.vertices.append(o.vertices[v1][:])
        obj.vertices.append(o.vertices[v2][:])
        obj.vertices.append(o.vertices[v3][:])
        obj.triangles.append([obj.nvert+1,obj.nvert+2,obj.nvert+3])
        obj.nvert += 3
        obj.ntri += 1

    # make any connections in new set of triangles

    obj.connections = connect(obj.nvert,obj.ntri,obj.triangles)
    obj.center()

  # --------------------------------------------------------------------
  # set binning parameters for a surface

  def bins(self,id,nx,ny):
    obj = self.objs[self.ids[id]]
    if obj.style != SURFACE:
      raise StandardError,"Can only set bins for surface"
    obj.nbinx = nx
    obj.nbiny = ny

  # --------------------------------------------------------------------
  # create a group object with npart particles
  # use inside and outside as constraints on particle coords

  def part(self,id,npart,in_id,out_id=None):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id

    obj = Group()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = GROUP
    obj.on_id = ""
    obj.npart = npart
    obj.xyz = []

    in_obj = self.objs[self.ids[in_id]]
    if out_id: out_obj = self.objs[self.ids[out_id]]

    # pre-process SURFACE objects to bin their triangles for faster searching
    
    if in_obj.style == SURFACE: in_obj.inside_prep()
    if out_id and out_obj.style == SURFACE: out_obj.inside_prep()

    # bounding box for generating points

    xlo,ylo,zlo,xhi,yhi,zhi = in_obj.bbox()
    xsize = xhi-xlo
    ysize = yhi-ylo
    zsize = zhi-zlo

    # generate particles until have enough that satisfy in/out constraints
    
    count = attempt = 0
    while count < npart:
      attempt += 1
      x = xlo + self.random() * xsize
      y = ylo + self.random() * ysize
      z = zlo + self.random() * zsize
      if not in_obj.inside(x,y,z): continue
      if out_id and out_obj.inside(x,y,z): continue
      obj.xyz.append([x,y,z])
      count += 1

    obj.center()
    print "Created %d particles in %d attempts" % (count,attempt)
    
  # --------------------------------------------------------------------
  # create a group object with npart 2d particles on surface of on_id object

  def part2d(self,id,npart,on_id):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id

    obj = Group()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = GROUP
    obj.on_id = on_id
    obj.npart = npart
    obj.xyz = []

    on_obj = self.objs[self.ids[on_id]]
    if on_obj.style != SURFACE and on_obj.style != REGION and \
           on_obj.style != UNION:
      raise StandardError,"Illegal ID to place particles on"
    totalarea = on_obj.area()
    
    for count in xrange(npart):
      area = self.random() * totalarea
      pt,norm = on_obj.loc2d(area,self.random)
      obj.xyz.append(pt)
    
    obj.center()
    print "Created %d particles on area of %g" % (npart,totalarea)

  # --------------------------------------------------------------------
  # create a 3d array of particles
  # array size = Nx by Ny by Nz
  # lower left corner of array = x,y,z
  # array spacing = dx,dy,dz
  
  def partarray(self,id,nx,ny,nz,x,y,z,dx,dy,dz):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id

    obj = Group()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = GROUP
    obj.on_id = ""
    obj.npart = nx * ny * nz
    obj.xyz = []

    for k in xrange(nz):
      znew = z + k*dz
      for j in xrange(ny):
        ynew = y + j*dy
        for i in xrange(nx):
          xnew = x + i*dx
          obj.xyz.append([xnew,ynew,znew])
          
    obj.center()
    print "Created %d particles" % (nx*ny*nz)

  # --------------------------------------------------------------------
  # create a ring of N particles
  # ring center and radius = x,y,z,r
  # ring axis = 'x' or 'y' or 'z'
  
  def partring(self,id,n,x,y,z,r,axis):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id

    obj = Group()
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = GROUP
    obj.on_id = ""
    obj.npart = n
    obj.xyz = []

    deltheta = 2.0*pi / n
    for i in xrange(n):
      if axis == 'x':
        xnew = x
        ynew = y + r * cos(i*deltheta)
        znew = z + r * sin(i*deltheta)
      elif axis == 'y':
        xnew = x + r * cos(i*deltheta)
        ynew = y
        znew = z + r * sin(i*deltheta)
      elif axis == 'z':
        xnew = x + r * cos(i*deltheta)
        ynew = y + r * sin(i*deltheta)
        znew = z
      obj.xyz.append([xnew,ynew,znew])
          
    obj.center()
    print "Created %d particles" % n

  # --------------------------------------------------------------------
  # change surface assignment for a 2d group of particles

  def partsurf(self,id,on_id):
    obj = self.objs[self.ids[id]]
    if obj.style != GROUP:
      raise StandardError,"Must use particle group with partsurf()"
    if obj.on_id == "":
      raise StandardError,"Must use partsurf() with 2d particles"
    obj.on_id = on_id

  # --------------------------------------------------------------------
  # set random number seed

  def seed(self,new_seed):
    self.random.seed = new_seed

  # --------------------------------------------------------------------
  # pick a random pt on surf of object and return pt and normal vec

  def random(self,id):
    obj = self.objs[self.ids[id]]
    if obj.style != SURFACE and obj.style != REGION:
      raise StandardError,"Must use surf or region with random()"

    totalarea = obj.area()
    area = self.random() * totalarea
    pt,norm = obj.loc2d(area,self.random)
    return pt,norm

  # --------------------------------------------------------------------
  # project particles in ID to surface of object ID2
  # direction to project depends on flag
  # if no flag, dir = (dx,dy,dz) from every particle
  # if yes flag, dir = particle coord to (dx,dy,dz)

  def project(self,id,id2,dx,dy,dz,EPS,flag=None):
    obj = self.objs[self.ids[id]]
    if obj.style != GROUP:
      raise StandardError,"Must use particle group as 1st obj of project()"
    obj_on = self.objs[self.ids[id2]]
    if obj_on.style != SURFACE and obj_on.style != REGION:
      raise StandardError,"Must use surf or region as 2nd obj of project()"

    # pre-process SURFACE to bin its triangles for faster searching

    if obj_on.style == SURFACE: obj_on.inside_prep()

    # for each particle, move it in dir from current location
    # move along dir until get within EPS of surf
    # factor = multiply bracketing distance by this amount each iteration
    # maxscale = max multiple of dir vector to bracket in each direction
    
    factor = 2
    maxscale = 10.0
    
    for i in xrange(obj.npart):
      x,y,z = obj.xyz[i]
      if flag: dir = [dx-x,dy-y,dz-z]
      else: dir = [dx,dy,dz]
      normalize(dir)
      
      # start = in/out at starting pt
      # stop = in/out at bracketing pt
      
      start = obj_on.inside(x,y,z)
      if start: stop = 0
      else: stop = 1

      # iterate to find bracketing point or until scale dist > maxdist
      # bracket pt = xyz +/- scale*dir
      # multiply scale by factor each iteration

      scale = EPS
      bracket = start
      while scale < maxscale:
        xnew,ynew,znew = x+scale*dir[0], y+scale*dir[1], z+scale*dir[2]
        bracket = obj_on.inside(xnew,ynew,znew)
        #print "BBB",xnew,ynew,znew,scale,bracket
        if bracket == stop: break
        xnew,ynew,znew = x-scale*dir[0], y-scale*dir[1], z-scale*dir[2]
        bracket = obj_on.inside(xnew,ynew,znew)
        #print "CCC",xnew,ynew,znew,scale,bracket
        if bracket == stop: break
        scale *= factor

      if bracket == start:
        raise StandardError,"Could not find bracket point for particle %d" % i

      # bisection search to zoom in to within EPS of surface
      # separation = distance between 2 points
      
      delx = xnew-x; dely = ynew-y; delz = znew-z
      separation = sqrt(delx*delx + dely*dely + delz*delz)
      while separation > EPS:
        xmid = 0.5 * (x + xnew)
        ymid = 0.5 * (y + ynew)
        zmid = 0.5 * (z + znew)
        value = obj_on.inside(xmid,ymid,zmid)
        if value == start:
          x = xmid; y = ymid; z = zmid
        else:
          xnew = xmid; ynew = ymid; znew = zmid
        delx = xnew-x; dely = ynew-y; delz = znew-z
        separation = sqrt(delx*delx + dely*dely + delz*delz)

      obj.xyz[i][0] = x
      obj.xyz[i][1] = y
      obj.xyz[i][2] = z
      
    obj.on_id = id2
    obj.center()

  # --------------------------------------------------------------------
  # set center pt of an object

  def center(self,id,x,y,z):
    obj = self.objs[self.ids[id]]
    if obj.style != SURFACE and obj.style != GROUP:
      raise StandardError,"Can only use center() on a surface or group object"
    obj.center(x,y,z)

  # --------------------------------------------------------------------
  # translate an object by dx,dy,dz displacement
  # add displacement to its vertices and center pt

  def trans(self,id,dx,dy,dz):
    obj = self.objs[self.ids[id]]
    if obj.style != SURFACE and obj.style != GROUP:
      raise StandardError,"Can only use trans() on a surface or group object"
    obj.xc += dx
    obj.yc += dy
    obj.zc += dz

    # apply translation to each vertex or part coord

    if obj.style == SURFACE:
      for i in xrange(obj.nvert):
        obj.vertices[i][0] += dx
        obj.vertices[i][1] += dy
        obj.vertices[i][2] += dz
    elif obj.style == GROUP:
      for i in xrange(obj.npart):
        obj.xyz[i][0] += dx
        obj.xyz[i][1] += dy
        obj.xyz[i][2] += dz

  # --------------------------------------------------------------------
  # rotate an object so current coord xyz axes align with new ones
  # rotate around center pt
  # old coord axes = 100, 010, 001
  # xyz new = new coord axes, orthonormalized
  # transformation: r_new = Lij r_old
  #   Lij = direction cosine of inew along jold where i,j = x,y,z
  #   direction cosine for 2 unit vecs is just cos(theta)

  def rotate(self,id,axis1,i1,j1,k1,axis2,i2,j2,k2):
    obj = self.objs[self.ids[id]]
    if obj.style != SURFACE and obj.style != GROUP:
      raise StandardError,"Can only use rotate() on a surface or group object"

    # create xyz old and new
    
    xnew = ynew = znew = None
    if axis1 == 'x': xnew = [i1,j1,k1]
    elif axis1 == 'y': ynew = [i1,j1,k1]
    elif axis1 == 'z': znew = [i1,j1,k1]
    else: raise StandardError,"Illegal rotate arguments"
    if axis2 == 'x': xnew = [i2,j2,k2]
    elif axis2 == 'y': ynew = [i2,j2,k2]
    elif axis2 == 'z': znew = [i2,j2,k2]
    else: raise StandardError,"Illegal rotate arguments"
    if not xnew: xnew = cross(ynew,znew)
    elif not ynew: ynew = cross(znew,xnew)
    elif not znew: znew = cross(xnew,ynew)
    else: raise StandardError,"Illegal rotate arguments"
    normalize(xnew)
    normalize(ynew)
    normalize(znew)

    # apply rotation matrix of direction cosines to each vertex or part coord

    if obj.style == SURFACE:
      for i in xrange(obj.nvert):
        x = obj.vertices[i][0] - obj.xc
        y = obj.vertices[i][1] - obj.yc
        z = obj.vertices[i][2] - obj.zc
        xn = xnew[0]*x + xnew[1]*y + xnew[2]*z
        yn = ynew[0]*x + ynew[1]*y + ynew[2]*z
        zn = znew[0]*x + znew[1]*y + znew[2]*z
        obj.vertices[i][0] = xn + obj.xc
        obj.vertices[i][1] = yn + obj.yc
        obj.vertices[i][2] = zn + obj.zc
    elif obj.style == GROUP:
      for i in xrange(obj.npart):
        x = obj.xyz[i][0] - obj.xc
        y = obj.xyz[i][1] - obj.yc
        z = obj.xyz[i][2] - obj.zc
        xn = xnew[0]*x + ynew[0]*y + znew[0]*z
        yn = xnew[1]*x + ynew[1]*y + znew[1]*z
        zn = xnew[2]*x + ynew[2]*y + znew[2]*z
        obj.xyz[i][0] = xn + obj.xc
        obj.xyz[i][1] = yn + obj.yc
        obj.xyz[i][2] = zn + obj.zc

  # --------------------------------------------------------------------
  # scale an object by sx,sy,sz factors
  # scale its vertices relative to center pt

  def scale(self,id,sx,sy,sz):
    obj = self.objs[self.ids[id]]
    if obj.style != SURFACE and obj.style != GROUP:
      raise StandardError,"Can only use scale() on a surface or group object"
    if obj.style == SURFACE:
      for i in xrange(obj.nvert):
        obj.vertices[i][0] = obj.xc + sx * (obj.vertices[i][0] - obj.xc)
        obj.vertices[i][1] = obj.yc + sy * (obj.vertices[i][1] - obj.yc)
        obj.vertices[i][2] = obj.zc + sz * (obj.vertices[i][2] - obj.zc)
    elif obj.style == GROUP:
      for i in xrange(obj.npart):
        obj.xyz[i][0] = obj.xc + sx * (obj.xyz[i][0] - obj.xc)
        obj.xyz[i][1] = obj.yc + sy * (obj.xyz[i][1] - obj.yc)
        obj.xyz[i][2] = obj.zc + sz * (obj.xyz[i][2] - obj.zc)
    
  # --------------------------------------------------------------------
  # create union object from list of other objects

  def union(self,id,*list):
    if self.ids.has_key(id):
      raise StandardError,"ID %s is already in use" % id
    obj = Union(self.ids,self.objs,*list)
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = UNION

  # --------------------------------------------------------------------
  # join objects in list to form a new object
  # only possible for group, surface, line objects
  
  def join(self,id,*list):
    style = self.objs[self.ids[list[0]]].style
    if style == GROUP: obj = Group()
    elif style == SURFACE: obj = Surface()
    elif style == LINE: obj = Line()
    else: raise StandardError,"Cannot perform join on these object styles"
    obj.select = 1
    self.ids[id] = len(self.objs)
    self.objs.append(obj)
    obj.id = id
    obj.style = style

    if style == GROUP:
      obj.on_id = self.objs[self.ids[list[0]]].on_id
      obj.npart = 0
      obj.xyz = []
    elif style == SURFACE:
      obj.nvert = obj.ntri = 0
      obj.vertices = []
      obj.triangles = []
      obj.connections = []
    elif style == LINE:
      obj.nline = 0
      obj.pairs = []
      
    for id in list:
      o = self.objs[self.ids[id]]
      if o.style != style:
        raise StandardError,"All joined objects must be of same style"

      # force deep copy of particle coords
      
      if style == GROUP:
        if o.on_id != obj.on_id:
          raise StandardError,"Particle group surfaces do not match"
        for i in xrange(o.npart):
          xyz = o.xyz[i][:]
          obj.xyz.append(xyz)
        obj.npart += o.npart
        obj.center()
        
      # force deep copy of triangle vertices and indices
      # increment vertex and triangle indices b/c now have previous surfaces
      
      elif style == SURFACE:
        for i in xrange(o.nvert):
          vert = o.vertices[i][:]
          obj.vertices.append(vert)
        for i in xrange(o.ntri):
          tri = [iv+obj.nvert for iv in o.triangles[i]]
          obj.triangles.append(tri)
          connect = o.connections[i][:]
          if connect[0]: connect[0] += obj.ntri
          if connect[2]: connect[2] += obj.ntri
          if connect[4]: connect[4] += obj.ntri
          obj.connections.append(connect)
        obj.nvert += o.nvert
        obj.ntri += o.ntri
        obj.center()

      # force deep copy of line pt pairs
      
      elif style == LINE:
        pairs = o.pairs[:]
        obj.pairs += (pairs)
        obj.nline += o.nline

  # --------------------------------------------------------------------
  # delete each object in list
  # reset values in ids since some indices are decremented
  
  def delete(self,*list):
    for id in list:
      i = self.ids[id]
      del self.ids[id]
      del self.objs[i]
      for key in self.ids.keys():
        j = self.ids[key]
        if j > i: self.ids[key] = j-1
        
  # --------------------------------------------------------------------
  # rename the ID of an object
  # check that new name doesn't already exist

  def rename(self,idold,idnew):
    i = self.ids[idold]
    if self.ids.has_key(idnew):
      raise StandardError,"ID %s is already in use" % id
    self.ids[idnew] = i
    self.objs[i].id = idnew
    del self.ids[idold]
    
  # --------------------------------------------------------------------
  # create a deep copy of an object and assign it a new ID
  # check that new name doesn't already exist

  def copy(self,idold,idnew):
    obj = deepcopy(self.objs[self.ids[idold]])
    obj.select = 1
    if self.ids.has_key(idnew):
      raise StandardError,"ID %s is already in use" % id
    self.ids[idnew] = len(self.objs)
    self.objs.append(obj)
    obj.id = idnew

  # --------------------------------------------------------------------
  # set selection flag for each object in list
  # if list is empty, select all
  
  def select(self,*list):
    if len(list) == 0: list = self.ids.keys()
    for id in list:
      obj = self.objs[self.ids[id]]
      obj.select = 1

  # --------------------------------------------------------------------
  # unset selection flag for each object in list
  # if list is empty, unselect all
  
  def unselect(self,*list):
    if len(list) == 0: list = self.ids.keys()
    for id in list:
      obj = self.objs[self.ids[id]]
      obj.select = 0

  # --------------------------------------------------------------------
  # write out list of objects to data file via filewrite()
  # open a new file
  
  def write(self,file,*list):
    if not len(list): vlist = range(len(self.objs))
    else:
      vlist = []
      for id in list: vlist.append(self.ids[id])

    fp = open(file,'w')
    self.filewrite(fp,vlist)
    fp.close()
    
  # --------------------------------------------------------------------
  # append list of objects to data file via filewrite()
  # open existing file for appending
  
  def append(self,file,*list):
    if not len(list): vlist = range(len(self.objs))
    else:
      vlist = []
      for id in list: vlist.append(self.ids[id])

    fp = open(file,'a')
    self.filewrite(fp,vlist)
    fp.close()

  # --------------------------------------------------------------------
  # write out list of objects to previously opened data file
  # for particles, write as 3d or 2d depending on on_id
  
  def filewrite(self,fp,vlist):
    for index in vlist:
      obj = self.objs[index]
      if not obj.select: continue
      if obj.style == GROUP:
        if not obj.on_id:
          print >>fp,"particles %s %d" % (obj.id,obj.npart)
        else:
          print >>fp,"particles %s %d %s" % (obj.id,obj.npart,obj.on_id)
        print >>fp
        for i in xrange(obj.npart):
          print >>fp,i+1,obj.xyz[i][0],obj.xyz[i][1],obj.xyz[i][2]
        print >>fp
      if obj.style == SURFACE:
        print >>fp,"triangles %s %d %d" % (obj.id,obj.nvert,obj.ntri)
        print >>fp
        for i in xrange(obj.nvert):
          print >>fp,i+1,obj.vertices[i][0],obj.vertices[i][1], \
                obj.vertices[i][2]
        print >>fp
        for i in xrange(obj.ntri):
          print >>fp,i+1,obj.triangles[i][0],obj.triangles[i][1], \
                obj.triangles[i][2]
        print >>fp
        for i in xrange(obj.ntri):
          print >>fp,i+1,obj.connections[i][0],obj.connections[i][1], \
                obj.connections[i][2],obj.connections[i][3], \
                obj.connections[i][4],obj.connections[i][5]
      if obj.style == REGION:
        print >>fp,"region %s" % obj.command()

  # --------------------------------------------------------------------
  # iterator called from other tools

  def iterator(self,flag):
    if flag == 0: return 0,0,1
    return 0,0,-1

  # --------------------------------------------------------------------
  # return list of atoms and triangles and lines to viz for cdata object

  def viz(self,isnap):
    if isnap:
      raise StandardError, "cannot call cdata.viz() with isnap != 0"
    
    # create atom list from sum of all particle groups
    # id = running count
    # type = running type of particle group

    id = itype = 0
    atoms = []
    for obj in self.objs:
      if obj.style != GROUP: continue
      if not obj.select: continue
      itype += 1
      for xyz in obj.xyz:
        id += 1
        atoms.append([id,itype,xyz[0],xyz[1],xyz[2]])

    # no bonds
    
    bonds = []

    # create triangle list from sum of all surfaces and regions
    # id = running count
    # type = type of set of tris

    id = itype = 0
    tris = []
    for obj in self.objs:
      if obj.style != SURFACE and obj.style != REGION: continue
      if not obj.select: continue
      if obj.style == REGION: obj.triangulate()
      itype += 1
      for tri in obj.triangles:
        v1 = obj.vertices[tri[0]-1]
        v2 = obj.vertices[tri[1]-1]
        v3 = obj.vertices[tri[2]-1]
        list = v1 + v2 + v3
        n = normal(list[0:3],list[3:6],list[6:9])
        id += 1
        tris.append([id,itype] + list + n)

    # create line list from sum of all line objects

    id = itype = 0
    lines = []
    for obj in self.objs:
      if obj.style != LINE: continue
      if not obj.select: continue
      itype += 1
      for pair in obj.pairs:
        id += 1
        lines.append([id,itype] + pair)
    
    return 0,self.bbox(),atoms,bonds,tris,lines

  # --------------------------------------------------------------------
  # time query from other tools

  def findtime(self,n):
    if n == 0: return 0
    raise StandardError, "no step %d exists" % (n)

  # --------------------------------------------------------------------
  # return box size

  def maxbox(self):
    return self.bbox()

  # --------------------------------------------------------------------
  # return box that bounds all selected objects

  def bbox(self):
    xlo = ylo = zlo = BIG
    xhi = yhi = zhi = -BIG
    for obj in self.objs:
      if not obj.select: continue
      xxlo,yylo,zzlo,xxhi,yyhi,zzhi = obj.bbox()
      if xxlo < xlo: xlo = xxlo
      if yylo < ylo: ylo = yylo
      if zzlo < zlo: zlo = zzlo
      if xxhi > xhi: xhi = xxhi
      if yyhi > yhi: yhi = yyhi
      if zzhi > zhi: zhi = zzhi

    return (xlo,ylo,zlo,xhi,yhi,zhi)

# --------------------------------------------------------------------
# object styles

EPSILON = 1.0e-6
BIG = 1.0e20

REGION =  1
SURFACE = 2
GROUP =   3
UNION =   4

BOX   =    5
SPHERE =   6
SHELL =    7
CYLINDER = 8
CAPPED =   9
LINE =     10

# --------------------------------------------------------------------
# random number generator class

IM = 2147483647
AM = 1.0/IM
IA = 16807
IQ = 127773
IR = 2836

class Random:
  def __init__(self,seed):
    self.seed = seed
    
  def __call__(self):
    k = self.seed/IQ
    self.seed = IA*(self.seed-k*IQ) - IR*k
    if self.seed < 0: self.seed += IM
    return AM*self.seed

# --------------------------------------------------------------------
# triangulated surface

class Surface:

  def __init__(self):
    self.nbinx = self.nbiny = 0
    
  # bounding box
  
  def bbox(self):
    list = [float(vert[0]) for vert in self.vertices]
    xlo = min(list)
    xhi = max(list)
    list = [float(vert[1]) for vert in self.vertices]
    ylo = min(list)
    yhi = max(list)
    list = [float(vert[2]) for vert in self.vertices]
    zlo = min(list)
    zhi = max(list)

    return (xlo,ylo,zlo,xhi,yhi,zhi)

  # set center point explicitly or set to middle of bounding box
  
  def center(self,*xyz):
    if len(xyz):
      self.xc = xyz[0]
      self.yc = xyz[1]
      self.zc = xyz[2]
    else:
      xlo,ylo,zlo,xhi,yhi,zhi = self.bbox()
      self.xc = 0.5 * (xlo + xhi)
      self.yc = 0.5 * (ylo + yhi)
      self.zc = 0.5 * (zlo + zhi)

  # bin 2d-projected triangles into xy plane for nbinx by nbiny bins
  # each bin overlayed by triangle's xy bounding box stores the triangle index
  # EPSILON insures that bins completely overlay surface bounding box
  # self.xlo, self.ylo, self.dxinv, self.dyinv, self.bin are used by inside()
  
  def inside_prep(self):
    self.xlo,self.ylo,self.zlo,self.xhi,self.yhi,self.zhi = self.bbox()
    xsize = self.xhi - self.xlo
    ysize = self.yhi - self.ylo
    if self.nbinx == 0 and self.nbiny == 0:
      binsize = sqrt(xsize*ysize / (self.ntri/2))
      self.nbinx = int(xsize/binsize)
      self.nbiny = int(ysize/binsize)
      if self.nbinx < 2: self.nbinx = 2
      if self.nbiny < 2: self.nbiny = 2
    self.dxinv = 1.0 / (xsize / self.nbinx + EPSILON)
    self.dyinv = 1.0 / (ysize / self.nbiny + EPSILON)

    print "Binning %d triangles into %d by %d bins ..." % \
          (self.ntri,self.nbinx,self.nbiny)
    
    self.bin = []
    for i in xrange(self.nbinx):
      self.bin.append(self.nbiny*[0])
      for j in xrange(self.nbiny):
        self.bin[i][j] = []
    
    for m in xrange(self.ntri):
      v1 = self.vertices[self.triangles[m][0]-1]
      v2 = self.vertices[self.triangles[m][1]-1]
      v3 = self.vertices[self.triangles[m][2]-1]
      xmin = min((v1[0],v2[0],v3[0]))
      xmax = max((v1[0],v2[0],v3[0]))
      ymin = min((v1[1],v2[1],v3[1]))
      ymax = max((v1[1],v2[1],v3[1]))
      ilo = int((xmin - self.xlo) * self.dxinv)
      ihi = int((xmax - self.xlo) * self.dxinv)
      jlo = int((ymin - self.ylo) * self.dyinv)
      jhi = int((ymax - self.ylo) * self.dyinv)
      for i in xrange(ilo,ihi+1):
        for j in xrange(jlo,jhi+1):
          self.bin[i][j].append(m)
          
    print "Done with binning"

  # check for inside assumes that surf is a closed set of triangles
  # consider a line segment from x,y,z to x,y,INF
  # project 3d triangles to xy plane
  # only need examine triangles in bin that point is in
  # if pt is outside bins, don't need to check
  # x,y,z is inside surf if line segment intersects an odd number of triangles
  # intersection test:
  #   is xy pt in bounding rectangle of tri in xy plane ?
  #   is xy pt inside tri in xy plane (including edges and vertices) ?
  #   compute ztri = z value of xy pt on plane of 3d tri
  #   if z of pt is <= ztri, then line segment intersects the tri
  
  def inside(self,x,y,z):
    ix = int((x - self.xlo) * self.dxinv)
    iy = int((y - self.ylo) * self.dyinv)
    if ix < 0 or ix >= self.nbinx or iy < 0 or iy >= self.nbiny: return 0
    n = len(self.bin[ix][iy])
    
    hit = 0
    for m in xrange(n):
      itri = self.bin[ix][iy][m]
      tri = self.triangles[itri]

      v1 = self.vertices[tri[0]-1]
      v2 = self.vertices[tri[1]-1]
      v3 = self.vertices[tri[2]-1]
      
      # is x,y in bounding box of 2d triangle ?
      
      if x < v1[0] and x < v2[0] and x < v3[0]: continue
      if x > v1[0] and x > v2[0] and x > v3[0]: continue
      if y < v1[1] and y < v2[1] and y < v3[1]: continue
      if y > v1[1] and y > v2[1] and y > v3[1]: continue

      # is x,y inside 2d triangle ?
      # cross product of each edge with vertex-to-point must have same sign
      # cross product = 0 is OK, means point is on edge or vertex
      
      c1 = (v2[0]-v1[0])*(y-v1[1]) - (v2[1]-v1[1])*(x-v1[0])
      c2 = (v3[0]-v2[0])*(y-v2[1]) - (v3[1]-v2[1])*(x-v2[0])
      c3 = (v1[0]-v3[0])*(y-v3[1]) - (v1[1]-v3[1])*(x-v3[0])
      if c1 < 0:
        if c2 > 0 or c3 > 0: continue
      elif c1 > 0:
        if c2 < 0 or c3 < 0: continue
      else:
        if c2 < 0 and c3 > 0: continue
        if c2 > 0 and c3 < 0: continue

      # represent x,y point in basis of 2 triangle edge vectors
      # p = (x,y) - v1, v = v2 - v1, w = v3 - v2
      # find alpha,beta such that p = alpha v + beta w
      #   by solving system of 2 linear eqs for their intersection
      # denom (vx - vy * wx/xy) cannot be 0 since would imply vx/vy = wx/wy
      #   and thus e1 would be parallel to e2
      # if wy = 0, vy cannot be 0, since e1 would be parallel to e2
      
      px = x - v1[0];
      py = y - v1[1]
      vx = v2[0] - v1[0]
      vy = v2[1] - v1[1]
      wx = v3[0] - v2[0]
      wy = v3[1] - v2[1]
      if wy: alpha = (px - py * wx/wy) / (vx - vy * wx/wy)
      else: alpha = py/vy
      if wx: beta = (px - alpha * vx) / wx
      else: beta = (py - alpha * vy) / wy
      pz = alpha * (v2[2] - v1[2]) + beta * (v3[2] - v2[2])
      if z <= pz + v1[2]: hit += 1

    return hit % 2

  # surface area
  # areas = cummulative total area of all triangles
  # triangle area = 1/2 of magnitude of cross product of 2 edge vectors
  
  def area(self):
    a = 3*[0]
    self.areas = []
    sum = 0.0
    for tri in self.triangles:
      v1 = self.vertices[tri[0]-1]
      v2 = self.vertices[tri[1]-1]
      v3 = self.vertices[tri[2]-1]
      a[0] = (v2[1]-v1[1])*(v3[2]-v1[2]) - (v2[2]-v1[2])*(v3[1]-v1[1])
      a[1] = (v2[2]-v1[2])*(v3[0]-v1[0]) - (v2[0]-v1[0])*(v3[2]-v1[2])
      a[2] = (v2[0]-v1[0])*(v3[1]-v1[1]) - (v2[1]-v1[1])*(v3[0]-v1[0])
      sum += 0.5 * sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
      self.areas.append(sum)
    return sum

  # return a random location on one of triangles
  
  def loc2d(self,area,random):
    for i in xrange(self.ntri):
      if area < self.areas[i]: break
    v1 = self.vertices[self.triangles[i][0]-1]
    v2 = self.vertices[self.triangles[i][1]-1]
    v3 = self.vertices[self.triangles[i][2]-1]
    r1 = random()
    r2 = random()
    if r2 > r1:
      tmp = r1
      r1 = r2
      r2 = tmp
    x = v1[0] + r1*(v2[0] - v1[0]) + r2*(v3[0] - v2[0])
    y = v1[1] + r1*(v2[1] - v1[1]) + r2*(v3[1] - v2[1])
    z = v1[2] + r1*(v2[2] - v1[2]) + r2*(v3[2] - v2[2])
    list = v1 + v2 + v3
    return [x,y,z],normal(list[0:3],list[3:6],list[6:9])

# --------------------------------------------------------------------
# group of particles

class Group:
  def bbox(self):
    list = [xyz[0] for xyz in self.xyz]
    xlo = min(list)
    xhi = max(list)
    list = [xyz[1] for xyz in self.xyz]
    ylo = min(list)
    yhi = max(list)
    list = [xyz[2] for xyz in self.xyz]
    zlo = min(list)
    zhi = max(list)

    return (xlo,ylo,zlo,xhi,yhi,zhi)

  # set center point explicitly or set to middle of bounding box
  
  def center(self,*xyz):
    if len(xyz):
      self.xc = xyz[0]
      self.yc = xyz[1]
      self.zc = xyz[2]
    else:
      xlo,ylo,zlo,xhi,yhi,zhi = self.bbox()
      self.xc = 0.5 * (xlo + xhi)
      self.yc = 0.5 * (ylo + yhi)
      self.zc = 0.5 * (zlo + zhi)

# --------------------------------------------------------------------
# box region

class Box:
  def __init__(self,*list):
    self.xlo = float(list[0])
    self.ylo = float(list[1])
    self.zlo = float(list[2])
    self.xhi = float(list[3])
    self.yhi = float(list[4])
    self.zhi = float(list[5])
    self.q1 = self.q2 = self.q3 = 1

  # bounding box around region
    
  def bbox(self):
    return (self.xlo,self.ylo,self.zlo,self.xhi,self.yhi,self.zhi)

  # return 1,0 for xyz inside/outside the region

  def inside(self,x,y,z):
    if x < self.xlo or x > self.xhi: return 0
    if y < self.ylo or y > self.yhi: return 0
    if z < self.zlo or z > self.zhi: return 0
    return 1

  # triangulate the region
  # set nvert,ntri,vertices,triangles,connections
  # convert vertices from unit box to lo/hi box
  # convert tuples returned by box_triangluate (used for dict lookup) to lists
  
  def triangulate(self):
    vertices,triangles = box_triangulate(self.q1,self.q2,self.q3)
    self.nvert = len(vertices)
    self.ntri = len(triangles)
    self.vertices = []
    for i in xrange(self.nvert):
      v1 = self.xlo + vertices[i][0]*(self.xhi-self.xlo)
      v2 = self.ylo + vertices[i][1]*(self.yhi-self.ylo)
      v3 = self.zlo + vertices[i][2]*(self.zhi-self.zlo)
      self.vertices.append([v1,v2,v3])
    self.triangles = []
    for i in xrange(self.ntri):
      self.triangles.append([triangles[i][0],triangles[i][1],triangles[i][2]])
    self.connections = connect(self.nvert,self.ntri,self.triangles)
    
  # surface area of region
  
  def area(self):
    xsize = self.xhi - self.xlo
    ysize = self.yhi - self.ylo
    zsize = self.zhi - self.zlo
    self.areas = []
    sum = ysize*zsize
    self.areas.append(sum)
    sum += ysize*zsize
    self.areas.append(sum)
    sum += xsize*zsize
    self.areas.append(sum)
    sum += xsize*zsize
    self.areas.append(sum)
    sum += xsize*ysize
    self.areas.append(sum)
    sum += xsize*ysize
    self.areas.append(sum)
    return sum

  # return a random location on surface of region
  
  def loc2d(self,area,random):
    xsize = self.xhi - self.xlo
    ysize = self.yhi - self.ylo
    zsize = self.zhi - self.zlo
    r1 = random()
    r2 = random()
    if area <= self.areas[0]:
      return [self.xlo,self.ylo + r1*ysize,self.zlo + r2*zsize],[-1,0,0]
    elif area <= self.areas[1]:
      return [self.xhi,self.ylo + r1*ysize,self.zlo + r2*zsize],[1,0,0]
    elif area <= self.areas[2]:
      return [self.xlo + r1*xsize,self.ylo,self.zlo + r2*zsize],[0,-1,0]
    elif area <= self.areas[3]:
      return [self.xlo + r1*xsize,self.yhi,self.zlo + r2*zsize],[0,1,0]
    elif area <= self.areas[4]:
      return [self.xlo + r1*xsize,self.ylo + r2*ysize,self.zlo],[0,0,-1]
    else:
      return [self.xlo + r1*xsize,self.ylo + r2*ysize,self.zhi],[0,0,1]

  # ChemCell text to create the region

  def command(self):
    return "%s box %g %g %g %g %g %g" % (self.id,self.xlo,self.ylo,self.zlo,
                                         self.xhi,self.yhi,self.zhi)

# --------------------------------------------------------------------
# sphere region

class Sphere:
  def __init__(self,*list):
    self.x = float(list[0])
    self.y = float(list[1])
    self.z = float(list[2])
    self.r = float(list[3])
    self.rsq = self.r*self.r
    self.q1 = 2

  # bounding box around region

  def bbox(self):
    return (self.x-self.r,self.y-self.r,self.z-self.r, \
            self.x+self.r,self.y+self.r,self.z+self.r)

  # return 1,0 for xyz inside/outside the region
    
  def inside(self,x,y,z):
    dx = x - self.x
    dy = y - self.y
    dz = z - self.z
    rsq = dx*dx + dy*dy + dz*dz
    if rsq > self.rsq: return 0
    return 1

  # triangulate the region
  # set nvert,ntri,vertices,triangles,connections
  # convert vertices from unit box to sphere at (x,y,z) with radius r
  # normalize() pushes vertices to surface of sphere
  # convert tuples returned by box_triangluate (used for dict lookup) to lists

  def triangulate(self):
    vertices,triangles = box_triangulate(self.q1,self.q1,self.q1)
    self.nvert = len(vertices)
    self.ntri = len(triangles)
    self.vertices = []
    for i in xrange(self.nvert):
      v1 = vertices[i][0] - 0.5
      v2 = vertices[i][1] - 0.5
      v3 = vertices[i][2] - 0.5
      c = [v1,v2,v3]
      normalize(c)
      c[0] = self.x + self.r*c[0]
      c[1] = self.y + self.r*c[1]
      c[2] = self.z + self.r*c[2]
      self.vertices.append(c)
    self.triangles = []
    for i in xrange(self.ntri):
      self.triangles.append([triangles[i][0],triangles[i][1],triangles[i][2]])
    self.connections = connect(self.nvert,self.ntri,self.triangles)
    
  # surface area of region
  
  def area(self):
    value = 4 * pi * self.r*self.r
    return value

  # return a random location on surface of region
  
  def loc2d(self,area,random):
    while 1:
      x = random() - 0.5
      y = random() - 0.5
      z = random() - 0.5
      if x*x + y*y + z*z <= 0.25: break
    c = [x,y,z]
    normalize(c)
    return [self.x + self.r*c[0], self.y + self.r*c[1], self.z + self.r*c[2]],c

  # ChemCell text to create the region
  
  def command(self):
    return "%s sphere %g %g %g %g" % (self.id,self.x,self.y,self.z,self.r)

# --------------------------------------------------------------------
# shell region
# inherits most functionality from Sphere, adds inner radius
# loc2d is only only outer surface, inherited from Sphere

class Shell(Sphere):
  def __init__(self,*list):
    self.rinner = float(list[4])
    self.innersq = self.rinner*self.rinner
    
  def inside(self,x,y,z):
    dx = x - self.x
    dy = y - self.y
    dz = z - self.z
    rsq = dx*dx + dy*dy + dz*dz
    if rsq > self.rsq or rsq < self.innersq: return 0
    return 1

  def command(self):
    return "%s shell %g %g %g %g %g" % (self.id,self.x,self.y,self.z,
                                        self.r,self.rinner)

# --------------------------------------------------------------------
# cylinder region

class Cylinder:
  def __init__(self,*list):
    self.axis = list[0]
    self.c1 = list[1]
    self.c2 = list[2]
    self.r = list[3]
    self.lo = list[4]
    self.hi = list[5]
    self.rsq = self.r*self.r
    self.q1 = 2
    self.q2 = 1

  # bounding box around region

  def bbox(self):
    if self.axis == 'x':
      return (self.lo,self.c1-self.r,self.c2-self.r, \
              self.hi,self.c1+self.r,self.c2+self.r)
    elif self.axis == 'y':
      return (self.c1-self.r,self.lo,self.c2-self.r, \
              self.c1+self.r,self.hi,self.c2+self.r)
    elif self.axis == 'z':
      return (self.c1-self.r,self.c2-self.r,self.lo, \
              self.c1+self.r,self.c2+self.r,self.hi)

  # return 1,0 for xyz inside/outside the region

  def inside(self,x,y,z):
    if self.axis == 'x':
      d1 = y - self.c1
      d2 = z - self.c2
      d3 = x
    elif self.axis == 'y':
      d1 = x - self.c1
      d2 = z - self.c2
      d3 = y
    elif self.axis == 'z':
      d1 = x - self.c1
      d2 = y - self.c2
      d3 = z
    rsq = d1*d1 + d2*d2
    if rsq > self.rsq: return 0
    if d3 < self.lo or d3 > self.hi: return 0
    return 1

  # triangulate the region
  # set nvert,ntri,vertices,triangles,connections
  # convert vertices from unit box to cylinder with correct axis
  # x,y have spacings for ends of cylinder, z has spacings for axis of cylinder
  # convert tuples returned by box_triangluate (used for dict lookup) to lists

  def triangulate(self):
    vertices,triangles = box_triangulate(self.q1,self.q1,self.q2)
    self.nvert = len(vertices)
    self.ntri = len(triangles)
    self.vertices = []
    for i in xrange(self.nvert):
      v1 = vertices[i][0] - 0.5
      v2 = vertices[i][1] - 0.5
      v3 = vertices[i][2]
      c = [v1,v2,0]
      normalize(c)
      if fabs(v1) >= fabs(v2): length = 2*fabs(v1)
      if fabs(v2) >= fabs(v1): length = 2*fabs(v2)
      c[0] *= length
      c[1] *= length
      p1 = self.c1 + self.r*c[0]
      p2 = self.c2 + self.r*c[1]
      p3 = self.lo + v3*(self.hi - self.lo)
      if self.axis == 'x': self.vertices.append([p3,p1,p2])
      elif self.axis == 'y': self.vertices.append([p1,p3,p2])
      elif self.axis == 'z': self.vertices.append([p1,p2,p3])
    self.triangles = []
    for i in xrange(self.ntri):
      self.triangles.append([triangles[i][0],triangles[i][1],triangles[i][2]])
    self.connections = connect(self.nvert,self.ntri,self.triangles)

  # surface area of region

  def area(self):
    self.areas = []
    sum = pi * self.rsq
    self.areas.append(sum)
    sum += pi * self.rsq
    self.areas.append(sum)
    sum += 2 * pi * self.r * (self.hi - self.lo)
    self.areas.append(sum)
    return sum

  # return a random location on surface of region

  def loc2d(self,area,random):
    if area <= self.areas[0]:
      while 1:
        r1 = random() - 0.5
        r2 = random() - 0.5
        if r1*r1 + r2*r2 < 0.25: break
      if self.axis == 'x':
        return [self.lo,self.c1 + 2*r1*self.r,self.c2 + 2*r2*self.r],[-1,0,0]
      elif self.axis == 'y':
        return [self.c1 + 2*r1*self.r,self.lo,self.c2 + 2*r2*self.r],[0,-1,0]
      elif self.axis == 'z':
        return [self.c1 + 2*r1*self.r,self.c2 + 2*r2*self.r,self.lo],[0,0,-1]
    elif area <= self.areas[1]:
      while 1:
        r1 = random() - 0.5
        r2 = random() - 0.5
        if r1*r1 + r2*r2 < 0.25: break
      if self.axis == 'x':
        return [self.hi,self.c1 + 2*r1*self.r,self.c2 + 2*r2*self.r],[1,0,0]
      elif self.axis == 'y':
        return [self.c1 + 2*r1*self.r,self.hi,self.c2 + 2*r2*self.r],[0,1,0]
      elif self.axis == 'z':
        return [self.c1 + 2*r1*self.r,self.c2 + 2*r2*self.r,self.hi],[0,0,1]
    else:
      r1 = 2*pi*random()
      r2 = random()
      if self.axis == 'x':
        n = [0,cos(r1),sin(r1)]
        return [self.lo + r2*(self.hi-self.lo),
                self.c1 + self.r*cos(r1),
                self.c2 + self.r*sin(r1)],n
      elif self.axis == 'y':
        n = [cos(r1),0,sin(r1)]
        return [self.c1 + self.r*cos(r1),
                self.lo + r2*(self.hi-self.lo),
                self.c2 + self.r*sin(r1)],n
      elif self.axis == 'z':
        n = [cos(r1),sin(r1),0]
        return [self.c1 + self.r*cos(r1),
                self.c2 + self.r*sin(r1),
                self.lo + r2*(self.hi-self.lo)],n

  # ChemCell text to create the region

  def command(self):
    return "%s cylinder %s %g %g %g %g %g" % (self.id,self.axis,
                                              self.c1,self.c2,self.r,
                                              self.lo,self.hi)

# --------------------------------------------------------------------
# capped-cylinder region

class Capped:
  def __init__(self,*list):
    self.axis = list[0]
    self.c1 = list[1]
    self.c2 = list[2]
    self.r = list[3]
    self.lo = list[4]
    self.hi = list[5]
    self.rsq = self.r*self.r
    self.q1 = 2
    self.q2 = 1

  # bounding box around region

  def bbox(self):
    if self.axis == 'x':
      return (self.lo-self.r,self.c1-self.r,self.c2-self.r, \
              self.hi+self.r,self.c1+self.r,self.c2+self.r)
    elif self.axis == 'y':
      return (self.c1-self.r,self.lo-self.r,self.c2-self.r, \
              self.c1+self.r,self.hi+self.r,self.c2+self.r)
    elif self.axis == 'z':
      return (self.c1-self.r,self.c2-self.r,self.lo-self.r, \
              self.c1+self.r,self.c2+self.r,self.hi+self.r)
  
  # return 1,0 for xyz inside/outside the region

  def inside(self,x,y,z):
    if self.axis == 'x':
      d1 = y - self.c1
      d2 = z - self.c2
      d3 = x
    elif sef.axis == 'y':
      d1 = x - self.c1
      d2 = z - self.c2
      d3 = y
    elif self.axis == 'z':
      d1 = x - self.c1
      d2 = y - self.c2
      d3 = z
    rsq = d1*d1 + d2*d2
    if d3 >= self.lo and d3 <= self.hi:
      if rsq > self.rsq: return 0
    elif d3 < self.lo:
      if d1*d1 + d2*d2 + (d3-self.lo)*(d3-self.lo) > self.rsq: return 0
    else:
      if d1*d1 + d2*d2 + (d3-self.hi)*(d3-self.hi) > self.rsq: return 0
    return 1

  # triangulate the region
  # set nvert,ntri,vertices,triangles,connections
  # convert vertices from unit box to cylinder with correct axis
  # x,y have spacings for ends of cylinder
  # z has spacings for axis of cylinder plus 1/2 of q1 on each end cap
  # convert tuples returned by box_triangluate (used for dict lookup) to lists

  def triangulate(self):
    if self.q1 % 2: raise StandardError,"Capped cylinder q1 must be even"
    vertices,triangles = box_triangulate(self.q1,self.q1,self.q2+self.q1)
    self.nvert = len(vertices)
    self.ntri = len(triangles)
    self.vertices = []
    cutlo = self.q1/2 * 1.0/(self.q2+self.q1) + EPSILON
    cuthi = 1.0 - cutlo
    for i in xrange(self.nvert):
      v1 = vertices[i][0]
      v2 = vertices[i][1]
      v3 = vertices[i][2]
      if v3 < cutlo:
        c = [v1-0.5,v2-0.5,(v3-cutlo)*(0.5/cutlo)]
        normalize(c)
        p1 = self.c1 + self.r*c[0]
        p2 = self.c2 + self.r*c[1]
        p3 = self.lo + self.r*c[2]
      elif v3 > cuthi:
        c = [v1-0.5,v2-0.5,(v3-cuthi)*(0.5/cutlo)]
        normalize(c)
        p1 = self.c1 + self.r*c[0]
        p2 = self.c2 + self.r*c[1]
        p3 = self.hi + self.r*c[2]
      else:
        c = [v1-0.5,v2-0.5,0.0]
        normalize(c)
        p1 = self.c1 + self.r*c[0]
        p2 = self.c2 + self.r*c[1]
        p3 = self.lo + (v3-cutlo)/(cuthi-cutlo) * (self.hi - self.lo)
      if self.axis == 'x': self.vertices.append([p3,p1,p2])
      elif self.axis == 'y': self.vertices.append([p1,p3,p2])
      elif self.axis == 'z': self.vertices.append([p1,p2,p3])
    self.triangles = []
    for i in xrange(self.ntri):
      self.triangles.append([triangles[i][0],triangles[i][1],triangles[i][2]])
    self.connections = connect(self.nvert,self.ntri,self.triangles)

  # surface area of region

  def area(self):
    self.areas = []
    sum = 2 * pi * self.rsq
    self.areas.append(sum)
    sum += 2 * pi * self.rsq
    self.areas.append(sum)
    sum += 2 * pi * self.r * (self.hi - self.lo)
    self.areas.append(sum)
    return sum

  # return a random location on surface of region

  def loc2d(self,area,random):
    if area <= self.areas[0]:
      while 1:
        r1 = random() - 0.5
        r2 = random() - 0.5
        r3 = random() - 0.5
        if r3 <= 0 and r1*r1 + r2*r2 + r3*r3 < 0.25: break
      c = [r1,r2,r3]
      normalize(c)
      if self.axis == 'x':
        return [self.lo + c[2]*self.r,
                self.c1 + c[0]*self.r,
                self.c2 + c[1]*self.r],[c[2],c[0],c[1]]
      elif self.axis == 'y':
        return [self.c1 + c[0]*self.r,
                self.lo + c[2]*self.r,
                self.c2 + c[1]*self.r],[c[0],c[2],c[1]]
      elif self.axis == 'z':
        return [self.c1 + c[0]*self.r,
                self.c2 + c[1]*self.r,
                self.lo + c[2]*self.r],[c[0],c[1],c[2]]
    elif area <= self.areas[1]:
      while 1:
        r1 = random() - 0.5
        r2 = random() - 0.5
        r3 = random() - 0.5
        if r3 >= 0 and r1*r1 + r2*r2 + r3*r3 < 0.25: break
      c = [r1,r2,r3]
      normalize(c)
      if self.axis == 'x':
        return [self.hi + c[2]*self.r,
                self.c1 + c[0]*self.r,
                self.c2 + c[1]*self.r],[c[2],c[0],c[1]]
      elif self.axis == 'y':
        return [self.c1 + c[0]*self.r,
                self.hi + c[2]*self.r,
                self.c2 + c[1]*self.r],[c[0],c[2],c[1]]
      elif self.axis == 'z':
        return [self.c1 + c[0]*self.r,
                self.c2 + c[1]*self.r,
                self.hi + c[2]*self.r],[c[0],c[1],c[2]]
    else:
      r1 = 2*pi*random()
      r2 = random()
      if self.axis == 'x':
        n = [0,cos(r1),sin(r1)]
        return [self.lo + r2*(self.hi-self.lo),
                self.c1 + self.r*cos(r1),
                self.c2 + self.r*sin(r1)],n
      elif self.axis == 'y':
        n = [cos(r1),0,sin(r1)]
        return [self.c1 + self.r*cos(r1),
                self.lo + r2*(self.hi-self.lo),
                self.c2 + self.r*sin(r1)],n
      elif self.axis == 'z':
        n = [cos(r1),sin(r1),0]
        return [self.c1 + self.r*cos(r1),
                self.c2 + self.r*sin(r1),
                self.lo + r2*(self.hi-self.lo)],n

  # ChemCell text to create the region

  def command(self):
    return "%s capped %s %g %g %g %g %g" % (self.id,self.axis,
                                            self.c1,self.c2,self.r,
                                            self.lo,self.hi)

# --------------------------------------------------------------------
# set of line segments

class Line:
  def bbox(self):
    list = [pair[0] for pair in self.pairs]
    list += [pair[3] for pair in self.pairs]
    xlo = min(list)
    xhi = max(list)
    list = [pair[1] for pair in self.pairs]
    list += [pair[4] for pair in self.pairs]
    ylo = min(list)
    yhi = max(list)
    list = [pair[2] for pair in self.pairs]
    list += [pair[5] for pair in self.pairs]
    zlo = min(list)
    zhi = max(list)
    return (xlo,ylo,zlo,xhi,yhi,zhi)

  def addline(self,coords):
    self.pairs.append(list(coords))
    
# --------------------------------------------------------------------
# union object that contains other objects
# pass in parent cdata IDs and object list to constructor
# self.objs = list of child objects

class Union:
  def __init__(self,ids,objs,*list):
    self.objs = []
    for id in list:
      obj = objs[ids[id]]
      if obj.style != SURFACE and obj.style != REGION and obj.style != UNION:
        raise StandardError,"Union child object is of invalid style"
      self.objs.append(obj)

  # bounding box around all child objects

  def bbox(self):
    xlo,ylo,zlo,xhi,yhi,zhi = self.objs[0].bbox()
    for obj in self.objs[1:]:
      xxlo,yylo,zzlo,xxhi,yyhi,zzhi = obj.bbox()
      if xxlo < xlo: xlo = xxlo
      if yylo < ylo: ylo = yylo
      if zzlo < zlo: zlo = zzlo
      if xxhi > xhi: xhi = xxhi
      if yyhi > yhi: yhi = yyhi
      if zzhi > zhi: zhi = zzhi
    return (xlo,ylo,zlo,xhi,yhi,zhi)

  # return 1,0 for xyz inside/outside the union
  # inside union if inside any of its child objects
  
  def inside(self,x,y,z):
    for obj in self.objs:
      if obj.inside(x,y,z): return 1
    return 0

  # surface area of union
  # areas = cummulative total area for child objects
  
  def area(self):
    self.areas = []
    sum = 0.0
    for obj in self.objs:
      sum += obj.area()
      self.areas.append(sum)
    return sum

  # return a random location on surface of one of child objects
  
  def loc2d(self,area,random):
    for i in xrange(len(self.objs)):
      if area < self.areas[i]: break
    if i > 0: area -= self.areas[i]
    return self.objs[i].loc2d(area,random)
    
# --------------------------------------------------------------------
# return c = a x b

def cross(a,b):
  c = 3*[0]
  c[0] = a[1]*b[2] - a[2]*b[1]
  c[1] = a[2]*b[0] - a[0]*b[2]
  c[2] = a[0]*b[1] - a[1]*b[0]
  return c

# --------------------------------------------------------------------
# normalize vector a to unit length

def normalize(a):
  length = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
  if length == 0.0: return
  a[0] /= length
  a[1] /= length
  a[2] /= length

# --------------------------------------------------------------------
# compute normal for a triangle with 3 vertices

def normal(x,y,z):
 v1 = 3*[0]
 v1[0] = y[0] - x[0]
 v1[1] = y[1] - x[1]
 v1[2] = y[2] - x[2]

 v2 = 3*[0]
 v2[0] = z[0] - y[0]
 v2[1] = z[1] - y[1]
 v2[2] = z[2] - y[2]

 n = 3*[0]
 n[0] = v1[1]*v2[2] - v1[2]*v2[1]
 n[1] = v1[2]*v2[0] - v1[0]*v2[2]
 n[2] = v1[0]*v2[1] - v1[1]*v2[0]

 length = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
 n[0] /= length
 n[1] /= length
 n[2] /= length

 return n

# --------------------------------------------------------------------
# create connections between triangles for a triangulated surface

def connect(nvert,ntri,triangles):

  # triangles have list of 3 vertices
  # create inverted v2tri where each vertex has list of triangles it is in

  v2tri = []
  for i in xrange(nvert): v2tri.append([]);
    
  for i in xrange(ntri):
    for vert in triangles[i]: v2tri[vert-1].append(i)
    
  # loop over triangles to reset connections

  connections = []
  for i in xrange(ntri):

    # zero connections

    connect = 6*[0]
    connections.append(connect)
      
    # tri123 = list of tris attached to each vertex of triangle i
      
    v = triangles[i]
    tri1 = v2tri[v[0]-1]
    tri2 = v2tri[v[1]-1]
    tri3 = v2tri[v[2]-1]
      
    # loop over attached tris and look for 2nd vertex in each of 3 edges
    # when find it, set connection triangle and edge values
    
    for itri in tri1:
      if itri == i: continue
      if v[1] in triangles[itri]:
        connect[0] = itri+1
        m = triangles[itri].index(v[0])
        n = triangles[itri].index(v[1])
        if (m == 0 and n == 1) or (m == 1 and n == 0): connect[1] = 1
        elif (m == 1 and n == 2) or (m == 2 and n == 1): connect[1] = 2
        else: connect[1] = 3
        break

    for itri in tri2:
      if itri == i: continue
      if v[2] in triangles[itri]:
        connect[2] = itri+1
        m = triangles[itri].index(v[1])
        n = triangles[itri].index(v[2])
        if (m == 0 and n == 1) or (m == 1 and n == 0): connect[3] = 1
        elif (m == 1 and n == 2) or (m == 2 and n == 1): connect[3] = 2
        else: connect[3] = 3
        break

    for itri in tri3:
      if itri == i: continue
      if v[0] in triangles[itri]:
        connect[4] = itri+1
        m = triangles[itri].index(v[2])
        n = triangles[itri].index(v[0])
        if (m == 0 and n == 1) or (m == 1 and n == 0): connect[5] = 1
        elif (m == 1 and n == 2) or (m == 2 and n == 1): connect[5] = 2
        else: connect[5] = 3
        break

  return connections

# --------------------------------------------------------------------
# add a vertex v to vertices list unless already exists in vdict dictionary
# return index of where v is in vertices list

def vertex(v,vertices,vdict):
  if vdict.has_key(v): return vdict[v]
  n = len(vertices)
  vertices.append(v)
  vdict[v] = n
  return n

# --------------------------------------------------------------------
# triangulate a unit box from (0,0,0) to (1,1,1) with spacings q1,q2,q3
# return lists of vertices and triangles
# insure right-hand rule for each tri points OUT of the box

def box_triangulate(q1,q2,q3):
  if q1: dx = 1.0 / q1
  if q2: dy = 1.0 / q2
  if q3: dz = 1.0 / q3
  vdict = {}
  vertices = []
  triangles = []
  for j in xrange(q2):
    for k in xrange(q3):
      v1 = (0, j*dy,     k*dz)
      v2 = (0, (j+1)*dy, k*dz)
      v3 = (0, (j+1)*dy, (k+1)*dz)
      v4 = (0, j*dy,     (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1+1,iv3+1,iv2+1])
      triangles.append([iv1+1,iv4+1,iv3+1])
      v1 = (1, j*dy,     k*dz)
      v2 = (1, (j+1)*dy, k*dz)
      v3 = (1, (j+1)*dy, (k+1)*dz)
      v4 = (1, j*dy,     (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1+1,iv2+1,iv3+1])
      triangles.append([iv1+1,iv3+1,iv4+1])
  for i in xrange(q1):
    for k in xrange(q3):
      v1 = (i*dx,     0, k*dz)
      v2 = ((i+1)*dx, 0, k*dz)
      v3 = ((i+1)*dx, 0, (k+1)*dz)
      v4 = (i*dx,     0, (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1+1,iv2+1,iv3+1])
      triangles.append([iv1+1,iv3+1,iv4+1])
      v1 = (i*dx,     1, k*dz)
      v2 = ((i+1)*dx, 1, k*dz)
      v3 = ((i+1)*dx, 1, (k+1)*dz)
      v4 = (i*dx,     1, (k+1)*dz)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1+1,iv3+1,iv2+1])
      triangles.append([iv1+1,iv4+1,iv3+1])
  for i in xrange(q1):
    for j in xrange(q2):
      v1 = (i*dx,     j*dy,     0)
      v2 = ((i+1)*dx, j*dy,     0)
      v3 = ((i+1)*dx, (j+1)*dy, 0)
      v4 = (i*dx,     (j+1)*dy, 0)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1+1,iv3+1,iv2+1])
      triangles.append([iv1+1,iv4+1,iv3+1])
      v1 = (i*dx,     j*dy,     1)
      v2 = ((i+1)*dx, j*dy,     1)
      v3 = ((i+1)*dx, (j+1)*dy, 1)
      v4 = (i*dx,     (j+1)*dy, 1)
      iv1 = vertex(v1,vertices,vdict)
      iv2 = vertex(v2,vertices,vdict)
      iv3 = vertex(v3,vertices,vdict)
      iv4 = vertex(v4,vertices,vdict)
      triangles.append([iv1+1,iv2+1,iv3+1])
      triangles.append([iv1+1,iv3+1,iv4+1])
  return vertices,triangles
