# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# raster tool

oneline = "3d visualization via Raster3d program"

docstr = """
r = raster(d)               create Raster3d wrapper for data in d

  d = atom snapshot object (dump, data)

r.bg("black")               set background color (def = "black")
r.size(N)		    set image size to NxN
r.size(N,M)		    set image size to NxM
r.rotate(60,135)            view from z theta and azimuthal phi (def = 60,30)
r.shift(x,y)                translate by x,y pixels in view window (def = 0,0)
r.zoom(0.5)                 scale image by factor (def = 1)
r.box(0/1/2)                0/1/2 = none/variable/fixed box
r.box(0/1/2,"green")        set box color
r.box(0/1/2,"red",4)        set box edge thickness
r.file = "image"            file prefix for created images (def = "image")

r.show(N)                   show image of snapshot at timestep N
  
r.all()                     make images of all selected snapshots
r.all(P)                    images of all, start file label at P
r.all(N,M,P)                make M images of snapshot N, start label at P

r.pan(60,135,1.0,40,135,1.5)    pan during all() operation
r.pan()                         no pan during all() (default)

  args = z theta, azimuthal phi, zoom factor at beginning and end
  values at each step are interpolated between beginning and end values

r.select = "$x > %g*3.0"    string to pass to d.aselect.test() during all()
r.select = ""               no extra aselect (default)
				
  %g varies from 0.0 to 1.0 from beginning to end of all()

r.label(x,y,"h",size,"red","This is a label")    add label to each image
r.nolabel()                                      delete all labels
  
  x,y coords = -0.5 to 0.5, "h" or "t" for Helvetica or Times font
  size = fontsize (e.g. 10), "red" = color of text
  
r.acol(2,"green")		   set atom colors by atom type (1-N)
r.acol([2,4],["red","blue"])	   1st arg = one type or list of types
r.acol(0,"blue")	           2nd arg = one color or list of colors
r.acol(range(20),["red","blue"])   if list lengths unequal, interpolate
r.acol(range(10),"loop")           assign colors in loop, randomly ordered

  if 1st arg is 0, set all types to 2nd arg
  if list of types has a 0 (e.g. range(10)), +1 is added to each value
  interpolate means colors blend smoothly from one value to the next

r.arad([1,2],[0.5,0.3])            set atom radii, same rules as acol()

r.bcol()			   set bond color, same args as acol()
r.brad()			   set bond thickness, same args as arad()

r.tcol()			   set triangle color, same args as acol()
r.tfill()			   set triangle fill, 0 fill, 1 line, 2 both

r.lcol()                           set line color, same args as acol()
r.lrad()                           set line thickness, same args as arad()

r.adef()                           set atom/bond/tri/line properties to default
r.bdef()			   default = "loop" for colors, 0.45 for radii
r.tdef()  			   default = 0.25 for bond/line thickness
r.ldef()  			   default = 0 fill

  by default 100 types are assigned
  if atom/bond/tri/line has type > # defined properties, is an error
  
from vizinfo import colors         access color list
print colors                       list defined color names and RGB values
colors["nickname"] = [R,G,B]       set new RGB values from 0 to 255

  140 pre-defined colors: red, green, blue, purple, yellow, black, white, etc
"""

# History
#   8/05, Steve Plimpton (SNL): original version
#   9/05, Steve Plimpton (SNL): adjusted box and label attributes

# ToDo list
#   when do aselect with select str while looping N times on same timestep
#     would not let you grow # of atoms selected
#   triangles are not drawn with fill type

# Variables
#   ztheta = vertical angle from z-azis to view from
#   azphi = azimuthal angle to view from
#   xshift,yshift = xy translation of scene (in pixels)
#   distance = size of simulation box (largest dim)
#   eye = viewpoint distance from center of scene
#   file = filename prefix to use for images produced
#   boxflag = 0/1/2 for drawing simulation box: none/variable/fixed
#   bxcol = color of box
#   bxthick = thickness of box
#   bgcol = color of background
#   vizinfo = scene attributes
#   xtrans,ytrans,ztrans = translation factors computed by Raster3d

# Imports and external programs

import sys, os, commands, re
from vizinfo import vizinfo
from math import fabs,atan,cos,sin

try: from DEFAULTS import PIZZA_RENDER
except: PIZZA_RENDER = "render"
try: from DEFAULTS import PIZZA_LABEL3D
except: PIZZA_LABEL3D = "label3d"
try: from DEFAULTS import PIZZA_DISPLAY
except: PIZZA_DISPLAY = "display"

# Class definition

class raster:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data
    self.xpixels = 512
    self.ypixels = 512
    self.ztheta = 60
    self.azphi = 30
    self.scale = 1.0
    self.xshift = self.yshift = 0
    self.eye = 50.0
    
    self.file = "image"
    self.boxflag = 0
    self.bxcol = [1,1,0]
    self.bxthick = 0.3
    self.bgcol = [0,0,0]
    self.labels = []
    self.panflag = 0
    self.select = ""

    self.vizinfo = vizinfo()
    self.adef()
    self.bdef()
    self.tdef()
    self.ldef()

  # --------------------------------------------------------------------

  def bg(self,color):
    from vizinfo import colors
    self.bgcol = [colors[color][0]/255.0,colors[color][1]/255.0,
                  colors[color][2]/255.0]
  
  # --------------------------------------------------------------------

  def size(self,xnew,ynew=None):
    self.xpixels = xnew
    if not ynew: self.ypixels = self.xpixels
    else: self.ypixels = ynew
        
  # --------------------------------------------------------------------

  def rotate(self,ztheta,azphi):
    self.ztheta = ztheta;
    self.azphi = azphi;

  # --------------------------------------------------------------------

  def shift(self,x,y):
    self.xshift = x;
    self.yshift = y;

  # --------------------------------------------------------------------

  def zoom(self,scale):
    self.scale = scale;
  
  # --------------------------------------------------------------------

  def box(self,*args):
    self.boxflag = args[0]
    if len(args) > 1:
      from vizinfo import colors
      self.bxcol = [colors[args[1]][0]/255.0,colors[args[1]][1]/255.0,
                    colors[args[1]][2]/255.0]
    if len(args) > 2: self.bxthick = args[2]

  # --------------------------------------------------------------------
  # scale down point-size by 3x
  
  def label(self,x,y,font,point,color,text):
    from vizinfo import colors
    scaledcolor = [colors[color][0]/255.0,colors[color][1]/255.0,
                   colors[color][2]/255.0]
    list = [x,y,fontlist[font],point/3.0,"Left",scaledcolor,text]
    self.labels.append(list)

  # --------------------------------------------------------------------

  def nolabel(self):
    self.labels = []
  
  # --------------------------------------------------------------------
  # show a single snapshot
  # distance from snapshot box or max box for all selected steps
  # always pre-call single() to re-center simulation data
  
  def show(self,ntime):
    data = self.data
    which = data.findtime(ntime)
    time,box,atoms,bonds,tris,lines = data.viz(which)
    if self.boxflag == 2: box = data.maxbox()
    self.distance = compute_distance(box)

    self.xtrans = self.ytrans = self.ztrans = 0.0
    output = self.single(1,self.file,box,atoms,bonds,tris,lines)
    print output
    nums = re.findall("translation to:\s*(\S*)\s*(\S*)\s*(\S*)\s",output)
    self.xtrans = float(nums[0][0])
    self.ytrans = float(nums[0][1])
    self.ztrans = float(nums[0][2])
    
    self.single(0,self.file,box,atoms,bonds,tris,lines)
    cmd = "%s %s.png" % (PIZZA_DISPLAY,self.file)
    commands.getoutput(cmd)

  # --------------------------------------------------------------------

  def pan(self,*list):
    if len(list) == 0: self.panflag = 0
    else:
      self.panflag = 1
      self.ztheta_start = list[0]
      self.azphi_start = list[1]
      self.scale_start = list[2]
      self.ztheta_stop = list[3]
      self.azphi_stop = list[4]
      self.scale_stop = list[5]
      
  # --------------------------------------------------------------------

  def all(self,*list):
    data = self.data
    if len(list) == 0:
      nstart = 0
      ncount = data.nselect
    elif len(list) == 1:
      nstart = list[0]
      ncount = data.nselect
    else:
      ntime = list[0]
      nstart = list[2]
      ncount = list[1]

    if self.boxflag == 2: box = data.maxbox()

    # loop over all selected steps
    # distance from 1st snapshot box or max box for all selected steps
    # pre-call single() to re-center simulation data on 1st step or if panning

    if len(list) <= 1:

      n = nstart
      i = flag = 0
      while 1:
        which,time,flag = data.iterator(flag)
        if flag == -1: break

        fraction = float(i) / (ncount-1)
        
        if self.select != "":
          newstr = self.select % fraction
          data.aselect.test(newstr,time)
        time,boxone,atoms,bonds,tris,lines = data.viz(which)

        if self.boxflag < 2: box = boxone
        if n == nstart: self.distance = compute_distance(box)

        if n < 10:     file = self.file + "000" + str(n)
        elif n < 100:  file = self.file + "00" + str(n)
        elif n < 1000: file = self.file + "0" + str(n)
        else:          file = self.file + str(n)

        if self.panflag:
          self.ztheta = self.ztheta_start + \
                        fraction*(self.ztheta_stop - self.ztheta_start)
          self.azphi = self.azphi_start + \
                       fraction*(self.azphi_stop - self.azphi_start)
          self.scale = self.scale_start + \
                          fraction*(self.scale_stop - self.scale_start)

	if n == nstart or self.panflag:
          self.xtrans = self.ytrans = self.ztrans = 0.0
	  output = self.single(1,file,box,atoms,bonds,tris,lines)
          nums = re.findall("translation to:\s*(\S*)\s*(\S*)\s*(\S*)\s",output)
          self.xtrans = float(nums[0][0])
          self.ytrans = float(nums[0][1])
          self.ztrans = float(nums[0][2])

        self.single(0,file,box,atoms,bonds,tris,lines)
        print time,
        sys.stdout.flush()
        i += 1
        n += 1

    # loop ncount times on same step
    # distance from 1st snapshot box or max box for all selected steps
    # pre-call single() to re-center simulation data on 1st step or if panning

    else:

      which = data.findtime(ntime)

      n = nstart
      for i in range(ncount):
        fraction = float(i) / (ncount-1)

        if self.select != "":
          newstr = self.select % fraction
          data.aselect.test(newstr,ntime)
        time,boxone,atoms,bonds,tris,lines = data.viz(which)

        if self.boxflag < 2: box = boxone
        if n == nstart: self.distance = compute_distance(box)

        if n < 10:     file = self.file + "000" + str(n)
        elif n < 100:  file = self.file + "00" + str(n)
        elif n < 1000: file = self.file + "0" + str(n)
        else:          file = self.file + str(n)

        if self.panflag:
          self.ztheta = self.ztheta_start + \
                        fraction*(self.ztheta_stop - self.ztheta_start)
          self.azphi = self.azphi_start + \
                       fraction*(self.azphi_stop - self.azphi_start)
          self.scale = self.scale_start + \
                          fraction*(self.scale_stop - self.scale_start)

	if n == nstart or self.panflag:
          self.xtrans = self.ytrans = self.ztrans = 0.0
	  output = self.single(1,file,box,atoms,bonds,tris,lines)
          nums = re.findall("translation to:\s*(\S*)\s*(\S*)\s*(\S*)\s",output)
          self.xtrans = float(nums[0][0])
          self.ytrans = float(nums[0][1])
          self.ztrans = float(nums[0][2])
        
        self.single(0,file,box,atoms,bonds,tris,lines)
        print n,
        sys.stdout.flush()
        n += 1

    print "\n%d images" % ncount

  # --------------------------------------------------------------------

  def single(self,flag,file,box,atoms,bonds,tris,lines):

    matrix = rotation_matrix('x',-self.ztheta,'z',270.0-self.azphi)

    f = open("tmp.r3d","w")

    color = self.bgcol
    xshift = 1.6*self.distance/self.scale * self.xshift/self.xpixels
    yshift = 1.6*self.distance/self.scale * self.yshift/self.xpixels
    header = template[1:] % \
             (self.xpixels,self.ypixels,color[0],color[1],color[2],
              self.eye,matrix,self.xtrans+xshift,self.ytrans+yshift,
              self.ztrans,1.6*self.distance/self.scale)
    print >>f,header,

    # draw box if boxflag or flag is set
    # flag = 1 is a pre-call of single to set frame size correctly
    #   this keeps the view fixed even if atoms move around
    
    if self.boxflag or flag: box_write(f,box,self.bxcol,self.bxthick)

    ncolor = self.vizinfo.nacolor
    for atom in atoms:
      itype = int(atom[1])
      if itype > ncolor: raise StandardError,"atom type too big"
      color = self.vizinfo.acolor[itype]
      rad = self.vizinfo.arad[itype]
      print >>f,2
      print >>f,atom[2],atom[3],atom[4],rad,color[0],color[1],color[2]

    # need to include vizinfo.tfill options

    ncolor = self.vizinfo.ntcolor
    for tri in tris:
      itype = int(tri[1])
      if itype > ncolor: raise StandardError,"tri type too big"
      color = self.vizinfo.tcolor[itype]
      print >>f,1
      print >>f,tri[2],tri[3],tri[4],tri[5],tri[6],tri[7], \
            tri[8],tri[9],tri[10],color[0],color[1],color[2]

    bound = 0.25 * self.distance
    ncolor = self.vizinfo.nbcolor
    for bond in bonds:
      itype = int(bond[1])
      if itype > ncolor: raise StandardError,"bond type too big"
      color = self.vizinfo.bcolor[itype]
      rad = self.vizinfo.brad[itype]
      if fabs(bond[2]-bond[5]) > bound or fabs(bond[3]-bond[6]) > bound:
        continue
      print >>f,5
      print >>f,bond[2],bond[3],bond[4],rad, \
            bond[5],bond[6],bond[7],0.0,color[0],color[1],color[2]

    ncolor = self.vizinfo.nlcolor
    for line in lines:
      itype = int(line[1])
      if itype > ncolor: raise StandardError,"line type too big"
      color = self.vizinfo.lcolor[itype]
      thick = self.vizinfo.lrad[itype]
      print >>f,3  
      print >>f,line[2],line[3],line[4],thick, \
                line[5],line[6],line[7],thick,color[0],color[1],color[2]

    for label in self.labels:
      print >>f,15
      print >>f,10
      print >>f,label[2],label[3],label[4]
      print >>f,11
      print >>f,label[0],label[1],0.0,label[5][0],label[5][1],label[5][2]
      print >>f,label[6]
      
    f.close()

    if len(self.labels) == 0:
      cmd = "cat tmp.r3d | %s -png %s.png" % (PIZZA_RENDER,file)
    else:
      cmd = "cat tmp.r3d | %s -png %s.png" % (PIZZA_LABEL3D,file)

    output = commands.getoutput(cmd)
    return output

  # --------------------------------------------------------------------

  def adef(self):
    self.vizinfo.setcolors("atom",range(100),"loop")
    self.vizinfo.setradii("atom",range(100),0.45)
 
  # --------------------------------------------------------------------

  def bdef(self):
    self.vizinfo.setcolors("bond",range(100),"loop")
    self.vizinfo.setradii("bond",range(100),0.25)
    
  # --------------------------------------------------------------------

  def tdef(self):
    self.vizinfo.setcolors("tri",range(100),"loop")  
    self.vizinfo.setfills("tri",range(100),0)  

  # --------------------------------------------------------------------

  def ldef(self):
    self.vizinfo.setcolors("line",range(100),"loop")  
    self.vizinfo.setradii("line",range(100),0.25)
   
  # --------------------------------------------------------------------

  def acol(self,atypes,colors):
    self.vizinfo.setcolors("atom",atypes,colors)
    
  # --------------------------------------------------------------------

  def arad(self,atypes,radii):
    self.vizinfo.setradii("atom",atypes,radii)  
    
  # --------------------------------------------------------------------

  def bcol(self,btypes,colors):
    self.vizinfo.setcolors("bond",btypes,colors)
  
  # --------------------------------------------------------------------

  def brad(self,btypes,radii):
    self.vizinfo.setradii("bond",btypes,radii)
  
  # --------------------------------------------------------------------

  def tcol(self,ttypes,colors):
    self.vizinfo.setcolors("tri",ttypes,colors)

  # --------------------------------------------------------------------

  def tfill(self,ttypes,flags):
    self.vizinfo.setfills("tri",ttypes,flags)

  # --------------------------------------------------------------------

  def lcol(self,ltypes,colors):
    self.vizinfo.setcolors("line",ltypes,colors)

  # --------------------------------------------------------------------

  def lrad(self,ltypes,radii):
    self.vizinfo.setradii("line",ltypes,radii)

# --------------------------------------------------------------------
# return characteristic distance of simulation domain = max dimension

def compute_distance(box):
  distance = box[3]-box[0]
  if box[4]-box[1] > distance: distance = box[4]-box[1]
  if box[5]-box[2] > distance: distance = box[5]-box[2]
  return distance

# --------------------------------------------------------------------
# draw a 12-edge box around simulation domain

def box_write(f,box,color,thick):
  xlo,ylo,zlo = box[0],box[1],box[2]
  xhi,yhi,zhi = box[3],box[4],box[5]
  red = color[0]
  green = color[1]
  blue = color[2]
  
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
  	      (xlo,ylo,zlo,thick,xhi,ylo,zlo,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xlo,yhi,zlo,thick,xhi,yhi,zlo,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xlo,ylo,zhi,thick,xhi,ylo,zhi,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xlo,yhi,zhi,thick,xhi,yhi,zhi,thick,red,green,blue)

  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xlo,ylo,zlo,thick,xlo,yhi,zlo,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xhi,ylo,zlo,thick,xhi,yhi,zlo,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xlo,ylo,zhi,thick,xlo,yhi,zhi,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xhi,ylo,zhi,thick,xhi,yhi,zhi,thick,red,green,blue)

  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xlo,ylo,zlo,thick,xlo,ylo,zhi,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xhi,ylo,zlo,thick,xhi,ylo,zhi,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xlo,yhi,zlo,thick,xlo,yhi,zhi,thick,red,green,blue)
  print >>f,"3\n%g %g %g %g %g %g %g %g %g %g %g" % \
	      (xhi,yhi,zlo,thick,xhi,yhi,zhi,thick,red,green,blue)
    
# --------------------------------------------------------------------
# compute 3x3 rotation matrix for viewing angle

# initially the scene is viewed:
#   (1) along the z-axis, looking towards the origin from (0,0,1)
#   (2) seeing xy plane, with +y up and +x to the right
# 1st rotation angle rotates the body
# 2nd rotation angle rotates the new body around the new rotated axis
# sign of rotation angles follow right-hand rule

# rotation_matrix(x/y/z,angle1,x/y/z,angle2)
#   x/y/z = 1st axis to rotate around
#   angle1 = 1st angle to rotate by
#   x/y/z = 2nd axis to rotate around
#   angle2 = 2nd angle to rotate by

def rotation_matrix(coord1,angle1,coord2,angle2):

  # convert rotation angles to radians

  pi = 4.0*atan(1.0)
  angle1 = angle1/180.0 * pi
  angle2 = angle2/180.0 * pi

  # sines/cosines of 2 angles

  c1 = cos(angle1)
  s1 = sin(angle1)
  c2 = cos(angle2)
  s2 = sin(angle2)

  # 1st rotation matrix

  a11 = a12 = a13 = a21 = a22 = a23 = a31 = a32 = a33 = 0.0
  if coord1 == 'x':
    a11 = 1.0
    a22 = a33 = c1
    a23 = s1; a32 = -s1
  elif coord1 == 'y':
    a22 = 1.0
    a11 = a33 = c1
    a13 = s1; a31 = -s1
  elif coord1 == 'z':
    a33 = 1.0
    a12 = a22 = c1
    a12 = s1; a21 = -s1
  
  # 2nd rotation matrix
  
  b11 = b12 = b13 = b21 = b22 = b23 = b31 = b32 = b33 = 0.0
  if coord2 == 'x':
    b11 = 1.0
    b22 = b33 = c2
    b23 = s2; b32 = -s2
  elif coord2 == 'y':
    b22 = 1.0
    b11 = b33 = c2
    b13 = s2; b31 = -s2
  elif coord2 == 'z':
    b33 = 1.0
    b11 = b22 = c2
    b12 = s2; b21 = -s2
  
  # full matrix c = b*a
  
  c11 = b11*a11 + b12*a21 + b13*a31
  c12 = b11*a12 + b12*a22 + b13*a32
  c13 = b11*a13 + b12*a23 + b13*a33
  
  c21 = b21*a11 + b22*a21 + b23*a31
  c22 = b21*a12 + b22*a22 + b23*a32
  c23 = b21*a13 + b22*a23 + b23*a33
  
  c31 = b31*a11 + b32*a21 + b33*a31
  c32 = b31*a12 + b32*a22 + b33*a32
  c33 = b31*a13 + b32*a23 + b33*a33
  
  # form rotation matrix
  # each line padded with 0.0 for 4x4 raster3d matrix

  matrix = "%g %g %g 0.0\n%g %g %g 0.0\n%g %g %g 0.0" % \
  	   (c11,c12,c13,c21,c22,c23,c31,c32,c33)
  return matrix

# --------------------------------------------------------------------
# template file for Raster3d

template = """
Raster3D header file
%g %g     tiles in x,y
0 0       pixels (x,y) per tile
4         anti-aliasing
%g %g %g  background color
F         no, shadowed rods look funny
25        Phong power
0.15      secondary light contribution
0.05      ambient light contribution
0.25      specular reflection component
%g        eye position
1 1 1     main light source postion
%s	  
%g %g %g %g
3         mixed object types
*
*
*
"""

# --------------------------------------------------------------------
# fontlist

fontlist = {}
fontlist["t"] = "Times-Roman"
fontlist["h"] = "Helvetica"
