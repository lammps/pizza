# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# svg tool

oneline = "3d visualization via SVG files"

docstr = """
s = svg(d)                  create SVG object for data in d

  d = atom snapshot object (dump, data)

s.bg("black")               set background color (def = "black")
s.size(N)		    set image size to NxN
s.size(N,M)		    set image size to NxM
s.rotate(60,135)            view from z theta and azimuthal phi (def = 60,30)
s.shift(x,y)                translate by x,y pixels in view window (def = 0,0)
s.zoom(0.5)                 scale image by factor (def = 1)
s.box(0/1/2)                0/1/2 = none/variable/fixed box
s.box(0/1/2,"green")        set box color
s.box(0/1/2,"red",4)        set box edge thickness
s.file = "image"            file prefix for created images (def = "image")

s.show(N)                   show image of snapshot at timestep N

s.all()                     make images of all selected snapshots
s.all(P)                    images of all, start file label at P
s.all(N,M,P)                make M images of snapshot N, start label at P

s.pan(60,135,1.0,40,135,1.5)    pan during all() operation
s.pan()                         no pan during all() (default)

  args = z theta, azimuthal phi, zoom factor at beginning and end
  values at each step are interpolated between beginning and end values

s.select = "$x > %g*3.0"    string to pass to d.aselect.test() during all()
s.select = ""               no extra aselect (default)
				
  %g varies from 0.0 to 1.0 from beginning to end of all()

s.label(x,y,"h",size,"red","This is a label")    add label to each image
s.nolabel()                                      delete all labels
  
  x,y coords = -0.5 to 0.5, "h" or "t" for Helvetica or Times font
  size = fontsize (e.g. 10), "red" = color of text
  
s.acol(2,"green")		   set atom colors by atom type (1-N)
s.acol([2,4],["red","blue"])	   1st arg = one type or list of types
s.acol(0,"blue")          	   2nd arg = one color or list of colors
s.acol(range(20),["red","blue"])   if list lengths unequal, interpolate
s.acol(range(10),"loop")           assign colors in loop, randomly ordered

  if 1st arg is 0, set all types to 2nd arg
  if list of types has a 0 (e.g. range(10)), +1 is added to each value
  interpolate means colors blend smoothly from one value to the next

s.arad([1,2],[0.5,0.3])            set atom radii, same rules as acol()

s.bcol()			   set bond color, same args as acol()
s.brad()			   set bond thickness, same args as arad()

s.tcol()			   set triangle color, same args as acol()
s.tfill()			   set triangle fill, 0 fill, 1 line, 2 both

s.lcol()                           set line color, same args as acol()
s.lrad()                           set line thickness, same args as arad()

s.adef()                           set atom/bond/tri/line properties to default
s.bdef()			   default = "loop" for colors, 0.45 for radii
s.tdef()  			   default = 0.25 for bond/line thickness
s.ldef()  			   default = 0 fill

  by default 100 types are assigned
  if atom/bond/tri/line has type > # defined properties, is an error

from vizinfo import colors         access color list
print colors                       list defined color names and RGB values
colors["nickname"] = [R,G,B]       set new RGB values from 0 to 255

  140 pre-defined colors: red, green, blue, purple, yellow, black, white, etc

Settings specific to svg tool:

s.thick = 2.0               pixel thickness of black atom border
"""

# History
#   8/05, Matt Jones (BYU): original version
#   9/05, Steve Plimpton: adjusted box and label attributes

# ToDo list
#   when do aselect with select str while looping N times on same timestep
#     would not let you grow # of atoms selected
#   triangles are not drawn with fill type

# Variables
#   ztheta = vertical angle from z-azis to view from
#   azphi = azimuthal angle to view from
#   xshift,yshift = xy translation of scene (in pixels)
#   distance = size of simulation box (largest dim)
#   file = filename prefix to use for images produced
#   boxflag = 0/1/2 for drawing simulation box: none/variable/fixed
#   bxcol = color of box
#   bxthick = thickness of box
#   bgcol = color of background
#   vizinfo = scene attributes

# Imports and external programs

import sys, os, commands, re
from vizinfo import vizinfo
from math import sqrt,atan,cos,sin,fabs

try: from DEFAULTS import PIZZA_DISPLAY
except: PIZZA_DISPLAY = "display"

# Class definition

class svg:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data
    self.xpixels = 512
    self.ypixels = 512
    self.ztheta = 60
    self.azphi = 30
    self.scale = 1.0
    self.xshift = self.yshift = 0

    self.file = "image"
    self.boxflag = 0
    self.bxcol = [1,1,0]
    self.bxthick = 0.3
    self.bgcol = [0,0,0]
    self.labels = []
    self.panflag = 0
    self.select = ""
    self.thick = 1.0

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

  def size(self,newx,newy=None):
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

  def box(self,*args):
    self.boxflag = args[0]
    if len(args) > 1:
      from vizinfo import colors
      self.bxcol = [colors[args[1]][0]/255.0,colors[args[1]][1]/255.0,
                    colors[args[1]][2]/255.0]
    if len(args) > 2: self.bxthick = args[2]
  
  # --------------------------------------------------------------------

  def zoom(self,factor):
    self.scale = factor
  
  # --------------------------------------------------------------------

  def show(self,ntime):
    data = self.data
    which = data.findtime(ntime)
    time,box,atoms,bonds,tris,lines = data.viz(which)
    if self.boxflag == 2: box = data.maxbox()
    self.distance = compute_distance(box)
    
    self.single(self.file,box,atoms,bonds,tris,lines,1)
    cmd = "%s %s.svg" % (PIZZA_DISPLAY,self.file)
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
    # call single() w/ scaling on 1st step or if panning

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
          
        scaleflag = 0
	if n == nstart or self.panflag: scaleflag = 1       

	self.single(file,box,atoms,bonds,tris,lines,scaleflag) 
        print time,
        sys.stdout.flush()
        i += 1
        n += 1

    # loop ncount times on same step
    # distance from 1st snapshot box or max box for all selected steps
    # call single() w/ scaling on 1st step or if panning

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

        scaleflag = 0
	if n == nstart or self.panflag: scaleflag = 1

	self.single(file,box,atoms,bonds,tris,lines,scaleflag) 
        print n,
        sys.stdout.flush()
        n += 1
        
    print "\n%d images" % ncount  
  
  # --------------------------------------------------------------------

  def label(self,x,y,font,point,color,text):
    from vizinfo import colors
    scaledcolor = [colors[color][0]/255.0,colors[color][1]/255.0,
                   colors[color][2]/255.0]
    list = [x,y,fontlist[font],point,scaledcolor,text]
    self.labels.append(list)

  # --------------------------------------------------------------------

  def nolabel(self):
    self.labels = []
    
  # --------------------------------------------------------------------

  def single(self,file,box,atoms,bonds,tris,lines,scaleflag):
    
    matrix = rotation_matrix('x',-self.ztheta,'z',270.0-self.azphi)
    if scaleflag:
      self.factor = self.xpixels*self.scale / (1.6*self.distance)
      xctr = 0.5 * (box[0]+box[3])
      yctr = 0.5 * (box[1]+box[4])
      zctr = 0.5 * (box[2]+box[5])
      self.offsetx = matrix[0]*xctr + matrix[3]*yctr + matrix[6]*zctr
      self.offsety = matrix[1]*xctr + matrix[4]*yctr + matrix[7]*zctr

    olist = []

    for atom in atoms:
      atom[0] = 0
      newatom = self.transform(atom,matrix)
      olist.append(newatom)

    for tri in tris:
      tri[0] = 1
      newtri = self.transform(tri,matrix)
      olist.append(newtri)
    
    bound = 0.25 * self.distance
    for bond in bonds:
      newbond = [2,bond[1]]
      dx = bond[5] - bond[2]
      dy = bond[6] - bond[3]
      dz = bond[7] - bond[4] 
      r = sqrt(dx*dx+dy*dy+dz*dz)
      if not r: r = 1
      rad = self.vizinfo.arad[int(bond[9])]
      newbond.append(bond[2] + (r/r - rad/r) * dx)
      newbond.append(bond[3] + (r/r - rad/r) * dy)
      newbond.append(bond[4] + (r/r - rad/r) * dz)
      
      # cut off second side of bond
      
      dx = bond[2] - bond[5]
      dy = bond[3] - bond[6]
      dz = bond[4] - bond[7] 
      r  = sqrt(dx*dx+dy*dy+dz*dz)
      if not r: r = 1
      rad = self.vizinfo.arad[int(bond[8])]
      newbond.append(bond[5] + (r/r - rad/r) * dx)
      newbond.append(bond[6] + (r/r - rad/r) * dy)
      newbond.append(bond[7] + (r/r - rad/r) * dz)     

      if fabs(newbond[2]-newbond[5]) > bound or \
             fabs(newbond[3]-newbond[6]) > bound: continue

      newbond = self.transform(newbond,matrix)
      if newbond[4] < newbond[7]: newbond[4] = newbond[7]
      olist.append(newbond)
      
    for line in lines:
      line[0] = 3
      newline = self.transform(line,matrix)
      olist.append(newline)

    if self.boxflag:
      x1,y1,z1 = box[0],box[1],box[2]
      x2,y2,z2 = box[3],box[4],box[5]
      blines = []
      blines.append([4,0,x1,y1,z1,x1,y1,z2])
      blines.append([4,0,x2,y1,z1,x2,y1,z2])
      blines.append([4,0,x2,y2,z1,x2,y2,z2])
      blines.append([4,0,x1,y2,z1,x1,y2,z2])
      blines.append([4,0,x1,y1,z1,x2,y1,z1])
      blines.append([4,0,x1,y2,z1,x2,y2,z1])
      blines.append([4,0,x1,y2,z2,x2,y2,z2])
      blines.append([4,0,x1,y1,z2,x2,y1,z2])
      blines.append([4,0,x1,y1,z1,x1,y2,z1])
      blines.append([4,0,x2,y1,z1,x2,y2,z1])
      blines.append([4,0,x2,y1,z2,x2,y2,z2])
      blines.append([4,0,x1,y1,z2,x1,y2,z2])
      for line in blines:
        newline = self.transform(line,matrix)
        olist.append(newline)

    # convert objects by factor/offset and sort by z-depth

    self.convert(olist)
    olist.sort(cmprz)

    # write SVG file

    file += ".svg"
    f = open(file,"w") 
    
    header = '<?xml version="1.0"?> <svg height="%s" width="%s" >' % \
             (self.ypixels,self.xpixels)
    header += '<g style="fill-opacity:1.0; stroke:black; stroke-width:0.001;">'
    print >>f,header

    color = '<rect x="0" y="0" height="%s" width="%s" ' % \
            (self.ypixels,self.xpixels)
    color += 'fill="rgb(%s,%s,%s)"/>' % \
             (self.bgcol[0]*255,self.bgcol[1]*255,self.bgcol[2]*255)
    print >>f,color
    
    for element in olist: self.write(f,0,element)
    for label in self.labels: self.write(f,1,label)

    footer = "</g></svg>"
    print >> f,footer 
    
    f.close()  
  
  # --------------------------------------------------------------------
  # rotate with matrix

  def transform(self,obj,matrix):

    onew = obj[0:2]

    if obj[0] == 0:                      # transform atom
      onew.append(matrix[0]*obj[2] + matrix[3]*obj[3] + matrix[6]*obj[4])
      onew.append(matrix[1]*obj[2] + matrix[4]*obj[3] + matrix[7]*obj[4])
      onew.append(matrix[2]*obj[2] + matrix[5]*obj[3] + matrix[8]*obj[4])

    elif obj[0] == 1:                    # transform triangle
      onew.append(matrix[0]*obj[2] + matrix[3]*obj[3] + matrix[6]*obj[4])
      onew.append(matrix[1]*obj[2] + matrix[4]*obj[3] + matrix[7]*obj[4])
      onew.append(matrix[2]*obj[2] + matrix[5]*obj[3] + matrix[8]*obj[4])
      onew.append(matrix[0]*obj[5] + matrix[3]*obj[6] + matrix[6]*obj[7])
      onew.append(matrix[1]*obj[5] + matrix[4]*obj[6] + matrix[7]*obj[7])
      onew.append(matrix[2]*obj[5] + matrix[5]*obj[6] + matrix[8]*obj[7])
      onew.append(matrix[0]*obj[8] + matrix[3]*obj[9] + matrix[6]*obj[10])
      onew.append(matrix[1]*obj[8] + matrix[4]*obj[9] + matrix[7]*obj[10])
      onew.append(matrix[2]*obj[8] + matrix[5]*obj[9] + matrix[8]*obj[10])

    else:                                # transform bond or line
      onew.append(matrix[0]*obj[2] + matrix[3]*obj[3] + matrix[6]*obj[4])
      onew.append(matrix[1]*obj[2] + matrix[4]*obj[3] + matrix[7]*obj[4])
      onew.append(matrix[2]*obj[2] + matrix[5]*obj[3] + matrix[8]*obj[4])
      onew.append(matrix[0]*obj[5] + matrix[3]*obj[6] + matrix[6]*obj[7])
      onew.append(matrix[1]*obj[5] + matrix[4]*obj[6] + matrix[7]*obj[7])
      onew.append(matrix[2]*obj[5] + matrix[5]*obj[6] + matrix[8]*obj[7])
    
    return onew
    
  # --------------------------------------------------------------------

  def convert(self,objlist):
    factor = self.factor
    offsetx = self.offsetx
    offsety = self.offsety
    xctr = 0.5 * self.xpixels + self.xshift
    yctr = 0.5 * self.ypixels - self.yshift

    for obj in objlist:
      if obj[0] == 0:                                # convert atom
        obj[2] = factor*(obj[2] - offsetx) + xctr
        obj[3] = yctr - factor*(obj[3] - offsety)
      elif obj[0] == 1:                              # convert triangle
        obj[2] = factor*(obj[2] - offsetx) + xctr
        obj[3] = yctr - factor*(obj[3] - offsety)
        obj[5] = factor*(obj[5] - offsetx) + xctr
        obj[6] = yctr - factor*(obj[6] - offsety)
        obj[8] = factor*(obj[8] - offsetx) + xctr
        obj[9] = yctr - factor*(obj[9] - offsety)
      else:                                          # convert bond or line
        obj[2] = factor*(obj[2] - offsetx) + xctr
        obj[3] = yctr - factor*(obj[3] - offsety)
        obj[5] = factor*(obj[5] - offsetx) + xctr
        obj[6] = yctr - factor*(obj[6] - offsety)
 
  # --------------------------------------------------------------------

  def write(self,f,flag,*args):
    if len(args): obj = args[0]
    
    if flag == 0:
      if obj[0] == 0:     # atom with its color and radius
	itype = int(obj[1])
	if itype > self.vizinfo.nacolor:
          raise StandardError,"atom type too big"
        color = self.vizinfo.acolor[itype]
        rad = self.vizinfo.arad[itype]
	print >>f,'<circle cx="%s" cy="%s" r="%s" fill="rgb(%s,%s,%s)" stroke-width="%s" />' % \
                   (obj[2],obj[3],rad*self.factor,
                    color[0]*255,color[1]*255,color[2]*255,self.thick)
        
      elif obj[0] == 1:    # tri with its color (need to add fill type)
        itype = int(obj[1]) 
	if itype > self.vizinfo.ntcolor:
          raise StandardError,"tri type too big"
        color = self.vizinfo.tcolor[itype]
        print >>f,'<polygon points= "%s,%s %s,%s %s,%s" fill="rgb(%s,%s,%s)" stroke="black" stroke-width="0.01" />' % \
              (obj[2],obj[3],obj[5],obj[6],obj[8],obj[9],
               color[0]*255,color[1]*255,color[2]*255)

      elif obj[0] == 2:    # bond with its color and thickness
        itype = int(obj[1])
	if itype > self.vizinfo.nbcolor:
          raise StandardError,"bond type too big"
        color = self.vizinfo.bcolor[itype]
        thick = self.vizinfo.brad[itype]
	print >>f,'<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="rgb(%s,%s,%s)" stroke-width="%s" />' % \
              (obj[2],obj[3],obj[5],obj[6],
               color[0]*255,color[1]*255,color[2]*255,thick*self.factor)
        
      elif obj[0] == 3:    # line with its color and thickness
        itype = int(obj[1])
	if itype > self.vizinfo.nlcolor:
          raise StandardError,"line type too big"
        color = self.vizinfo.lcolor[itype]
        thick = self.vizinfo.lrad[itype]
        print >>f,'<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="rgb(%s,%s,%s)" stroke-width="%s" />' % \
              (obj[2],obj[3],obj[5],obj[6],
               color[0]*255,color[1]*255,color[2]*255,thick*self.factor)

      elif obj[0] == 4:    # box line with built-in color and thickness
        color = self.bxcol
        thick = self.bxthick
        print >>f,'<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="rgb(%s,%s,%s)" stroke-width="%s" />' % \
              (obj[2],obj[3],obj[5],obj[6],
               color[0]*255,color[1]*255,color[2]*255,thick*self.factor)

    elif flag == 1:
      x = (obj[0]*self.xpixels) + (self.xpixels/2.0)
      y = (self.ypixels/2.0) - (obj[1]*self.ypixels)
      color = obj[4]
      print >>f,'<text x="%s" y="%s" font-size="%s" font-family="%s" stroke="rgb(%s,%s,%s)" fill="rgb(%s,%s,%s"> %s </text>' % \
              (x,y,obj[3],obj[2],color[0]*255,color[1]*255,color[2]*255,
               color[0]*255,color[1]*255,color[2]*255,obj[5])
  
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
# compare function for the sort method, orders according to z coordinate

def cmprz(a,b):
  if a[4] > b[4]:
    return 1
  elif a[4] < b[4]:
    return -1
  elif a[4] == b[4]:
    return 0

# --------------------------------------------------------------------
# return characteristic distance of simulation domain = max dimension

def compute_distance(box):
  distance = box[3]-box[0]
  if box[4]-box[1] > distance: distance = box[4]-box[1]
  if box[5]-box[2] > distance: distance = box[5]-box[2]
  return distance

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
#   returns the rotation matrix as a string for now

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
  
  matrix = (c11,c12,c13,c21,c22,c23,c31,c32,c33)
  return matrix

# --------------------------------------------------------------------
# fontlist

fontlist = {}
fontlist["t"] = "Times"
fontlist["h"] = "Helvetica"
