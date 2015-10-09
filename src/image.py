# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# image tool

oneline = "View and manipulate images"

docstr = """
i = image("my1.gif my2.gif")    display thumbnails of matching images
i = image("*.png *.gif")        wildcards allowed
i = image("*.png *.gif",0)      2nd arg = sort filenames, 0 = no sort, def = 1
i = image("")	   		blank string matches all image suffixes
i = image()			no display window opened if no arg

  image suffixes for blank string = *.png, *.bmp, *.gif, *.tiff, *.tif
  click on a thumbnail to view it full-size
  click on thumbnail again to remove full-sized version

i.view("*.png *.gif")	        display thumbnails of matching images

  view arg is same as constructor arg

i.convert("image*.svg","new*.png")                      each SVG file to PNG
i.convert("image*.svg","new*.jpg","-quality 50")        3rd arg is switch
i.convert("image*.png","movie.mpg")                     all PNGs to MPG movie
i.convert("image*.png","movie.mpg","-resize 128x128")   3rd arg is switch
i.montage("","image*.png","plot*.png","two*.png")       image + plot = two
i.montage("-geometry 512x512","i*.png","new.png")       1st arg is switch

  convert with 2 wildcard args loops over 1st set of files to make 2nd set
  convert with not all wildcard args will issue single convert command
  montage with all wildcard args loops over 1st set of files,
    combines with one file from other sets, to make last set of files
  montage with not all wildcard args will issue single montage command
"""

# History
#   8/05, Matt Jones (BYU): original version
#   9/05, Steve Plimpton: added convert() and montage() methods

# ToDo list

# Variables

# Imports and external programs

import sys, os, commands, re, glob
from math import *
from Tkinter import *
import Pmw
import Image,ImageTk

try: from DEFAULTS import PIZZA_CONVERT
except: PIZZA_CONVERT = "convert"
try: from DEFAULTS import PIZZA_MONTAGE
except: PIZZA_MONTAGE = "montage"

# Class definition

class image:
  
  # --------------------------------------------------------------------

  def __init__(self,filestr=None,sortflag=1):
    if filestr == None: return
    self.view(filestr,sortflag)

  def view(self,filestr,sortflag):

    # convert filestr into full list of files

    if filestr == "": filestr = ' '.join(extensions)
    list = str.split(filestr)
    files = []
    for file in list: files += glob.glob(file)
    if len(files) == 0: raise StandardError, "no image files to load"
    if sortflag: files.sort()

    # grab Tk instance from main

    from __main__ import tkroot

    # GUI control window
    
    gui = Toplevel(tkroot)
    gui.title('Pizza.py image tool')
    
    scroll = \
      Pmw.ScrolledFrame(gui,usehullsize=1,hull_width=420,hull_height=500)
    pane = scroll.interior()
    
    ncolumns = 4
    for i in xrange(len(files)):
      
      # create new row frame if 1st in column
      
      if i % ncolumns == 0: rowframe = Frame(pane)
      oneframe = Frame(rowframe)
      
      # create a thumbnail of image
      
      im = Image.open(files[i])
      imt = im.copy()
      imt.thumbnail((60,60),Image.ANTIALIAS)
      basename = os.path.basename(files[i])
      imt.save("tmp." + basename)
      thumbnail = ImageTk.PhotoImage(file = "tmp." + basename)
      os.remove("tmp." + basename)
      
      # read in full size image
      # create a thumbnail object that links to it
      # create button that calls the thumbnail, label with filename
      # buttton needs to store thumbnail else it is garbage collected

      big = ImageTk.PhotoImage(file=files[i])
      obj = thumbnails(gui,files[i],big,thumbnail)
      Button(oneframe,image=thumbnail,command=obj.display).pack(side=TOP)
      Label(oneframe,text=basename).pack(side=BOTTOM)
      
      # pack into row frame
      
      oneframe.pack(side=LEFT)
      if (i+1) % ncolumns == 0: rowframe.pack(side=TOP)
      
    if len(files) % ncolumns != 0: rowframe.pack(side=TOP)
    scroll.pack(side=LEFT)     

  # --------------------------------------------------------------------
  # wrapper on ImageMagick convert command
  
  def convert(self,file1,file2,switch=""):
    if file1.find('*') < 0 or file2.find('*') < 0:
      cmd = "%s %s %s %s" % (PIZZA_CONVERT,switch,file1,file2)
      commands.getoutput(cmd)
      return

    index = file1.index('*')
    pre1 = file1[:index]
    post1 = file1[index+1:]
    index = file2.index('*')
    pre2 = file2[:index]
    post2 = file2[index+1:]
    expr = "%s(.*)%s" % (pre1,post1)

    filelist = glob.glob(file1)
    for file1 in filelist:
      middle = re.search(expr,file1).group(1)
      file2 = "%s%s%s" % (pre2,middle,post2)
      cmd = "%s %s %s %s" % (PIZZA_CONVERT,switch,file1,file2)
      print middle,
      sys.stdout.flush()
      commands.getoutput(cmd)
    print

  # --------------------------------------------------------------------
  # wrapper on ImageMagick montage command

  def montage(self,switch,*fileargs):
    nsets = len(fileargs)
    if nsets < 2: raise StandardError,"montage requires 2 or more file args"

    for i in range(nsets):
      if fileargs[i].find('*') < 0:
        cmd = "%s %s" % (PIZZA_MONTAGE,switch)
        for j in range(nsets): cmd += " %s" % fileargs[j]
        commands.getoutput(cmd)
        return
    
    nfiles = len(glob.glob(fileargs[0]))
    filesets = []
    for i in range(nsets-1):
      filesets.append(glob.glob(fileargs[i]))
      if len(filesets[-1]) != nfiles:
        raise StandardError,"each montage arg must represent equal # of files"

    index = fileargs[0].index('*')
    pre1 = fileargs[0][:index]
    post1 = fileargs[0][index+1:]
    index = fileargs[-1].index('*')
    preN = fileargs[-1][:index]
    postN = fileargs[-1][index+1:]
    expr = "%s(.*)%s" % (pre1,post1)

    for i in range(nfiles):
      cmd = "%s %s" % (PIZZA_MONTAGE,switch)
      for j in range(nsets-1): cmd += " %s" % filesets[j][i]
      middle = re.search(expr,filesets[0][i]).group(1)
      fileN = "%s%s%s" % (preN,middle,postN)
      cmd += " %s" % fileN
      commands.getoutput(cmd)
      print middle,
      sys.stdout.flush()
    print

# --------------------------------------------------------------------
# thumbnail class
     
class thumbnails:

  def __init__(self,root,name,bigimage,thumbimage):
    self.root = root
    self.big = bigimage
    self.thumb = thumbimage
    self.name = name
    self.bigexist = 0
    self.window = None
    
  def display(self):

    # destroy the big image window
    
    if self.bigexist:
      self.bigexist = 0
      if self.window:
        self.window.destroy()
	self.window = None
    
    # create a new window with the big image
    
    else:
      self.bigexist = 1
      self.window = Toplevel(self.root)
      Label(self.window,text=self.name).pack(side=TOP)
      Label(self.window,image=self.big).pack(side=BOTTOM)

# --------------------------------------------------------------------
# list of file extensions to test for
# could add any extensions that PIL recognizes

extensions = ["*.png", "*.bmp", "*.gif", "*.tiff", "*.tif"] 	
