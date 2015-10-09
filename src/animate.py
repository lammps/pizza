# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# animate tool

oneline = "Animate a series of image files"

docstr = """
a = animate("image*.png")     create GUI to animate set of image files
a = animate("image*.png",1)   2nd arg = sort filenames, 0 = no sort, def = 1
    
Actions (same as GUI widgets):

a.first()		      go to first frame
a.prev()                      go to previous frame
a.back()                      play backwards from current frame to start
a.stop()         	      stop on current frame		
a.play()                      play from current frame to end
a.next()		      go to next frame
a.last()		      go to last frame

a.frame(31)    	              set frame slider
a.delay(0.4)     	      set delay slider
"""

# History
#   8/05, Matt Jones (BYU): original version

# ToDo list
#   make image window non-resizable while displaying an image
#   allow for file conversion, e.g. animate("png","image*jpg")

# Variables
#   tkroot = root of entire Tk
#   files = list of file names
#   images = list of Tkimage objects
#   nframes = number of images
#   index = current frame (0 to N-1)
#   loop_flag = set to 1 or -1 when play or back pushed
#               set to 0 when stop is pushed
#   delay_value = delay between frames (secs)
#   delay_msec = delay in millisec

# Imports and external programs

import sys, os, commands, re, glob
from Tkinter import *
from ImageTk import PhotoImage

# Class definition

class animate:

  # --------------------------------------------------------------------

  def __init__(self,filestr,sortflag=1):
    self.loop_flag = 0
    self.delay_value = 0.0
    self.delay_msec = 0

    # convert filestr into full list of files
    
    list = str.split(filestr)
    self.files = []
    for file in list: self.files += glob.glob(file)
    self.nframes = len(self.files)
    if self.nframes == 0: raise StandardError, "No files to load"
    if sortflag: self.files.sort()
    
    # load all images
    
    self.images = []
    for i in xrange(self.nframes):
      self.images.append(PhotoImage(file=self.files[i]))

    # grab Tk instance from main
    
    from __main__ import tkroot
    self.tkroot = tkroot

    # GUI control window
    
    win1 = Toplevel(tkroot)
    win1.title("Pizza.py animate tool")

    holder1 = Frame(win1)
    button1 = Button(holder1,text="<<",command=self.first).pack(side=LEFT)
    button2 = Button(holder1,text="<",command=self.previous).pack(side=LEFT)
    button3 = Button(holder1,text="Back",command=self.back).pack(side=LEFT)
    button4 = Button(holder1,text="Stop",command=self.stop).pack(side=LEFT)
    button5 = Button(holder1,text="Play",command=self.play).pack(side=LEFT)
    button6 = Button(holder1,text=">",command=self.next).pack(side=LEFT)
    button7 = Button(holder1,text=">>",command=self.last).pack(side=LEFT)
    holder1.pack(side=TOP)
    
    holder2 = Frame(win1)
    self.slider_frame = Scale(holder2,from_=0,to=self.nframes-1,
                              command=self.frame,orient=HORIZONTAL,
                              label="       Frame")
    self.slider_delay = Scale(holder2,from_=0.0,to=1.0,resolution=0.1,
                              command=self.delay,orient=HORIZONTAL,
                              label="       Delay")
    self.slider_frame.pack(side=LEFT)
    self.slider_delay.pack(side=LEFT)
    holder2.pack(side=TOP)
    
    holder3 = Frame(win1)
    self.label_frame = Label(holder3)
    self.label_frame.pack(side=LEFT)
    holder3.pack(side=TOP)
    
    # image window

    win2 = Toplevel(tkroot)
    self.image_pane = Label(win2,image=self.images[0])
    self.image_pane.pack(side=BOTTOM)
    tkroot.update_idletasks()              # force window to appear

    # display 1st image
    
    self.index = 0
    self.display(self.index)

  # --------------------------------------------------------------------

  def first(self):
    self.index = 0
    self.display(self.index)

  # --------------------------------------------------------------------

  def last(self):
    self.index = self.nframes - 1
    self.display(self.index)
  
  # --------------------------------------------------------------------

  def previous(self):
    if self.index > 0: self.index -= 1
    self.display(self.index)
  
  # --------------------------------------------------------------------

  def next(self):
    if self.index < self.nframes - 1: self.index += 1
    self.display(self.index)
  
  # --------------------------------------------------------------------

  def back(self):
    if self.loop_flag != 0: return
    self.loop_flag = -1
    if self.index == 0:
      self.index = self.nframes - 1
      self.display(self.index)
    self.loop()
  
  # --------------------------------------------------------------------

  def play(self):
    if self.loop_flag != 0: return
    self.loop_flag = 1
    if self.index == self.nframes - 1:
      self.index = 0
      self.display(self.index)
    self.loop()
  
  # --------------------------------------------------------------------

  def stop(self):
    self.loop_flag = 0
  
  # --------------------------------------------------------------------
  # loop forward or back until end of animation
  
  def loop(self):
    if self.loop_flag == 1 and self.index == self.nframes - 1:
      self.loop_flag = 0
    if self.loop_flag == -1 and self.index == 0:
      self.loop_flag = 0
    if self.loop_flag == 0: return

    self.index += self.loop_flag
    self.display(self.index)
    self.tkroot.update_idletasks()
    self.tkroot.after(self.delay_msec,self.loop)

  # --------------------------------------------------------------------
  # display a frame corresponding to iframe
  
  def display(self,iframe):
    self.image_pane.configure(image=self.images[iframe])
    self.slider_frame.set(iframe)
    textstr = "Frame: %d    File: %s" % (iframe,self.files[iframe])
    self.label_frame.configure(text=textstr)

  # --------------------------------------------------------------------

  def frame(self,value):
    self.index = int(value)
    self.display(self.index)

  # --------------------------------------------------------------------

  def delay(self,value):
    self.delay_value = float(value)
    self.slider_delay.set(self.delay_value)
    self.delay_msec = int(1000*self.delay_value)
  
