# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# vcr tool

oneline = "VCR-style GUI for 3d interactive OpenGL visualization"

docstr = """
v = vcr(gl1,gl2,...)       start vcr GUI with one or more gl windows
v.add(gl)                  add a gl window to vcr GUI

Actions (same as GUI widgets):

v.first()		   go to first frame
v.prev()                   go to previous frame
v.back()                   play backwards from current frame to start
v.stop()         	   stop on current frame		
v.play()                   play from current frame to end
v.next()		   go to next frame
v.last()		   go to last frame

v.frame(31)    	           set frame slider
v.delay(0.4)     	   set delay slider
v.q(5)          	   set quality slider

v.xaxis()		   view scene from x axis
v.yaxis()  		   view scene from y axis
v.zaxis()		   view scene from z axis
v.box()		 	   toggle bounding box
v.axis() 	           toggle display of xyz axes
v.norm()	           recenter and resize the view
v.ortho()	           toggle ortho/perspective button
v.reload()           	   reload all frames from gl viewer data files

v.clipxlo(0.2)             clip scene at x lo fraction of box
v.clipxhi(1.0)             clip at x hi
v.clipylo(0.2)             clip in y
v.clipyhi(1.0)
v.clipzlo(0.2)             clip in z
v.clipzhi(1.0)

v.save()		   save current scene to file.png
v.file("image")		   set filename
v.saveall()		   toggle save-all checkbox
"""

# History
#   8/05, Matt Jones (BYU): original version
#   9/05, Steve Plimpton: modified for GL viewer

# ToDo list

# Variables
#   tkroot = root of entire Tk
#   index = current frame (0 to N-1)
#   nframes = number of frames in data
#   loop_flag = set to 1 or -1 when play or back pushed
#               set to 0 when stop is pushed
#   delay_value = delay between frames (secs)
#   delay_msec = delay in millisec

# Imports and external programs

from Tkinter import *
import types

# Class definition

class vcr:

  # --------------------------------------------------------------------

  def __init__(self,*views):
    self.boxflag = 0
    self.axisflag = 0
    self.orthoflag = 1
    self.saveflag = 0

    self.loop_flag = 0
    self.delay_value = 0.0
    self.delay_msec = 0

    # load data for each viewer
    # if each viewer has different data set, nframes is for 1st viewer

    if len(views) == 0: raise StandardError,"must have at least one GL viewer"
    self.viewlist = []
    for view in views: self.viewlist.append(view)
    for view in self.viewlist: view.reload()
    self.nframes = self.viewlist[0].nframes

    from __main__ import tkroot
    self.tkroot = tkroot
    root = Toplevel(tkroot)
    root.title('Pizza.py vcr tool')

    frame1 = Frame(root)
    Button(frame1,text="<<",command=self.first).pack(side=LEFT)
    Button(frame1,text="<",command=self.previous).pack(side=LEFT)
    Button(frame1,text="Back",command=self.back).pack(side=LEFT)
    Button(frame1,text="Stop",command=self.stop).pack(side=LEFT)
    Button(frame1,text="Play",command=self.play).pack(side=LEFT)
    Button(frame1,text=">",command=self.next).pack(side=LEFT)
    Button(frame1,text=">>",command=self.last).pack(side=LEFT)
    frame1.pack()
    
    frame2 = Frame(root)
    self.slider_frame = Scale(frame2,from_=0,to=self.nframes-1,
                              command=self.frame,orient=HORIZONTAL,
                              label="       Frame")
    self.slider_quality = Scale(frame2,from_=3,to=12,orient=HORIZONTAL,
                                label="       Quality")
    self.slider_quality.bind('<ButtonRelease-1>',self.q)
    self.slider_quality.set(5)
    self.slider_delay = Scale(frame2,from_=0.0,to=1.0,resolution=0.1,
                              command=self.delay,orient=HORIZONTAL,
                              label="       Delay")
    self.slider_frame.pack(side=LEFT)
    self.slider_quality.pack(side=LEFT) 
    self.slider_delay.pack(side=LEFT)
    frame2.pack()

    frame3 = Frame(root)
    self.slider_xlo = Scale(frame3,from_=0.0,to=1.0,resolution=0.01,
                            command=self.clipxlo,orient=HORIZONTAL,
                            showvalue=0,label="       xclip")
    self.slider_ylo = Scale(frame3,from_=0.0,to=1.0,resolution=0.01,
                            command=self.clipylo,orient=HORIZONTAL,
                            showvalue=0,label="       yclip")
    self.slider_zlo = Scale(frame3,from_=0.0,to=1.0,resolution=0.01,
                            command=self.clipzlo,orient=HORIZONTAL,
                            showvalue=0,label="       zclip")
    self.slider_xlo.set(0.0)
    self.slider_ylo.set(0.0)
    self.slider_zlo.set(0.0)
    self.slider_xlo.pack(side=LEFT)
    self.slider_ylo.pack(side=LEFT)
    self.slider_zlo.pack(side=LEFT)
    frame3.pack()

    frame4 = Frame(root)
    self.slider_xhi = Scale(frame4,from_=0.0,to=1.0,resolution=0.01,
                            command=self.clipxhi,orient=HORIZONTAL,
                            showvalue=0)
    self.slider_yhi = Scale(frame4,from_=0.0,to=1.0,resolution=0.01,
                            command=self.clipyhi,orient=HORIZONTAL,
                            showvalue=0)
    self.slider_zhi = Scale(frame4,from_=0.0,to=1.0,resolution=0.01,
                            command=self.clipzhi,orient=HORIZONTAL,
                            showvalue=0)
    self.slider_xhi.set(1.0)
    self.slider_yhi.set(1.0)
    self.slider_zhi.set(1.0)
    self.slider_xhi.pack(side=LEFT)
    self.slider_yhi.pack(side=LEFT)
    self.slider_zhi.pack(side=LEFT)
    frame4.pack()

    frame5 = Frame(root)
    Button(frame5,text="X",command=self.xaxis).pack(side=LEFT)
    Button(frame5,text="Y",command=self.yaxis).pack(side=LEFT)
    Button(frame5,text="Z",command=self.zaxis).pack(side=LEFT)
    Button(frame5,text="Box",command=self.box).pack(side=LEFT)
    Button(frame5,text="Axes",command=self.axis).pack(side=LEFT)
    Button(frame5,text="Recenter",command=self.recenter).pack(side=LEFT)
    self.button_ortho = Button(frame5,text="Persp",command=self.ortho)
    self.button_ortho.pack(side=LEFT)
    Button(frame5,text="Reload",command=self.reload).pack(side=LEFT)
    frame5.pack()

    frame6 = Frame(root)
    Button(frame6,text="Save As:",command=self.save).pack(side=LEFT)
    self.entry_file = Entry(frame6,width = 16)
    self.entry_file.insert(0,"image")
    self.entry_file.pack(side=LEFT) 
    self.button_save = Checkbutton(frame6,text="SaveAll",command=self.saveall)
    self.button_save.pack(side=LEFT)
    frame6.pack()

    frame7 = Frame(root)
    Label(frame7,text="Mouse LMR = Shift/Rotate/Zoom").pack(padx=5,side=LEFT)
    Label(frame7,text="Axes XYZ = Red/Green/Blue").pack(padx=5,side=LEFT)
    frame7.pack()

    frame8 = Frame(root)
    self.label_frame = Label(frame8,text="Frame 0")
    self.label_frame.pack(padx=5,side=LEFT)
    self.label_time = Label(frame8,text="Time 0")
    self.label_time.pack(padx=5,side=LEFT)
    self.label_atoms = Label(frame8,text="Atoms 0")
    self.label_atoms.pack(padx=5,side=LEFT)
    frame8.pack()

    # display 1st image
    
    self.index = 0
    self.display()

  # --------------------------------------------------------------------

  def add(self,view):
    self.viewlist.append(view)
    view.reload()
    view.display(self.index)

  # --------------------------------------------------------------------

  def first(self):
    self.index = 0
    self.slider_frame.set(self.index)

  # --------------------------------------------------------------------

  def last(self):
    self.index = self.nframes - 1
    self.slider_frame.set(self.index)
  
  # --------------------------------------------------------------------

  def next(self):
    if self.index < self.nframes - 1:
      self.index += 1
      self.slider_frame.set(self.index)
  
  # --------------------------------------------------------------------

  def previous(self):
    if self.index > 0:
      self.index -= 1
      self.slider_frame.set(self.index)
  
  # --------------------------------------------------------------------
  # play backward loop
  # disable GL caching while animating
  
  def back(self):
    if self.loop_flag != 0: return
    self.loop_flag = -1
    for view in self.viewlist: view.cache = 0
    if self.index == 0:
      self.index = self.nframes - 1
      self.slider_frame.set(self.index)
      self.tkroot.update()
    if self.saveflag: self.saveloop(0)
    self.tkroot.after(self.delay_msec,self.loop)
  
  # --------------------------------------------------------------------
  # play forward loop
  # disable GL caching while animating

  def play(self):
    if self.loop_flag != 0: return
    self.loop_flag = 1
    for view in self.viewlist: view.cache = 0
    if self.index == self.nframes - 1:
      self.index = 0
      self.slider_frame.set(self.index)
      self.tkroot.update()
    if self.saveflag: self.saveloop(0)
    self.tkroot.after(self.delay_msec,self.loop)

  # --------------------------------------------------------------------
  # stop looping
  # re-enable GL caching when not animating

  def stop(self):
    self.loop_flag = 0
    for view in self.viewlist: view.cache = 1
  
  # --------------------------------------------------------------------
  # loop forward or back until end of animation
  # if save flag is set, change file name and save snapshots
  
  def loop(self):
    if self.loop_flag == 1 and self.index == self.nframes - 1:
      self.loop_flag = 0
    if self.loop_flag == -1 and self.index == 0:
      self.loop_flag = 0
    if self.loop_flag == 0:
      if self.saveflag: self.saveloop(-1)
      return

    self.index += self.loop_flag
    self.slider_frame.set(self.index)

    # call to view.display() should not be necessary
    # since slider_frame already called it
    # but seems to be necessary before GL save() saves window to file
    # else get previous image saved into file at each frame
    
    if self.saveflag:
      for view in self.viewlist: time,natoms = view.display(self.index)
      self.saveloop(1)

    self.tkroot.update_idletasks()
    self.tkroot.after(self.delay_msec,self.loop)

  # --------------------------------------------------------------------
  # display a new frame, called in one of 3 ways:
  #   when frame slider is set by user in GUI
  #   when user calls v.frame(N)
  #   when other vcr method calls slider_frame.set(N)

  def frame(self,value):
    self.index = int(value)
    self.display()
  
  # --------------------------------------------------------------------

  def delay(self,value):
    self.delay_value = float(value)
    self.slider_delay.set(self.delay_value)
    self.delay_msec = int(1000*self.delay_value)

  # --------------------------------------------------------------------

  def q(self,value):
    if type(value) is not types.IntType: value = self.slider_quality.get()
    self.slider_quality.set(value)
    for view in self.viewlist: view.q(value)

  # --------------------------------------------------------------------

  def xaxis(self):
    for view in self.viewlist:
      view.zoom(1)
      view.shift(0,0)
      view.rotate(90,0)

  # --------------------------------------------------------------------

  def yaxis(self):
    for view in self.viewlist:
      view.zoom(1)
      view.shift(0,0)
      view.rotate(90,270)

  # --------------------------------------------------------------------

  def zaxis(self):
    for view in self.viewlist:
      view.zoom(1)
      view.shift(0,0)
      view.rotate(0,270)

  # --------------------------------------------------------------------

  def box(self):
    if self.boxflag:
      self.boxflag = 0
      for view in self.viewlist: view.box(0)
    else:
      self.boxflag = 1
      for view in self.viewlist: view.box(1)
  
  # --------------------------------------------------------------------

  def axis(self):
    if self.axisflag:
      self.axisflag = 0
      for view in self.viewlist: view.axis(0)
    else:
      self.axisflag = 1
      for view in self.viewlist: view.axis(1)
      
  # --------------------------------------------------------------------

  def recenter(self):
    for view in self.viewlist:
      view.shift(0,0)
      view.zoom(1)

  # --------------------------------------------------------------------

  def ortho(self):
    if self.orthoflag:
      self.orthoflag = 0
      self.button_ortho.config(text="Ortho")
      for view in self.viewlist: view.ortho(0)
    else:
      self.orthoflag = 1
      self.button_ortho.config(text="Persp")
      for view in self.viewlist: view.ortho(1)
  
  # --------------------------------------------------------------------

  def reload(self):
    for view in self.viewlist: view.reload()
    self.display()

  # --------------------------------------------------------------------

  def clipxlo(self,value):
    value = float(value)
    hi = float(self.slider_xhi.get())
    if value > hi: value = hi
    self.slider_xlo.set(value)
    for view in self.viewlist: view.clip('xlo',value)

  # --------------------------------------------------------------------

  def clipxhi(self,value):
    value = float(value)
    lo = float(self.slider_xlo.get())
    if value < lo: value = lo
    self.slider_xhi.set(value)
    for view in self.viewlist: view.clip('xhi',value)

  # --------------------------------------------------------------------

  def clipylo(self,value):
    value = float(value)
    hi = float(self.slider_yhi.get())
    if value > hi: value = hi
    self.slider_ylo.set(value)
    for view in self.viewlist: view.clip('ylo',value)

  # --------------------------------------------------------------------

  def clipyhi(self,value):
    value = float(value)
    lo = float(self.slider_ylo.get())
    if value < lo: value = lo
    self.slider_yhi.set(value)
    for view in self.viewlist: view.clip('yhi',value)

  # --------------------------------------------------------------------

  def clipzlo(self,value):
    value = float(value)
    hi = float(self.slider_zhi.get())
    if value > hi: value = hi
    self.slider_zlo.set(value)
    for view in self.viewlist: view.clip('zlo',value)

  # --------------------------------------------------------------------

  def clipzhi(self,value):
    value = float(value)
    lo = float(self.slider_zlo.get())
    if value < lo: value = lo
    self.slider_zhi.set(value)
    for view in self.viewlist: view.clip('zhi',value)

  # --------------------------------------------------------------------
  # set filename for saving
  
  def file(self,newtext):
    oldtext = self.entry_file.get()
    self.entry_file.delete(0,len(oldtext))
    self.entry_file.insert(0,newtext)

  # --------------------------------------------------------------------
  # toggle save all checkbox
  
  def saveall(self):
    if self.saveflag:
      self.saveflag = 0
      self.button_save.deselect()
    else:
      self.saveflag = 1
      self.button_save.select()
  
  # --------------------------------------------------------------------
  # save current image to file
  # if multiple windows change filenames to file.0.png, file.1.png, etc
  
  def save(self):
    file = self.entry_file.get()
    if len(self.viewlist) == 1:
      self.viewlist[0].file = file
      self.viewlist[0].save()
    else:
      n = 0
      for view in self.viewlist:
        view.file = "%s.%d" % (file,n)
        view.save()
        n += 1

  # --------------------------------------------------------------------
  # save images when in a play/back loop
  # flag 0 = first save, flag 1 = continuing save, flag -1 = stop
  
  def saveloop(self,flag):
    if flag == -1:
      self.file(self.fileroot)
      return
    if flag == 0: self.fileroot = self.entry_file.get()

    if self.index < 10:     file = self.fileroot + "000" + str(self.index)
    elif self.index < 100:  file = self.fileroot + "00" + str(self.index)
    elif self.index < 1000: file = self.fileroot + "0" + str(self.index)
    else:                   file = self.fileroot + str(self.index)

    self.file(file)
    self.save()

  # --------------------------------------------------------------------
  # display index frame and set status strings
  
  def display(self):
    for view in self.viewlist: time,natoms = view.display(self.index)
    self.label_frame.config(text="Frame: %d" % self.index)
    self.label_time.config(text="Time: %d" % time)
    self.label_atoms.config(text="Atoms: %d" % natoms)
