# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# matlab tool

oneline = "Create plots via MatLab numerical analysis program"

docstr = """
m = matlab()		       start up MatLab
m.stop()		       shut down MatLab process
    
m.plot(a)                      plot vector A against linear index
m.plot(a,b)	 	       plot B against A
m.plot(a,b,c,d,...)	       plot B against A, D against C, etc
m.mplot(M,N,S,"file",a,b,...)  multiple plots saved to file0000.eps, etc

  each plot argument can be a tuple, list, or Numeric/NumPy vector
  mplot loops over range(M,N,S) and create one plot per iteration
    last args are same as list of vectors for plot(), e.g. 1, 2, 4 vectors
    each plot is made from a portion of the vectors, depending on loop index i
      Ith plot is of b[0:i] vs a[0:i], etc
    series of plots saved as file0000.eps, file0001.eps, etc
    if use xrange(),yrange() then plot axes will be same for all plots

m("c = a + b")                 execute string in MatLab

m.enter()	   	       enter MatLab shell
matlab> c = a + b              type commands directly to MatLab
matlab> exit, quit	       exit MatLab shell
    
m.export("data",range(100),a,...)       create file with columns of numbers

  all vectors must be of equal length
  could plot from file with MatLab commands:
    cols = importdata('data')
    plot(cols(:,1),cols(:,2))

m.select(N)  	               figure N becomes the current plot
  
  subsequent commands apply to this plot

m.hide(N)  	               delete window for figure N
m.save("file")	               save current plot as file.eps

Set attributes for current plot:

m.erase()                      reset all attributes to default values
m.aspect(1.3)                  aspect ratio
m.xtitle("Time")               x axis text
m.ytitle("Energy")             y axis text
m.title("My Plot")             title text
m.title("title","x","y")       title, x axis, y axis text
m.xrange(xmin,xmax)            x axis range
m.xrange()                     default x axis range
m.yrange(ymin,ymax)            y axis range
m.yrange()                     default y axis range
m.xlog()                       toggle x axis between linear and log
m.ylog()                       toggle y axis between linear and log
m.label(x,y,"text")            place label at x,y coords
m.curve(N,'r')                 set color of curve N
m.curve(N,'g','--')            set color and line style of curve N
m.curve(N,'b','-','v')         set color, line style, symbol of curve N

  colors:  'k' = black, 'r' = red, 'g' = green, 'b' = blue
           'm' = magenta, 'c' = cyan, 'y' = yellow
  styles:  '-' = solid, '--' = dashed, ':' = dotted, '-.' = dash-dot
  symbols: '+' = plus, 'o' = circle, '*' = asterik, 'x' = X,
           's' = square, 'd' = diamond, '^' = up triangle,
           'v' = down triangle, '>' = right triangle,
           '<' = left triangle, 'p' = pentagram, 'h' = hexagram
"""

# History
#   8/05, Matt Jones (BYU): original version

# ToDo list
#   allow choice of JPG or PNG or GIF when saving via "saveas" command
#   have vec/array import/export functions to pass variables between
#     Pizza.py and MatLab and have them named in MatLab

# Variables
#   current = index of current figure (1-N)
#   figures = list of figure objects with each plot's attributes
#             so they aren't lost between replots
#   import command to yank MatLab variables back to Python
#   allow for alternate export command to name the variables from Python
#     in MatLab and vice versa

# Imports and external programs

import types, os

try: from DEFAULTS import PIZZA_MATLAB
except: PIZZA_MATLAB = "matlab -nosplash -nodesktop -nojvm"

# Class definition

class matlab:
  
  # --------------------------------------------------------------------

  def __init__(self):
    self.MATLAB = os.popen(PIZZA_MATLAB,'w')
    self.file = "tmp.matlab"
    self.figures = []
    self.select(1)
              
  # --------------------------------------------------------------------

  def stop(self):
    self.__call__("quit")
    del self.MATLAB

  # --------------------------------------------------------------------

  def __call__(self,command):
    self.MATLAB.write(command + '\n')
    self.MATLAB.flush()
    
  # --------------------------------------------------------------------

  def enter(self):
    while 1:
      command = raw_input("matlab> ")
      if command == "quit" or command == "exit": return
      self.__call__(command)

  # --------------------------------------------------------------------
  # write plot vectors to files and plot them

  def plot(self,*vectors):
    if len(vectors) == 1:
      file = self.file + ".%d.1" % self.current
      linear = range(len(vectors[0]))
      self.export(file,linear,vectors[0])
      self.figures[self.current-1].ncurves = 1
    else:
      if len(vectors) % 2: raise StandardError,"vectors must come in pairs"
      for i in range(0,len(vectors),2):
        file = self.file + ".%d.%d" % (self.current,i/2+1)
        self.export(file,vectors[i],vectors[i+1])
      self.figures[self.current-1].ncurves = len(vectors)/2
    self.draw()
    
  # --------------------------------------------------------------------
  # create multiple plots from growing vectors, save to numbered files
  # don't plot empty vector, create a [0] instead
  
  def mplot(self,start,stop,skip,file,*vectors):
    n = 0
    for i in range(start,stop,skip):
      partial_vecs = []
      for vec in vectors:
        if i: partial_vecs.append(vec[:i])
        else: partial_vecs.append([0])
      self.plot(*partial_vecs)

      if n < 10:     newfile = file + "000" + str(n)
      elif n < 100:  newfile = file + "00" + str(n)
      elif n < 1000: newfile = file + "0" + str(n)
      else:          newfile = file + str(n)

      self.save(newfile)
      n += 1

  # --------------------------------------------------------------------
  # write list of equal-length vectors to filename

  def export(self,filename,*vectors):
    n = len(vectors[0])
    for vector in vectors:
      if len(vector) != n: raise StandardError,"vectors must be same length"
    f = open(filename,'w')
    nvec = len(vectors)
    for i in xrange(n):
      for j in xrange(nvec):
        print >>f,vectors[j][i],
      print >>f
    f.close()

  # --------------------------------------------------------------------
  # select plot N as current plot

  def select(self,n):
    self.current = n
    if len(self.figures) < n:
      for i in range(n - len(self.figures)):
        self.figures.append(figure())
    self.draw()

  # --------------------------------------------------------------------
  # delete window for plot N

  def hide(self,n):
    cmd = "set(figure(%d),'Visible','off')" % n
    self.__call__(cmd)

  # --------------------------------------------------------------------
  # save plot to file.eps
  # do not continue until plot file is written out
  #   else script could go forward and change data file
  #   use tmp.done as semaphore to indicate plot is finished

  def save(self,file):
    if os.path.exists("tmp.done"): os.remove("tmp.done")
    cmd = "saveas(gcf,'%s.eps','psc2')" % file
    self.__call__(cmd)
    self.__call__("!touch tmp.done")
    while not os.path.exists("tmp.done"): continue
  
  # --------------------------------------------------------------------
  # restore default attributes by creating a new fig object

  def erase(self):
    fig = figure()
    fig.ncurves = self.figures[self.current-1].ncurves
    self.figures[self.current-1] = fig
    self.draw()
  
  # --------------------------------------------------------------------

  def aspect(self,value):
    self.figures[self.current-1].aspect = value
    self.draw()

  # --------------------------------------------------------------------

  def xrange(self,*values):
    if len(values) == 0:
      self.figures[self.current-1].xlimit = 0
    else:
      self.figures[self.current-1].xlimit = (values[0],values[1])
    self.draw()

  # --------------------------------------------------------------------

  def yrange(self,*values):
    if len(values) == 0:
      self.figures[self.current-1].ylimit = 0
    else:
      self.figures[self.current-1].ylimit = (values[0],values[1])
    self.draw()
    	    
  # --------------------------------------------------------------------

  def label(self,x,y,text):
    self.figures[self.current-1].labels.append((x,y,text))
    self.figures[self.current-1].nlabels += 1  	    
    self.draw()

  # --------------------------------------------------------------------

  def nolabels(self):
    self.figures[self.current-1].nlabel = 0
    self.figures[self.current-1].labels = []
    self.draw()
      
  # --------------------------------------------------------------------

  def title(self,*strings):
    if len(strings) == 1:
      self.figures[self.current-1].title = strings[0]
    else:
      self.figures[self.current-1].title = strings[0]
      self.figures[self.current-1].xtitle = strings[1]
      self.figures[self.current-1].ytitle = strings[2]
    self.draw()

  # --------------------------------------------------------------------

  def xtitle(self,label):
    self.figures[self.current-1].xtitle = label
    self.draw()
    
  # --------------------------------------------------------------------

  def ytitle(self,label):
    self.figures[self.current-1].ytitle = label
    self.draw()
    
  # --------------------------------------------------------------------

  def xlog(self):
    if self.figures[self.current-1].xlog:
      self.figures[self.current-1].xlog = 0
    else:
      self.figures[self.current-1].xlog = 1
    self.draw()
 
  # --------------------------------------------------------------------

  def ylog(self):
    if self.figures[self.current-1].ylog:
      self.figures[self.current-1].ylog = 0
    else:
      self.figures[self.current-1].ylog = 1
    self.draw()
  
  # --------------------------------------------------------------------

  def curve(self,num,*settings):
    fig = self.figures[self.current-1]
    while len(fig.curves) < num: fig.curves.append('')
    value = settings[0]
    if len(settings) >= 2: value += settings[1]
    if len(settings) >= 3: value += settings[2]
    fig.curves[num-1] = value
    self.draw()

  # --------------------------------------------------------------------
  # draw a plot with all its settings
  # just return if no files of vectors defined yet

  def draw(self):
    fig = self.figures[self.current-1]
    if not fig.ncurves: return
    self.__call__("figure(%d)" % self.current)

    if not fig.xlog and not fig.ylog: cmd = "plot"
    elif fig.xlog and not fig.ylog: cmd = "semilogx"
    elif not fig.xlog and fig.ylog: cmd = "semilogy"
    elif fig.xlog and fig.ylog: cmd = "loglog"

    cmd += '('
    for i in range(fig.ncurves):
      file = self.file + ".%d.%d" % (self.current,i+1)
      readcmd = "pizza%d = importdata('%s');" % (i+1,file)
      self.__call__(readcmd)
      if len(fig.curves) > i:
        cmd += "pizza%d(:,1),pizza%d(:,2),'%s'," % (i+1,i+1,fig.curves[i])
      else:
        cmd += "pizza%d(:,1),pizza%d(:,2),''," % (i+1,i+1)
    cmd = cmd[:-1] + ",'LineWidth',1.5)"  # kludge on line width for now
                                          # problem is it applies to all curves
    self.__call__(cmd)                    # should allow other attributes set

    cmd = "set(gca,'PlotBoxAspectRatio',[%g,1,1])" % float(fig.aspect)
    self.__call__(cmd)

                                          # kludge to set title sizes
    self.__call__("title('%s','FontSize',16)" % fig.title)
    self.__call__("xlabel('%s','FontSize',16)" % fig.xtitle)
    self.__call__("ylabel('%s','FontSize',16)" % fig.ytitle)

    if fig.xlimit == 0 or fig.ylimit == 0: self.__call__("axis auto")
    if fig.xlimit: 
      self.__call__("xlim([%g,%g])" % (fig.xlimit[0],fig.xlimit[1]))
    if fig.ylimit: 
      self.__call__("ylim([%g,%g])" % (fig.ylimit[0],fig.ylimit[1]))

    for i in range(fig.nlabels):
      x = fig.labels[i][0]
      y = fig.labels[i][1]
      text = fig.labels[i][2]           # kludge to set label font size
      self.__call__("text(%g,%g,'%s','FontSize',16)" % (x,y,text))      

# --------------------------------------------------------------------
# class to store settings for a single plot

class figure:

  def __init__(self):
    self.ncurves = 0
    self.curves  = []
    self.title   = ""
    self.xtitle  = ""
    self.ytitle  = ""
    self.aspect  = 1.3
    self.xlimit  = 0
    self.ylimit  = 0
    self.xlog    = 0
    self.ylog    = 0
    self.nlabels = 0
    self.labels  = []
