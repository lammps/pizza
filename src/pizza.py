#!/usr/bin/env python -i

# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# Change log:
#   8/05, Steve Plimpton (SNL): original version
#  12/09, David Hart (SNL): except hook for Tkinter no-display error
#   5/11, David Hart (SNL): began list of excludes for no-display machines

# ToDo list:

# Help strings:

version = "9 Oct 2015"

intro = """
Pizza.py (%s), a toolkit written in Python
type ? for help, CTRL-D to quit
"""

help = """
pizza.py switch arg(s) switch arg(s) ...
  -s	 	              silent (else print start-up help)
  -t log dump raster	      load only these tools
  -x raster rasmol	      load all tools except these
  -f mine.py arg1 arg2        run script file with args
  -c "vec = range(100)"	      run Python command
  -q 	   		      quit (else interactive)

Everything typed at the ">" prompt is a Python command

Additional commands available at ">" prompt:
  ?                           print help message
  ??                          one-line for each tool and script
  ? raster                    list tool commands or script syntax
  ?? energy.py                full documentation of tool or script 
  !ls -l                      shell command
  @cd ..                      cd to a new directory
  @log tmp.log                log all commands typed so far to file
  @run block.py arg1 arg2     run script file with args
  @time d = dump("*.dump")    time a command
  
Tools:
"""

# -------------------------------------------------------------------------
# modules needed by pizza.py

import sys, commands, os, string, exceptions, glob, re
from time import clock

# readline not available in all Pythons

try:
  import readline
  readline_flag = 1
except ImportError, exception:
  print "readline option not available"
  readline_flag = 0

# create global Tk root if Tkinter is loaded
# used by all tools that do GUIs via Tkinter

nodisplay = False
try:
  import Tkinter
  tkroot = Tkinter.Tk()
  tkroot.withdraw()
except ImportError, exception:
  nodisplay = True
  pass
except Exception, exception:
  nodisplay = True
  pass

# -------------------------------------------------------------------------
# error trap that enables special commands at interactive prompt

def trap(type,value,tback):
   global argv

   # only check SyntaxErrors
   
   if not isinstance(value,exceptions.SyntaxError):
     sys.__excepthook__(type,value,tback)
     return
        
   # special commands at top level only, not in indented text entry
   
   if value.text[0].isspace():
     sys.__excepthook__(type,value,tback)
     return

   # ? = top-level help
   # ?? = one-line description of each tool and script
   # ? name = one-line for each tool command or script purpose/syntax
   # ?? name = entire documentation for tool or script
   # name with no .py suffix is tool, name with .py suffix is script

   if value.text[0] == "?":
     words = value.text.split()
     
     if len(words) == 1 and words[0] == "?":
       print intro[1:] % version
       print help[1:]," ",
       for tool in tools: print tool,
       print
       
     elif len(words) == 1 and words[0] == "??":
       for tool in tools:
         exec "oneline = oneline_%s" % tool
         print "%-11s%s" % (tool,oneline)
       print
       
       scripts = []
       for dir in PIZZA_SCRIPTS[1:]:
         list = glob.glob("%s/*.py" % dir)
         list.sort()
         scripts += list
       for script in scripts:
         filename = os.path.basename(script)
         lines = open(script,'r').readlines()
         flag = 0
         for line in lines:
           if line.find("Purpose:") >= 0:
             flag = 1
             break
         if flag: doc = line[line.find("Purpose:")+8:]
         else: doc = " not available\n"
         print "%-20s%s" % (filename,doc),
           
     elif len(words) == 2 and words[0] == "?":
       if words[1][-3:] == ".py":
         fileflag = 0
         for dir in PIZZA_SCRIPTS:
           filename = "%s/%s" % (dir,words[1])
           if os.path.isfile(filename):
             fileflag = 1
             lineflag = 0
             lines = open(filename,'r').readlines()
             for line in lines:
               if line.find("# Purpose:") >= 0: print line[2:],
               if line.find("# Syntax:") >= 0:
                 lineflag = 1
                 break
             if not lineflag: print "%s has no Syntax line" % words[1]
             else: print line[2:],
             break
         if not fileflag:
           print "%s is not a recognized script" % words[1]
           
       else:
         if words[1] in tools:
           exec "txt = docstr_%s" % words[1]
           txt = re.sub("\n\s*\n","\n",txt)
           txt = re.sub("\n .*","",txt)
           exec "print oneline_%s" % words[1]
           print txt
         else:
           print "%s is not a recognized tool" % words[1]
           
     elif len(words) == 2 and words[0] == "??":
       if words[1][-3:] == ".py":
         fileflag = 0
         for dir in PIZZA_SCRIPTS:
           filename = "%s/%s" % (dir,words[1])
           if os.path.isfile(filename):
             fileflag = 1
             lines = open(filename,'r').readlines()
             for line in lines:
               if len(line.strip()) == 0: continue
               if line[0] == '#': print line,
               else: break
             break
         if not fileflag:
           print "%s is not a recognized script" % words[1]

       else:
         if words[1] in tools:
           exec "print oneline_%s" % words[1]
           exec "print docstr_%s" % words[1]
         else:
           print "%s is not a recognized class" % words[1]
         
     return

   # shell command like !ls, !ls -l

   if value.text[0] == "!":
     os.system(value.text[1:])
     return

   # @ commands = @cd, @log, @run, @time
   # for run and time, use namespace in execfile and exec commands
   #   else variables defined in script/command
   #   won't be set in top-level Pizza.py
   
   if value.text[0] == "@":
     words = value.text.split()
     if words[0][1:] == "cd":
       os.chdir(words[1])
       return
     elif words[0][1:] == "log":
       if readline_flag == 0:
         print "cannot use @log without readline module"
         return
       f = open(words[1],"w")
       print >>f,"# pizza.py log file\n"
       nlines = readline.get_current_history_length()
       for i in xrange(1,nlines):
         print >>f,readline.get_history_item(i)
       f.close()
       return
     elif words[0][1:] == "run":
       argv = words[1:]
       file = argv[0]
       flag = 0
       for dir in PIZZA_SCRIPTS:
         fullfile = dir + '/' + file
         if os.path.exists(fullfile):
           flag = 1
           print "Executing file:",fullfile
           execfile(fullfile,namespace)
           break
       if not flag: print "Could not find file",file
       return
     elif words[0][1:] == "time":
       cmd = string.join(words[1:])
       t1 = clock()
       exec cmd in namespace
       t2 = clock()
       print "CPU time = ",t2-t1
       return
     
   # unrecognized command, let system handle error
   
   sys.__excepthook__(type,value,tback)

# -------------------------------------------------------------------------
# process command-line switches
# store scripts and commands in tasks list

silent = 0
yes_tools = []
no_tools = []
tasks = []
quitflag = 0

iarg = 1
while (iarg < len(sys.argv)):
  if (sys.argv[iarg][0] != '-'):
    print "ERROR: arg is not a switch: %s" % (sys.argv[iarg])
    sys.exit()
  if (sys.argv[iarg] == "-s"):
    silent = 1
    iarg += 1
  elif (sys.argv[iarg] == "-t"):
    jarg = iarg + 1
    while (jarg < len(sys.argv) and sys.argv[jarg][0] != '-'):
      yes_tools.append(sys.argv[jarg])
      jarg += 1
    iarg = jarg
  elif (sys.argv[iarg] == "-x"):
    jarg = iarg + 1
    while (jarg < len(sys.argv) and sys.argv[jarg][0] != '-'):
      no_tools.append(sys.argv[jarg])
      jarg += 1
    iarg = jarg

  # allow for "--" as arg to script and not Pizza.py arg
    
  elif (sys.argv[iarg] == "-f"):
    jarg = iarg + 1
    list = []
    while (jarg < len(sys.argv) and
           (sys.argv[jarg][0] != '-' or
            (len(sys.argv[jarg]) >= 3 and sys.argv[jarg][0:2] == "--"))):
      list.append(sys.argv[jarg])
      jarg += 1
    task = ("script",list)
    tasks.append(task)
    iarg = jarg
    
  elif (sys.argv[iarg] == "-c"):
    jarg = iarg + 1
    list = []
    while (jarg < len(sys.argv) and sys.argv[jarg][0] != '-'):
      list.append(sys.argv[jarg])
      jarg += 1
    task = ("command",list)
    tasks.append(task)
    iarg = jarg
  elif (sys.argv[iarg] == "-q"):
    quitflag = 1
    iarg += 1
  else:
    print "ERROR: unknown switch: %s" % (sys.argv[iarg])
    sys.exit()

# print intro message

if not silent: print intro[1:] % version,

# error test on m,x command-line switches

if len(yes_tools) > 0 and len(no_tools) > 0:
  print "ERROR: cannot use -t and -x switches together"
  sys.exit()
  
# -------------------------------------------------------------------------
# tools = list of tool names to import
# if -t switch was used, tools = just those files
# else scan for *.py files in all dirs in PIZZA_TOOLS list
#   and then Pizza.py src dir (sys.path[0])

if not silent: print "Loading tools ..."
if not silent and nodisplay: print "Display not available ... no GUIs"

try: from DEFAULTS import PIZZA_TOOLS
except: PIZZA_TOOLS = []
PIZZA_TOOLS = map(os.path.expanduser,PIZZA_TOOLS)
PIZZA_TOOLS.append(sys.path[0])

if len(yes_tools) > 0: tools = yes_tools
else:
  tools = []
  for dir in PIZZA_TOOLS:
    tools += glob.glob(dir + "/*.py")
  for i in range(len(tools)):
    tools[i] = os.path.basename(tools[i])
    tools[i] = tools[i][:-3]

# remove duplicate entries, reverse enables removing all but first entry

tools.reverse()
for tool in tools:
  while tools.count(tool) > 1: tools.remove(tool)
tools.reverse()

# remove tools in EXCLUDE list and command-line -x list

try: from DEFAULTS import PIZZA_EXCLUDE
except: PIZZA_EXCLUDE = []
for tool in PIZZA_EXCLUDE:
  if tool in tools: tools.remove(tool)
for tool in no_tools:
  if tool in tools: tools.remove(tool)

# add PIZZA_TOOLS dirs to front of module search path (sys.path)
# import each tool as a Python module and its documentation strings
# restore sys.path

sys.path = PIZZA_TOOLS + sys.path

failed = []
for tool in tools:
  #print "loading tool '%s'"%tool
  if nodisplay and tool in ['gl']:
    failed.append(tool)
    continue
  try:
    exec "from %s import %s" % (tool,tool)
    exec "from %s import oneline as oneline_%s" % (tool,tool)
    exec "from %s import docstr as docstr_%s" % (tool,tool)
  except Exception, exception:
    print "%s tool did not load:" % tool
    print " ",exception
    failed.append(tool)

for dir in PIZZA_TOOLS: sys.path = sys.path[1:]

# final list of tools: remove tools where import failed, sort them

for tool in failed: tools.remove(tool)
tools.sort()

# add current working dir to sys.path so user can import own modules
# cwd isn't in sys.path when Pizza.py is launched

sys.path.insert(0,'')

# -------------------------------------------------------------------------
# PIZZA_SCRIPTS = list of dirs to look in to find scripts

try: from DEFAULTS import PIZZA_SCRIPTS
except: PIZZA_SCRIPTS = []
PIZZA_SCRIPTS = map(os.path.expanduser,PIZZA_SCRIPTS)
PIZZA_SCRIPTS.insert(0,'.')
PIZZA_SCRIPTS.append(sys.path[1][:-3] + "scripts")   # path for pizza.py

# run specified script files and commands in order specified
# put arguments in argv so script can access them
# check list of PIZZA_SCRIPTS dirs to find script file
# catch errors so pizza.py will continue even if script is bad
# traceback logic prints where in the script the error occurred

for task in tasks:
  if task[0] == "script":
    argv = task[1]
    file = argv[0]
    try:
      flag = 0
      for dir in PIZZA_SCRIPTS:
        fullfile = dir + '/' + file
        if os.path.exists(fullfile):
          print "Executing file:",fullfile
          execfile(fullfile)
          flag = 1
          break
      if not flag: print "Could not find file",file
    except StandardError, exception:
      (type,value,tback) = sys.exc_info()
      print type,value,tback
      type = str(type)
      type = type[type.find(".")+1:]
      print "%s with value: %s" % (type,value)
      tback = tback.tb_next
      while tback:
        print "error on line %d of file %s" % \
              (tback.tb_lineno,tback.tb_frame.f_code.co_filename)
        tback = tback.tb_next
  elif task[0] == "command":
    argv = task[1]
    cmd = ""
    for arg in argv: cmd += arg + " "
    exec cmd
    
# -------------------------------------------------------------------------
# store global namespace
# swap in a new exception handler
# change interactive prompts

namespace = sys.modules['__main__'].__dict__
sys.excepthook = trap
sys.ps1 = "> "
sys.ps2 = ". "

# should now go interactive if launched with "python -i"
# unless -q switch used

if quitflag > 0: sys.exit()
