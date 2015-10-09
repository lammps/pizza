# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# ensight tool

oneline = "Convert LAMMPS snapshots or meshes to Ensight format"

docstr = """
e = ensight(d)	     d = object with atoms or elements (dump,data,mdump)
e.change = 1         set to 1 if element nodal xyz change with time (def = 0)
e.maxtype = 10       max particle type, set if query to data will be bad

e.one()
e.one("new")
e.one("cns","Centro","eng","Energy")
e.one("new","cns","Centro","eng","Energy")
                     write all snapshots as an Ensight data set
                     Ensight header file = tmp.case (no 1st arg) or new.case
                     Ensight coord file = tmp.xyz or new.xyz
                     additional pairs of args create auxiliary files:
                       tmp.cns, tmp.eng or new.cns, new.eng
                     cns,eng = column name in dump file and file name suffix
                     Centro,Energy = Ensight name for the variable

e.increment()        same args as one(), but process dump out-of-core

e.many()             same args as one(), but create multiple Ensight files
                     tmp0000.xyz, tmp0001.xyz, etc
                     new0000.cns, new0001.cns, etc
                     new0000.eng, new0001.eng, etc

e.single(N)          same args as one() prepended by N, but write a single snap
"""

# History
#   10/06, Steve Plimpton (SNL): original version

# ToDo list
#   binary files
#   create vector or tensor variable files, not just scalar
#     via pair of args like ["vx","vy","vz"],"vel"

# Variables
#   data = data file to read from
#   which = 0 for particles, 1 for elements
#   change = 0 for unchanging mesh coords, 1 for changing mesh coords (def = 0)

# Imports and external programs

import sys, types

# Class definition

class ensight:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.change = 0
    self.maxtype = 0
    self.data = data
    if type(data) is types.InstanceType and ".dump" in str(data.__class__):
      self.which = 0
    elif type(data) is types.InstanceType and ".data" in str(data.__class__):
      self.which = 0
    elif type(data) is types.InstanceType and ".mdump" in str(data.__class__):
      self.which = 1
    elif type(data) is types.InstanceType and ".cdata" in str(data.__class__):
      self.which = 1
    else:
      raise StandardError,"unrecognized object passed to ensight"
    
  # --------------------------------------------------------------------

  def one(self,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()
    
    # write Ensight *.case header file

    f = open("%s.case" % root,"w")
    times = self.data.time()
    self.case_file(f,root,pairs,0,len(times),times)
    f.close()

    # open additional files

    f = open(root + ".xyz","w")
    vfiles = []
    for pair in pairs: vfiles.append(open(root + "." + pair[0],"w"))

    # loop over snapshots
    # write coords into xyz file, variables into their files

    first = 1
    n = flag = etype = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break

      if self.which == 0:
        print >>f,"BEGIN TIME STEP"
        time,box,atoms,bonds,tris,lines = self.data.viz(which)
        self.coord_file_atoms(f,box,atoms)
        print >>f,"END TIME STEP"
      elif self.change == 0 and first:
        print >>f,"BEGIN TIME STEP"
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        first = 0
        print >>f,"END TIME STEP"
      elif self.change:
        print >>f,"BEGIN TIME STEP"
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        print >>f,"END TIME STEP"

      for i in range(len(pairs)):
        print >>vfiles[i],"BEGIN TIME STEP"
        values = self.data.vecs(time,pairs[i][0])
        if self.which == 0:
          self.variable_file_atoms(vfiles[i],pairs[i][1],atoms,values)
        else:
          self.variable_file_elements(vfiles[i],pairs[i][1],etype,values)
        print >>vfiles[i],"END TIME STEP"
        
      print time,
      sys.stdout.flush()
      n += 1

    # close additional files

    f.close()
    for f in vfiles: f.close()

    print "\nwrote %s snapshots in Ensight format" % n

  # --------------------------------------------------------------------

  def increment(self,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()

    # open additional files

    f = open(root + ".xyz","w")
    vfiles = []
    for pair in pairs: vfiles.append(open(root + "." + pair[0],"w"))

    # loop over snapshots
    # write coords into xyz file, variables into their files

    times = []
    first = 1
    n = etype = 0
    while 1:
      time = self.data.next()
      if time == -1: break
      times.append(time)
      self.data.tselect.one(time)
      self.data.delete()

      if self.which == 0:
        print >>f,"BEGIN TIME STEP"
        time,box,atoms,bonds,tris,lines = self.data.viz(0)
        self.coord_file_atoms(f,box,atoms)
        print >>f,"END TIME STEP"
      elif self.change == 0 and first:
        print >>f,"BEGIN TIME STEP"
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(0)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        first = 0
        print >>f,"END TIME STEP"
      elif self.change:
        print >>f,"BEGIN TIME STEP"
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(0)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        print >>f,"END TIME STEP"

      for i in range(len(pairs)):
        print >>vfiles[i],"BEGIN TIME STEP"
        values = self.data.vecs(time,pairs[i][0])
        if self.which == 0:
          self.variable_file_atoms(vfiles[i],pairs[i][1],atoms,values)
        else:
          self.variable_file_elements(vfiles[i],pairs[i][1],etype,values)
        print >>vfiles[i],"END TIME STEP"
        
      print time,
      sys.stdout.flush()
      n += 1

    # close additional files

    f.close()
    for f in vfiles: f.close()

    # write Ensight *.case header file now that know all timesteps

    f = open("%s.case" % root,"w")
    self.case_file(f,root,pairs,0,len(times),times)
    f.close()

    print "\nwrote %s snapshots in Ensight format" % n

  # --------------------------------------------------------------------

  def many(self,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()

    # write Ensight *.case header file

    f = open("%s.case" % root,"w")
    times = self.data.time()
    self.case_file(f,root,pairs,1,len(times),times)
    f.close()

    # loop over snapshots
    # generate unique filenames
    # write coords into one xyz file per snapshot, variables into their files

    first = 1
    n = flag = etype = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break

      files = []
      if n < 10:
        file = root + "000" + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + "000" + str(n) + "." + pair[0])
      elif n < 100:
        file = root + "00" + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + "00" + str(n) + "." + pair[0])
      elif n < 1000:
        file = root + "0" + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + "0" + str(n) + "." + pair[0])
      else:
        file = root + str(n) + ".xyz"
	for pair in pairs:
	  files.append(root + str(n) + "." + pair[0])

      if self.which == 0:
        f = open(file,"w")
        time,box,atoms,bonds,tris,lines = self.data.viz(which)
        self.coord_file_atoms(f,box,atoms)
        f.close()
      elif self.change == 0 and first:
        f = open(root + ".xyz","w")
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        first = 0
        f.close()
      elif self.change:
        f = open(file,"w")
        time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
        self.coord_file_elements(f,box,nodes,elements)
        etype = len(elements[0])
        f.close()

      for i in range(len(pairs)):
        values = self.data.vecs(time,pairs[i][0])
        f = open(files[i],"w")
        if self.which == 0:
          self.variable_file_atoms(f,pairs[i][1],atoms,values)
        else:
          self.variable_file_elements(f,pairs[i][1],etype,values)
	f.close()

      print time,
      sys.stdout.flush()
      n += 1
  
    print "\nwrote %s snapshots in Ensight format" % n

  # --------------------------------------------------------------------

  def single(self,time,*args):
    if len(args) % 2 == 0: root = "tmp"
    else: 
      root = args[0]
      args = args[1:]

    pairs = []
    for i in range(0,len(args),2): pairs.append([args[i],args[i+1]])

    # max # of types for all steps in Ensight files

    if self.which == 0 and self.maxtype == 0:
      self.maxtype = self.data.maxtype()

    # write Ensight *.case header file

    f = open("%s.case" % root,"w")
    self.case_file(f,root,pairs,0,1,[time])
    f.close()

    # write coords into xyz file, variables into their files

    which = self.data.findtime(time)
    etype = 0
    
    f = open(root + ".xyz","w")
    if self.which == 0:
      time,box,atoms,bonds,tris,lines = self.data.viz(which)
      self.coord_file_atoms(f,box,atoms)
    else:
      time,box,nodes,elements,nvalues,evalues = self.data.mviz(which)
      self.coord_file_elements(f,box,nodes,elements)
      etype = len(elements[0])
    f.close()
    
    for i in range(len(pairs)):
      values = self.data.vecs(time,pairs[i][0])
      f = open(root + "." + pairs[i][0],"w")
      if self.which == 0:
        self.variable_file_atoms(f,pairs[i][1],atoms,values)
      else:
        self.variable_file_elements(f,pairs[i][1],etype,values)
      f.close()
      
  # --------------------------------------------------------------------
  # write Ensight case file

  def case_file(self,f,root,pairs,multifile,nsnaps,times):
    print >>f,"# Ensight case file\n"
    print >>f,"FORMAT\ntype: ensight gold\n"
    
    if self.which == 0:
      if multifile:
#        print >>f,"GEOMETRY\nmodel: %s****.xyz change_coords_only\n" % root
        print >>f,"GEOMETRY\nmodel: %s****.xyz\n" % root
      else:
#        print >>f,"GEOMETRY\nmodel: 1 1 %s.xyz change_coords_only\n" % root
        print >>f,"GEOMETRY\nmodel: 1 1 %s.xyz\n" % root
    else:
      if self.change == 0:
        print >>f,"GEOMETRY\nmodel: %s.xyz\n" % root
      elif multifile:
        print >>f,"GEOMETRY\nmodel: %s****.xyz\n" % root
      else:
        print >>f,"GEOMETRY\nmodel: 1 1 %s.xyz\n" % root

    if len(pairs):
      print >>f,"VARIABLE"
      for pair in pairs:
        if self.which == 0:
          if multifile:
            print >>f,"scalar per node: %s %s****.%s" % (pair[1],root,pair[0])
          else:
            print >>f,"scalar per node: 1 1 %s %s.%s" % (pair[1],root,pair[0])
        else:
          if multifile:
            print >>f,"scalar per element: %s %s****.%s" % (pair[1],root,pair[0])
          else:
            print >>f,"scalar per element: 1 1 %s %s.%s" % (pair[1],root,pair[0])
      print >>f

    print >>f,"TIME"
    print >>f,"time set: 1"
    print >>f,"number of steps:",nsnaps
    print >>f,"filename start number: 0"
    print >>f,"filename increment: 1"
    print >>f,"time values:"
    for i in range(nsnaps):
      print >>f,times[i],
      if i % 10 == 9: print >>f
    print >>f
    print >>f
    
    if not multifile:
      print >>f,"FILE"
      print >>f,"file set: 1"
      print >>f,"number of steps:",nsnaps

  # --------------------------------------------------------------------
  # write Ensight coordinates for atoms
  # partition into "parts"
  # one part = coords for all atoms of a single type

  def coord_file_atoms(self,f,box,atoms):
    print >>f,"Particle geometry\nfor a collection of atoms"
    print >>f,"node id given"
    print >>f,"element id off"
    print >>f,"extents"
    print >>f,"%12.5e%12.5e" % (box[0],box[3])
    print >>f,"%12.5e%12.5e" % (box[1],box[4])
    print >>f,"%12.5e%12.5e" % (box[2],box[5])

    for type in xrange(1,self.maxtype+1):
      print >>f,"part"
      print >>f,"%10d" % type
      print >>f,"type",type
      print >>f,"coordinates"
      group = [atom for atom in atoms if int(atom[1]) == type]
      print >>f,"%10d" % len(group)
      for atom in group: print >>f,"%10d" % int(atom[0])
      for atom in group: print >>f,"%12.5e" % atom[2]
      for atom in group: print >>f,"%12.5e" % atom[3]
      for atom in group: print >>f,"%12.5e" % atom[4]
      print >>f,"point"
      print >>f,"%10d" % len(group)
      for i in xrange(1,len(group)+1): print >>f,"%10d" % i

  # --------------------------------------------------------------------
  # write Ensight coordinates for elements

  def coord_file_elements(self,f,box,nodes,elements):
    print >>f,"Element geometry\nfor a collection of elements"
    print >>f,"node id given"
    print >>f,"element id given"
    print >>f,"extents"
    print >>f,"%12.5e%12.5e" % (box[0],box[3])
    print >>f,"%12.5e%12.5e" % (box[1],box[4])
    print >>f,"%12.5e%12.5e" % (box[2],box[5])

    print >>f,"part"
    print >>f,"%10d" % 1
    print >>f,"all elements"
    print >>f,"coordinates"
    print >>f,"%10d" % len(nodes)
    for node in nodes: print >>f,"%10d" % int(node[0])
    for node in nodes: print >>f,"%12.5e" % node[2]
    for node in nodes: print >>f,"%12.5e" % node[3]
    for node in nodes: print >>f,"%12.5e" % node[4]

    if len(elements[0]) == 5: print >>f,"tria3"
    elif len(elements[0]) == 6: print >>f,"tetra4"
    else: raise StandardError,"unrecognized element type"
    print >>f,"%10d" % len(elements)

    for element in elements: print >>f,"%10d" % int(element[0])
    if len(elements[0]) == 5:
      for element in elements:
        print >>f,"%10d%10d%10d" % \
              (int(element[2]),int(element[3]),int(element[4]))
    elif len(elements[0]) == 6:
      for element in elements:
        print >>f,"%10d%10d%10d%10d" % \
              (int(element[2]),int(element[3]),int(element[4]),int(element[5]))

  # --------------------------------------------------------------------
  # write Ensight variable values for atoms
  # partition into "parts"
  # one part = values for all atoms of a single type

  def variable_file_atoms(self,f,name,atoms,values):
    print >>f,"Particle %s" % name
    for type in xrange(1,self.maxtype+1):
      print >>f,"part"
      print >>f,"%10d" % type
      print >>f,"coordinates"
      group = [values[i] for i in xrange(len(atoms))
               if int(atoms[i][1]) == type]
      for value in group: print >>f,"%12.5e" % value

  # --------------------------------------------------------------------
  # write Ensight variable values for elements

  def variable_file_elements(self,f,name,etype,values):
    print >>f,"Element %s" % name
    print >>f,"part"
    print >>f,"%10d" % 1
    if etype == 5: print >>f,"tria3"
    elif etype == 6: print >>f,"tetra4"
    for value in values: print >>f,"%12.5e" % value
