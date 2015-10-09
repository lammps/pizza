# compute pairwise energy between 2 groups of atoms

# Syntax: group_energy.py data.file dump.file1 dump.file2 ...
# Author: Paul Crozier (Sandia)

# return distance sq between 2 atoms with PBC

def distance(box,x1,y1,z1,x2,y2,z2):
  
  delx = x2 - x1
  dely = y2 - y1
  delz = z2 - z1
  
  xprd = box[3] - box[0]
  yprd = box[4] - box[1]
  zprd = box[5] - box[2]

  if abs(delx) > 0.5*xprd:
    if delx < 0.0:
      delx += xprd
    else:
      delx -= xprd
  if abs(dely) > 0.5*yprd:
    if dely < 0.0:
      dely += yprd
    else:
      dely -= yprd
  if abs(delz) > 0.5*zprd:
    if delz < 0.0:
      delz += zprd
    else:
      delz -= zprd
      
  distsq = delx*delx + dely*dely + delz*delz
  return distsq

# main script

if len(argv) < 3:
  raise StandardError,"group_energy.py data.file dump.file1 dump.file2 ..."

dt = data(argv[1])				# data file
q = dt.get("Atoms",4)

files = ' '.join(argv[2:])		        # dump files
d = dump(files,0)
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

p = pair("lj/charmm/coul/charmm")
p.coeff(dt)

cut1 = 8.0                                      # potential cutoffs
cut2 = 10.0
cut3 = 8.0
cut4 = 10.0
p.init(cut1,cut2,cut3,cut4)

maxcut = cut1
if cut2 > maxcut: maxcut = cut2
if cut3 > maxcut: maxcut = cut3
if cut4 > maxcut: maxcut = cut4
maxcut_sq = maxcut*maxcut

while 1:
  time = d.next()
  if time < 0: break
  d.unscale(time)

  box = (d.snaps[0].xlo,d.snaps[0].ylo,d.snaps[0].zlo,
         d.snaps[0].xhi,d.snaps[0].yhi,d.snaps[0].zhi)
  d.aselect.all(time)
  d.aselect.test("$id >= 14306 and $id <= 14516",time)          # 1st group
  id1,type1,x1,y1,z1 = d.vecs(time,"id","type","x","y","z")
  d.aselect.all(time)                                           # 2nd group
  d.aselect.test("$id >= 1 and $id <= 7243 or $id >= 7274 and $id <= 14283",
                 time)
  id2,type2,x2,y2,z2 = d.vecs(time,"id","type","x","y","z")
  id1 = map(int,id1)
  id2 = map(int,id2)
  type1 = map(int,type1)
  type2 = map(int,type2)
  n1 = len(type1)
  n2 = len(type2)
  for i in xrange(n1): 
    id1[i] -= 1
    type1[i] -= 1
  for i in xrange(n2): 
    id2[i] -= 1
    type2[i] -= 1

  e_coul_sum = 0.0
  e_vdwl_sum = 0.0
  for i in xrange(n1):
    typei = type1[i]
    qi = q[id1[i]]
    for j in xrange(n2):
      rsq = distance(box,x1[i],y1[i],z1[i],x2[j],y2[j],z2[j])
      if rsq < maxcut_sq:
        eng_coul,eng_vdwl = p.single(rsq,typei,type2[j],qi,q[id2[j]])
        e_coul_sum += eng_coul     
        e_vdwl_sum += eng_vdwl    
  print "eng_coul = %g at timestep %d" % (e_coul_sum,time)
  print "eng_vdwl = %g at timestep %d" % (e_vdwl_sum,time)

  d.tselect.none()
  d.tselect.one(time)
  d.delete()
