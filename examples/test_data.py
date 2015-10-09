# simple test of data tool
# requires files/data.micelle and dump.micelle
# creates tmp.data

d = data("files/data.micelle")
d.map(1,"id",3,"type",4,"x",5,"y",6,"z")
coeffs = d.get("Masses")
print "Masses",coeffs
x = d.get("Atoms",4)
print "X of 1st atom",x[0]

d.title = "New LAMMPS data file"

natoms = d.headers["atoms"]
vec = range(1,natoms+1)
vec.reverse()
d.replace("Atoms",1,vec)

dm = dump("files/dump.micelle")
d.newxyz(dm,1000)

flag = 0
while 1:
  index,time,flag = d.iterator(flag)
  if flag == -1: break
  time,box,atoms,bonds,tris,lines= d.viz(index)
  
d.write("tmp.data")

print "all done ... type CTRL-D to exit Pizza.py"
