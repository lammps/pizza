# simple test of dump tool
# requires files/dump.peptide and files/dump.peptide.*
# creates tmp.* files

d = dump("files/dump.peptide")

d.tselect.none()
d.tselect.one(2000)
d.tselect.all()
d.tselect.skip(2)
d.tselect.all()
d.tselect.test("$t >= 4000 and $t <= 6000")
d.delete()

d.aselect.all()
d.aselect.test("$id <= 10")

d.write("tmp.dump")
d.scatter("tmp")

d = dump("files/dump.peptide")
d.scale()
d.unscale()
d.unwrap()
d.wrap()
d.set("$xyz = $x + $y + $z")
d.spread("x",100,"color")
d.clone(0,"color")

print "Time",d.time()
color = d.atom(1,"color")
print "Color of atom 1",color
d.aselect.test("$id <= 10")
color = d.vecs(1000,"color")
print "Color of 1st 10 atoms in step 1000",color

d.atype = "color"

flag = 0
while 1:
  index,time,flag = d.iterator(flag)
  if flag == -1: break
  time,box,atoms,bonds,tris,lines = d.viz(index)
  colors = [atom[1] for atom in atoms]
  print time,colors

d = dump("files/dump.peptide.*",0)
while 1:
  time = d.next()
  if time < 0: break

print "Incrementally read snaps =",d.nsnaps

print "all done ... type CTRL-D to exit Pizza.py"
