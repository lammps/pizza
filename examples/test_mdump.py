# simple test of mdump tool
# requires files/mesh.grain
# creates tmp.* files

m = mdump("files/mesh.grain")

m.tselect.none()
m.tselect.one(200)
m.tselect.all()
m.tselect.skip(2)
m.tselect.all()
m.tselect.test("$t >= 100 and $t <= 200")
m.delete()

print "Time",m.time()

m.map(2,"spin")
m.etype = "spin"

flag = 0
while 1:
  index,time,flag = m.iterator(flag)
  if flag == -1: break
  time,box,atoms,bonds,tris,lines = m.viz(index)
  colors = [tri[1] for tri in tris]
  print time,colors

m = dump("files/mesh.grain",0)
while 1:
  time = m.next()
  if time < 0: break

print "Incrementally read snaps =",m.nsnaps

print "all done ... type CTRL-D to exit Pizza.py"
