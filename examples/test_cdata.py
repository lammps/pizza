# simple test of cdata tool
# creates tmp.cdata

c = cdata()

c.box("box",0,0,0,10,5,5)
c.sphere("nucleus",7,2,2,1)
c.cap("organelle",'x',2.5,2.5,1,1.5,3.5)
c.q("nucleus",5)
c.union("interior","nucleus","organelle")

c.surf("nuc","nucleus")
c.surfselect("nuchalf","nuc","$z < 2.0")

c.part("A",100,"box","interior")
c.part("B",100,"nucleus")
c.part2d("C",100,"organelle")

c.write("tmp.cdata","box","nucleus","organelle","A","B","C")

c.lbox("linebox",0,0,0,10,5,5)
c.unselect()
c.select("A","B","C","linebox","organelle","nuchalf")

#s = svg(c)
#s.rotate(0,0)
#s.lrad(0,0.1)
#s.lcol(0,"yellow")
#s.show(0)

g = gl(c)
v = vcr(g)
 
print "all done ... type CTRL-D to exit Pizza.py"
