# simple test of raster tool
# requires files/dump.kinase
# creates tmp*.png

d = dump("files/dump.kinase")
r = raster(d)

r.bg("white")
r.rotate(60,130)
r.box(1)
r.file = "tmp"

print "kill image window when ready to contine ..."
r.show(0)
r.all()
a1 = animate("tmp0*png")

from vizinfo import colors

r.acol([1,4,6,8,9],["gray","red","blue","green","yellow"])
r.arad(range(9),0.3)
r.label(0.2,0.4,'h',15,"red","test label #1")
r.label(-0.2,-0.4,'h',15,"yellow","test label #2")

print "kill image window when ready to contine ..."
r.show(0)
r.pan(60,130,1,60,30,0.5)
r.all(0,10,0)
a2 = animate("tmp0*png")

print "all done ... type CTRL-D to exit Pizza.py"
