# simple test of gl tool
# requires files/dump.kinase
# creates tmp*.png

d = dump("files/dump.kinase")
g = gl(d)

g.bg("white")
g.rotate(60,130)
g.box(1)
g.q(10)
g.file = "tmp"

g.show(0)
g.all()

from vizinfo import colors

g.acol([1,4,6,8,9],["gray","red","blue","green","yellow"])
g.arad(range(9),0.3)
#g.label(0.2,0.4,'h',15,"red","test label #1")
#g.label(-0.2,-0.4,'h',15,"yellow","test label #2")

g.show(0)
g.pan(60,130,1,60,30,0.5)
g.all(0,10,0)

print "all done ... type CTRL-D to exit Pizza.py"
