# simple test of vcr tool

d = dump("files/dump.micelle")
dt = data("files/data.micelle")
d.extra(dt)
g = gl(d)
g.rotate(0,270)
v = vcr(g)
v.q(10)
v.box()
v.axis()
v.clipxlo(0.2)
v.clipxhi(0.5)
v.play()

print "all done ... type CTRL-D to exit Pizza.py"
