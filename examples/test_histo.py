# simple test of histo tool
# requires files/dump.kinase

d = dump("files/dump.kinase")
h = histo(d)
x,y = h.compute('x',25)
g = gnu()
g.xtitle("X position")
g.ytitle("Particle Count")
g.title("Histogram of Particle Density")
g.plot(x,y)

print "all done ... type CTRL-D to exit Pizza.py"
