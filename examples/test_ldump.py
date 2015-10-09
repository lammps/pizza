# simple test of dump tool
# requires files/dump.lines
# uses gl and vcr tools to visualize line segments

d = dump("files/dump.lines")
ld = ldump("files/dump.lines")
ld.map(1,"id",2,"type",6,"end1x",7,"end1y",8,"end2x",9,"end2y")
d.extra(ld)
g = gl(d)
g.arad(0,0.2)
g.lrad(1,5)
g.lcol(1,"green")
v = vcr(g)

print "all done ... type CTRL-D to exit Pizza.py"
