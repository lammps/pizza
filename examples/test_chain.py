# simple test of chain tool
# creates tmp.data.chain file

c = chain(500,0.7,1,1,2)
c.seed = 54321
c.build(25,10)

c.mtype = 2
c.btype = 2
c.blen = 1.5
c.dmin = 1.2
c.id = "end1"
c.build(10,25)

c.write("tmp.data.chain")

print "all done ... type CTRL-D to exit Pizza.py"
