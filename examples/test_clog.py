# simple test of clog tool
# requires files/log.ccell
# creates tmp.clog and tmp.clog.two

c = log("files/log.ccell")

print "# of vectors =",c.nvec
print "length of vectors =",c.nlen
print "names of vectors =",c.names

time,a,b = c.get("Step","prey","predator")
print a,b
c.write("tmp.clog")
c.write("tmp.clog.two","Step","prey")

print "all done ... type CTRL-D to exit Pizza.py"
