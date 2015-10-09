# simple test of log tool
# requires files/log.obstacle
# creates tmp.log and tmp.log.two

lg = log("files/log.obstacle")

print "# of vectors =",lg.nvec
print "length of vectors =",lg.nlen
print "names of vectors =",lg.names

time,temp,press = lg.get("Step","Temp","Press")
print temp,press
lg.write("tmp.log")
lg.write("tmp.log.two","Step","E_pair")

print "all done ... type CTRL-D to exit Pizza.py"
