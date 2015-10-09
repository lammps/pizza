# simple test of vec tool
# requires files/vec.txt
# creates tmp.vec and tmp.vec.two

v = vec("files/vec.txt")

print "# of vectors =",v.nvec
print "length of vectors =",v.nlen
print "names of vectors =",v.names

time,temp,press = v.get(1,"col2",6)
print temp,press
v.write("tmp.vec")
v.write("tmp.vec.two","col1",3)

print "all done ... type CTRL-D to exit Pizza.py"
