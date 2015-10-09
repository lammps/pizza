# simple test of patch tool
# creates tmp.data.patch file

p = patch(0.5)
p.seed = 54321
p.build(100,"hex2",1,2,3)
p.build(50,"tri5",4,5)
p.write("tmp.data.patch")

print "all done ... type CTRL-D to exit Pizza.py"
