# simple test of vtk tool
# requires files/dump.peptide.*
# creates tmp*.vtk

d = dump("files/dump.peptide.*")
v = vtk(d)
v.one()
v.many()
v.single(0,"tmp.single")

print "all done ... type CTRL-D to exit Pizza.py"
