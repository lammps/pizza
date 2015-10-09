# simple test of cfg tool
# requires files/dump.peptide.*
# creates tmp*.cfg

d = dump("files/dump.peptide.*")
d.sort()
x = cfg(d)
x.one()
x.many()
x.single(0,"tmp.single")

print "all done ... type CTRL-D to exit Pizza.py"
