# simple test of xyz tool
# requires files/dump.micelle.*
# creates tmp*.case, etc

d = dump("files/dump.micelle")
e = ensight(d)
e.one()
e.many()
e.single(0)

print "all done ... type CTRL-D to exit Pizza.py"
