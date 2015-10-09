# simple test of rasmol tool
# requires files/dump.peptide and files/peptide.pdb
# creates tmp*.pdb and tmp*.gif

d = dump("files/dump.peptide")
d.unwrap()
p = pdbfile("files/peptide",d)
r = rasmol(p)
r.file = "tmp"

print "kill image window when ready to contine ..."
r.show(0)
r.all()

# change RasMol settings, run all() with those settings

print "change RasMol settings as desired, then type 'quit' ..."
r.run(0)
r.all("tmp.rasmol")

print "all done ... type CTRL-D to exit Pizza.py"
