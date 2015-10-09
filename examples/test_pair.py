# simple test of pair tool
# requires files/data.rhodo

p = pair("lj/charmm/coul/charmm")
d = data("files/data.rhodo")
p.coeff(d)
p.init(8.0,10.0)
ev,ec = p.single(5.0,1,2,0.5,-0.5)
print "Energies",ev,ec

print "all done ... type CTRL-D to exit Pizza.py"
