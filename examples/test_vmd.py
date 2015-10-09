# simple test of vmd tool
# requires files/dump.peptide.only and files/dump.bond
# uses gl and vcr tools to visualize peptide molecule with bonds

v = vmd()
v('menu main off')
v.rep('VDW')
v.new('files/peptide.pdb','pdb')
v.flush()

print "all done ... type CTRL-D to exit Pizza.py"
