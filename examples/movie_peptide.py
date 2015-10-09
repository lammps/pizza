# movie of solvated peptide data

d = dump("dump.peptide")
d.unwrap()
p = pdb("peptide",d)
r.file = "peptide"
r = rasmol(p)
r.all()
